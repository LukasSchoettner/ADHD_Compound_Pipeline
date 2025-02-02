#!/usr/bin/env python
import os
import shutil
import subprocess
import argparse
import pandas as pd
import traceback
import psycopg2
import re
import concurrent.futures
import logging

from Bio.PDB import PDBParser, PDBIO, Select
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("md_simulation_pipeline.log"),
        logging.StreamHandler()
    ]
)

class KeepAll(Select):
    def accept_residue(self, residue):
        return True

###############################################################################
# 1. Database Queries
###############################################################################


def insert_md_record(db_connection_string: str, record: dict):
    """
    Insert a single MD simulation record into the 'molecular_dynamics' table using psycopg2.
    """
    insert_sql = """
        INSERT INTO molecular_dynamics (
            docking_id, ligand_cid, deg_id, binding_energy,
            simulation_time, temperature, pressure, simulation_type,
            status, created_at, updated_at
        ) VALUES (
            %(docking_id)s, %(ligand_cid)s, %(deg_id)s, %(binding_energy)s,
            %(simulation_time)s, %(temperature)s, %(pressure)s, %(simulation_type)s,
            %(status)s, NOW(), NOW()
        )
    """

    try:
        with psycopg2.connect(db_connection_string) as conn:
            with conn.cursor() as cur:
                cur.execute(insert_sql, record)
        logging.info(f"[DB] Inserted MD record for docking_id={record['docking_id']}.")
    except Exception as e:
        traceback.print_exc()
        logging.error(f"[DB] Failed to insert MD record for docking_id={record.get('docking_id', 'unknown')}: {e}")

###############################################################################
# 2. Atom Counting and Validation
###############################################################################


def count_atoms_in_itp(itp_file: Path) -> int:
    """
    Count the number of atoms in an .itp or .top file by parsing the [ atoms ] section.

    Args:
        itp_file (Path): Path to the .itp or .top file.

    Returns:
        int: Number of atoms.
    """
    count = 0
    with itp_file.open("r") as f:
        in_atoms_section = False
        for line in f:
            line = line.strip()
            if line.startswith("[ atoms ]"):
                in_atoms_section = True
                continue
            if in_atoms_section:
                if line.startswith("[") and line.endswith("]"):
                    break
                if line and not line.startswith(";"):
                    count += 1
    logging.debug(f"[AtomCount] {itp_file} has {count} atoms.")
    return count


def generate_protein_topology(protein_pdb: Path, force_field: str, water_model: str, output_dir: Path):
    """
    Generate protein topology using pdb2gmx.

    Args:
        protein_pdb (Path): Path to the cleaned protein PDB file.
        force_field (str): Force field to use (e.g., 'amber99sb-ildn').
        water_model (str): Water model to use (e.g., 'tip3p').
        output_dir (Path): Directory to save the generated topology and GRO files.

    Returns:
        tuple: Paths to the processed protein GRO and topology files.
    """
    protein_gro = output_dir / f"{protein_pdb.stem}_processed.gro"
    protein_top = output_dir / f"{protein_pdb.stem}.top"

    cmd = [
        "gmx", "pdb2gmx",
        "-f", str(protein_pdb),
        "-o", str(protein_gro),
        "-p", str(protein_top),
        "-ff", force_field,
        "-water", water_model
    ]

    logging.info(f"[Topology] Running pdb2gmx for {protein_pdb.name} with force field {force_field} and water model {water_model}.")
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        logging.info(f"[Topology] Generated protein GRO: {protein_gro}")
        logging.info(f"[Topology] Generated protein TOP: {protein_top}")
    except subprocess.CalledProcessError as e:
        logging.error(f"[Topology] pdb2gmx failed: {e.stderr.decode()}")
        raise RuntimeError(f"pdb2gmx failed: {e.stderr.decode()}")

    return protein_gro, protein_top


def validate_topology_gro_consistency(gro_file: Path, top_file: Path, molecule_atom_counts: dict):
    try:
        # Read the .gro file
        with gro_file.open("r") as f:
            lines = f.readlines()
            num_atoms_gro = int(lines[1].strip())
            logging.debug(f"[Validate] .gro file atom count: {num_atoms_gro}")

        # Read the [ molecules ] section in the .top file
        with top_file.open("r") as f:
            top_lines = f.readlines()

        molecule_section = False
        molecule_counts = {}
        total_atoms_top = 0
        for line in top_lines:
            line = line.strip()
            if line.startswith("[ molecules ]"):
                molecule_section = True
                continue
            if molecule_section:
                if not line or line.startswith(";"):
                    continue
                if line.startswith("[") and line.endswith("]"):
                    break
                try:
                    molname, count_str = re.split(r'\s+', line)
                    count = int(count_str)
                    molecule_counts[molname.strip()] = count
                except ValueError:
                    raise ValueError(f"Invalid line in topology file: '{line}'")

        # Log the molecules found and their counts
        logging.debug(f"[Validate] Molecule counts from [ molecules ] section: {molecule_counts}")
        logging.debug(f"[Validate] Molecule atom counts dictionary: {molecule_atom_counts}")

        # Calculate total atoms based on molecule counts
        for mol, cnt in molecule_counts.items():
            mol_key = mol.strip()
            if mol_key not in molecule_atom_counts:
                logging.warning(f"[Validate] Molecule '{mol_key}' not found in molecule_atom_counts. Skipping.")
                continue
            total_atoms_top += molecule_atom_counts[mol_key] * cnt

        logging.debug(f"[Validate] Topology file total atom count: {total_atoms_top}")
        logging.debug(f"[Validate] Molecule counts: {molecule_counts}")

        if num_atoms_gro != total_atoms_top:
            raise ValueError(f"Mismatch in .gro ({num_atoms_gro}) vs .top ({total_atoms_top}) atom counts.")
        else:
            logging.info("[Validate] Validation successful: .gro and .top atom counts match.")
    except Exception as e:
        logging.error(f"[Validate] Validation failed: {e}")
        raise

###############################################################################
# 3. System Preparation Functions
###############################################################################


def combine_protein_ligand(protein_gro: Path, protein_top: Path,
                           ligand_gro: Path, ligand_itp: Path,
                           out_dir: Path,
                           protein_name: str,
                           ligand_name: str,
                           system_name: str) -> tuple:
    """
    Combine protein and ligand GRO files and update the topology accordingly.

    Args:
        protein_gro (Path): Path to the protein GRO file.
        protein_top (Path): Path to the protein topology (.top) file.
        ligand_gro (Path): Path to the ligand GRO file.
        ligand_itp (Path): Path to the ligand ITP file.
        out_dir (Path): Output directory for combined files.
        protein_name (str): Molecule name for the protein in topology.
        ligand_name (str): Molecule name for the ligand in topology.
        system_name (str): Description of the system.

    Returns:
        tuple: Paths to the combined GRO and topology files.
    """
    try:
        # Ensure output directory exists
        out_dir.mkdir(parents=True, exist_ok=True)
        logging.info(f"[Combine] Starting to combine protein and ligand in {out_dir}")

        # Check if required files exist
        for file in [protein_gro, protein_top, ligand_gro, ligand_itp]:
            if not file.exists():
                logging.error(f"[Combine] Required file {file} does not exist.")
                raise FileNotFoundError(f"Required file {file} does not exist.")

        # Read protein GRO
        with protein_gro.open("r") as pf:
            protein_lines = pf.readlines()
        num_protein_atoms = int(protein_lines[1].strip())
        protein_atom_lines = protein_lines[2:-1]  # Exclude header and 'END'

        # Read ligand GRO
        with ligand_gro.open("r") as lf:
            ligand_lines = lf.readlines()
        num_ligand_atoms = int(ligand_lines[1].strip())
        ligand_atom_lines = ligand_lines[2:-1]  # Exclude header and 'END'

        # Total atoms
        total_atoms = num_protein_atoms + num_ligand_atoms

        # Log atom counts
        logging.info(f"[Combine] Protein atoms: {num_protein_atoms}")
        logging.info(f"[Combine] Ligand atoms: {num_ligand_atoms}")
        logging.info(f"[Combine] Total atoms: {total_atoms}")

        # Combine GRO files
        combined_gro = out_dir / "combined.gro"
        with combined_gro.open("w") as fout:
            fout.write(f"{system_name}\n")
            fout.write(f"{total_atoms}\n")
            for line in protein_atom_lines:
                fout.write(line)
            for line in ligand_atom_lines:
                fout.write(line)
            # Write dummy box dimensions; they will be redefined by 'editconf'
            fout.write(" 0.000 0.000 0.000\n")
        logging.info(f"[Combine] Combined GRO file created at {combined_gro}")

        # Combine topology files
        combined_top = out_dir / "combined.top"
        with combined_top.open("w") as fout:
            # 1. Include the forcefield from protein.top
            with protein_top.open("r") as pt:
                for line in pt:
                    fout.write(line)

            # 2. Extract and include [ atomtypes ] from ligand.itp
            with ligand_itp.open("r") as lt:
                lines = lt.readlines()

            # Find [ atomtypes ] section in ligand.itp
            atomtypes_start = None
            atomtypes_end = None
            for i, line in enumerate(lines):
                if line.strip().startswith("[ atomtypes ]"):
                    atomtypes_start = i
                    break

            if atomtypes_start is not None:
                # Find the end of [ atomtypes ] section
                for j in range(atomtypes_start + 1, len(lines)):
                    if lines[j].startswith("[") and not lines[j].strip().startswith("[ atomtypes ]"):
                        atomtypes_end = j
                        break
                if atomtypes_end is None:
                    atomtypes_end = len(lines)
                # Write [ atomtypes ] section to combined.top
                fout.write("\n")
                fout.write("\n".join(lines[atomtypes_start:atomtypes_end]))
                fout.write("\n")
                logging.info(f"[Combine] Included [ atomtypes ] from {ligand_itp}")

            # 3. Include the rest of ligand.itp excluding [ atomtypes ]
            with ligand_itp.open("r") as lt:
                in_atomtypes = False
                for line in lt:
                    if line.strip().startswith("[ atomtypes ]"):
                        in_atomtypes = True
                        continue  # Skip writing [ atomtypes ] section
                    elif line.strip().startswith("[ moleculetype ]"):
                        in_atomtypes = False
                    elif line.strip().startswith("[") and not line.strip().startswith("[ moleculetype ]"):
                        in_atomtypes = False
                    if not in_atomtypes:
                        fout.write(line)
            logging.info(f"[Combine] Included [ moleculetype ] and [ atoms ] from {ligand_itp}")

            # 4. Write [ system ] and [ molecules ] sections
            fout.write(f"\n[ system ]\n{system_name}\n\n")
            fout.write(f"[ molecules ]\n{protein_name} 1\n{ligand_name} 1\n")
        logging.info(f"[Combine] Combined topology file created at {combined_top}")

        # 5. Copy ligand.itp to out_dir to ensure it's accessible (if needed)
        destination_itp = out_dir / ligand_itp.name
        shutil.copy(ligand_itp, destination_itp)
        logging.info(f"[Combine] Copied {ligand_itp} to {destination_itp}")

        # Verify combined TOP
        with combined_top.open("r") as ft:
            combined_top_contents = ft.read()
        logging.debug(f"[Combine] Combined TOP Contents:\n{combined_top_contents}")

        return combined_gro, combined_top
    except Exception as e:
        logging.error(f"[Combine] Error combining protein and ligand: {e}")
        raise

###############################################################################
# 4. GROMACS Pipeline Functions
###############################################################################


def run_gromacs_pipeline(combined_gro: Path, combined_top: Path, out_dir: Path, data_dir: Path,
                         simulation_time: int, temperature: float, pressure: float,
                         simulation_type: str, use_gpu: bool = False):
    """
    Execute the GROMACS simulation pipeline:
      1. Define box
      2. Solvate
      3. Add ions
      4. Energy minimization
      5. Production MD

    Args:
        combined_gro (Path): Combined GRO file.
        combined_top (Path): Combined topology file.
        out_dir (Path): Output directory for simulation files.
        data_dir (Path): Directory containing MDP files.
        simulation_time (int): Simulation time in ns.
        temperature (float): Simulation temperature in Kelvin.
        pressure (float): Simulation pressure in bar.
        simulation_type (str): Type of simulation ('em' or 'md').
        use_gpu (bool): Whether to utilize GPU acceleration.
    """
    try:
        # Ensure output directory exists
        out_dir.mkdir(parents=True, exist_ok=True)
        logging.info(f"[GROMACS] Starting GROMACS pipeline in {out_dir}")

        # 1. Define Box
        box_gro = out_dir / "box.gro"
        cmd = [
            "gmx", "editconf",
            "-f", str(combined_gro),
            "-o", str(box_gro),
            "-c", "-d", "1.2", "-bt", "cubic"
        ]
        subprocess.run(cmd, check=True)
        logging.info(f"[GROMACS] Defined simulation box: {box_gro}")

        # 2. Solvate
        solvated_gro = out_dir / "solvated.gro"
        topol_combined = combined_top

        # Ensure the solvent model matches the force field's solvent model
        solvent_itp = "tip3p.gro"  # Adjust as per force field

        # Check if the solvent_itp file exists, else raise an error
        solvent_itp_path = data_dir / "solvents" / solvent_itp
        if not solvent_itp_path.exists():
            logging.error(f"[GROMACS] Solvent topology file {solvent_itp_path} does not exist.")
            raise FileNotFoundError(f"Solvent topology file {solvent_itp_path} does not exist.")

        cmd = [
            "gmx", "solvate",
            "-cp", str(box_gro),
            "-cs", str(solvent_itp_path),
            "-o", str(solvated_gro),
            "-p", str(topol_combined)
        ]
        subprocess.run(cmd, check=True)
        logging.info(f"[GROMACS] Solvated system: {solvated_gro}")

        # 3. Add Ions
        ions_tpr = out_dir / "ions.tpr"
        ions_mdp = data_dir / "molecular_dynamics" / "ions.mdp"

        if not ions_mdp.exists():
            logging.error(f"[GROMACS] Ions MDP file {ions_mdp} does not exist.")
            raise FileNotFoundError(f"Ions MDP file {ions_mdp} does not exist.")

        cmd = [
            "gmx", "grompp",
            "-f", str(ions_mdp),
            "-c", str(solvated_gro),
            "-p", str(topol_combined),
            "-o", str(ions_tpr),
            "-maxwarn", "1"  # danger zone: todo: fix this.
        ]
        subprocess.run(cmd, check=True)
        logging.info(f"[GROMACS] Prepared ions topology: {ions_tpr}")

        # Run genion to add ions
        solvated_ions_gro = out_dir / "solvated_ions.gro"
        cmd = [
            "gmx", "genion",
            "-s", str(ions_tpr),
            "-o", str(solvated_ions_gro),
            "-p", str(topol_combined),
            "-pname", "NA", "-nname", "CL", "-neutral"
        ]
        try:
            # Automate input by echoing the group to replace
            # Using subprocess with input parameter
            process = subprocess.Popen(["gmx", "genion", "-s", str(ions_tpr),
                                        "-o", str(solvated_ions_gro),
                                        "-p", str(topol_combined),
                                        "-pname", "NA", "-nname", "CL", "-neutral"],
                                       stdin=subprocess.PIPE,
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE,
                                       text=True)
            stdout, stderr = process.communicate(input="SOL\n")
            if process.returncode != 0:
                logging.error(f"[GROMACS] genion failed:\n{stderr}")
                raise RuntimeError(f"genion failed: {stderr}")
            logging.info(f"[GROMACS] Added ions: {solvated_ions_gro}")
        except subprocess.CalledProcessError as e:
            logging.error(f"[GROMACS] genion failed: {e}")
            raise

        # 4. Energy Minimization
        em_tpr = out_dir / "em.tpr"
        minim_mdp = data_dir / "molecular_dynamics" / "minim.mdp"

        if not minim_mdp.exists():
            logging.error(f"[GROMACS] Minimization MDP file {minim_mdp} does not exist.")
            raise FileNotFoundError(f"Minimization MDP file {minim_mdp} does not exist.")

        cmd = [
            "gmx", "grompp",
            "-f", str(minim_mdp),
            "-c", str(solvated_ions_gro),
            "-p", str(topol_combined),
            "-o", str(em_tpr)
        ]
        subprocess.run(cmd, check=True)
        logging.info(f"[GROMACS] Prepared energy minimization TPR: {em_tpr}")

        # Run energy minimization
        em_defnm = out_dir / "em"
        cmd = [
            "gmx", "mdrun",
            "-v", "-deffnm", str(em_defnm)
        ]
        if use_gpu:
            cmd.extend(["-gpu_id", "0"])  # Adjust if multiple GPUs
        subprocess.run(cmd, check=True)
        logging.info(f"[GROMACS] Completed energy minimization: {em_defnm}.gro")

        # 5. Production MD Simulation
        if simulation_type.lower() == "md":
            md_tpr = out_dir / "md.tpr"
            md_mdp = data_dir / "molecular_dynamics" / "md.mdp"

            if not md_mdp.exists():
                logging.error(f"[GROMACS] Production MD MDP file {md_mdp} does not exist.")
                raise FileNotFoundError(f"Production MD MDP file {md_mdp} does not exist.")

            cmd = [
                "gmx", "grompp",
                "-f", str(md_mdp),
                "-c", f"{em_defnm}.gro",  # Start from minimized structure
                "-p", str(topol_combined),
                "-o", str(md_tpr)
            ]
            subprocess.run(cmd, check=True)
            logging.info(f"[GROMACS] Prepared production MD TPR: {md_tpr}")

            # Run production MD
            md_defnm = out_dir / "md"
            cmd = [
                "gmx", "mdrun",
                "-v", "-deffnm", str(md_defnm)
            ]
            if use_gpu:
                cmd.extend(["-gpu_id", "0"])  # Adjust if multiple GPUs
            subprocess.run(cmd, check=True)
            logging.info(f"[GROMACS] Completed production MD simulation: {md_defnm}.gro")

    except Exception as e:
        logging.error(f"[GROMACS] GROMACS pipeline failed: {e}")
        logging.error(traceback.format_exc())
        raise

###############################################################################
# 5. System Validation Function
###############################################################################

def validate_system(combined_gro: Path, combined_top: Path, molecule_atom_counts: dict):
    """
    Validate the combined GRO and topology files for atom count consistency.
    """
    try:
        validate_topology_gro_consistency(combined_gro, combined_top, molecule_atom_counts)
    except Exception as e:
        logging.error(f"[Validation] System validation failed: {e}")
        raise

###############################################################################
# 6. Process MD Simulation
###############################################################################


def process_md_simulation(row: dict, docking_dir: Path, md_exp_results_dir: Path,
                          db_connection_string: str, experiment_id: str,
                          data_dir: Path, simulation_type: str,
                          force_field: str, water_model: str):
    try:
        docking_id = row["docking_id"]
        uniprot_id = row["uniprot_id"]
        ligand_cid = row["ligand_cid"]
        deg_id = row["deg_id"]
        binding_energy = row["binding_energy"]
        gene_name = row["gene_name"]

        sim_time = row.get("generic_sim_time", 100)
        temp = row.get("generic_temp", 300)
        pressure = row.get("generic_pressure", 1)
        solvent = row.get("generic_solvent", "water")

        # Locate the clean protein PDB file
        protein_pdb = docking_dir / str(experiment_id) / f"{uniprot_id}_clean.pdb"
        if not protein_pdb.exists():
            logging.error(f"[ProcessMD] Clean protein PDB {protein_pdb} does not exist.")
            return

        # The final MD results folder
        md_exp_sim_results_dir = md_exp_results_dir / f"docking_{gene_name}_{docking_id}"
        md_exp_sim_results_dir.mkdir(parents=True, exist_ok=True)

        # Define file paths for ligand
        ligand_gro = docking_dir / str(experiment_id) / "ligand.gro"
        ligand_itp = docking_dir / str(experiment_id) / "ligand.itp"

        if not ligand_gro.exists() or not ligand_itp.exists():
            logging.error(f"[ProcessMD] Ligand files {ligand_gro} or {ligand_itp} do not exist.")
            return

        # Generate protein topology
        protein_processed_gro, protein_top = generate_protein_topology(
            protein_pdb=protein_pdb,
            force_field=force_field,
            water_model=water_model,
            output_dir=md_exp_sim_results_dir
        )

        # Combine protein and ligand
        combined_gro, combined_top = combine_protein_ligand(
            protein_gro=protein_processed_gro,
            protein_top=protein_top,
            ligand_gro=ligand_gro,
            ligand_itp=ligand_itp,
            out_dir=md_exp_sim_results_dir,
            protein_name="Protein_chain_A",
            ligand_name="LIG",
            system_name=f"Protein-Ligand System: {gene_name}"
        )

        # Define molecule atom counts (Only include molecules present in [ molecules ] section)
        molecule_atom_counts = {
            "Protein_chain_A": count_atoms_in_itp(protein_top),
            "LIG": count_atoms_in_itp(ligand_itp),
        }

        # Validate atom counts
        validate_system(
            combined_gro=combined_gro,
            combined_top=combined_top,
            molecule_atom_counts=molecule_atom_counts
        )

        # Run GROMACS pipeline
        run_gromacs_pipeline(
            combined_gro=combined_gro,
            combined_top=combined_top,
            out_dir=md_exp_sim_results_dir,
            data_dir=data_dir,
            simulation_time=sim_time,
            temperature=temp,
            pressure=pressure,
            simulation_type=simulation_type,
            use_gpu=True  # Set based on configuration
        )

        # Insert DB record
        record = {
            "docking_id": docking_id,
            "ligand_cid": ligand_cid,
            "deg_id": deg_id,
            "binding_energy": binding_energy,
            "simulation_time": sim_time,
            "temperature": temp,
            "pressure": pressure,
            "simulation_type": simulation_type,
            "status": "Completed"
        }
        insert_md_record(db_connection_string, record)
        logging.info(f"[ProcessMD] MD simulation completed for docking_id={docking_id}")

    except Exception as e:
        traceback.print_exc()
        logging.error(f"[ProcessMD] Simulation failed for docking_id={row.get('docking_id', 'unknown')}: {e}")

def parse_experiment_id(exp_id: str) -> str:
    """
    If exp_id is a file, read its content as the real experiment ID.
    Otherwise, just return it.
    """
    path_obj = Path(exp_id)
    if path_obj.is_file():
        experiment_id = path_obj.read_text().strip()
        logging.info(f"[ExperimentID] Read experiment ID {experiment_id} from file {exp_id}")
        return experiment_id
    logging.info(f"[ExperimentID] Using experiment ID: {exp_id}")
    return exp_id

###############################################################################
# 7. Main Function
###############################################################################


def main():
    # Argument parsing and setup
    parser = argparse.ArgumentParser(description="Perform MD simulations using parameters from md_parameters.csv")
    parser.add_argument("--db_connection_string", required=True, help="Database connection string.")
    parser.add_argument("--experiment_id", required=True, help="Experiment ID or path to a file containing it.")
    parser.add_argument("--md_param_file", required=True, help="Path to md_parameters.csv.")
    parser.add_argument("--force_field", default="amber99sb-ildn", help="Force field to use for pdb2gmx.")
    parser.add_argument("--water_model", default="tip3p", help="Water model to use for pdb2gmx.")
    parser.add_argument("--output_dir", required=True, help="Output directory for MD simulations.")
    parser.add_argument("--data_dir", required=True, help="Data directory containing MDP files.")
    parser.add_argument("--num_workers", type=int, default=1, help="Number of parallel workers.")
    parser.add_argument("--simulation_type", choices=["em", "md"], default="md",
                        help="Choose simulation type: 'em' for energy minimization or 'md' for production MD.")
    parser.add_argument("--project_dir", required=True, help="Project directory path.")
    args = parser.parse_args()

    simulation_type = args.simulation_type

    # Paths
    project_dir = Path(args.project_dir)
    src = project_dir / "results"
    final_dest = project_dir / "results_out"

    output_dir = Path(args.output_dir)
    data_dir = Path(args.data_dir)
    experiment_id = parse_experiment_id(args.experiment_id)

    md_exp_results_dir = output_dir / "molecular_dynamics" / experiment_id
    md_exp_results_dir.mkdir(parents=True, exist_ok=True)

    # Load parameters
    try:
        md_params = pd.read_csv(args.md_param_file)
        logging.info(f"[Main] Loaded MD parameters from {args.md_param_file}")
    except Exception as e:
        logging.error(f"[Main] Failed to load MD parameters: {e}")
        return

    docking_dir = output_dir / "molecular_docking"

    logging.info(f"[Main] Docking directory: {docking_dir}")
    logging.info(f"[Main] Experiment ID: {experiment_id}")

    # Run simulations
    if args.num_workers > 1:
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.num_workers) as executor:
            futures = []
            for _, row in md_params.iterrows():
                futures.append(executor.submit(
                    process_md_simulation,
                    row,
                    docking_dir,
                    md_exp_results_dir,
                    args.db_connection_string,
                    experiment_id,
                    data_dir,
                    simulation_type,
                    args.force_field,
                    args.water_model
                ))

            for future in concurrent.futures.as_completed(futures):
                try:
                    future.result()
                except Exception as e:
                    logging.error(f"[Main] A simulation task encountered an exception: {e}")
    else:
        for _, row in md_params.iterrows():
            process_md_simulation(
                row=row,
                docking_dir=docking_dir,
                md_exp_results_dir=md_exp_results_dir,
                db_connection_string=args.db_connection_string,
                experiment_id=experiment_id,
                data_dir=data_dir,
                simulation_type=simulation_type,
                force_field=args.force_field,
                water_model=args.water_model
            )

    # Copy results, skipping symlinks
    try:
        shutil.copytree(src, final_dest, dirs_exist_ok=True, ignore=shutil.ignore_patterns('*', '*/.*'))
        logging.info(f"[Main] Copied simulation results from {src} to {final_dest}")
    except Exception as e:
        logging.error(f"[Main] Failed to copy simulation results: {e}")

    logging.info("[Pipeline] All simulations completed.")


if __name__ == "__main__":
    main()
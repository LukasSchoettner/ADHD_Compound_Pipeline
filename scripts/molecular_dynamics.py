#!/usr/bin/env python
import shutil
import time
import os
import argparse
from os import mkdir

import pandas as pd
from pathlib import Path, PurePath
import subprocess
import traceback
from sqlalchemy import create_engine, text

def wait_for_file(file_path, max_retries=5, delay=1):
    """
    Wait for a file to be created and ensure it's not empty.
    """
    for _ in range(max_retries):
        if file_path.exists() and file_path.stat().st_size > 0:
            return
        time.sleep(delay)
    raise FileNotFoundError(f"[Error] File {file_path} was not created or is empty after {max_retries * delay} seconds.")

def extract_moleculetype(itp_file):
    with open(itp_file, 'r') as f:
        for line in f:
            if line.strip().startswith("[ moleculetype ]"):
                next(f)  # Skip header line
                first_line = next(f).strip()
                return first_line.split()[0]
    raise ValueError(f"Couldn't find moleculetype in {itp_file}")

def parse_md_parameters(md_param_file):
    """
    Parse the md_parameters.csv file into a DataFrame.
    """
    try:
        df = pd.read_csv(md_param_file)
        return df
    except Exception as e:
        print(f"[Error] Failed to parse MD parameters: {e}")
        raise

def run_command(cmd, cwd=None):
    """
    Helper function to run a shell command with error checking.
    """
    print("[CMD]", " ".join(str(x) for x in cmd))
    result = subprocess.run(cmd, check=False, cwd=cwd,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    print(result.stdout)
    if result.returncode != 0:
        print(result.stderr)
        raise subprocess.CalledProcessError(result.returncode, cmd)
    return result

def split_protein_ligand(cleaned_complex_pdb, protein_only_pdb, ligand_only_pdb,
                         ligand_resname="LIG"):
    """
    Split the 'cleaned_complex_pdb' into two files:
      - protein_only_pdb (ATOM lines not matching ligand_resname)
      - ligand_only_pdb  (ATOM/HETATM lines with residue name == ligand_resname)

    Adjust the residue numbering if needed. This function assumes the ligand is labeled with
    'resname == ligand_resname' in the PDB file.
    """
    protein_lines = []
    ligand_lines = []

    with open(cleaned_complex_pdb, "r") as fin:
        for line in fin:
            if not line.startswith(("ATOM", "HETATM")):
                # Copy non-atom lines (e.g. REMARK, TER, etc.) to the protein file or ignore
                # Typically you might keep them in the protein PDB
                protein_lines.append(line)
                continue

            # Parse columns: see PDB format, resName is columns 17-20 (1-based)
            # or line[17:20] in zero-based indexing
            residue_name = line[17:20].strip()
            if residue_name.upper() == ligand_resname.upper():
                ligand_lines.append(line)
            else:
                protein_lines.append(line)

    # Write out protein
    with open(protein_only_pdb, "w") as fprot:
        fprot.writelines(protein_lines)

    # Write out ligand
    with open(ligand_only_pdb, "w") as flig:
        flig.writelines(ligand_lines)

    print(f"[DEBUG] Wrote protein to {protein_only_pdb}, ligand to {ligand_only_pdb}")

def parameterize_protein(protein_pdb, out_dir, ff="amber99sb", water="tip3p"):
    """
    Run pdb2gmx on the protein-only PDB to get protein.gro and topol.top.
    """
    protein_gro = out_dir / "protein.gro"
    protein_top = out_dir / "topol.top"

    cmd = [
        "gmx", "pdb2gmx",
        "-f", protein_pdb,
        "-o", protein_gro,
        "-p", protein_top,
        "-ff", ff,
        "-water", water
    ]
    run_command(cmd)
    return protein_gro, protein_top

def parameterize_ligand(ligand_pdb, md_exp_sim_results_dir):
    """
    Use acpype to generate ligand topology (ligand.acpype/*) with GAFF/Amber.
    """
    ligand_pdb = ligand_pdb.resolve()  # Convert to absolute path
    if not ligand_pdb.exists() or ligand_pdb.stat().st_size == 0:
        raise FileNotFoundError(f"[Error] Ligand file {ligand_pdb} does not exist or is empty!")

    cmd = [
        "acpype",
        "-i", str(ligand_pdb),
        "-c", "bcc",  # partial charge method
        "-n", "0"     # net charge (adjust if your ligand is charged)
    ]
    run_command(cmd, cwd=md_exp_sim_results_dir)

    base = ligand_pdb.stem  # e.g. "ligand_only"
    acpype_dir = md_exp_sim_results_dir / f"{base}.acpype"
    ligand_itp = acpype_dir / f"{base}_GMX.itp"
    ligand_gro = acpype_dir / f"{base}_GMX.gro"

    if not ligand_itp.exists():
        raise FileNotFoundError(f"[ACPYPE Error] Could not find {ligand_itp}")
    if not ligand_gro.exists():
        raise FileNotFoundError(f"[ACPYPE Error] Could not find {ligand_gro}")

    return ligand_itp, ligand_gro

def create_combined_system(protein_gro, protein_top,
                           ligand_gro, ligand_itp,
                           out_dir,
                           protein_name="Protein_chain_A",
                           ligand_name="LIG",
                           system_name="Protein-Ligand System"):
    """
    Combine protein and ligand coordinates into one .gro,
    and produce a single 'topol_combined.top' with [system] and [molecules].
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    combined_gro = out_dir / "combined.gro"
    combined_top = out_dir / "topol_combined.top"

    # Step 1) Read original topol and remove existing [ system ] / [ molecules ].
    with open(protein_top, "r") as fin:
        original_lines = fin.readlines()

    cleaned_lines = []
    skip_block = False
    for line in original_lines:
        # Detect old blocks
        if line.strip().startswith("[ system ]"):
            skip_block = True
        elif line.strip().startswith("[ molecules ]"):
            skip_block = True
        elif line.strip().startswith("[") and line.strip().endswith("]"):
            skip_block = False
        if not skip_block:
            cleaned_lines.append(line)

    # Step 2) Write out adjusted lines, insert ligand .itp
    with open(combined_top, "w") as fout:
        for line in cleaned_lines:
            fout.write(line)
            if "forcefield.itp" in line:
                # Insert ligand .itp right after the forcefield includes
                ligand_relative_path = ligand_itp.relative_to(out_dir)
                fout.write(f'#include "{ligand_relative_path}"\n')

        fout.write("\n[ system ]\n")
        fout.write(f"{system_name}\n\n")

        fout.write("[ molecules ]\n")
        fout.write("; Compound        #mols\n")
        fout.write(f"{protein_name:<20} 1\n")
        fout.write(f"{ligand_name:<20} 1\n")

    # Step 3) Combine GROs
    with open(protein_gro, "r") as f_prot, open(ligand_gro, "r") as f_lig, open(combined_gro, "w") as f_out:
        prot_lines = f_prot.readlines()
        lig_lines = f_lig.readlines()

        # Protein
        num_atoms_prot = int(prot_lines[1].strip())
        atom_lines_prot = prot_lines[2:-1]
        box_line_prot   = prot_lines[-1]

        # Ligand
        num_atoms_lig = int(lig_lines[1].strip())
        atom_lines_lig = lig_lines[2:-1]
        box_line_lig   = lig_lines[-1]

        combined_num_atoms = num_atoms_prot + num_atoms_lig

        # Write combined
        # Title line: use the protein's first line or a custom line
        f_out.write("Protein + Ligand combined\n")
        f_out.write(f"{combined_num_atoms}\n")
        f_out.writelines(atom_lines_prot)
        f_out.writelines(atom_lines_lig)

        # Usually we want the box line from the *protein*, ignoring ligand box
        # But if the ligand box is different, choose whichever is correct
        f_out.write(box_line_prot)

    # Debug info
    print(f"[DEBUG] Combined system in {combined_gro}: protein = {num_atoms_prot}, ligand = {num_atoms_lig}, total = {combined_num_atoms}")

    return combined_gro, combined_top

def validate_topology_gro_consistency(gro_file, top_file):
    """
    Check the .gro and .top have the same total atom count (including water, ions, etc.).
    """
    # Read the .gro
    with open(gro_file, 'r') as f:
        lines = f.readlines()
        num_atoms_gro = int(lines[1].strip())

    # Sum up from the [ molecules ] block in the .top
    with open(top_file, 'r') as f:
        top_lines = f.readlines()

    molecule_section = False
    total_atoms_top = 0
    molecule_counts = {}
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
                molname, count_str = line.split()
                count = int(count_str)
                molecule_counts[molname] = count

                # Hard-coded example: you must customize the actual
                # per-molecule atom counts for your system
                if molname == "Protein_chain_A":
                    total_atoms_top += 7195 * count
                elif molname == "LIG":
                    total_atoms_top += 38 * count
                elif molname == "SOL":
                    total_atoms_top += 3 * count
                elif molname in ["NA", "CL"]:
                    total_atoms_top += 1 * count  # if ions appear as separate molecule lines
                else:
                    raise ValueError(f"Unknown molecule type: {molname}")
            except ValueError:
                raise ValueError(f"Invalid line in topology file: '{line}'")

    print(f".gro file atom count: {num_atoms_gro}")
    print(f"Topology file total atom count: {total_atoms_top}")
    print(f"Breakdown: {molecule_counts}")

    if num_atoms_gro != total_atoms_top:
        raise ValueError(f"Mismatch in .gro ({num_atoms_gro}) vs .top ({total_atoms_top})")

def run_gromacs_pipeline(combined_gro, combined_top, out_dir, data_dir,
                         simulation_time, temperature, pressure, simulation_type, use_gpu=False):
    """
    The usual steps: define box, solvate, add ions, energy minimize, etc.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1. editconf -> box
    box_gro = out_dir / "box.gro"
    cmd = [
        "gmx", "editconf",
        "-f", combined_gro,
        "-o", box_gro,
        "-c", "-d", "1.2", "-bt", "cubic"
    ]
    run_command(cmd)

    # 2. solvate
    solvated_gro = out_dir / "solvated.gro"
    topol_combined = out_dir / "topol_combined.top"

    cmd = [
        "gmx", "solvate",
        "-cp", box_gro,
        "-cs", "spc216.gro",
        "-o", solvated_gro,
        "-p", topol_combined
    ]
    run_command(cmd)

    # 3. Ion param: grompp -> genion
    ions_tpr = out_dir / "ions.tpr"
    ions_mdp = Path(data_dir) / "molecular_dynamics" / "ions.mdp"

    cmd = [
        "gmx", "grompp",
        "-f", ions_mdp,
        "-c", solvated_gro,
        "-p", topol_combined,
        "-o", ions_tpr,
        "-maxwarn", "1"
    ]
    run_command(cmd)

    solvated_ions_gro = out_dir / "solvated_ions.gro"
    cmd = [
        "gmx", "genion",
        "-s", ions_tpr,
        "-o", solvated_ions_gro,
        "-p", topol_combined,
        "-pname", "NA", "-nname", "CL", "-neutral"
    ]
    result = subprocess.run(cmd, input="SOL\n", text=True, check=True)

    # validate_topology_gro_consistency(solvated_ions_gro, topol_combined)

    # 4. Minimization
    em_tpr = out_dir / "em.tpr"
    minim_mdp = Path(data_dir) / "molecular_dynamics" / "minim.mdp"
    cmd = [
        "gmx", "grompp",
        "-f", minim_mdp,
        "-c", solvated_ions_gro,
        "-p", topol_combined,
        "-o", em_tpr
    ]
    run_command(cmd)

    # Use GPU for energy minimization
    em_defnm = str(out_dir / "em")

    if use_gpu:
        cmd = [
            "gmx", "mdrun",
            "-v", "-deffnm", em_defnm,
            "-gpu_id", "0"  # Specify GPU device index; adjust if needed
        ]
    else:
        cmd = [
            "gmx", "mdrun",
            "-v", "-deffnm", em_defnm,
        ]
    run_command(cmd)

    print(f"[MD] Completed EM. Results in {out_dir}")

    # If simulation type is MD, proceed with production MD simulation
    if simulation_type == "md":
        # 4. Production MD simulation setup
        md_tpr = out_dir / "md.tpr"
        md_mdp = Path(data_dir) / "molecular_dynamics" / "md.mdp"  # Ensure md.mdp exists here

        cmd = [
            "gmx", "grompp",
            "-f", str(md_mdp),
            "-c", str(out_dir / "em.gro"),  # start from the minimized structure
            "-p", str(topol_combined),
            "-o", str(md_tpr)
        ]
        run_command(cmd)

        # Use GPU for production MD
        md_defnm = str(out_dir / "md")
        if use_gpu:
            cmd = [
                "gmx", "mdrun",
                "-v", "-deffnm", md_defnm,
                "-gpu_id", "0"  # Adjust if needed
            ]
        else:
            cmd = [
                "gmx", "mdrun",
                "-v", "-deffnm", md_defnm
            ]
        run_command(cmd)

        print(f"[MD] Completed production MD. Results in {out_dir}")

def insert_md_record(db_connection_string, record):
    """
    Insert a row into the molecular_dynamic table with the MD simulation info.
    """
    engine = create_engine(db_connection_string)
    insert_sql = text("""
        INSERT INTO molecular_dynamic (
            docking_id, ligand_cid, deg_id, binding_energy, gene_name,
            simulation_time, temperature, pressure, solvent_type, status
        )
        VALUES (
            :docking_id, :ligand_cid, :deg_id, :binding_energy, :gene_name,
            :simulation_time, :temperature, :pressure, :solvent_type, :status
        )
    """)

    try:
        with engine.begin() as conn:
            conn.execute(insert_sql, record)
        print(f"[DB] Record inserted for docking_id={record['docking_id']}")
    except Exception as e:
        traceback.print_exc()
        print(f"[DB] Failed to insert record for docking_id={record['docking_id']}")
    finally:
        engine.dispose()

def process_md_simulation(row, docking_dir, md_exp_results_dir,
                          db_connection_string, experiment_id, data_dir, md_exp_data_dir, simulation_type):
    """
    Process a single MD simulation:
      1) Extract first model (if multiple).
      2) Split into protein_only.pdb + ligand_only.pdb.
      3) Parameterize protein (pdb2gmx).
      4) Parameterize ligand (acpype).
      5) Create combined system.
      6) Run GROMACS pipeline.
      7) Insert DB record.
    """
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

    # Input PDB from docking
    docked_pdb = Path(docking_dir) / str(experiment_id) / f"{uniprot_id}_docked.pdb"
    cleaned_complex_pdb = Path(docking_dir) / str(experiment_id) / f"{uniprot_id}_docked_complex.pdb"

    # The final MD results folder
    md_exp_sim_results_dir = md_exp_results_dir / f"docking_{gene_name}_{docking_id}"
    md_exp_sim_results_dir.mkdir(parents=True, exist_ok=True)

    # 2. Split
    md_exp_sim_results_dir.mkdir(parents=True, exist_ok=True)
    protein_only_pdb = md_exp_sim_results_dir / "protein_only.pdb"
    ligand_only_pdb  = md_exp_sim_results_dir / "ligand_only.pdb"

    # Adjust "LIG" below if your actual ligand residue name differs
    split_protein_ligand(
        cleaned_complex_pdb,
        protein_only_pdb,
        ligand_only_pdb,
        ligand_resname="LIG"
    )

    # 3. Parameterize protein
    protein_gro, protein_top = parameterize_protein(protein_only_pdb, md_exp_sim_results_dir)

    # 4. Parameterize ligand
    ligand_itp, ligand_gro = parameterize_ligand(ligand_only_pdb, md_exp_sim_results_dir)

    ligand_moleculetype = extract_moleculetype(ligand_itp)

    # 5. Combine
    combined_gro, combined_top = create_combined_system(
        protein_gro, protein_top,
        ligand_gro, ligand_itp,
        out_dir=md_exp_sim_results_dir,
        protein_name="Protein_chain_A",
        ligand_name=ligand_moleculetype,
        system_name=f"Protein-Ligand System: {gene_name}"
    )

    # 6. Run GROMACS
    run_gromacs_pipeline(
        combined_gro, combined_top,
        md_exp_sim_results_dir,
        data_dir,
        sim_time, temp, pressure,
        simulation_type
    )

    # 7. Insert DB record
    record = {
        "docking_id": docking_id,
        "ligand_cid": ligand_cid,
        "deg_id": deg_id,
        "binding_energy": binding_energy,
        "gene_name": gene_name,
        "simulation_time": sim_time,
        "temperature": temp,
        "pressure": pressure,
        "solvent_type": solvent,
        "status": "Completed"
    }
    insert_md_record(db_connection_string, record)
    print(f"[MD] Done for docking_id={docking_id}")

def get_experiment_id(experiment_id):
    if os.path.isfile(experiment_id):
        with open(experiment_id, "r") as f:
            return int(f.read().strip())
    return int(experiment_id)

def main():
    parser = argparse.ArgumentParser(description="Perform MD simulations using parameters from md_parameters.csv")
    parser.add_argument("--db_connection_string", required=True)
    parser.add_argument("--experiment_id", required=True)
    parser.add_argument("--md_param_file", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--data_dir", required=True)
    parser.add_argument("--num_workers", type=int, default=1)
    parser.add_argument("--simulation_type", choices=["em", "md"], default="md",
                        help="Choose simulation type: 'em' for energy minimization or 'md' for production MD")
    args = parser.parse_args()

    simulation_type = args.simulation_type

    # Paths
    output_dir = Path(args.output_dir)
    data_dir = Path(args.data_dir)
    experiment_id = str(get_experiment_id(args.experiment_id))

    md_exp_results_dir = output_dir / "molecular_dynamics" / experiment_id
    md_exp_results_dir.mkdir(parents=True, exist_ok=True)

    md_exp_data_dir = data_dir / "molecular_dynamics" / experiment_id
    md_exp_data_dir.mkdir(parents=True, exist_ok=True)

    # Load parameters
    md_params = parse_md_parameters(args.md_param_file)
    docking_dir = output_dir / "molecular_docking"

    print(md_params)
    print("docking_dir =", docking_dir)
    print("experiment_id =", experiment_id)

    # Run each row
    for idx, row in md_params.iterrows():
        process_md_simulation(
            row,
            docking_dir,
            md_exp_results_dir,
            args.db_connection_string,
            experiment_id,
            data_dir,
            md_exp_data_dir,
            simulation_type
        )

    print("[Pipeline] All simulations completed.")

if __name__ == "__main__":
    main()

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


def parameterize_protein(protein_pdb, out_dir, ff="amber99sb", water="tip3p"):
    """
    Run pdb2gmx on the protein alone to get protein.gro and topol.top.
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

    # Ensure the file exists and is not empty
    ligand_pdb = ligand_pdb.resolve()  # Convert to absolute path
    if not ligand_pdb.exists() or ligand_pdb.stat().st_size == 0:
        raise FileNotFoundError(f"[Error] Ligand file {ligand_pdb} does not exist or is empty!")

    # Run ACPYPE with the absolute path
    cmd = [
        "acpype",
        "-i", str(ligand_pdb),
        "-c", "bcc",  # partial charge method
        "-n", "0"     # net charge; change if your ligand is not neutral
    ]
    run_command(cmd, cwd=md_exp_sim_results_dir)

    # Locate output files
    base = ligand_pdb.stem  # e.g., "P17752_docked_clean"
    acpype_dir = md_exp_sim_results_dir / f"{base}.acpype"
    ligand_itp = acpype_dir / f"{base}_GMX.itp"
    ligand_gro = acpype_dir / f"{base}_GMX.gro"

    if not ligand_itp.exists():
        raise FileNotFoundError(f"[ACPYPE Error] Could not find {ligand_itp}")
    if not ligand_gro.exists():
        raise FileNotFoundError(f"[ACPYPE Error] Could not find {ligand_gro}")

    return ligand_itp, ligand_gro


def extract_first_model(input_pdb, output_pdb):
    """
    Extract the first complete MODEL block from a PDB file, ignoring duplicate MODEL lines.
    """
    inside_model = False
    model_found = False
    with open(input_pdb, "r") as infile, open(output_pdb, "w") as outfile:
        for line in infile:
            if line.startswith("MODEL") and not model_found:
                inside_model = True
                model_found = True
                outfile.write(line)
                continue
            if inside_model:
                outfile.write(line)
            if line.startswith("ENDMDL") and inside_model:
                break
    if not model_found:
        raise ValueError(f"No valid MODEL block found in {input_pdb}")

def create_combined_system(protein_gro, protein_top, ligand_gro, ligand_itp, out_dir, ligand_name, protein_name="Protein_chain_A", system_name="Protein-Ligand System"):
    """
    Combine protein and ligand coordinates into one .gro,
    and update topol.top to include the ligand .itp and a single [ system ] and [ molecules ] section.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    combined_gro = out_dir / "combined.gro"
    combined_top = out_dir / "topol_combined.top"

    # 1) Read the original topol
    with open(protein_top, "r") as fin:
        original_lines = fin.readlines()

    # 2) Remove any existing [ system ] / [ molecules ] blocks
    cleaned_lines = []
    skip_block = False
    for line in original_lines:
        if line.strip().startswith("[ system ]"):
            skip_block = True
        elif line.strip().startswith("[ molecules ]"):
            skip_block = True
        elif line.strip().startswith("[") and line.strip().endswith("]"):
            # Found a new block, so we stop skipping
            skip_block = False

        if not skip_block:
            cleaned_lines.append(line)

    # 3) Write out the adjusted lines and insert the ligand itp
    with open(combined_top, "w") as fout:
        for line in cleaned_lines:
            fout.write(line)
            if "forcefield.itp" in line:
                # Insert the ligand .itp line right after the forcefield includes
                ligand_relative_path = ligand_itp.relative_to(out_dir)
                fout.write(f'#include "{ligand_relative_path}"\n')

        # Now add your single [ system ] block
        fout.write("\n[ system ]\n")
        fout.write(f"{system_name}\n\n")

        # And a single [ molecules ] block
        fout.write("[ molecules ]\n")
        fout.write("; Compound     #mols\n")
        fout.write(f"{protein_name:<20} 1\n")
        fout.write(f"{ligand_name:<20} 1\n")

    # 4) Combine the GRO files (same as you already do)
    with open(protein_gro, "r") as f_prot, open(ligand_gro, "r") as f_lig, open(combined_gro, "w") as f_out:
        prot_lines = f_prot.readlines()
        lig_lines = f_lig.readlines()

        header = prot_lines[0]
        num_atoms_prot = int(prot_lines[1].strip())
        atom_lines_prot = prot_lines[2:-1]
        box_line = prot_lines[-1]

        num_atoms_lig = int(lig_lines[1].strip())
        atom_lines_lig = lig_lines[2:-1]

        combined_num_atoms = num_atoms_prot + num_atoms_lig

        f_out.write(header)
        f_out.write(f"{combined_num_atoms}\n")
        f_out.writelines(atom_lines_prot)
        f_out.writelines(atom_lines_lig)
        f_out.write(box_line)

    return combined_gro, combined_top


def validate_topology_gro_consistency(gro_file, top_file):
    """
    Validate that the number of atoms in the .gro file matches the total
    atoms specified in the topology file.
    """
    # Read the number of atoms in the .gro file
    with open(gro_file, 'r') as f:
        lines = f.readlines()
        num_atoms_gro = int(lines[1].strip())

    # Calculate the total number of atoms from the topology file
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
                # End parsing if we hit an empty line or a comment
                continue

            if line.startswith("[") and line.endswith("]"):
                # End parsing when we hit the next section
                break

            try:
                molecule, count = line.split()
                count = int(count)
                molecule_counts[molecule] = count

                # Atom counts per molecule type (adjust for your system)
                if molecule == "Protein_chain_A":
                    total_atoms_top += 7195 * count  # Replace with actual atom count for protein
                elif molecule == "P17752_docked_clean":
                    total_atoms_top += 38 * count  # Replace with actual atom count for ligand
                elif molecule == "SOL":
                    total_atoms_top += 3 * count  # 3 atoms per water molecule
                else:
                    raise ValueError(f"Unknown molecule type in topology: {molecule}")

            except ValueError:
                raise ValueError(f"Invalid line in topology file: {line.strip()}")

    print(f".gro file atom count: {num_atoms_gro}")
    print(f"Topology file total atom count: {total_atoms_top}")
    print(f"Breakdown by molecule type: {molecule_counts}")

    if num_atoms_gro != total_atoms_top:
        raise ValueError(f"Mismatch in atom counts: .gro file ({num_atoms_gro}) vs topology file ({total_atoms_top})")


def run_gromacs_pipeline(combined_gro, combined_top, out_dir, data_dir,
                         simulation_time, temperature, pressure):
    """
    Once we have combined.gro and topol_combined.top, do the usual:
    1) editconf -> box.gro
    2) solvate -> solvated.gro
    3) grompp -> ions.tpr
    4) genion -> solvated_ions.gro
    5) grompp -> em.tpr
    6) mdrun -> em
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1. Define the box
    box_gro = out_dir / "box.gro"
    cmd = [
        "gmx", "editconf",
        "-f", combined_gro,
        "-o", box_gro,
        "-c", "-d", "1.0", "-bt", "cubic"
    ]
    run_command(cmd)

    # 2. Solvate
    solvated_gro = out_dir / "solvated.gro"
    topol_combined = out_dir / "topol_combined.top"  # Destination for GROMACS topology

    # Ensure we only copy if the source and destination are different
    if combined_top != topol_combined:
        shutil.copy2(combined_top, topol_combined)

    cmd = [
        "gmx", "solvate",
        "-cp", box_gro,
        "-cs", "spc216.gro",
        "-o", solvated_gro,
        "-p", topol_combined
    ]
    run_command(cmd)

    validate_topology_gro_consistency(solvated_gro, topol_combined)

    # 3. Ion addition
    ions_tpr = out_dir / "ions.tpr"
    ions_mdp = Path(data_dir) / "molecular_dynamics" / "ions.mdp"
    cmd = [
        "gmx", "grompp",
        "-f", ions_mdp,  # Ensure ions.mdp is available in your working directory
        "-c", solvated_gro,
        "-p", topol_combined,
        "-o", ions_tpr
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
    run_command(cmd)

    validate_topology_gro_consistency(solvated_ions_gro, topol_combined)

    # 4. Minimization
    em_tpr = out_dir / "em.tpr"
    cmd = [
        "gmx", "grompp",
        "-f", "minim.mdp",  # Ensure minim.mdp is in your working directory
        "-c", solvated_ions_gro,
        "-p", topol_combined,
        "-o", em_tpr
    ]
    run_command(cmd)

    cmd = [
        "gmx", "mdrun",
        "-v", "-deffnm", str(out_dir / "em")
    ]
    run_command(cmd)

    print(f"[MD] Molecular dynamics prep (EM) completed. Results in {out_dir}")



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


def process_md_simulation(row, docking_dir, md_exp_results_dir, db_connection_string, experiment_id, data_dir, md_exp_data_dir):
    """
    Process a single MD simulation entry:
      1. Parameterize protein (pdb2gmx).
      2. Clean docked file (extract first model).
      3. Parameterize ligand (acpype).
      4. Merge topologies -> combined.gro/top.
      5. Solvate, add ions, minimize.
      6. Insert record in DB.
    """
    docking_id = row["docking_id"]
    uniprot_id = row["uniprot_id"]    # Must exist in your CSV
    ligand_cid = row["ligand_cid"]
    deg_id = row["deg_id"]
    binding_energy = row["binding_energy"]
    gene_name = row["gene_name"]

    sim_time = row.get("generic_sim_time", 100)
    temp = row.get("generic_temp", 300)
    pressure = row.get("generic_pressure", 1)
    solvent = row.get("generic_solvent", "water")

    ligand_name = f"{uniprot_id}_docked_clean"

    # Paths
    protein_pdb = Path(docking_dir) / str(experiment_id) / f"{uniprot_id}_clean.pdb"
    docked_complex_pdb = Path(docking_dir) / str(experiment_id) / f"{uniprot_id}_docked.pdb"
    cleaned_complex_pdb = Path(docking_dir) / str(experiment_id) / f"{uniprot_id}_docked_clean.pdb"


    # The final MD results folder
    md_exp_sim_results_dir = md_exp_results_dir / f"docking_{gene_name}_{docking_id}"
    md_exp_sim_results_dir.mkdir(parents=True, exist_ok=True)

    # 1. Parameterize protein
    protein_gro, protein_top = parameterize_protein(protein_pdb, md_exp_sim_results_dir)

    # Extract the first model from the docked protein-ligand complex
    extract_first_model(docked_complex_pdb, cleaned_complex_pdb)

    # Wait for the cleaned file to be ready
    wait_for_file(cleaned_complex_pdb)

    # Parameterize ligand using the cleaned protein-ligand complex
    ligand_itp, ligand_gro = parameterize_ligand(cleaned_complex_pdb, md_exp_sim_results_dir)

    # 4. Combine them
    combined_gro, combined_top = create_combined_system(
        protein_gro, protein_top,
        ligand_gro, ligand_itp,
        md_exp_sim_results_dir,
        ligand_name
    )

    # 5. Run GROMACS pipeline: box -> solvate -> ions -> minim
    run_gromacs_pipeline(combined_gro, combined_top, md_exp_sim_results_dir, data_dir, sim_time, temp, pressure)

    # 6. Insert DB record
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
    """
    Parse experiment_id from file if it is a file path.
    """
    if os.path.isfile(experiment_id):
        with open(experiment_id, "r") as f:
            return int(f.read().strip())
    return int(experiment_id)


def main():
    parser = argparse.ArgumentParser(description="Perform MD simulations using parameters from md_parameters.csv")
    parser.add_argument("--db_connection_string", required=True, help="Database connection string (SQLAlchemy format).")
    parser.add_argument("--experiment_id", required=True, help="ID or file containing experiment ID.")
    parser.add_argument("--md_param_file", required=True, help="Path to the MD parameters CSV file.")
    parser.add_argument("--output_dir", required=True, help="Directory to store MD results.")
    parser.add_argument("--data_dir", required=True, help="Directory to store data.")
    parser.add_argument("--num_workers", type=int, default=1, help="Number of parallel workers (optional).")
    args = parser.parse_args()

    # Convert output_dir and data_dir to Path objects
    output_dir = Path(args.output_dir)
    data_dir = Path(args.data_dir)

    # Convert experiment_id to a string if it's a number
    experiment_id = str(get_experiment_id(args.experiment_id))

    # Define MD results and data directories
    md_exp_results_dir = output_dir / "molecular_dynamics" / experiment_id
    md_exp_results_dir.mkdir(parents=True, exist_ok=True)

    md_exp_data_dir = data_dir / "molecular_dynamics" / experiment_id
    md_exp_data_dir.mkdir(parents=True, exist_ok=True)

    # Load MD parameters
    md_params = parse_md_parameters(args.md_param_file)
    docking_dir = output_dir / "molecular_docking"

    print(md_params)
    print("docking_dir =", docking_dir)
    print("experiment_id =", experiment_id)

    # Loop over each row
    for idx, row in md_params.iterrows():
        process_md_simulation(row, docking_dir, md_exp_results_dir, args.db_connection_string, experiment_id, args.data_dir, md_exp_data_dir)

    print("[Pipeline] All simulations completed.")



if __name__ == "__main__":
    main()

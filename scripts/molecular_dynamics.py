import os
import argparse
import pandas as pd
from pathlib import Path
import subprocess
import traceback
from concurrent.futures import ProcessPoolExecutor
from sqlalchemy import create_engine, text


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


def combine_protein_ligand(protein_pdb, ligand_pdb, output_pdb):
    """
    Combine protein and ligand into a single PDB file for MD simulation.
    """
    with open(output_pdb, "w") as outfile:
        with open(protein_pdb, "r") as f1:
            outfile.write(f1.read())
        with open(ligand_pdb, "r") as f2:
            for line in f2:
                if line.startswith(("ATOM", "HETATM")):
                    outfile.write(line)
    print(f"[MD] Combined protein and ligand into {output_pdb}")


def perform_md_simulation_with_gromacs(protein_ligand_pdb, output_dir, simulation_time, temperature, pressure):
    """
    Prepare and run MD simulation using GROMACS.
    """
    os.makedirs(output_dir, exist_ok=True)

    # Step 1: Generate topology
    cmd = [
        "gmx", "pdb2gmx",
        "-f", protein_ligand_pdb,
        "-o", os.path.join(output_dir, "complex.gro"),
        "-p", os.path.join(output_dir, "topol.top"),
        "-ff", "amber99sb",
        "-water", "tip3p"
    ]
    subprocess.run(cmd, check=True)

    # Step 2: Define simulation box
    cmd = [
        "gmx", "editconf",
        "-f", os.path.join(output_dir, "complex.gro"),
        "-o", os.path.join(output_dir, "box.gro"),
        "-c", "-d", "1.0", "-bt", "cubic"
    ]
    subprocess.run(cmd, check=True)

    # Step 3: Solvate
    cmd = [
        "gmx", "solvate",
        "-cp", os.path.join(output_dir, "box.gro"),
        "-cs", "spc216.gro",
        "-o", os.path.join(output_dir, "solvated.gro"),
        "-p", os.path.join(output_dir, "topol.top")
    ]
    subprocess.run(cmd, check=True)

    # Step 4: Add ions
    cmd = [
        "gmx", "grompp",
        "-f", "ions.mdp",
        "-c", os.path.join(output_dir, "solvated.gro"),
        "-p", os.path.join(output_dir, "topol.top"),
        "-o", os.path.join(output_dir, "ions.tpr")
    ]
    subprocess.run(cmd, check=True)

    cmd = [
        "gmx", "genion",
        "-s", os.path.join(output_dir, "ions.tpr"),
        "-o", os.path.join(output_dir, "solvated_ions.gro"),
        "-p", os.path.join(output_dir, "topol.top"),
        "-pname", "NA", "-nname", "CL", "-neutral"
    ]
    subprocess.run(cmd, check=True)

    # Step 5: Energy minimization
    cmd = [
        "gmx", "grompp",
        "-f", "minim.mdp",
        "-c", os.path.join(output_dir, "solvated_ions.gro"),
        "-p", os.path.join(output_dir, "topol.top"),
        "-o", os.path.join(output_dir, "em.tpr")
    ]
    subprocess.run(cmd, check=True)

    cmd = [
        "gmx", "mdrun",
        "-v", "-deffnm", os.path.join(output_dir, "em")
    ]
    subprocess.run(cmd, check=True)

    print(f"[MD] Molecular dynamics simulation completed. Results in {output_dir}")


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


def process_md_simulation(row, docking_dir, output_dir, db_connection_string):
    """
    Process a single MD simulation entry.
    """
    docking_id = row["docking_id"]
    ligand_cid = row["ligand_cid"]
    deg_id = row["deg_id"]
    binding_energy = row["binding_energy"]
    gene_name = row["gene_name"]
    sim_time = row.get("generic_sim_time", 100)
    temp = row.get("generic_temp", 300)
    pressure = row.get("generic_pressure", 1)
    solvent = row.get("generic_solvent", "water")

    # Paths
    protein_pdb = Path(docking_dir) / f"{docking_id}_protein_clean.pdb"
    ligand_pdb = Path(docking_dir) / f"{docking_id}_docked.pdb"
    output_pdb = Path(output_dir) / f"{docking_id}_complex.pdb"
    sim_output_dir = Path(output_dir) / f"docking_{gene_name}_{docking_id}"

    # Prepare protein-ligand complex
    combine_protein_ligand(protein_pdb, ligand_pdb, output_pdb)

    # Run MD simulation
    perform_md_simulation_with_gromacs(output_pdb, sim_output_dir, sim_time, temp, pressure)

    # Insert results into the database
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


def main():
    parser = argparse.ArgumentParser(description="Perform MD simulations using parameters from md_parameters.csv")
    parser.add_argument("--db_connection_string", required=True, help="Database connection string (SQLAlchemy format).")
    parser.add_argument("--md_param_file", required=True, help="Path to the MD parameters CSV file.")
    parser.add_argument("--docking_dir", required=True, help="Directory containing docking results.")
    parser.add_argument("--output_dir", required=True, help="Directory to store MD results.")
    parser.add_argument("--num_workers", type=int, default=8, help="Number of parallel workers.")
    args = parser.parse_args()

    # Load MD parameters
    md_params = parse_md_parameters(args.md_param_file)

    # Parallelize MD simulations
    tasks = [
        (row, args.docking_dir, args.output_dir, args.db_connection_string)
        for _, row in md_params.iterrows()
    ]

    with ProcessPoolExecutor(max_workers=args.num_workers) as executor:
        results = executor.map(lambda t: process_md_simulation(*t), tasks)

    print("[Pipeline] All simulations completed.")


if __name__ == "__main__":
    main()

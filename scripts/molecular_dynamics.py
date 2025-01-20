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


def process_md_simulation(row, docking_dir, output_dir, db_connection_string, experiment_id):
    """
    Process a single MD simulation entry.
    """

    # -------------------------------------------------------------------------
    # 1) EXTRACT COLUMNS FROM THE ROW
    # Ensure your CSV or DB result has uniprot_id, docking_id, gene_name, etc.
    # -------------------------------------------------------------------------
    docking_id = row["docking_id"]
    uniprot_id = row["uniprot_id"]  # <-- Make sure your CSV / DB has this column
    ligand_cid = row["ligand_cid"]
    deg_id = row["deg_id"]
    binding_energy = row["binding_energy"]
    gene_name = row["gene_name"]

    sim_time = row.get("generic_sim_time", 100)
    temp = row.get("generic_temp", 300)
    pressure = row.get("generic_pressure", 1)
    solvent = row.get("generic_solvent", "water")

    # -------------------------------------------------------------------------
    # 2) BUILD THE PATHS BASED ON uniprot_id (or docking_id if that matches)
    # -------------------------------------------------------------------------
    # Example: results/molecular_docking/<experiment_id>/<uniprot_id>_clean.pdb
    # Make sure these files exist or your script will fail.
    protein_pdb = Path(docking_dir) / str(experiment_id) / f"{uniprot_id}_clean.pdb"
    ligand_pdb  = Path(docking_dir) / str(experiment_id) / f"{uniprot_id}_docked.pdb"

    # The combined output PDB:
    output_pdb  = Path(docking_dir) / str(experiment_id) / f"{uniprot_id}_complex.pdb"

    # The final MD results folder.  E.g.:
    # results/molecular_dynamics/<experiment_id>/docking_<gene_name>_<docking_id>
    sim_output_dir = Path(output_dir) / "molecular_dynamics" / str(experiment_id) / f"docking_{gene_name}/{docking_id}"

    # -------------------------------------------------------------------------
    # 3) COMBINE AND RUN MD
    # -------------------------------------------------------------------------
    combine_protein_ligand(protein_pdb, ligand_pdb, output_pdb)
    perform_md_simulation_with_gromacs(output_pdb, sim_output_dir, sim_time, temp, pressure)

    # -------------------------------------------------------------------------
    # 4) INSERT RESULTS INTO DB
    # -------------------------------------------------------------------------
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
    parser.add_argument("--num_workers", type=int, default=8, help="Number of parallel workers.")
    args = parser.parse_args()

    # 1. Load MD parameters
    md_params = parse_md_parameters(args.md_param_file)
    docking_dir = os.path.join(args.output_dir, "molecular_docking")
    experiment_id = get_experiment_id(args.experiment_id)

    # 3. Create the tasks
    tasks = [
        (row, docking_dir, args.output_dir, args.db_connection_string, experiment_id)
        for _, row in md_params.iterrows()
    ]

    # ---------------------------------------------------------------------
    # 4. For clarity, let's run them in a simple for-loop
    #    (You can re-enable concurrency later if you want.)
    # ---------------------------------------------------------------------
    for t in tasks:
        process_md_simulation(*t)

    print("[Pipeline] All simulations completed.")


if __name__ == "__main__":
    main()

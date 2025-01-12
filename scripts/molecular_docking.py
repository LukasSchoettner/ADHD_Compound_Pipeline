import os
import subprocess
import requests
import argparse
import psycopg2
import pandas as pd


def fetch_therapeutic_targets(db_connection_string, experiment_id):
    """
    Fetch therapeutic targets for the given experiment ID, including disease gene IDs.

    Args:
        db_connection_string (str): Database connection string.
        experiment_id (int): Experiment ID.

    Returns:
        list: List of disease gene IDs from therapeutic targets.
    """
    query = """
    SELECT DISTINCT tt.disease_gene_id
    FROM therapeutic_targets tt
    WHERE tt.experiment_id = %s;
    """
    try:
        conn = psycopg2.connect(db_connection_string)
        therapeutic_targets = pd.read_sql_query(query, conn, params=(experiment_id,))
        conn.close()
        return therapeutic_targets['disease_gene_id'].tolist()
    except Exception as e:
        raise Exception(f"Error fetching therapeutic targets: {e}")


def fetch_pdb_ids_from_ensembl(ensembl_id):
    """
    Query RCSB PDB for PDB IDs corresponding to a given Ensembl ID.

    Args:
        ensembl_id (str): Ensembl ID.

    Returns:
        list: List of PDB IDs.
    """
    url = f"https://search.rcsb.org/rcsbsearch/v1/query"
    query_payload = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {"value": ensembl_id}
        },
        "return_type": "entry",
        "request_options": {"return_all_hits": True}
    }
    try:
        response = requests.post(url, json=query_payload)
        if response.status_code == 200:
            data = response.json()
            return [entry["identifier"] for entry in data.get("result_set", [])]
        else:
            print(f"Failed to fetch PDB IDs for Ensembl ID {ensembl_id}: {response.status_code}")
    except requests.RequestException as e:
        print(f"Error querying RCSB PDB for Ensembl ID {ensembl_id}: {e}")
    return []


def identify_docking_site(protein_pdb, output_dir):
    """
    Identify docking site automatically using AutoDockTools.

    Args:
        protein_pdb (str): Path to the protein structure file in PDB format.
        output_dir (str): Directory to store output grid parameters.

    Returns:
        tuple: Center (x, y, z) and size (x, y, z) of the docking grid box.
    """
    try:
        print(f"Identifying docking site for {protein_pdb}...")
        grid_file = os.path.join(output_dir, "grid_params.txt")
        subprocess.run([
            "pythonsh", "prepare_gpf4.py",
            "-r", protein_pdb,
            "-o", grid_file
        ], check=True)
        with open(grid_file, "r") as f:
            lines = f.readlines()
            center = None
            size = None
            for line in lines:
                if line.startswith("GRID_CENTER"):
                    center = tuple(map(float, line.split()[1:4]))
                if line.startswith("GRID_SIZE"):
                    size = tuple(map(float, line.split()[1:4]))
            if center and size:
                print(f"Docking site identified: Center={center}, Size={size}")
                return center, size
        raise Exception("Failed to identify docking site parameters.")
    except subprocess.CalledProcessError as e:
        raise Exception(f"AutoDockTools docking site identification failed: {e}")


def download_ligand(pubchem_id, output_file):
    """
    Download ligand 3D structure from PubChem.

    Args:
        pubchem_id (str): PubChem Compound ID (CID).
        output_file (str): Output file path.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{pubchem_id}/record/SDF/?record_type=3d"
    response = requests.get(url)
    if response.status_code == 200:
        with open(output_file, "wb") as f:
            f.write(response.content)
        print(f"Ligand file saved as {output_file}")
    else:
        raise Exception(f"Failed to download ligand. Status code: {response.status_code}")


def download_protein(pdb_id, output_file):
    """
    Download protein structure from RCSB PDB.

    Args:
        pdb_id (str): Protein Data Bank ID.
        output_file (str): File path to save the protein structure.
    """
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url)
    if response.status_code == 200:
        with open(output_file, "wb") as f:
            f.write(response.content)
        print(f"Protein structure saved as {output_file}")
    else:
        raise Exception(f"Failed to download protein structure for {pdb_id}. HTTP status: {response.status_code}")


def prepare_ligand(input_file, output_file):
    """
    Prepare ligand for docking.

    Args:
        input_file (str): Path to input ligand file in SDF format.
        output_file (str): Path to output ligand file in PDBQT format.
    """
    print("Preparing ligand...")
    subprocess.run(["obabel", input_file, "-O", "temp.pdb", "--gen3d"], check=True)
    subprocess.run(["obabel", "temp.pdb", "-O", "temp_minimized.pdb", "--minimize"], check=True)
    subprocess.run(["prepare_ligand", "-l", "temp_minimized.pdb", "-o", output_file, "-A", "hydrogens"], check=True)
    os.remove("temp.pdb")
    os.remove("temp_minimized.pdb")
    print(f"Ligand prepared: {output_file}")


def prepare_protein(input_file, output_file):
    """
    Prepare protein for docking.

    Args:
        input_file (str): Path to input protein file in PDB format.
        output_file (str): Path to output protein file in PDBQT format.
    """
    print(f"Preparing protein: {input_file} -> {output_file}")
    subprocess.run(["prepare_receptor", "-r", input_file, "-o", output_file, "-A", "hydrogens"], check=True)
    print(f"Protein prepared: {output_file}")


def perform_docking(vina_executable, protein_pdbqt, ligand_pdbqt, output_dir, center, size):
    """
    Perform molecular docking using AutoDock Vina.

    Args:
        vina_executable (str): Path to the AutoDock Vina executable.
        protein_pdbqt (str): Path to prepared protein file in PDBQT format.
        ligand_pdbqt (str): Path to prepared ligand file in PDBQT format.
        output_dir (str): Directory to store docking results.
        center (tuple): Coordinates for docking box center (x, y, z).
        size (tuple): Dimensions for docking box size (x, y, z).
    """
    print("Running docking simulation...")
    docking_command = [
        vina_executable,
        "--receptor", protein_pdbqt,
        "--ligand", ligand_pdbqt,
        "--center_x", str(center[0]),
        "--center_y", str(center[1]),
        "--center_z", str(center[2]),
        "--size_x", str(size[0]),
        "--size_y", str(size[1]),
        "--size_z", str(size[2]),
        "--out", os.path.join(output_dir, f"{os.path.basename(protein_pdbqt)}_docked.pdbqt"),
        "--log", os.path.join(output_dir, f"{os.path.basename(protein_pdbqt)}_docking.log"),
    ]
    result = subprocess.run(docking_command, capture_output=True, text=True)
    if result.returncode == 0:
        print("Docking completed successfully!")
    else:
        raise Exception(f"Docking failed: {result.stderr}")


def get_experiment_id(experiment_id):
    """
    Parse experiment_id from file if it is a file path.

    Args:
        experiment_id (str): Either the experiment ID as a string or a file path.

    Returns:
        int: Parsed experiment ID.
    """
    if os.path.isfile(experiment_id):
        with open(experiment_id, "r") as f:
            return int(f.read().strip())
    return int(experiment_id)


if __name__ == "__main__":
    # Parse arguments from Nextflow
    parser = argparse.ArgumentParser(description="Perform molecular docking")
    parser.add_argument("--ligand_cid", required=True, help="PubChem CID of the ligand")
    parser.add_argument("--db_connection_string", required=True, help="Database connection string")
    parser.add_argument("--experiment_id", type=str, required=True, help="Experiment ID or file containing the ID")
    parser.add_argument("--vina_executable", default="vina", help="Path to AutoDock Vina executable")
    parser.add_argument("--output_dir", default="docking_results", help="Directory to save docking results")
    parser.add_argument("--docking_params", required=True, help="Path to docking parameters file")
    args = parser.parse_args()

    # Prepare paths and directories
    os.makedirs(args.output_dir, exist_ok=True)
    ligand_sdf = os.path.join(args.output_dir, "ligand.sdf")
    ligand_pdbqt = os.path.join(args.output_dir, "ligand.pdbqt")

    experiment_id = get_experiment_id(args.experiment_id)

    try:
        # Fetch therapeutic targets
        ensembl_ids = fetch_therapeutic_targets(args.db_connection_string, experiment_id)

        # Map Ensembl IDs to PDB IDs
        pdb_ids = []
        for ensembl_id in ensembl_ids:
            pdb_ids.extend(fetch_pdb_ids_from_ensembl(ensembl_id))

        if not pdb_ids:
            raise Exception("No PDB IDs could be resolved from therapeutic targets.")

        # Download and prepare ligand
        download_ligand(args.ligand_cid, ligand_sdf)
        prepare_ligand(ligand_sdf, ligand_pdbqt)

        # Process each protein
        docking_params = pd.read_csv(args.docking_params)
        for _, row in docking_params.iterrows():
            pdb_id = row["Protein"]
            protein_pdb = os.path.join(args.output_dir, f"{pdb_id}.pdb")
            protein_pdbqt = os.path.join(args.output_dir, f"{pdb_id}.pdbqt")

            # Download and prepare protein
            download_protein(pdb_id, protein_pdb)
            prepare_protein(protein_pdb, protein_pdbqt)

            # Determine docking site
            if pd.notna(row["Docking Site Center"]) and pd.notna(row["Docking Site Size"]):
                center = tuple(map(float, row["Docking Site Center"].strip("()").split(",")))
                size = tuple(map(float, row["Docking Site Size"].strip("()").split(",")))
            else:
                center, size = identify_docking_site(protein_pdb, args.output_dir)

            # Perform docking
            perform_docking(
                vina_executable=args.vina_executable,
                protein_pdbqt=protein_pdbqt,
                ligand_pdbqt=ligand_pdbqt,
                output_dir=args.output_dir,
                center=center,
                size=size
            )

        print("Molecular docking completed successfully.")

    except Exception as e:
        print(f"Error during molecular docking: {e}")

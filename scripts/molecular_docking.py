import os
import subprocess
import requests
import argparse
import pandas as pd
from sqlalchemy import create_engine

###############################################################################
# 1. Database: fetch targets, including uniprot_id if present
###############################################################################

def fetch_therapeutic_targets(db_connection_string, experiment_id):
    """
    Fetch therapeutic targets for a given experiment_id, including disease_gene_id,
    uniprot_id (if in DB), and gene_name for reference.

    Returns a DataFrame with columns:
      - disease_gene_id
      - uniprot_id (might be NULL in DB)
      - gene_name
    """
    engine = create_engine(db_connection_string)
    query = """
        SELECT
            tt.disease_gene_id,
            tt.uniprot_id,
            dg.gene_name
        FROM therapeutic_targets tt
        INNER JOIN disease_genes dg
            ON tt.disease_gene_id = dg.disease_gene_id
        WHERE tt.experiment_id = %s
    """
    df = pd.read_sql_query(query, engine, params=(experiment_id,))
    engine.dispose()
    return df


###############################################################################
# 2. Ensembl → UniProt if uniprot_id is missing
###############################################################################

def fetch_uniprot_ids_from_ensembl(ensembl_id):
    """
    Given an Ensembl gene ID, fetch any UniProt primaryAccessions from UniProt's REST API.
    Returns a list of UniProt IDs (can be multiple isoforms).
    """
    url = f"https://rest.uniprot.org/uniprotkb/stream?query=ensembl:{ensembl_id}&format=json"
    try:
        resp = requests.get(url, timeout=15)
        if resp.status_code == 200:
            data = resp.json()
            return [entry['primaryAccession'] for entry in data.get('results', [])]
        else:
            print(f"[Ensembl→UniProt] Failed for {ensembl_id}: HTTP {resp.status_code}")
    except requests.RequestException as e:
        print(f"[Ensembl→UniProt] Request error for {ensembl_id}: {e}")
    return []


def select_canonical_uniprot(uniprot_ids):
    """
    Given a list of possible UniProt IDs (isoforms), try to pick the 'canonical' one.
    If no canonical is found, just return the first in the list.
    """
    for uid in uniprot_ids:
        meta_url = f"https://rest.uniprot.org/uniprotkb/{uid}.json"
        try:
            resp = requests.get(meta_url, timeout=15)
            if resp.status_code == 200:
                meta = resp.json()
                # If the 'isCanonical' flag is True, pick this isoform
                if meta.get("isCanonical", False):
                    return uid
        except requests.RequestException:
            pass
    # fallback: first if none is labeled canonical
    return uniprot_ids[0] if uniprot_ids else None


###############################################################################
# 3. UniProt → PDB
###############################################################################

def fetch_pdb_ids_for_uniprot(uniprot_id):
    """
    Fetch all PDB IDs for a UniProt ID from RCSB.
    Returns a list of PDB IDs (strings).
    """
    url = "https://search.rcsb.org/rcsbsearch/v1/query"
    payload = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {"value": uniprot_id}
        },
        "return_type": "entry",
        "request_options": {"return_all_hits": True}
    }
    try:
        resp = requests.post(url, json=payload, timeout=15)
        if resp.status_code == 200:
            data = resp.json()
            return [entry["identifier"] for entry in data.get("result_set", [])]
        else:
            print(f"[UniProt→PDB] Failed for {uniprot_id}: HTTP {resp.status_code}")
    except requests.RequestException as e:
        print(f"[UniProt→PDB] Request error for {uniprot_id}: {e}")
    return []


###############################################################################
# 4. File Downloads & Preparation
###############################################################################

def download_ligand(pubchem_id, output_file):
    """
    Download a 3D SDF of the ligand from PubChem by CID.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{pubchem_id}/record/SDF/?record_type=3d"
    resp = requests.get(url, timeout=15)
    if resp.status_code == 200:
        with open(output_file, "wb") as f:
            f.write(resp.content)
        print(f"[Ligand] Downloaded: {output_file}")
    else:
        raise RuntimeError(f"[Ligand] Download error for CID {pubchem_id}: HTTP {resp.status_code}")


def prepare_ligand(input_sdf, output_pdbqt):
    """
    Convert ligand to 3D coordinates, minimize, output as PDBQT with hydrogens.
    """
    print("[Ligand] Preparing ligand...")
    subprocess.run(["obabel", input_sdf, "-O", "temp.pdb", "--gen3d"], check=True)
    subprocess.run(["obabel", "temp.pdb", "-O", "temp_minimized.pdb", "--minimize"], check=True)
    subprocess.run([
        "prepare_ligand",
        "-l", "temp_minimized.pdb",
        "-o", output_pdbqt,
        "-A", "hydrogens"
    ], check=True)
    os.remove("temp.pdb")
    os.remove("temp_minimized.pdb")
    print(f"[Ligand] Prepared: {output_pdbqt}")


def download_protein(pdb_id, output_file):
    """
    Download a PDB file from RCSB.
    """
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    resp = requests.get(url, timeout=15)
    if resp.status_code == 200:
        with open(output_file, "wb") as f:
            f.write(resp.content)
        print(f"[Protein] Downloaded: {output_file}")
    else:
        raise RuntimeError(f"[Protein] Download error for PDB {pdb_id}: HTTP {resp.status_code}")


def prepare_protein(input_pdb, output_pdbqt):
    """
    Run prepare_receptor to add hydrogens & convert to PDBQT for docking.
    """
    print("[Protein] Preparing receptor...")
    subprocess.run([
        "prepare_receptor",
        "-r", input_pdb,
        "-o", output_pdbqt,
        "-A", "hydrogens"
    ], check=True)
    print(f"[Protein] Prepared: {output_pdbqt}")


###############################################################################
# 5. Docking
###############################################################################

def perform_docking(vina_exec, protein_pdbqt, ligand_pdbqt, output_dir, center, size):
    """
    Run AutoDock Vina with the specified bounding box (center, size).
    """
    out_pdbqt = os.path.join(output_dir, f"{os.path.basename(protein_pdbqt)}_docked.pdbqt")
    log_file = os.path.join(output_dir, f"{os.path.basename(protein_pdbqt)}_docking.log")

    cmd = [
        vina_exec,
        "--receptor", protein_pdbqt,
        "--ligand", ligand_pdbqt,
        "--center_x", str(center[0]),
        "--center_y", str(center[1]),
        "--center_z", str(center[2]),
        "--size_x", str(size[0]),
        "--size_y", str(size[1]),
        "--size_z", str(size[2]),
        "--out", out_pdbqt,
        "--log", log_file
    ]
    print(f"[Docking] Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"[Docking] Failed: {result.stderr}")
    print("[Docking] Completed successfully!")


###############################################################################
# 6. Possibly parse experiment_id from file
###############################################################################

def get_experiment_id(experiment_id):
    """
    If experiment_id is a path, read from file. Otherwise parse as int directly.
    """
    if os.path.isfile(experiment_id):
        with open(experiment_id, "r") as f:
            return int(f.read().strip())
    return int(experiment_id)


###############################################################################
# 7. Main: bring it all together
###############################################################################

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Automatic ADHD drug discovery pipeline")
    parser.add_argument("--db_connection_string", required=True,
                        help="SQLAlchemy-style DB URI (e.g. postgresql+psycopg2://user:pw@host/db)")
    parser.add_argument("--experiment_id", required=True,
                        help="Experiment ID integer or file containing the ID")
    parser.add_argument("--ligand_cid", required=True,
                        help="PubChem CID for the ligand to dock")
    parser.add_argument("--vina_executable", default="vina",
                        help="Path to the 'vina' executable")
    parser.add_argument("--output_dir", default="docking_results",
                        help="Directory to store output PDB/PDBQT files")
    parser.add_argument("--docking_params", required=False,
                        help="CSV with 'Protein','Docking Site Center','Docking Site Size'. If omitted, use fallback box.")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Convert experiment_id if it's a file
    exp_id = get_experiment_id(args.experiment_id)

    try:
        # 1) Fetch therapeutic targets: disease_gene_id + uniprot_id
        targets_df = fetch_therapeutic_targets(args.db_connection_string, exp_id)
        if targets_df.empty:
            raise RuntimeError("No therapeutic targets found for this experiment.")

        # 2) For each target, if uniprot_id is missing, do Ensembl → UniProt → canonical pick
        for idx, row in targets_df.iterrows():
            if not row["uniprot_id"] or pd.isna(row["uniprot_id"]):
                ensembl_id = row["disease_gene_id"]
                possible_uniprots = fetch_uniprot_ids_from_ensembl(ensembl_id)
                chosen_uniprot = select_canonical_uniprot(possible_uniprots)
                targets_df.at[idx, "uniprot_id"] = chosen_uniprot  # store in memory
                print(f"[Auto] Resolved {ensembl_id} → {chosen_uniprot}")

        # 3) For each uniprot, fetch PDB IDs, pick first one automatically
        #    (or you could store them all and loop further)
        all_pdbs = []
        for idx, row in targets_df.iterrows():
            uni_id = row["uniprot_id"]
            if not uni_id or pd.isna(uni_id):
                print(f"[Warning] Could not resolve UniProt for row:\n{row}")
                continue
            pdb_list = fetch_pdb_ids_for_uniprot(uni_id)
            if pdb_list:
                chosen_pdb = pdb_list[0]
                all_pdbs.append(chosen_pdb)
                print(f"[Auto] {uni_id} → picking PDB {chosen_pdb}")
            else:
                print(f"[Warning] No PDBs found for {uni_id}")

        if not all_pdbs:
            raise RuntimeError("No PDB IDs resolved for any targets.")

        # 4) Download and prepare the ligand
        ligand_sdf = os.path.join(args.output_dir, "ligand.sdf")
        ligand_pdbqt = os.path.jo

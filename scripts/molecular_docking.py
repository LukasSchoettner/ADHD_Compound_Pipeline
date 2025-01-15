import os
import subprocess
import requests
import argparse
import pandas as pd
from sqlalchemy import create_engine

from Bio.PDB import PDBParser, PDBIO
import numpy as np


###############################################################################
# 1. Database Queries
###############################################################################

def fetch_therapeutic_targets(db_connection_string, experiment_id):
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
# 2. Ligand Handling with Open Babel
###############################################################################

def download_ligand(pubchem_cid, output_sdf):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{pubchem_cid}/record/SDF/?record_type=3d"
    resp = requests.get(url, timeout=15)
    if resp.status_code == 200:
        with open(output_sdf, "wb") as f:
            f.write(resp.content)
        print(f"[Ligand] Downloaded to {output_sdf}")
    else:
        raise RuntimeError(f"[Ligand] Download error for CID={pubchem_cid}. HTTP {resp.status_code}")

def prepare_ligand_with_openbabel(input_sdf, output_pdbqt):
    print("[Ligand] Generating PDB + Minimization with Open Babel...")
    subprocess.run(["obabel", input_sdf, "-O", "temp_ligand.pdb",
                    "--gen3d", "--minimize"], check=True)

    print("[Ligand] Converting minimized PDB -> PDBQT (with partial charges)...")
    subprocess.run([
        "obabel", "temp_ligand.pdb",
        "-O", output_pdbqt,
        "--partialcharge", "gasteiger",
        "-xh"
    ], check=True)

    os.remove("temp_ligand.pdb")
    print(f"[Ligand] Prepared: {output_pdbqt}")


###############################################################################
# 3. Protein Fetching & Fallback
###############################################################################

def fetch_alphafold_pdb(uniprot_id, output_pdb):
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
    resp = requests.get(url, timeout=15)
    if resp.status_code == 200:
        with open(output_pdb, "w") as f:
            f.write(resp.text)
        print(f"[AlphaFold] Downloaded {uniprot_id} -> {output_pdb}")
        return True
    else:
        print(f"[AlphaFold] No structure or error for {uniprot_id} (HTTP {resp.status_code}).")
        return False

def fetch_rcsb_structure(uniprot_id, output_pdb):
    """
    Attempt to find a PDB referencing this UniProt and download it from RCSB.
    """
    search_payload = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.accession",
                "operator": "exact_match",
                "value": uniprot_id
            }
        },
        "return_type": "entry",
        "request_options": {"return_all_hits": True}
    }
    url_search = "https://search.rcsb.org/rcsbsearch/v1/query"

    try:
        search_resp = requests.post(url_search, json=search_payload, timeout=15)
        if search_resp.status_code == 200:
            data = search_resp.json()
            results = data.get("result_set", [])
            if not results:
                print(f"[RCSB] No PDB found for {uniprot_id}.")
                return False
            pdb_id = results[0]["identifier"]
            print(f"[RCSB] Found PDB ID={pdb_id} for {uniprot_id}. Downloading...")

            url_download = f"https://files.rcsb.org/download/{pdb_id}.pdb"
            dresp = requests.get(url_download, timeout=15)
            if dresp.status_code == 200:
                with open(output_pdb, "wb") as f:
                    f.write(dresp.content)
                print(f"[RCSB] Downloaded {pdb_id} -> {output_pdb}")
                return True
            else:
                print(f"[RCSB] Could not download {pdb_id} (HTTP {dresp.status_code}).")
                return False
        else:
            print(f"[RCSB] Search request failed (HTTP {search_resp.status_code}).")
            return False
    except Exception as e:
        print(f"[RCSB] error: {e}")
        return False

def quick_pdb_check(pdb_file):
    from Bio.PDB import PDBParser
    parser = PDBParser(QUIET=True)
    try:
        _ = parser.get_structure("protein", pdb_file)
        return True
    except Exception as e:
        print(f"[ParseCheck] BioPython parse error: {e}")
        return False

###############################################################################
# 4. Cleaning & Converting Protein
###############################################################################

def clean_pdb_biopython(input_pdb, output_pdb):
    from Bio.PDB import PDBParser, PDBIO
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", input_pdb)
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb)

def prepare_protein_with_openbabel(input_pdb, output_pdbqt):
    print("[Protein] Converting PDB -> PDBQT with Open Babel...")
    subprocess.run([
        "obabel", input_pdb,
        "-O", output_pdbqt,
        "--partialcharge", "gasteiger",
        "-xh"
    ], check=True)
    print(f"[Protein] Prepared: {output_pdbqt}")


def remove_ligand_tags_in_receptor(receptor_pdbqt):
    """
    Some tools accidentally add ligand-style 'ROOT', 'BRANCH', 'ENDBRANCH', 'ENDROOT' lines
    in a receptor file. This function removes them so Vina doesn't parse them as errors.
    """
    cleaned_lines = []
    with open(receptor_pdbqt, "r") as infile:
        for line in infile:
            # skip any line that starts with these ligand-only tags
            if line.startswith(("ROOT", "ENDROOT", "BRANCH", "ENDBRANCH", "TORSDOF")):
                continue
            cleaned_lines.append(line)
    # rewrite the file
    with open(receptor_pdbqt, "w") as outfile:
        outfile.writelines(cleaned_lines)


###############################################################################
# 5. Auto-Detect Pocket
###############################################################################

def autodetect_pocket(protein_pdb):
    base_name = os.path.splitext(os.path.basename(protein_pdb))[0]
    out_dir = f"{base_name}_out"
    cmd = ["fpocket", "-f", protein_pdb]
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"[AutoDetect] fpocket failed: {e}")
        return (0.0,0.0,0.0), (20.0,20.0,20.0)

    pockets_dir = os.path.join(f"{base_name}_out", "pockets")
    pocket0 = os.path.join(pockets_dir, "pocket0_vert.pdb")
    if not os.path.isfile(pocket0):
        print("[AutoDetect] No pocket0 found. Fallback used.")
        return (0.0,0.0,0.0), (20.0,20.0,20.0)

    coords = []
    with open(pocket0, "r") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM", "TORSDOF")):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append([x, y, z])
    if not coords:
        print("[AutoDetect] pocket0 empty. Fallback used.")
        return (0.0,0.0,0.0), (20.0,20.0,20.0)

    arr = np.array(coords)
    min_xyz = arr.min(axis=0)
    max_xyz = arr.max(axis=0)
    center = tuple((min_xyz + max_xyz)/2.0)
    size   = tuple(max_xyz - min_xyz)
    print(f"[AutoDetect] center={center}, size={size}")
    return center, size

def get_docking_box(row, protein_pdb,
                    fallback_center=(0.0,0.0,0.0),
                    fallback_size=(20.0,20.0,20.0)):
    center_str = row.get("Docking Site Center", "") or ""
    size_str   = row.get("Docking Site Size", "") or ""

    if center_str.strip() and size_str.strip():
        try:
            cx, cy, cz = [float(x.strip()) for x in center_str.split(",")]
            sx, sy, sz = [float(x.strip()) for x in size_str.split(",")]
            return (cx, cy, cz), (sx, sy, sz)
        except ValueError:
            print("[Warning] parse error for center/size. Using fallback.")
            return fallback_center, fallback_size
    else:
        center, size = autodetect_pocket(protein_pdb)
        return center, size


###############################################################################
# 6. Docking with AutoDock Vina
###############################################################################

def perform_docking(vina_exec, protein_pdbqt, ligand_pdbqt, output_dir, center, size):
    out_pdbqt = os.path.join(output_dir, f"{os.path.basename(protein_pdbqt)}_docked.pdbqt")
    log_file  = os.path.join(output_dir, f"{os.path.basename(protein_pdbqt)}_docking.log")

    # remove any ROOT/BRANCH lines from the receptor just in case
    remove_ligand_tags_in_receptor(protein_pdbqt)

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
    subprocess.run(cmd, check=True)
    print("[Docking] Completed successfully!")


###############################################################################
# 7. Main
###############################################################################

def parse_experiment_id(exp_id):
    if os.path.isfile(exp_id):
        with open(exp_id, "r") as f:
            return f.read().strip()
    return exp_id

def fetch_structure_fallback(uniprot_id, output_pdb):
    """
    1) Attempt alphaFold, parse check
    2) If fail, attempt RCSB
    """
    alpha_ok = fetch_alphafold_pdb(uniprot_id, output_pdb)
    parse_ok = False
    if alpha_ok:
        if quick_pdb_check(output_pdb):
            return True
        else:
            print(f"[Alpha] {uniprot_id} parse fail. Trying RCSB fallback.")
    rcsb_ok = fetch_rcsb_structure(uniprot_id, output_pdb)
    if rcsb_ok and quick_pdb_check(output_pdb):
        return True
    return False

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="ADHD pipeline fallback & remove ROOT from receptor.")
    parser.add_argument("--db_connection_string", required=True)
    parser.add_argument("--experiment_id", required=True)
    parser.add_argument("--ligand_cid", required=True)
    parser.add_argument("--vina_executable", default="vina")
    parser.add_argument("--output_dir", default="/home/scmbag/Desktop/ADHD_Compound_Pipeline/results/molecular_docking")
    parser.add_argument("--docking_params", required=False,
                        help="Path to a CSV file with docking parameters.")
    args = parser.parse_args()

    output_dir = os.path.join(args.output_dir, f"{args.experiment_id}")

    os.makedirs(output_dir, exist_ok=True)

    exp_id = parse_experiment_id(args.experiment_id)
    targets = fetch_therapeutic_targets(args.db_connection_string, exp_id)
    if targets.empty:
        raise RuntimeError("No therapeutic targets found for this experiment.")

    # 1) Prepare ligand
    ligand_sdf   = os.path.join(output_dir, "ligand.sdf")
    ligand_pdbqt = os.path.join(output_dir, "ligand.pdbqt")
    download_ligand(args.ligand_cid, ligand_sdf)
    prepare_ligand_with_openbabel(ligand_sdf, ligand_pdbqt)

    # 2) Docking params if any
    if args.docking_params and os.path.exists(args.docking_params):
        df_params = pd.read_csv(args.docking_params)
    else:
        df_params = None

    # 3) For each target
    for idx, row in targets.iterrows():
        uniprot_id = row["uniprot_id"]
        gene_id    = row["disease_gene_id"]
        if not uniprot_id or pd.isna(uniprot_id):
            print(f"[Warning] Missing UniProt ID for row={row}, skipping.")
            continue

        protein_pdb    = os.path.join(output_dir, f"{uniprot_id}.pdb")
        protein_clean  = os.path.join(output_dir, f"{uniprot_id}_clean.pdb")
        protein_pdbqt  = os.path.join(output_dir, f"{uniprot_id}.pdbqt")

        # fallback approach
        if not fetch_structure_fallback(uniprot_id, protein_pdb):
            print(f"[Error] No structure found for {uniprot_id} on AlphaFold or RCSB. Skipping.")
            continue

        # Clean
        clean_pdb_biopython(protein_pdb, protein_clean)
        # Convert
        prepare_protein_with_openbabel(protein_clean, protein_pdbqt)

        # parse docking box
        if df_params is not None:
            match = df_params[df_params["disease_gene_id"] == gene_id]
            if not match.empty:
                row_params = match.iloc[0]
                center, size = get_docking_box(row_params, protein_clean)
            else:
                center, size = autodetect_pocket(protein_clean)
        else:
            center, size = autodetect_pocket(protein_clean)

        # docking (this calls remove_ligand_tags_in_receptor internally)
        perform_docking(
            vina_exec=args.vina_executable,
            protein_pdbqt=protein_pdbqt,
            ligand_pdbqt=ligand_pdbqt,
            output_dir=output_dir,
            center=center,
            size=size
        )

    print("[Pipeline] Molecular docking completed successfully.")

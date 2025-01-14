import os
import subprocess
import requests
import argparse
import pandas as pd
from sqlalchemy import create_engine
import numpy as np  # if you want to handle arrays for autodetect

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
    """
    Download a 3D SDF from PubChem for the given CID.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{pubchem_cid}/record/SDF/?record_type=3d"
    resp = requests.get(url, timeout=15)
    if resp.status_code == 200:
        with open(output_sdf, "wb") as f:
            f.write(resp.content)
        print(f"[Ligand] Downloaded to {output_sdf}")
    else:
        raise RuntimeError(f"[Ligand] Download error for CID={pubchem_cid}. HTTP {resp.status_code}")

def prepare_ligand_with_openbabel(input_sdf, output_pdbqt):
    """
    Use Open Babel to:
      1) SDF -> PDB (generate 3D, minimize)
      2) PDB -> PDBQT (partial charges, add hydrogens)
    """
    print("[Ligand] Generating PDB + Minimization with Open Babel...")
    subprocess.run(["obabel", input_sdf, "-O", "temp_ligand.pdb",
                    "--gen3d", "--minimize"], check=True)

    print("[Ligand] Converting minimized PDB -> PDBQT (with partial charges)...")
    subprocess.run([
        "obabel", "temp_ligand.pdb",
        "-O", output_pdbqt,
        "--addcharges",
        "-xh"
    ], check=True)

    os.remove("temp_ligand.pdb")
    print(f"[Ligand] Prepared: {output_pdbqt}")

###############################################################################
# 3. Protein Fetching & Preparation
###############################################################################

def fetch_alphafold_pdb(uniprot_id, output_pdb):
    """
    Download the AlphaFold-predicted structure for the given UniProt ID.
    """
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
    resp = requests.get(url, timeout=15)
    if resp.status_code == 200:
        with open(output_pdb, "w") as f:
            f.write(resp.text)
        print(f"[AlphaFold] Downloaded {uniprot_id} -> {output_pdb}")
    else:
        raise RuntimeError(f"[AlphaFold] No structure found for {uniprot_id}. HTTP {resp.status_code}")

def prepare_protein_with_openbabel(input_pdb, output_pdbqt):
    """
    Use Open Babel to convert a PDB to PDBQT, adding partial charges & hydrogens.
    """
    print("[Protein] Converting PDB -> PDBQT with Open Babel...")
    subprocess.run([
        "obabel", input_pdb,
        "-O", output_pdbqt,
        "--addcharges",
        "-xh"
    ], check=True)

    print(f"[Protein] Prepared: {output_pdbqt}")

###############################################################################
# 4. Auto-Detect Pocket (Placeholder Example)
###############################################################################

def autodetect_pocket(protein_pdb):
    """
    Example using 'fpocket' to detect the main pocket and compute bounding box.
    Adjust as needed for your tool (P2Rank, AutoSite, etc.).
    """
    base_name = os.path.splitext(os.path.basename(protein_pdb))[0]
    out_dir = f"{base_name}_out"

    # 1) Run fpocket
    cmd = ["fpocket", "-f", protein_pdb]
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"[AutoDetect] fpocket failed: {e}")
        return (0.0,0.0,0.0), (20.0,20.0,20.0)

    # 2) Check pocket0
    pockets_dir = os.path.join(f"{base_name}_out", "pockets")
    pocket0 = os.path.join(pockets_dir, "pocket0_vert.pdb")
    if not os.path.isfile(pocket0):
        print(f"[AutoDetect] No pocket0 found. Fallback used.")
        return (0.0,0.0,0.0), (20.0,20.0,20.0)

    coords = []
    with open(pocket0, "r") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append([x, y, z])

    if not coords:
        print("[AutoDetect] pocket0 empty. Fallback used.")
        return (0.0,0.0,0.0), (20.0,20.0,20.0)

    arr = np.array(coords)
    min_xyz = arr.min(axis=0)  # [minx, miny, minz]
    max_xyz = arr.max(axis=0)  # [maxx, maxy, maxz]
    center = tuple((min_xyz + max_xyz)/2.0)
    size   = tuple(max_xyz - min_xyz)
    print(f"[AutoDetect] center={center}, size={size}")
    return center, size

###############################################################################
# 5. Determine Docking Box
###############################################################################

def get_docking_box(row, protein_pdb,
                    fallback_center=(0.0,0.0,0.0),
                    fallback_size=(20.0,20.0,20.0)):
    """
    If 'Docking Site Center' and 'Docking Site Size' are non-empty,
    parse them. Otherwise call autodetect_pocket or fallback.
    """
    center_str = row.get("Docking Site Center", "") or ""
    size_str   = row.get("Docking Site Size", "") or ""

    if center_str.strip() and size_str.strip():
        # The user entered something like "10, 10, 10"
        try:
            cx, cy, cz = [float(x.strip()) for x in center_str.split(",")]
            sx, sy, sz = [float(x.strip()) for x in size_str.split(",")]
            return (cx, cy, cz), (sx, sy, sz)
        except ValueError:
            print("[Warning] parse error for center/size. Using fallback.")
            return fallback_center, fallback_size
    else:
        # auto-detect
        center, size = autodetect_pocket(protein_pdb)
        return center, size

###############################################################################
# 6. Docking with AutoDock Vina
###############################################################################

def perform_docking(vina_exec, protein_pdbqt, ligand_pdbqt, output_dir, center, size):
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

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Pipeline using Open Babel for ligand/protein prep & Vina for docking.")
    parser.add_argument("--db_connection_string", required=True)
    parser.add_argument("--experiment_id", required=True)
    parser.add_argument("--ligand_cid", required=True)
    parser.add_argument("--vina_executable", default="vina")
    parser.add_argument("--output_dir", default="docking_results")
    parser.add_argument("--docking_params", required=False,
                        help="Path to a CSV file with docking parameters.")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # 1) Convert experiment_id if file
    exp_id = parse_experiment_id(args.experiment_id)

    # 2) Fetch targets
    targets = fetch_therapeutic_targets(args.db_connection_string, exp_id)
    if targets.empty:
        raise RuntimeError("No therapeutic targets found for this experiment.")

    # 3) Prepare ligand with openbabel
    ligand_sdf   = os.path.join(args.output_dir, "ligand.sdf")
    ligand_pdbqt = os.path.join(args.output_dir, "ligand.pdbqt")
    download_ligand(args.ligand_cid, ligand_sdf)
    prepare_ligand_with_openbabel(ligand_sdf, ligand_pdbqt)

    # 4) Read docking parameters if provided
    if args.docking_params and os.path.exists(args.docking_params):
        df_params = pd.read_csv(args.docking_params)
    else:
        df_params = None

    # 5) For each target: fetch protein, prepare, auto-detect or parse box, run docking
    for idx, row in targets.iterrows():
        uniprot_id = row["uniprot_id"]
        gene_id    = row["disease_gene_id"]

        if not uniprot_id or pd.isna(uniprot_id):
            print(f"[Warning] Missing UniProt ID for gene_id={gene_id}, skipping.")
            continue

        protein_pdb    = os.path.join(args.output_dir, f"{uniprot_id}.pdb")
        protein_pdbqt  = os.path.join(args.output_dir, f"{uniprot_id}.pdbqt")

        # fetch from alphafold
        fetch_alphafold_pdb(uniprot_id, protein_pdb)

        # prepare
        prepare_protein_with_openbabel(protein_pdb, protein_pdbqt)

        # parse docking box from CSV or autodetect
        if df_params is not None:
            match = df_params[df_params["disease_gene_id"] == gene_id]
            if not match.empty:
                row_params = match.iloc[0]
                center, size = get_docking_box(row_params, protein_pdb)
            else:
                # fallback or autodetect
                center, size = autodetect_pocket(protein_pdb)
        else:
            # no CSV at all => autodetect
            center, size = autodetect_pocket(protein_pdb)

        # docking
        perform_docking(
            vina_exec=args.vina_executable,
            protein_pdbqt=protein_pdbqt,
            ligand_pdbqt=ligand_pdbqt,
            output_dir=args.output_dir,
            center=center,
            size=size
        )

    print("[Pipeline] Molecular docking completed successfully.")

import os
import shutil
import subprocess
import requests
import argparse
import pandas as pd
import traceback
import re
import concurrent

import yaml
from sqlalchemy import create_engine, text

from Bio.PDB import PDBParser, PDBIO, Select
import numpy as np
from pathlib import Path  # <-- added

class KeepAll(Select):
    def accept_residue(self, residue):
        return True
###############################################################################
# 1. Database Queries
###############################################################################

def fetch_therapeutic_targets(db_connection_string, experiment_id):
    engine = create_engine(db_connection_string)
    query = """
        SELECT
            tt.disease_gene_id,
            tt.uniprot_id,
            tt.deg_id,
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

def download_ligand(pubchem_cid, output_sdf: Path):
    """
    Download the ligand SDF from PubChem, saving to 'output_sdf' (Path).
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{pubchem_cid}/record/SDF/?record_type=3d"
    resp = requests.get(url, timeout=15)
    if resp.status_code == 200:
        output_sdf.write_bytes(resp.content)
        print(f"[Ligand] Downloaded to {output_sdf}")
    else:
        raise RuntimeError(f"[Ligand] Download error for CID={pubchem_cid}. HTTP {resp.status_code}")

def fetch_compound_properties(pubchem_cid):
    """
    Fetch detailed compound properties from PubChem for a given CID.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{pubchem_cid}/property/MolecularWeight,MolecularFormula,IUPACName/JSON"
    resp = requests.get(url, timeout=15)
    if resp.status_code == 200:
        data = resp.json()
        if "PropertyTable" in data and "Properties" in data["PropertyTable"]:
            properties = data["PropertyTable"]["Properties"][0]
            return {
                "pubchem_id": pubchem_cid,
                "compound_name": properties.get("IUPACName", "Unknown"),
                "molecular_weight": properties.get("MolecularWeight", None),
                "chemical_structure": None  # Placeholder if needed
            }
    raise RuntimeError(f"[PubChem] Error fetching properties for CID={pubchem_cid}. HTTP {resp.status_code}")

def prepare_ligand_with_openbabel(input_sdf: Path, output_pdbqt: Path):
    """
    Generate a minimized PDB from the input SDF, then convert to PDBQT.
    """
    print("[Ligand] Generating PDB + Minimization with Open Babel...")

    # We'll create a temporary PDB in the same directory as 'input_sdf'
    temp_ligand = input_sdf.with_name("temp_ligand.pdb")

    subprocess.run([
        "obabel",
        str(input_sdf),
        "-O", str(temp_ligand),
        "--gen3d",
        "--minimize"
    ], check=True)

    print("[Ligand] Converting minimized PDB -> PDBQT (with partial charges)...")
    subprocess.run([
        "obabel",
        str(temp_ligand),
        "-O", str(output_pdbqt),
        "--partialcharge", "gasteiger",
        "-xh"
    ], check=True)

    # Clean up the temp ligand file
    if temp_ligand.exists():
        temp_ligand.unlink()
    print(f"[Ligand] Prepared: {output_pdbqt}")

###############################################################################
# 3. Protein Fetching & Fallback
###############################################################################

def fetch_alphafold_pdb(uniprot_id: str, output_pdb: Path):
    """
    Attempt to download the AlphaFold model for 'uniprot_id' and save to output_pdb.
    Returns True on success, False otherwise.
    """
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
    resp = requests.get(url, timeout=15)
    if resp.status_code == 200:
        output_pdb.write_text(resp.text)
        print(f"[AlphaFold] Downloaded {uniprot_id} -> {output_pdb}")
        return True
    else:
        print(f"[AlphaFold] No structure or error for {uniprot_id} (HTTP {resp.status_code}).")
        return False

def fetch_rcsb_structure(uniprot_id: str, output_pdb: Path):
    """
    Attempt to find a PDB referencing this UniProt and download it from RCSB.
    """
    import json
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
                output_pdb.write_bytes(dresp.content)
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

def quick_pdb_check(pdb_file: Path):
    """
    Quick parse check on 'pdb_file' to see if BioPython can read it.
    """
    parser = PDBParser(QUIET=True)
    try:
        _ = parser.get_structure("protein", str(pdb_file))
        return True
    except Exception as e:
        print(f"[ParseCheck] BioPython parse error: {e}")
        return False

###############################################################################
# 4. Cleaning & Converting Protein
###############################################################################

def clean_pdb_biopython(input_pdb: Path, output_pdb: Path):
    """
    Load 'input_pdb' with BioPython, then save a 'clean' version to 'output_pdb'.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", str(input_pdb))
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(output_pdb), select=KeepAll())

def prepare_protein_with_openbabel(input_pdb: Path, output_pdbqt: Path):
    """
    Use obabel to convert the protein from PDB to PDBQT.
    """
    print("[Protein] Converting PDB -> PDBQT with Open Babel...")
    subprocess.run([
        "obabel", str(input_pdb),
        "-O", str(output_pdbqt),
        "--partialcharge", "gasteiger",
        "-xh"
    ], check=True)
    print(f"[Protein] Prepared: {output_pdbqt}")

def remove_ligand_tags_in_receptor(receptor_pdbqt: Path):
    """
    Some tools accidentally add ligand-style 'ROOT', 'BRANCH', 'ENDBRANCH', 'ENDROOT' lines
    in a receptor file. This function removes them so Vina doesn't parse them as errors.
    """
    cleaned_lines = []
    with receptor_pdbqt.open("r") as infile:
        for line in infile:
            if line.startswith(("ROOT", "ENDROOT", "BRANCH", "ENDBRANCH", "TORSDOF")):
                continue
            cleaned_lines.append(line)
    # rewrite the file
    with receptor_pdbqt.open("w") as outfile:
        outfile.writelines(cleaned_lines)

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
        # If there's no MODEL, maybe it's already a single-model PDB. We can just copy it.
        print(f"[Warning] No MODEL/ENDMDL found in {input_pdb}. Copying whole PDB.")
        shutil.copy2(input_pdb, output_pdb)

def merge_protein_ligand(protein_pdb: Path, ligand_pdb: Path, output_pdb: Path):
    """
    Merge the receptor coordinates from protein_pdb with ligand coordinates from ligand_pdb.
    The result is saved to output_pdb.
    """
    with protein_pdb.open("r") as pf, output_pdb.open("w") as of:
        for line in pf:
            # Skip any END or similar termination lines in the protein file
            if line.startswith("END"):
                continue
            of.write(line)
        with ligand_pdb.open("r") as lf:
            for line in lf:
                # Write only ATOM/HETATM lines from the ligand file
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    of.write(line)
        of.write("END\n")
    print(f"[Merge] Created merged protein-ligand PDB: {output_pdb}")


###############################################################################
# 5. Auto-Detect Pocket
###############################################################################
import glob

def autodetect_pocket(protein_pdb: Path):
    """
    Runs fpocket on 'protein_pdb', parses the 'base_name_info.txt' to find
    the best (highest Druggability Score) pocketN, then loads 'pocketN_*.pdb'
    to compute bounding box center & size.
    """
    base_name = protein_pdb.stem  # e.g. 'P04798_clean'
    out_dir = protein_pdb.parent / f"{base_name}_out"

    # 1) Run fpocket
    cmd = ["fpocket", "-f", str(protein_pdb)]
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"[AutoDetect] fpocket failed: {e}")
        return (0.0, 0.0, 0.0), (20.0, 20.0, 20.0)

    # 2) Parse the info file to find the best pocket number
    info_file = out_dir / f"{base_name}_info.txt"  # e.g. P04798_clean_out/P04798_clean_info.txt
    best_pocket_num = choose_best_pocket_number(info_file)

    # 3) Find the actual pocket file, e.g. "pocket6_*.pdb"
    pockets_dir = out_dir / "pockets"
    pattern = str(pockets_dir / f"pocket{best_pocket_num}_*.pdb")
    candidates = glob.glob(pattern)
    if not candidates:
        print(f"[AutoDetect] No pocket{best_pocket_num}_*.pdb found. Fallback used.")
        return (0.0, 0.0, 0.0), (20.0, 20.0, 20.0)

    pocket_file = Path(candidates[0])
    print(f"[AutoDetect] Best pocket = {best_pocket_num}, file = {pocket_file}")

    # 4) Parse that file’s coordinates
    coords = []
    with pocket_file.open("r") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append([x, y, z])

    if not coords:
        print(f"[AutoDetect] pocket{best_pocket_num} file empty. Fallback used.")
        return (0.0, 0.0, 0.0), (20.0, 20.0, 20.0)

    arr = np.array(coords)
    min_xyz = arr.min(axis=0)
    max_xyz = arr.max(axis=0)
    center = tuple((min_xyz + max_xyz) / 2.0)
    size   = tuple(max_xyz - min_xyz)
    print(f"[AutoDetect] center={center}, size={size}")
    return center, size


def get_docking_box(row: dict, protein_pdb: Path,
                    fallback_center=(0.0,0.0,0.0),
                    fallback_size=(20.0,20.0,20.0)):
    """
    If row has manual center/size, parse them. Otherwise call autodetect_pocket().
    """
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

def choose_best_pocket_number(info_file):
    """
    Parse an fpocket info file (like 'P04798_clean_info.txt') to find
    the pocket with the highest 'Druggability Score'.
    Returns an integer pocket_number (e.g. 6 if 'Pocket 6' is best).
    If none found, returns 1 by default.
    """
    if not info_file.is_file():
        print(f"[choose_best_pocket_number] Info file not found: {info_file}")
        return 1

    # Regex to match lines like "Pocket 6 :" or "Pocket 12 :"
    pocket_header_re = re.compile(r"^Pocket\s+(\d+)\s*:\s*$")
    # Regex to match lines like "Druggability Score :    0.988"
    drugg_score_re   = re.compile(r"^\s*Druggability Score\s*:\s*([\d\.]+)")

    best_pocket_num = 1
    best_score      = -999.0

    current_pocket_num = None

    with info_file.open("r") as f:
        for line in f:
            line = line.strip()

            # If it's a pocket header, e.g. "Pocket 6 :"
            match_header = pocket_header_re.match(line)
            if match_header:
                current_pocket_num = int(match_header.group(1))

            # If we see the Druggability Score
            match_score = drugg_score_re.match(line)
            if match_score and current_pocket_num is not None:
                score_value = float(match_score.group(1))
                if score_value > best_score:
                    best_score = score_value
                    best_pocket_num = current_pocket_num

    print(f"[choose_best_pocket_number] Best pocket = {best_pocket_num} (score={best_score})")
    return best_pocket_num

###############################################################################
# 6. Docking with AutoDock Vina
###############################################################################

def perform_docking(vina_exec: str,
                    protein_pdbqt: Path,
                    ligand_pdbqt: Path,
                    output_dir: Path,
                    center: tuple,
                    size: tuple):
    """
    Run vina with the given center/size, saving output to <output_dir>.
    """
    out_pdbqt = output_dir / f"{protein_pdbqt.stem}_docked.pdbqt"
    log_file  = output_dir / f"{protein_pdbqt.stem}_docking.log"

    remove_ligand_tags_in_receptor(protein_pdbqt)

    cmd = [
        vina_exec,
        "--receptor", str(protein_pdbqt),
        "--ligand", str(ligand_pdbqt),
        "--center_x", str(center[0]),
        "--center_y", str(center[1]),
        "--center_z", str(center[2]),
        "--size_x", str(size[0]),
        "--size_y", str(size[1]),
        "--size_z", str(size[2]),
        "--out", str(out_pdbqt),
        "--log", str(log_file)
    ]
    print(f"[Docking] Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    print("[Docking] Completed successfully!")

###############################################################################
# 7. Visualization
###############################################################################

def visualize_docking(receptor_pdb: Path, docked_ligand_pdb: Path, output_image: Path):
    """
    Visualize docking results using PyMOL (via pymol2) with better handling for visibility and file locations.

    Args:
        receptor_pdb (Path): Path to the receptor PDB file.
        docked_ligand_pdb (Path): Path to the docked ligand PDB file.
        output_image (Path): Path to save the PNG visualization.
    """
    import pymol2

    pse_output = output_image.with_suffix(".pse")  # Save .pse file in the same directory as the image

    with pymol2.PyMOL() as pymol:
        # Load receptor and ligand
        pymol.cmd.load(str(receptor_pdb), "receptor")
        pymol.cmd.load(str(docked_ligand_pdb), "ligand")

        # Display receptor as surface
        pymol.cmd.show("surface", "receptor")
        pymol.cmd.color("cyan", "receptor")

        # Display ligand as sticks
        pymol.cmd.show("sticks", "ligand")
        pymol.cmd.color("red", "ligand")

        # Highlight binding site (optional)
        pymol.cmd.select("binding_site", "receptor within 5.0 of ligand")
        pymol.cmd.show("sticks", "binding_site")
        pymol.cmd.color("yellow", "binding_site")

        # Adjust view to ensure ligand visibility
        pymol.cmd.orient("all")  # Ensures the entire system is visible
        pymol.cmd.zoom("all", buffer=15)  # Adds extra space around the objects
        pymol.cmd.turn("x", 90)  # Rotate along x-axis for better view
        pymol.cmd.turn("y", 45)  # Slight tilt for a 3D perspective

        # Save image
        pymol.cmd.bg_color("white")  # Set background color to white
        pymol.cmd.ray(1200, 900)  # High-resolution rendering
        pymol.cmd.png(str(output_image))
        print(f"[Visualization] Image saved as {output_image}")

        # Save PyMOL session
        pymol.cmd.save(str(pse_output))
        print(f"[Visualization] PyMOL session saved as {pse_output}")


###############################################################################
# 8. Save results to DB
###############################################################################

def insert_docking_results(db_connection_string: str, docking_results: list):
    """
    Insert each docking result row into the 'docking_results' table.
    """
    print(docking_results)
    engine = create_engine(db_connection_string)
    try:
        with engine.begin() as conn:
            query = text("""
                INSERT INTO docking_results (
                    experiment_id, uniprot_id, ligand_cid, deg_id, binding_energy,
                    rmsd_lower_bound, rmsd_upper_bound, docking_pose_rank,
                    center_x, center_y, center_z,
                    size_x, size_y, size_z,
                    created_at, updated_at
                ) VALUES (
                    :experiment_id, :uniprot_id, :ligand_cid, :deg_id, :binding_energy,
                    :rmsd_lower_bound, :rmsd_upper_bound, :docking_pose_rank,
                    :center_x, :center_y, :center_z,
                    :size_x, :size_y, :size_z,
                    NOW(), NOW()
                )
            """)

            for result in docking_results:
                conn.execute(query, result)

    except Exception as e:
        traceback.print_exc()
        print(f"[DB] Insertion failed: {e}")
    finally:
        engine.dispose()

def insert_natural_compound(db_connection_string, compound_data):
    """
    Insert or update the natural_compounds table with compound data.
    """
    engine = create_engine(db_connection_string)
    try:
        with engine.begin() as conn:
            query = text("""
                INSERT INTO natural_compounds (
                    pubchem_id, compound_name, molecular_weight, chemical_structure
                ) VALUES (
                    :pubchem_id, :compound_name, :molecular_weight, :chemical_structure
                )
                ON CONFLICT (pubchem_id) DO UPDATE SET
                    compound_name = EXCLUDED.compound_name,
                    molecular_weight = EXCLUDED.molecular_weight,
                    chemical_structure = EXCLUDED.chemical_structure
            """)
            conn.execute(query, compound_data)
    except Exception as e:
        print(f"[DB] Failed to insert natural compound data: {e}")
    finally:
        engine.dispose()


def parse_docking_results(log_file_path: Path,
                          experiment_id: str,
                          deg_id: int,
                          ligand_cid: str,
                          center: tuple,
                          size: tuple):
    """
    Example parser for a docking log: if it has lines with "mode | affinity | ...", etc.
    (But currently not used in your final code.)
    """
    results = []
    if not log_file_path.is_file():
        return results

    lines = log_file_path.read_text().splitlines()
    for line in lines:
        if line.startswith("mode"):
            continue  # Skip header
        if not line.strip():
            break  # End of results
        parts = line.split()
        docking_pose_rank = int(parts[0])
        binding_energy = float(parts[1])
        rmsd_lower_bound = float(parts[2])
        rmsd_upper_bound = float(parts[3])
        result = {
            "experiment_id": experiment_id,
            "ligand_cid": ligand_cid,
            "deg_id": deg_id,
            "binding_energy": binding_energy,
            "rmsd_lower_bound": rmsd_lower_bound,
            "rmsd_upper_bound": rmsd_upper_bound,
            "docking_pose_rank": docking_pose_rank,
            "center_x": center[0],
            "center_y": center[1],
            "center_z": center[2],
            "size_x": size[0],
            "size_y": size[1],
            "size_z": size[2],
        }
        results.append(result)
    return results

def parse_docked_pdbqt(pdbqt_file_path: Path,
                       experiment_id: str,
                       deg_id: int,
                       ligand_cid: str,
                       center: tuple,
                       size: tuple,
                       uniprot_id: str):
    """
    Parse lines from the final .pdbqt that look like:
      REMARK VINA RESULT:   -5.3  0.000  0.000
    Returns a list of dicts for 'docking_results'.
    """
    results = []
    docking_pose_rank = 1

    if not pdbqt_file_path.is_file():
        return results

    lines = pdbqt_file_path.read_text().splitlines()
    for line in lines:
        if line.startswith("REMARK VINA RESULT:"):
            parts = line.split()
            # parts might look like:
            # ["REMARK", "VINA", "RESULT:", "-5.3", "0.0", "0.0"]
            binding_energy = float(parts[3])
            rmsd_lower_bound = float(parts[4])
            rmsd_upper_bound = float(parts[5])

            result = {
                "experiment_id": experiment_id,
                "deg_id": deg_id,
                "ligand_cid": ligand_cid,
                "binding_energy": binding_energy,
                "rmsd_lower_bound": rmsd_lower_bound,
                "rmsd_upper_bound": rmsd_upper_bound,
                "docking_pose_rank": docking_pose_rank,
                "center_x": center[0],
                "center_y": center[1],
                "center_z": center[2],
                "size_x": size[0],
                "size_y": size[1],
                "size_z": size[2],
                "uniprot_id": uniprot_id
            }
            results.append(result)
            docking_pose_rank += 1
    return results

###############################################################################
# 9. Main
###############################################################################

def parse_experiment_id(project_dir, exp_id: str):
    """
    If exp_id is a file, read its content as the real experiment ID.
    Otherwise just return it.
    """
    path_obj = Path(exp_id)
    if path_obj.is_file():
        return path_obj.read_text().strip()
    return exp_id

def fetch_structure_fallback(uniprot_id: str, output_pdb: Path):
    """
    1) Attempt alphaFold
    2) If fail, attempt RCSB
    3) Return True if we end with a valid PDB
    """
    alpha_ok = fetch_alphafold_pdb(uniprot_id, output_pdb)
    if alpha_ok:
        if quick_pdb_check(output_pdb):
            return True
        else:
            print(f"[Alpha] {uniprot_id} parse fail. Trying RCSB fallback.")
    rcsb_ok = fetch_rcsb_structure(uniprot_id, output_pdb)
    if rcsb_ok and quick_pdb_check(output_pdb):
        return True
    return False

def process_target(args):
    """
    Worker function to process a single target for molecular docking.
    """
    row, args, experiment_dir, ligand_pdbqt, df_params = args
    uniprot_id = row["uniprot_id"]
    gene_id = row["disease_gene_id"]
    if not uniprot_id or pd.isna(uniprot_id):
        print(f"[Warning] Missing UniProt ID for row={row}, skipping.")
        return None

    protein_pdb = experiment_dir / f"{uniprot_id}.pdb"
    protein_clean = experiment_dir / f"{uniprot_id}_clean.pdb"
    protein_pdbqt = experiment_dir / f"{uniprot_id}.pdbqt"

    docked_pdbqt = experiment_dir / f"{uniprot_id}_docked.pdbqt"
    docked_pdb = experiment_dir / f"{uniprot_id}_docked.pdb"
    docked_cleaned_pdb = experiment_dir / f"{uniprot_id}_docked_clean.pdb"
    log_file = experiment_dir / f"{uniprot_id}_docking.log"

    try:
        # Fetch protein structure
        if not fetch_structure_fallback(uniprot_id, protein_pdb):
            print(f"[Error] No structure found for {uniprot_id} on AlphaFold or RCSB. Skipping.")
            return None

        # Clean and prepare protein
        clean_pdb_biopython(protein_pdb, protein_clean)
        prepare_protein_with_openbabel(protein_clean, protein_pdbqt)

        # Determine docking box
        if df_params is not None:
            match = df_params[df_params["disease_gene_id"] == gene_id]
            if not match.empty:
                row_params = match.iloc[0]
                auto_detect = row_params["auto_detect"]
                if str(auto_detect).lower() == "true":
                    center, size = autodetect_pocket(protein_clean)
                else:
                    center_str = str(row_params["center"])
                    size_str = str(row_params["size"])
                    if center_str.strip() and size_str.strip():
                        center, size = get_docking_box(
                            {"Docking Site Center": center_str,
                             "Docking Site Size": size_str},
                            protein_clean
                        )
                    else:
                        center, size = autodetect_pocket(protein_clean)
            else:
                center, size = autodetect_pocket(protein_clean)
        else:
            center, size = autodetect_pocket(protein_clean)

        # Perform docking
        perform_docking(
            vina_exec=args.vina_executable,
            protein_pdbqt=protein_pdbqt,
            ligand_pdbqt=ligand_pdbqt,
            output_dir=experiment_dir,
            center=center,
            size=size
        )

        # Parse docking results
        results = parse_docked_pdbqt(
            pdbqt_file_path=docked_pdbqt,
            experiment_id=args.experiment_id,
            deg_id=row["deg_id"],
            ligand_cid=args.ligand_cid,
            center=center,
            size=size,
            uniprot_id=uniprot_id
        )
        insert_docking_results(args.db_connection_string, results)

        # Convert docked PDBQT to PDB for visualization
        subprocess.run(["obabel", "-ipdbqt", str(docked_pdbqt), "-opdb", "-O", str(docked_pdb)], check=True)

        extract_first_model(docked_pdb, docked_cleaned_pdb)
        pdb_content = docked_cleaned_pdb.read_text()
        pdb_content = pdb_content.replace(" UNL ", " LIG ")
        docked_cleaned_pdb.write_text(pdb_content)

        # Merge protein and docked ligand to create a complex
        docked_complex = experiment_dir / f"{uniprot_id}_docked_complex.pdb"
        merge_protein_ligand(protein_clean, docked_cleaned_pdb, docked_complex)


        # Visualize results with PyMOL if enabled
        if args.visualize:
            output_image_path = experiment_dir / f"{uniprot_id}_visualization.png"
            visualize_docking(protein_clean, docked_pdb, output_image=output_image_path)


        return f"[Success] Docking completed for {uniprot_id}"
    except Exception as e:
        print(f"[Error] Docking failed for {uniprot_id}: {e}")
        return None

def main():
    # Argument parsing and setup
    parser = argparse.ArgumentParser(description="ADHD pipeline fallback & docking.")
    parser.add_argument("--db_connection_string", required=True)
    parser.add_argument("--experiment_id", required=True)
    parser.add_argument("--ligand_cid", required=True)
    parser.add_argument("--vina_executable", default="vina")
    parser.add_argument("--output_dir", default="/home/scmbag/Desktop/ADHD_Compound_Pipeline/results/molecular_docking")
    parser.add_argument("--docking_params", required=False, help="Path to a CSV file with docking parameters.")
    parser.add_argument("--visualize", default=True)
    parser.add_argument("--project_dir", required=True)
    args = parser.parse_args()

    experiment_dir = Path(args.output_dir) / str(args.experiment_id)
    experiment_dir.mkdir(parents=True, exist_ok=True)

    # Fetch and save ligand properties
    try:
        compound_data = fetch_compound_properties(args.ligand_cid)
        insert_natural_compound(args.db_connection_string, compound_data)
        print(f"[NaturalCompounds] Saved data for CID={args.ligand_cid}")
    except Exception as e:
        print(f"[Error] Failed to fetch or save compound data: {e}")
        return

    # Fetch therapeutic targets
    targets = fetch_therapeutic_targets(args.db_connection_string, args.experiment_id)
    if targets.empty:
        raise RuntimeError("No therapeutic targets found for this experiment.")

    # Download and prepare ligand
    ligand_sdf = experiment_dir / "ligand.sdf"
    ligand_pdbqt = experiment_dir / "ligand.pdbqt"
    download_ligand(args.ligand_cid, ligand_sdf)
    prepare_ligand_with_openbabel(ligand_sdf, ligand_pdbqt)

    ligand_pdb = experiment_dir / "ligand.pdb"  # Define where to save ligand.pdb
    subprocess.run([
        "obabel",
        str(ligand_pdbqt),
        "-opdb",
        "-O",
        str(ligand_pdb)
    ], check=True)
    print(f"[Ligand] Converted ligand.pdbqt to {ligand_pdb}")

    df_params = pd.read_csv(args.docking_params) if args.docking_params else None

    # Parallel processing
    tasks = [(row, args, experiment_dir, ligand_pdbqt, df_params) for _, row in targets.iterrows()]
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = list(executor.map(process_target, tasks))

    print("[Pipeline] Molecular docking completed successfully.")


if __name__ == "__main__":
    main()

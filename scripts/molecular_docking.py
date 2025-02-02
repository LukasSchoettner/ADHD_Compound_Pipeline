#!/usr/bin/env python
import os
import shutil
import subprocess
import requests
import argparse
import pandas as pd
import traceback
import psycopg2
import re
import concurrent.futures
import logging
import glob
from rdkit import Chem
from Bio.PDB import PDBParser, PDBIO, Select
import numpy as np
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("docking_pipeline.log"),
        logging.StreamHandler()
    ]
)

class KeepAll(Select):
    def accept_residue(self, residue):
        return True

###############################################################################
# 1. Database Queries
###############################################################################

def fetch_therapeutic_targets(db_connection_string, experiment_id):
    """
    Fetch therapeutic targets for a given experiment ID using psycopg2.

    Args:
        db_connection_string (str): Database connection string.
        experiment_id (int): Experiment ID.

    Returns:
        pd.DataFrame: DataFrame containing therapeutic targets.
    """
    try:
        # Connect to the PostgreSQL database
        with psycopg2.connect(db_connection_string) as connection:
            with connection.cursor() as cursor:
                # Define the query
                query = """
                    SELECT
                        tt.disease_gene_id,
                        tt.uniprot_id,
                        tt.deg_id,
                        dg.gene_name
                    FROM therapeutic_target tt
                    INNER JOIN disease_genes dg
                        ON tt.disease_gene_id = dg.disease_gene_id
                    WHERE tt.experiment_id = %s
                """
                # Execute the query
                cursor.execute(query, (experiment_id,))

                # Fetch results and convert them into a DataFrame
                columns = [desc[0] for desc in cursor.description]
                data = cursor.fetchall()
                df = pd.DataFrame(data, columns=columns)

        logging.info(f"Fetched {len(df)} therapeutic targets for experiment {experiment_id}.")
        return df

    except Exception as e:
        logging.error(f"Error fetching therapeutic targets for experiment {experiment_id}: {e}")
        return pd.DataFrame()

###############################################################################
# 2. Ligand Handling with Open Babel
###############################################################################

def download_ligand(pubchem_cid, output_sdf: Path):
    """
    Download the ligand SDF from PubChem, saving to 'output_sdf' (Path).
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{pubchem_cid}/record/SDF/?record_type=3d"
    try:
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        output_sdf.write_bytes(resp.content)
        logging.info(f"[Ligand] Downloaded CID={pubchem_cid} to {output_sdf}")
    except requests.HTTPError as http_err:
        logging.error(f"[Ligand] HTTP error occurred while downloading CID={pubchem_cid}: {http_err}")
        raise
    except Exception as err:
        logging.error(f"[Ligand] Unexpected error occurred while downloading CID={pubchem_cid}: {err}")
        raise


def fetch_compound_properties(pubchem_cid):
    """
    Fetch detailed compound properties from PubChem for a given CID.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{pubchem_cid}/property/MolecularWeight,MolecularFormula,IUPACName,CanonicalSMILES/JSON"
    try:
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        data = resp.json()
        if "PropertyTable" in data and "Properties" in data["PropertyTable"]:
            properties = data["PropertyTable"]["Properties"][0]
            compound_info = {
                "pubchem_id": pubchem_cid,
                "compound_name": properties.get("IUPACName", "Unknown"),
                "molecular_weight": float(properties.get("MolecularWeight", 0.0)),
                "molecular_formula": properties.get("MolecularFormula", ""),
                "canonical_smiles": properties.get("CanonicalSMILES", ""),
                "chemical_structure": None  # Placeholder if needed
            }
            logging.info(f"[PubChem] Fetched properties for CID={pubchem_cid}")
            return compound_info
        else:
            raise ValueError(f"[PubChem] No properties found for CID={pubchem_cid}.")
    except requests.HTTPError as http_err:
        logging.error(f"[PubChem] HTTP error occurred for CID={pubchem_cid}: {http_err}")
        raise
    except Exception as e:
        logging.error(f"[PubChem] Error fetching properties for CID={pubchem_cid}: {e}")
        raise


def add_hydrogens_with_rdkit(input_pdb: Path, output_pdb: Path):
    """
    Add explicit hydrogens to a ligand PDB file using RDKit and save the updated PDB.

    Args:
        input_pdb (Path): Path to the input PDB file without sufficient hydrogens.
        output_pdb (Path): Path to save the hydrogen-added PDB file.

    Raises:
        ValueError: If RDKit fails to parse the PDB file or add hydrogens.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    try:
        # Read the PDB file into an RDKit molecule (without removing existing hydrogens)
        mol = Chem.MolFromPDBFile(str(input_pdb), removeHs=False)
        if mol is None:
            raise ValueError(f"RDKit failed to parse PDB file: {input_pdb}")

        # Remove existing hydrogens to avoid duplicates
        mol = Chem.RemoveHs(mol)

        # Add hydrogens explicitly
        mol_with_h = Chem.AddHs(mol)

        # Ensure the molecule remains a single fragment
        frags = Chem.GetMolFrags(mol_with_h, asMols=True)
        if len(frags) > 1:
            raise ValueError("Ligand consists of multiple disconnected fragments after adding hydrogens.")

        # Embed the molecule in 3D space
        result = AllChem.EmbedMolecule(mol_with_h, randomSeed=0xf00d)
        if result != 0:
            raise ValueError(f"RDKit failed to embed molecule: {input_pdb}")

        # Optimize the geometry using UFF force field
        optimize_result = AllChem.UFFOptimizeMolecule(mol_with_h)
        if optimize_result != 0:
            logging.warning(f"[RDKit] UFF optimization may not have fully converged for {input_pdb}")

        # Write the molecule back to a PDB file with hydrogens
        Chem.MolToPDBFile(mol_with_h, str(output_pdb))
        logging.info(f"[RDKit] Added hydrogens and saved updated PDB to {output_pdb}")

    except Exception as e:
        logging.error(f"[RDKit] Error adding hydrogens with RDKit: {e}")
        raise


def enforce_single_ligand_residue(pdb_file: Path):
    """
    Ensure that all atoms in the PDB file belong to a single 'LIG' residue and the same chain.

    Args:
        pdb_file (Path): Path to the PDB file to be modified.

    Raises:
        ValueError: If the PDB file contains multiple chains or residues after processing.
    """
    from Bio.PDB import PDBParser, PDBIO, Structure, Model, Chain, Residue, Atom

    class SingleLIGSelect(Select):
        def accept_atom(self, atom):
            return True

    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("ligand", str(pdb_file))

        # Create a new structure
        new_structure = Structure.Structure("ligand")
        new_model = Model.Model(0)
        new_structure.add(new_model)

        new_chain = Chain.Chain('A')
        new_model.add(new_chain)

        # Create a single residue named 'LIG'
        new_residue = Residue.Residue(('H', 1, ' '), 'LIG', '')
        new_chain.add(new_residue)

        # Add all atoms to the new residue
        for atom in structure.get_atoms():
            new_atom = Atom.Atom(
                name=atom.get_name(),
                coord=atom.get_coord(),
                bfactor=atom.get_bfactor(),
                occupancy=atom.get_occupancy(),
                altloc=atom.get_altloc(),
                fullname=atom.get_fullname(),
                serial_number=atom.get_serial_number(),
                element=atom.element
            )
            new_residue.add(new_atom)

        # Save the new structure to the original file
        io = PDBIO()
        io.set_structure(new_structure)
        io.save(str(pdb_file), select=SingleLIGSelect())
        logging.info(f"[BioPython] Enforced single 'LIG' residue and single chain in {pdb_file}")

    except Exception as e:
        logging.error(f"[BioPython] Error enforcing single 'LIG' residue: {e}")
        raise


def prepare_ligand_with_openbabel(input_sdf: Path, output_pdbqt: Path, ligand_resname="LIG", ligand_charge=0):
    """
    Generate a minimized PDB from the input SDF, ensure residue name is 'LIG',
    then convert to PDBQT.
    """
    logging.info("[Ligand] Generating PDB + Minimization with Open Babel...")
    temp_pdb = input_sdf.with_suffix(".pdb")
    try:
        subprocess.run([
            "obabel",
            str(input_sdf),
            "-O", str(temp_pdb),
            "--gen3d",
            "--minimize",
            "--strict"       # Enforce strict atom naming and geometry
        ], check=True)
        logging.info(f"[Ligand] Generated minimized PDB: {temp_pdb}")
    except subprocess.CalledProcessError as e:
        logging.error(f"[Ligand] Open Babel failed during PDB generation: {e}")
        raise

    # Modify residue names to 'LIG'
    temp_clean_pdb = input_sdf.parent / "temp_ligand_clean.pdb"
    try:
        with temp_pdb.open("r") as fin, temp_clean_pdb.open("w") as fout:
            for line in fin:
                if line.startswith(("ATOM", "HETATM")):
                    # Replace residue name (columns 17-20) with 'LIG'
                    line = line[:17] + "LIG" + line[20:]
                fout.write(line)
        logging.info(f"[Ligand] Cleaned PDB with residue name 'LIG': {temp_clean_pdb}")
    except Exception as e:
        logging.error(f"[Ligand] Error cleaning PDB residue names: {e}")
        raise

    # Convert to PDBQT without assigning partial charges
    # Let ACPyPE handle charge assignment and hydrogen addition
    logging.info("[Ligand] Converting cleaned PDB -> PDBQT without partial charges...")
    try:
        subprocess.run([
            "obabel",
            str(temp_clean_pdb),
            "-O", str(output_pdbqt),
            "-xh"  # Ensure hydrogens are explicit in PDBQT
        ], check=True)
        logging.info(f"[Ligand] Prepared PDBQT: {output_pdbqt}")
    except subprocess.CalledProcessError as e:
        logging.error(f"[Ligand] Open Babel failed during PDBQT conversion: {e}")
        raise

    # Note: Do not delete temp_clean_pdb here. It is needed for ACPyPE.
    # Cleanup will be handled after ACPyPE successfully generates ligand.itp
    logging.debug("[Ligand] Prepared ligand PDBQT without cleaning up temporary files.")


def remove_hydrogens(pdb_file: Path, output_pdb: Path):
    """
    Remove all hydrogen atoms from a PDB file and save to output_pdb.
    """
    class NoHydrogenSelect(Select):
        def accept_atom(self, atom):
            return atom.element != 'H'

    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("ligand", str(pdb_file))
        io = PDBIO()
        io.set_structure(structure)
        io.save(str(output_pdb), select=NoHydrogenSelect())
        logging.info(f"[Ligand] Removed hydrogens and saved to {output_pdb}")
    except Exception as e:
        logging.error(f"[Ligand] Failed to remove hydrogens from {pdb_file}: {e}")
        raise


def convert_pdbqt_to_gro(input_pdbqt: Path, output_gro: Path):
    """
    Convert a PDBQT file to GRO format using Open Babel.

    Args:
        input_pdbqt (Path): Path to the input PDBQT file.
        output_gro (Path): Path to the output GRO file.
    """
    try:
        subprocess.run([
            "obabel",
            "-ipdbqt", str(input_pdbqt),
            "-ogro",
            "-O", str(output_gro)
        ], check=True)
        logging.info(f"[Ligand] Converted {input_pdbqt} to GRO: {output_gro}")
    except subprocess.CalledProcessError as e:
        logging.error(f"[Ligand] Open Babel failed during PDBQT to GRO conversion: {e}")
        raise


def generate_ligand_itp_from_smiles(smiles: str, output_dir: Path, force_field='gaff2'):
    """
    Generate ligand.itp using ACPyPE from a SMILES string.

    Args:
        smiles (str): Canonical SMILES string of the ligand.
        output_dir (Path): Directory to save ligand.itp and related files.
        force_field (str): Force field to use (default: 'gaff2').
    """
    try:
        # Define basename
        basename = "LIG"
        acpype_dir = output_dir / f"{basename}.acpype"  # Path to ACPyPE output directory
        generated_itp = acpype_dir / f"{basename}_GMX.itp"  # Path to generated itp file

        # Define the ACPyPE command with correct force field flag
        cmd = [
            "acpype",
            "-i", smiles,         # Passing SMILES directly
            "-c", "bcc",          # Using 'bcc' for charge assignment
            "-b", basename,
            "-o", "gmx",
            "-a", force_field     # Correct flag for atom type
        ]
        logging.info(f"[Ligand] Running ACPyPE to generate ligand.itp from SMILES with force field {force_field}")

        # Execute ACPyPE with specified working directory
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=output_dir)
        logging.debug(f"ACPyPE STDOUT: {result.stdout}")
        logging.debug(f"ACPyPE STDERR: {result.stderr}")
        if result.returncode != 0:
            logging.error(f"[Ligand] ACPyPE failed:\n{result.stderr}")
            raise RuntimeError(f"ACPyPE failed: {result.stderr}")
        else:
            logging.info(f"[Ligand] ACPyPE completed successfully.")

        # Verify the ACPyPE output directory and itp file
        if not acpype_dir.exists():
            logging.error(f"[Ligand] ACPyPE output directory not found: {acpype_dir}")
            raise FileNotFoundError(f"ACPyPE output directory not found: {acpype_dir}")

        if not generated_itp.exists():
            logging.error(f"[Ligand] Expected itp file not found: {generated_itp}")
            raise FileNotFoundError(f"Expected itp file not found: {generated_itp}")

        # Define the target itp path
        target_itp = output_dir / "ligand.itp"

        # Move the generated itp to the output directory with the desired name
        shutil.move(str(generated_itp), str(target_itp))
        logging.info(f"[Ligand] ligand.itp generated at {target_itp}")

        # Verify that ligand.itp exists
        if target_itp.exists():
            logging.info(f"[Ligand] ligand.itp successfully moved to {target_itp}")
            # log the contents of the itp file for verification
            try:
                with target_itp.open("r") as itp_file:
                    itp_contents = itp_file.read()
                    logging.debug(f"[Ligand] Contents of ligand.itp:\n{itp_contents}")
            except Exception as e:
                logging.error(f"[Ligand] Failed to read ligand.itp: {e}")
        else:
            logging.error(f"[Ligand] ligand.itp was not found after moving to {target_itp}")
            raise FileNotFoundError(f"ligand.itp was not found after moving to {target_itp}")

    except Exception as e:
        logging.error(f"[Ligand] Failed to generate ligand.itp from SMILES: {e}")
        raise


def generate_ligand_itp(ligand_pdb: Path, output_dir: Path, force_field='gaff2'):
    """
    Generate ligand.itp using ACPyPE from a PDB file.

    Args:
        ligand_pdb (Path): Path to the ligand PDB file.
        output_dir (Path): Directory to save ligand.itp.
        force_field (str): Force field to use (default: 'gaff2').
    """
    try:
        # Define basename
        basename = ligand_pdb.stem  # 'temp_ligand_clean_noh'
        acpype_dir = output_dir / f"{basename}.acpype"  # Path to ACPyPE output directory
        generated_itp = acpype_dir / f"{basename}_GMX.itp"  # Path to generated itp file

        # Define the ACPyPE command with force field specification
        cmd = [
            "acpype",
            "-i", str(ligand_pdb),
            "-c", "bcc",        # Using 'bcc' for charge assignment
            "-b", basename,
            "-o", "gmx",
            "-at", force_field  # Specify the atom type / force field
        ]
        logging.info(f"[Ligand] Running ACPyPE to generate ligand.itp for {ligand_pdb} with force field {force_field}")

        # Execute ACPyPE with specified working directory
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=output_dir)
        logging.debug(f"ACPyPE STDOUT: {result.stdout}")
        logging.debug(f"ACPyPE STDERR: {result.stderr}")
        if result.returncode != 0:
            logging.error(f"[Ligand] ACPyPE failed:\n{result.stderr}")
            raise RuntimeError(f"ACPyPE failed: {result.stderr}")
        else:
            logging.info(f"[Ligand] ACPyPE completed successfully.")

        # Verify the ACPyPE output directory and itp file
        if not acpype_dir.exists():
            logging.error(f"[Ligand] ACPyPE output directory not found: {acpype_dir}")
            raise FileNotFoundError(f"ACPyPE output directory not found: {acpype_dir}")

        if not generated_itp.exists():
            logging.error(f"[Ligand] Expected itp file not found: {generated_itp}")
            raise FileNotFoundError(f"Expected itp file not found: {generated_itp}")

        # Define the target itp path
        target_itp = output_dir / "ligand.itp"

        # Move the generated itp to the output directory with the desired name
        shutil.move(str(generated_itp), str(target_itp))
        logging.info(f"[Ligand] ligand.itp generated at {target_itp}")

        # Verify that ligand.itp exists
        if target_itp.exists():
            logging.info(f"[Ligand] ligand.itp successfully moved to {target_itp}")
        else:
            logging.error(f"[Ligand] ligand.itp was not found after moving to {target_itp}")
            raise FileNotFoundError(f"ligand.itp was not found after moving to {target_itp}")

    except Exception as e:
        logging.error(f"[Ligand] Failed to generate ligand.itp: {e}")
        raise


def verify_mol2_atom_types(mol2_file: Path):
    """
    Verify that the MOL2 file does not contain undefined atom types.

    Args:
        mol2_file (Path): Path to the MOL2 file.

    Raises:
        ValueError: If undefined atom types are found.
    """
    undefined_types = set()
    try:
        with mol2_file.open("r") as f:
            for line in f:
                if line.startswith("@<TRIPOS>ATOM"):
                    continue
                if line.startswith("@<TRIPOS>BOND"):
                    break
                parts = line.split()
                if len(parts) >= 6:
                    atom_type = parts[5]
                    if atom_type == "DU":
                        undefined_types.add(atom_type)
        if undefined_types:
            raise ValueError(f"Undefined atom types found in {mol2_file}: {undefined_types}")
        logging.info(f"[MOL2 Verify] All atom types in {mol2_file} are defined.")
    except Exception as e:
        logging.error(f"[MOL2 Verify] Error verifying MOL2 file {mol2_file}: {e}")
        raise


###############################################################################
# 3. Protein Fetching & Fallback
###############################################################################

def fetch_alphafold_pdb(uniprot_id: str, output_pdb: Path):
    """
    Attempt to download the AlphaFold model for 'uniprot_id' and save to output_pdb.
    Returns True on success, False otherwise.
    """
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
    try:
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        output_pdb.write_text(resp.text)
        logging.info(f"[AlphaFold] Downloaded {uniprot_id} -> {output_pdb}")
        return True
    except requests.HTTPError as http_err:
        logging.warning(f"[AlphaFold] No structure or error for {uniprot_id} (HTTP {resp.status_code}).")
        return False
    except Exception as e:
        logging.error(f"[AlphaFold] Unexpected error for {uniprot_id}: {e}")
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
        search_resp = requests.post(url_search, json=search_payload, timeout=30)
        search_resp.raise_for_status()
        data = search_resp.json()
        results = data.get("result_set", [])
        if not results:
            logging.warning(f"[RCSB] No PDB found for {uniprot_id}.")
            return False
        pdb_id = results[0]["identifier"]
        logging.info(f"[RCSB] Found PDB ID={pdb_id} for {uniprot_id}. Downloading...")

        url_download = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        dresp = requests.get(url_download, timeout=30)
        dresp.raise_for_status()
        output_pdb.write_bytes(dresp.content)
        logging.info(f"[RCSB] Downloaded {pdb_id} -> {output_pdb}")
        return True
    except requests.HTTPError as http_err:
        logging.error(f"[RCSB] HTTP error occurred while fetching PDB for {uniprot_id}: {http_err}")
        return False
    except Exception as e:
        logging.error(f"[RCSB] Error fetching PDB for {uniprot_id}: {e}")
        return False

def quick_pdb_check(pdb_file: Path):
    """
    Quick parse check on 'pdb_file' to see if BioPython can read it.
    """
    parser = PDBParser(QUIET=True)
    try:
        _ = parser.get_structure("protein", str(pdb_file))
        logging.info(f"[ParseCheck] Successfully parsed {pdb_file}.")
        return True
    except Exception as e:
        logging.error(f"[ParseCheck] BioPython parse error for {pdb_file}: {e}")
        return False

def fetch_structure_fallback(uniprot_id: str, output_pdb: Path):
    """
    1) Attempt AlphaFold
    2) If fail, attempt RCSB
    3) Return True if we end with a valid PDB
    """
    alpha_ok = fetch_alphafold_pdb(uniprot_id, output_pdb)
    if alpha_ok:
        if quick_pdb_check(output_pdb):
            return True
        else:
            logging.warning(f"[AlphaFold] {uniprot_id} parse fail. Trying RCSB fallback.")
    rcsb_ok = fetch_rcsb_structure(uniprot_id, output_pdb)
    if rcsb_ok and quick_pdb_check(output_pdb):
        return True
    logging.error(f"[StructureFetch] Failed to fetch valid PDB for {uniprot_id} from both AlphaFold and RCSB.")
    return False

###############################################################################
# 4. Cleaning & Converting Protein
###############################################################################

def clean_pdb_biopython(input_pdb: Path, output_pdb: Path):
    """
    Load 'input_pdb' with BioPython, then save a 'clean' version to 'output_pdb'.
    """
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("protein", str(input_pdb))
        io = PDBIO()
        io.set_structure(structure)
        io.save(str(output_pdb), select=KeepAll())
        logging.info(f"[Protein] Cleaned PDB saved to {output_pdb}")
    except Exception as e:
        logging.error(f"[Protein] Error cleaning PDB {input_pdb}: {e}")
        raise

def prepare_protein_with_openbabel(input_pdb: Path, output_pdbqt: Path):
    """
    Use obabel to convert the protein from PDB to PDBQT.
    """
    logging.info("[Protein] Converting PDB -> PDBQT with Open Babel...")
    try:
        subprocess.run([
            "obabel",
            str(input_pdb),
            "-O", str(output_pdbqt),
            "--partialcharge", "gasteiger",
            "-xh"
        ], check=True)
        logging.info(f"[Protein] Prepared PDBQT: {output_pdbqt}")
    except subprocess.CalledProcessError as e:
        logging.error(f"[Protein] Open Babel failed during PDBQT conversion: {e}")
        raise

def remove_ligand_tags_in_receptor(receptor_pdbqt: Path):
    """
    Remove ligand-specific tags from the receptor PDBQT to prevent parsing errors.
    """
    cleaned_lines = []
    try:
        with receptor_pdbqt.open("r") as infile:
            for line in infile:
                if line.startswith(("ROOT", "ENDROOT", "BRANCH", "ENDBRANCH", "TORSDOF")):
                    continue
                cleaned_lines.append(line)
        with receptor_pdbqt.open("w") as outfile:
            outfile.writelines(cleaned_lines)
        logging.debug(f"[Protein] Removed ligand-specific tags from {receptor_pdbqt}")
    except Exception as e:
        logging.error(f"[Protein] Error cleaning receptor PDBQT {receptor_pdbqt}: {e}")
        raise

def extract_first_model(input_pdb: Path, output_pdb: Path):
    """
    Extract the first complete MODEL block from a PDB file, ignoring duplicate MODEL lines.
    """
    inside_model = False
    model_found = False
    try:
        with input_pdb.open("r") as infile, output_pdb.open("w") as outfile:
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
            # If there's no MODEL, copy the entire file
            shutil.copy2(input_pdb, output_pdb)
            logging.warning(f"[ModelExtract] No MODEL/ENDMDL found in {input_pdb}. Copied entire PDB.")
        else:
            logging.info(f"[ModelExtract] Extracted first model to {output_pdb}.")
    except Exception as e:
        logging.error(f"[ModelExtract] Error extracting first model from {input_pdb}: {e}")
        raise

def merge_protein_ligand(protein_pdb: Path, ligand_pdb: Path, output_pdb: Path):
    """
    Merge the receptor coordinates from protein_pdb with ligand coordinates from ligand_pdb.
    The result is saved to output_pdb.
    """
    try:
        with output_pdb.open("w") as of:
            # Merge protein structure
            with protein_pdb.open("r") as pf:
                for line in pf:
                    if line.startswith("END"):
                        continue  # Skip termination lines
                    of.write(line)

            # Merge ligand structure
            with ligand_pdb.open("r") as lf:
                for line in lf:
                    if line.startswith(("ATOM", "HETATM")):
                        # Ensure residue name is 'LIG'
                        line = line[:17] + "LIG" + line[20:]
                        of.write(line)

            # Terminate the PDB file
            of.write("END\n")

        logging.info(f"[Merge] Created merged protein-ligand PDB: {output_pdb}")
    except Exception as e:
        logging.error(f"[Merge] Error merging protein and ligand PDBs: {e}")
        raise

###############################################################################
# 5. Auto-Detect Pocket
###############################################################################

def autodetect_pocket(protein_pdb: Path):
    """
    Runs fpocket on 'protein_pdb', parses the 'base_name_info.txt' to find
    the best (highest Druggability Score) pocketN, then loads 'pocketN_*.pdb'
    to compute bounding box center & size.
    """
    base_name = protein_pdb.stem  # e.g., 'P04798_clean'
    out_dir = protein_pdb.parent / f"{base_name}_out"

    # 1) Run fpocket
    cmd = ["fpocket", "-f", str(protein_pdb)]
    try:
        subprocess.run(cmd, check=True)
        logging.info(f"[AutoDetect] Ran fpocket on {protein_pdb}")
    except subprocess.CalledProcessError as e:
        logging.error(f"[AutoDetect] fpocket failed: {e}")
        return (0.0, 0.0, 0.0), (20.0, 20.0, 20.0)

    # 2) Parse the info file to find the best pocket number
    info_file = out_dir / f"{base_name}_info.txt"  # e.g., P04798_clean_out/P04798_clean_info.txt
    best_pocket_num = choose_best_pocket_number(info_file)

    # 3) Find the actual pocket file, e.g., "pocket6_*.pdb"
    pockets_dir = out_dir / "pockets"
    pattern = str(pockets_dir / f"pocket{best_pocket_num}_*.pdb")
    candidates = glob.glob(pattern)
    if not candidates:
        logging.warning(f"[AutoDetect] No pocket{best_pocket_num}_*.pdb found. Using fallback center and size.")
        return (0.0, 0.0, 0.0), (20.0, 20.0, 20.0)

    pocket_file = Path(candidates[0])
    logging.info(f"[AutoDetect] Best pocket = {best_pocket_num}, file = {pocket_file}")

    # 4) Parse that file’s coordinates
    coords = []
    try:
        with pocket_file.open("r") as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x, y, z])
    except Exception as e:
        logging.error(f"[AutoDetect] Error parsing pocket file {pocket_file}: {e}")
        return (0.0, 0.0, 0.0), (20.0, 20.0, 20.0)

    if not coords:
        logging.warning(f"[AutoDetect] pocket{best_pocket_num} file empty. Using fallback center and size.")
        return (0.0, 0.0, 0.0), (20.0, 20.0, 20.0)

    arr = np.array(coords)
    min_xyz = arr.min(axis=0)
    max_xyz = arr.max(axis=0)
    center = tuple((min_xyz + max_xyz) / 2.0)
    size = tuple(max_xyz - min_xyz)
    logging.info(f"[AutoDetect] Detected pocket center={center}, size={size}")
    return center, size

def get_docking_box(row: dict, protein_pdb: Path,
                    fallback_center=(0.0,0.0,0.0),
                    fallback_size=(20.0,20.0,20.0)):
    """
    If row has manual center/size, parse them. Otherwise call autodetect_pocket().
    """
    center_str = row.get("Docking Site Center", "").strip()
    size_str   = row.get("Docking Site Size", "").strip()

    if center_str and size_str:
        try:
            cx, cy, cz = [float(x.strip()) for x in center_str.split(",")]
            sx, sy, sz = [float(x.strip()) for x in size_str.split(",")]
            logging.info(f"[DockingBox] Using manual center and size for docking.")
            return (cx, cy, cz), (sx, sy, sz)
        except ValueError:
            logging.warning("[DockingBox] Parse error for center/size. Using fallback.")
            return fallback_center, fallback_size
    else:
        center, size = autodetect_pocket(protein_pdb)
        return center, size

def choose_best_pocket_number(info_file: Path):
    """
    Parse an fpocket info file (like 'P04798_clean_info.txt') to find
    the pocket with the highest 'Druggability Score'.
    Returns an integer pocket_number (e.g., 6 if 'Pocket 6' is best).
    If none found, returns 1 by default.
    """
    if not info_file.is_file():
        logging.warning(f"[choose_best_pocket_number] Info file not found: {info_file}")
        return 1

    # Regex to match lines like "Pocket 6 :" or "Pocket 12 :"
    pocket_header_re = re.compile(r"^Pocket\s+(\d+)\s*:\s*$")
    # Regex to match lines like "Druggability Score :    0.988"
    drugg_score_re   = re.compile(r"^\s*Druggability Score\s*:\s*([\d\.]+)")

    best_pocket_num = 1
    best_score      = -999.0

    current_pocket_num = None

    try:
        with info_file.open("r") as f:
            for line in f:
                line = line.strip()

                # If it's a pocket header, e.g., "Pocket 6 :"
                match_header = pocket_header_re.match(line)
                if match_header:
                    current_pocket_num = int(match_header.group(1))
                    continue

                # If we see the Druggability Score
                match_score = drugg_score_re.match(line)
                if match_score and current_pocket_num is not None:
                    score_value = float(match_score.group(1))
                    if score_value > best_score:
                        best_score = score_value
                        best_pocket_num = current_pocket_num

        logging.info(f"[choose_best_pocket_number] Best pocket = {best_pocket_num} (score={best_score})")
    except Exception as e:
        logging.error(f"[choose_best_pocket_number] Error parsing {info_file}: {e}")
        return 1

    return best_pocket_num

###############################################################################
# 5. Docking with AutoDock Vina
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
    logging.info(f"[Docking] Running: {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True)
        logging.info("[Docking] Completed successfully!")
    except subprocess.CalledProcessError as e:
        logging.error(f"[Docking] AutoDock Vina failed: {e}")
        raise

###############################################################################
# 6. Visualization
###############################################################################


def calculate_box_size_from_pdb(pdb_file: Path, padding=1.0):
    """
    Calculate the simulation box size based on the PDB file's coordinates with padding.

    Args:
        pdb_file (Path): Path to the PDB file.
        padding (float): Padding to add around the molecule in nm.

    Returns:
        tuple: Box dimensions (x, y, z) in nm.
    """
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("protein", str(pdb_file))
    except Exception as e:
        logging.error(f"[BoxSize] Failed to parse PDB file {pdb_file}: {e}")
        raise

    coords = [atom.get_coord() for atom in structure.get_atoms()]
    if not coords:
        logging.error(f"[BoxSize] No atomic coordinates found in {pdb_file}.")
        raise ValueError(f"No atomic coordinates found in {pdb_file}.")

    coords = np.array(coords)
    min_xyz = coords.min(axis=0)
    max_xyz = coords.max(axis=0)
    box_dimensions = (max_xyz - min_xyz) / 10.0  # Convert Å to nm
    box_size = tuple(box_dimensions + 2 * padding)  # Adding padding on all sides

    # Enforce maximum box size to prevent impractical simulations
    MAX_BOX_SIZE = (20.0, 20.0, 20.0)  # Maximum dimensions in nm
    box_size = tuple(min(bs, max_bs) for bs, max_bs in zip(box_size, MAX_BOX_SIZE))

    logging.info(f"[BoxSize] Calculated box size from PDB with limits: {box_size} nm")
    return box_size


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

    try:
        with pymol2.PyMOL() as pymol:
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

            # Adjust view to focus on the ligand
            pymol.cmd.center("ligand")  # Center view on ligand
            pymol.cmd.zoom("ligand", buffer=3)  # Zoom in closely on the ligand
            pymol.cmd.turn("x", 30)  # Small rotation for better perspective
            pymol.cmd.turn("y", 30)  # Small tilt to adjust the angle


            # Save image
            pymol.cmd.bg_color("white")  # Set background color to white
            pymol.cmd.ray(1200, 900)  # High-resolution rendering
            pymol.cmd.png(str(output_image))
            logging.info(f"[Visualization] Image saved as {output_image}")

            # Save PyMOL session
            pymol.cmd.save(str(pse_output))
            logging.info(f"[Visualization] PyMOL session saved as {pse_output}")

    except Exception as e:
        logging.error(f"[Visualization] PyMOL failed: {e}")
        raise

###############################################################################
# 7. Save results to DB
###############################################################################

def insert_docking_results(db_connection_string: str, docking_results: list):
    """
    Insert each docking result row into the 'docking_results' table using psycopg2.
    """
    insert_sql = """
        INSERT INTO docking_results (
            experiment_id, uniprot_id, ligand_cid, deg_id, binding_energy,
            rmsd_lower_bound, rmsd_upper_bound, docking_pose_rank,
            center_x, center_y, center_z,
            size_x, size_y, size_z,
            created_at, updated_at
        ) VALUES (
            %(experiment_id)s, %(uniprot_id)s, %(ligand_cid)s, %(deg_id)s, %(binding_energy)s,
            %(rmsd_lower_bound)s, %(rmsd_upper_bound)s, %(docking_pose_rank)s,
            %(center_x)s, %(center_y)s, %(center_z)s,
            %(size_x)s, %(size_y)s, %(size_z)s,
            NOW(), NOW()
        )
    """

    try:
        with psycopg2.connect(db_connection_string) as conn:
            with conn.cursor() as cur:
                for result in docking_results:
                    cur.execute(insert_sql, result)
        logging.info(f"[DB] Successfully inserted {len(docking_results)} docking records.")
    except Exception as e:
        traceback.print_exc()
        logging.error(f"[DB] Insertion failed: {e}")

def insert_natural_compound(db_connection_string, compound_data):
    """
    Insert or update the natural_compounds table with compound data using psycopg2.
    """
    insert_sql = """
        INSERT INTO natural_compound (
            pubchem_id, compound_name, molecular_weight, chemical_structure
        ) VALUES (
            %(pubchem_id)s, %(compound_name)s, %(molecular_weight)s, %(chemical_structure)s
        )
        ON CONFLICT (pubchem_id) DO UPDATE SET
            compound_name = EXCLUDED.compound_name,
            molecular_weight = EXCLUDED.molecular_weight,
            chemical_structure = EXCLUDED.chemical_structure
    """

    try:
        with psycopg2.connect(db_connection_string) as conn:
            with conn.cursor() as cur:
                cur.execute(insert_sql, compound_data)
        logging.info(f"[DB] Successfully inserted/updated natural compound CID={compound_data['pubchem_id']}.")
    except Exception as e:
        traceback.print_exc()
        logging.error(f"[DB] Failed to insert natural compound data: {e}")

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
        logging.warning(f"[DockingResults] Docked PDBQT file not found: {pdbqt_file_path}")
        return results

    try:
        with pdbqt_file_path.open("r") as f:
            for line in f:
                if line.startswith("REMARK VINA RESULT:"):
                    parts = line.strip().split()
                    # Expected format: "REMARK VINA RESULT: -5.3 0.0 0.0"
                    if len(parts) >= 5:
                        binding_energy = float(parts[3])
                        rmsd_lower_bound = float(parts[4])
                        rmsd_upper_bound = float(parts[5]) if len(parts) > 5 else 0.0

                        result = {
                            "experiment_id": experiment_id,
                            "uniprot_id": uniprot_id,
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
                        docking_pose_rank += 1
    except Exception as e:
        logging.error(f"[DockingResults] Error parsing docking results from {pdbqt_file_path}: {e}")

    logging.info(f"[DockingResults] Parsed {len(results)} docking results from {pdbqt_file_path}")
    return results

###############################################################################
# 8. Conversion and Box Calculation
###############################################################################

def convert_pdb_to_gro(input_pdb: Path, output_gro: Path, box_size=(10.0, 10.0, 10.0), gmx_exec="gmx"):
    """
    Convert a PDB file to GRO format using GROMACS editconf.

    Args:
        input_pdb (Path): Path to the input PDB file.
        output_gro (Path): Path to the output GRO file.
        box_size (tuple): Dimensions of the simulation box in nm (x, y, z).
        gmx_exec (str): Path to the GROMACS executable.
    """
    cmd = [
        gmx_exec, "editconf",
        "-f", str(input_pdb),
        "-o", str(output_gro),
        "-box",
        str(box_size[0]),
        str(box_size[1]),
        str(box_size[2]),
        "-c"  # Center the molecule in the box
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logging.error(f"[Conversion] GROMACS editconf failed:\n{result.stderr}")
            raise RuntimeError(f"editconf failed: {result.stderr}")
        else:
            logging.info(f"[Conversion] Successfully converted {input_pdb} to {output_gro} with box size {box_size}")
    except subprocess.CalledProcessError as e:
        logging.error(f"[Conversion] GROMACS editconf failed: {e}")
        logging.error(traceback.format_exc())
        raise

def calculate_box_size(gro_file: Path, padding=1.0):
    """
    Calculate the simulation box size based on the GRO file dimensions with added padding.

    Args:
        gro_file (Path): Path to the GRO file.
        padding (float): Padding to add around the molecule in nm.

    Returns:
        tuple: Box dimensions (x, y, z) in nm.
    """
    with gro_file.open("r") as f:
        lines = f.readlines()
    box_line = lines[-1]
    box_dims = [float(dim) for dim in box_line.strip().split()]
    box_size = tuple(dim + 2 * padding for dim in box_dims[:3])
    logging.info(f"[BoxSize] Calculated box size with padding: {box_size}")
    return box_size


def count_hydrogens(smiles: str):
    """
    Count the number of hydrogen atoms in a molecule given its SMILES string.

    Args:
        smiles (str): Canonical SMILES string of the molecule.

    Returns:
        int: Number of hydrogen atoms.

    Raises:
        ValueError: If RDKit fails to parse the SMILES string.
    """
    if not smiles:
        raise ValueError("Empty SMILES string provided.")

    # Parse the SMILES string to create an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit failed to parse SMILES string: {smiles}")

    # Add explicit hydrogens
    mol_with_h = Chem.AddHs(mol)

    # Count the number of hydrogen atoms
    hydrogen_count = sum(1 for atom in mol_with_h.GetAtoms() if atom.GetSymbol() == 'H')

    logging.info(f"[RDKit] Calculated hydrogen count: {hydrogen_count} for SMILES={smiles}")
    return hydrogen_count



###############################################################################
# 9. Main
###############################################################################

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


def validate_ligand_pdb(pdb_file: Path, compound_info: dict, sdf_file: Path = None):
    """
    Validate the ligand PDB file to ensure it contains only one residue named 'LIG'
    and has the correct number of hydrogens. If validation fails due to insufficient
    hydrogens, add hydrogens using RDKit and re-validate.

    Args:
        pdb_file (Path): Path to the ligand PDB file.
        compound_info (dict): Dictionary containing compound properties, including SMILES.
        sdf_file (Path, optional): Path to the ligand SDF file for additional validation.

    Raises:
        ValueError: If validation fails even after adding hydrogens.
    """

    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("ligand", str(pdb_file))
        residues = list(structure.get_residues())
        logging.debug(f"[Validation] Number of residues: {len(residues)}")
        for residue in residues:
            logging.debug(f"[Validation] Residue: {residue.get_resname()} in chain {residue.get_parent().id}")
        if len(residues) != 1:
            raise ValueError("Ligand PDB must contain exactly one residue.")
        if residues[0].get_resname() != "LIG":
            raise ValueError("Ligand residue name must be 'LIG'.")

        # Check for hydrogens
        atoms = list(structure.get_atoms())
        hydrogen_atoms = [atom for atom in atoms if atom.element == 'H' or atom.get_name().startswith('H')]
        actual_hydrogen_count = len(hydrogen_atoms)
        logging.debug(f"[Validation] Number of atoms: {len(atoms)}")
        logging.debug(f"[Validation] Number of hydrogens: {actual_hydrogen_count}")

        # Calculate expected hydrogens using RDKit
        expected_hydrogen_count = count_hydrogens(compound_info.get("canonical_smiles", ""))
        logging.debug(f"[Validation] Expected number of hydrogens: {expected_hydrogen_count}")

        if actual_hydrogen_count < expected_hydrogen_count:
            logging.warning(f"Ligand PDB has insufficient hydrogens: expected {expected_hydrogen_count}, found {actual_hydrogen_count}. Adding missing hydrogens.")
            # Add hydrogens with RDKit
            hydrogen_added_pdb = pdb_file.parent / f"{pdb_file.stem}_with_h.pdb"
            add_hydrogens_with_rdkit(pdb_file, hydrogen_added_pdb)

            # Enforce single residue and chain
            enforce_single_ligand_residue(hydrogen_added_pdb)

            # Re-validate the updated PDB
            structure = parser.get_structure("ligand", str(hydrogen_added_pdb))
            residues = list(structure.get_residues())
            if len(residues) != 1 or residues[0].get_resname() != "LIG":
                raise ValueError("Ligand PDB must contain exactly one residue named 'LIG' after adding hydrogens.")
            atoms = list(structure.get_atoms())
            hydrogen_atoms = [atom for atom in atoms if atom.element == 'H' or atom.get_name().startswith('H')]
            actual_hydrogen_count = len(hydrogen_atoms)
            logging.debug(f"[Validation] After adding hydrogens - Number of hydrogens: {actual_hydrogen_count}")
            if actual_hydrogen_count < expected_hydrogen_count:
                raise ValueError(f"Ligand PDB still has insufficient hydrogens after addition: expected {expected_hydrogen_count}, found {actual_hydrogen_count}.")
            else:
                # Replace the original PDB with the hydrogen-added PDB
                pdb_file.unlink()  # Remove the original PDB
                hydrogen_added_pdb.rename(pdb_file)  # Rename the hydrogen-added PDB to original name
                logging.info(f"[Validation] {pdb_file} now has the correct number of hydrogens and a single 'LIG' residue.")

        else:
            logging.info(f"[Validation] {pdb_file} passed validation with {actual_hydrogen_count} hydrogen atoms.")
    except Exception as e:
        logging.error(f"[Validation] {pdb_file} failed validation: {e}")
        raise


def process_target(params):
    """
    Worker function to process a single target for molecular docking.
    """
    row = params["row"]
    args = params["args"]
    experiment_dir = params["experiment_dir"]
    ligand_itp = params["ligand_itp"]  # Use ligand.itp generated in main
    df_params = params["df_params"]
    gpu_id = params.get("gpu_id", None)

    uniprot_id = row["uniprot_id"]
    gene_id = row["disease_gene_id"]
    deg_id = row["deg_id"]

    if not uniprot_id or pd.isna(uniprot_id):
        logging.warning(f"[ProcessTarget] Missing UniProt ID for row={row}, skipping.")
        return None

    protein_pdb = experiment_dir / f"{uniprot_id}.pdb"
    protein_clean = experiment_dir / f"{uniprot_id}_clean.pdb"
    protein_pdbqt = experiment_dir / f"{uniprot_id}.pdbqt"

    docked_pdbqt = experiment_dir / f"{uniprot_id}_docked.pdbqt"
    docked_pdb = experiment_dir / f"{uniprot_id}_docked.pdb"
    docked_cleaned_pdb = experiment_dir / f"{uniprot_id}_docked_clean.pdb"
    docked_complex_pdb = experiment_dir / f"{uniprot_id}_docked_complex.pdb"  # Added
    docked_complex_gro = experiment_dir / f"{uniprot_id}_docked_complex.gro"  # Added
    log_file = experiment_dir / f"{uniprot_id}_docking.log"

    try:
        # Assign GPU if specified
        if gpu_id is not None:
            os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
            logging.info(f"[ProcessTarget] Assigned GPU {gpu_id} to docking of UniProt ID {uniprot_id}.")
        else:
            logging.info(f"[ProcessTarget] No GPU assigned. Docking of UniProt ID {uniprot_id} will use CPU.")

        # Fetch protein structure
        if not fetch_structure_fallback(uniprot_id, protein_pdb):
            logging.error(f"[ProcessTarget] No structure found for {uniprot_id} on AlphaFold or RCSB. Skipping.")
            return None

        # Clean and prepare protein
        clean_pdb_biopython(protein_pdb, protein_clean)
        prepare_protein_with_openbabel(protein_clean, protein_pdbqt)

        # Determine docking box
        if df_params is not None:
            match = df_params[df_params["disease_gene_id"] == gene_id]
            if not match.empty:
                row_params = match.iloc[0]
                auto_detect = str(row_params.get("auto_detect", "True")).lower()
                if auto_detect == "true":
                    center, size = autodetect_pocket(protein_clean)
                else:
                    center, size = get_docking_box(
                        {"Docking Site Center": row_params.get("Docking Site Center", ""),
                         "Docking Site Size": row_params.get("Docking Site Size", "")},
                        protein_clean
                    )
            else:
                logging.warning(f"[ProcessTarget] No docking parameters found for gene_id={gene_id}. Using auto-detection.")
                center, size = autodetect_pocket(protein_clean)
        else:
            center, size = autodetect_pocket(protein_clean)

        # Perform docking
        perform_docking(
            vina_exec=args.vina_executable,
            protein_pdbqt=protein_pdbqt,
            ligand_pdbqt=experiment_dir / "ligand.pdbqt",  # Use ligand.pdbqt from main
            output_dir=experiment_dir,
            center=center,
            size=size
        )

        # Parse docking results
        docking_results = parse_docked_pdbqt(
            pdbqt_file_path=docked_pdbqt,
            experiment_id=args.experiment_id,
            deg_id=deg_id,
            ligand_cid=args.ligand_cid,
            center=center,
            size=size,
            uniprot_id=uniprot_id
        )
        if docking_results:
            insert_docking_results(args.db_connection_string, docking_results)
        else:
            logging.warning(f"[ProcessTarget] No docking results parsed for UniProt ID {uniprot_id}.")

        # Convert docked PDBQT to PDB for visualization
        try:
            subprocess.run(["obabel", "-ipdbqt", str(docked_pdbqt), "-opdb", "-O", str(docked_pdb)], check=True)
            logging.info(f"[ProcessTarget] Converted {docked_pdbqt} to PDB: {docked_pdb}")
        except subprocess.CalledProcessError as e:
            logging.error(f"[ProcessTarget] Open Babel failed during PDB conversion: {e}")
            raise

        # Extract first model
        extract_first_model(docked_pdb, docked_cleaned_pdb)

        # Rename ligand residue from "UNL" to "LIG" if necessary
        pdb_content = docked_cleaned_pdb.read_text()
        pdb_content = pdb_content.replace(" UNL ", " LIG ")
        docked_cleaned_pdb.write_text(pdb_content)

        # Merge protein and docked ligand to create a complex
        merge_protein_ligand(protein_clean, docked_cleaned_pdb, docked_complex_pdb)

        # Calculate box size based on the merged PDB
        box_size = calculate_box_size_from_pdb(docked_complex_pdb, padding=1.0)

        # Convert merged complex PDB to GRO
        try:
            convert_pdb_to_gro(
                input_pdb=docked_complex_pdb,
                output_gro=docked_complex_gro,
                box_size=box_size
            )
            logging.info(f"[ProcessTarget] Generated GRO file: {docked_complex_gro}")
        except Exception as e:
            logging.error(f"[ProcessTarget] Failed to convert merged PDB to GRO: {e}")
            raise

        output_image_path = experiment_dir / f"{uniprot_id}_visualization.png"
        visualize_docking(protein_clean, docked_pdb, output_image=output_image_path)

    except Exception as e:
        logging.error(f"[ProcessTarget] An error occurred while processing UniProt ID {uniprot_id}: {e}")
        return None

    return f"[Success] Docking completed for UniProt ID={uniprot_id}"


def main():
    # Argument parsing and setup
    parser = argparse.ArgumentParser(description="ADHD pipeline fallback & docking.")
    parser.add_argument("--db_connection_string", required=True, help="Database connection string.")
    parser.add_argument("--experiment_id", required=True, help="Experiment ID or path to a file containing it.")
    parser.add_argument("--ligand_cid", required=True, help="PubChem CID of the ligand.")
    parser.add_argument("--vina_executable", default="vina", help="Path to the AutoDock Vina executable.")
    parser.add_argument("--output_dir", default="/home/scmbag/Desktop/ADHD_Compound_Pipeline/results/molecular_docking",
                        help="Output directory for molecular docking results.")
    parser.add_argument("--docking_params", required=False, help="Path to a CSV file with docking parameters.")
    parser.add_argument("--visualize", action="store_true", help="Enable visualization of docking results.")
    parser.add_argument("--project_dir", required=True, help="Project directory path.")
    parser.add_argument("--ligand_resname", default="LIG", help="Residue name for the ligand.")
    parser.add_argument("--num_workers", type=int, default=10, help="Number of parallel workers.")
    parser.add_argument("--gpu_ids", nargs='*', type=int, default=[0], help="List of GPU IDs to assign to workers.")
    args = parser.parse_args()

    experiment_dir = Path(args.output_dir) / str(parse_experiment_id(args.experiment_id))
    experiment_dir.mkdir(parents=True, exist_ok=True)

    # Fetch and save ligand properties
    try:
        compound_data = fetch_compound_properties(args.ligand_cid)
        insert_natural_compound(args.db_connection_string, compound_data)
        logging.info(f"[NaturalCompounds] Saved data for CID={args.ligand_cid}")
    except Exception as e:
        logging.error(f"[Main] Failed to fetch or save compound data: {e}")
        return

    # Fetch therapeutic targets
    targets = fetch_therapeutic_targets(args.db_connection_string, args.experiment_id)
    if targets.empty:
        logging.error("No therapeutic targets found for this experiment. Exiting.")
        return

    # Download and prepare ligand
    ligand_sdf = experiment_dir / "ligand.sdf"
    ligand_pdbqt = experiment_dir / "ligand.pdbqt"
    ligand_gro = experiment_dir / "ligand.gro"  # Define the path for ligand.gro
    ligand_itp = experiment_dir / "ligand.itp"  # Define the path for ligand.itp

    try:
        download_ligand(args.ligand_cid, ligand_sdf)
        prepare_ligand_with_openbabel(ligand_sdf, ligand_pdbqt, ligand_resname=args.ligand_resname)
        convert_pdbqt_to_gro(ligand_pdbqt, ligand_gro)  # Convert PDBQT to GRO

        # Extract SMILES from SDF
        smiles = compound_data["canonical_smiles"]

        # Validate the ligand SMILES could be implemented here in future iterations

        # Generate ligand.itp from SMILES using ACPyPE
        generate_ligand_itp_from_smiles(smiles=smiles, output_dir=experiment_dir, force_field='gaff2')
        logging.info(f"[Ligand] Generated ligand.itp successfully.")
    except Exception as e:
        logging.error(f"[Main] Failed to download or prepare ligand: {e}")
        return

    # Load docking parameters if provided
    df_params = pd.read_csv(args.docking_params) if args.docking_params else None
    if df_params is not None:
        logging.info(f"[Main] Loaded docking parameters from {args.docking_params}")

    # Prepare list of GPU IDs
    gpu_ids = args.gpu_ids
    num_workers = args.num_workers
    if len(gpu_ids) < num_workers:
        # Repeat GPU IDs if not enough
        gpu_ids = (gpu_ids * (num_workers // len(gpu_ids) + 1))[:num_workers]
    else:
        gpu_ids = gpu_ids[:num_workers]

    # Prepare tasks for parallel processing
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = []
        for i, row in targets.iterrows():
            gpu_id = gpu_ids[i % len(gpu_ids)]
            task_args = {
                "row": row,
                "args": args,
                "experiment_dir": experiment_dir,
                "ligand_itp": ligand_itp,  # Use ligand.itp generated in main
                "df_params": df_params,
                "gpu_id": gpu_id
            }
            futures.append(executor.submit(process_target, task_args))

        # Collect results
        for future in concurrent.futures.as_completed(futures):
            try:
                result = future.result()
                if result:
                    logging.info(result)
                else:
                    logging.warning("A docking task failed or was skipped.")
            except Exception as e:
                logging.error(f"[Main] A docking task encountered an exception: {e}")

    logging.info("[Pipeline] Molecular docking completed successfully.")


if __name__ == "__main__":
    main()
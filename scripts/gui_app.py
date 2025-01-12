import streamlit as st
import yaml
import subprocess
import os
import requests
import pandas as pd


def fetch_ligand_cid(compound_name):
    """
    Fetch PubChem CID for a given compound name.

    Args:
        compound_name (str): Name of the compound.

    Returns:
        str: PubChem CID or None if not found.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{compound_name}/cids/JSON"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        cids = data.get("IdentifierList", {}).get("CID", [])
        if cids:
            return str(cids[0])  # Return the first CID as a string
    return None


# Load configuration
with open("config.yaml", "r") as file:
    config = yaml.safe_load(file)

# GUI defaults from config
geo_id = st.text_input("GEO ID:", config["pipeline"]["default_geo_id"])
samples = st.text_input("Samples:", config["pipeline"]["default_samples"])
compound_name = st.text_input("Natural Compound (Name):", config["pipeline"]["default_compound"])
ligand_cid = st.text_input("Ligand PubChem CID (Optional):", config["pipeline"]["default_cid"])
description = st.text_input("Description:", config["pipeline"]["default_description"])

# Checkbox for optional analyses
include_ppi = st.checkbox("Include PPI Network Analysis", value=True)
include_pathway = st.checkbox("Perform Pathway Enrichment Analysis", value=False)

# Button to run the pipeline
run_button = st.button("Run Pipeline")

if run_button:
    # Resolve PubChem CID if not provided
    if not ligand_cid:
        ligand_cid = fetch_ligand_cid(compound_name)
        if ligand_cid:
            st.success(f"PubChem CID for '{compound_name}' resolved: {ligand_cid}")
        else:
            st.error(f"Unable to resolve PubChem CID for compound: {compound_name}. Please provide a valid name or CID.")
            st.stop()

    # Confirm the pipeline execution details
    st.write(f"Running pipeline with the following parameters:")
    st.write(f"- GEO ID: {geo_id}")
    st.write(f"- Samples: {samples}")
    st.write(f"- Compound Name: {compound_name}")
    st.write(f"- PubChem CID: {ligand_cid}")
    st.write(f"- Experiment Description: {description}")

    # Build the command
    cmd = [
        "nextflow", "run", "main.nf",
        f"--geo_id={geo_id}",
        f"--samples={samples}",
        f"--compound_name={compound_name}",
        f"--ligand_cid={ligand_cid}",
        f"--description={description}"
    ]

    # Add optional parameters
    if include_ppi:
        cmd.append("--ppi_analysis")
    if include_pathway:
        cmd.append("--pathway_analysis")

    # Show the command being executed
    st.write("Executing the following command:")
    st.code(" ".join(cmd))

    # Run the pipeline and capture output
    try:
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        for line in process.stdout:
            st.text(line)  # Stream output to the GUI
        process.wait()
        if process.returncode == 0:
            st.success("Pre-Docking Analysis completed successfully!")
        else:
            st.error(f"Workflow 1 failed with return code {process.returncode}")
            st.error(process.stderr.read())
    except Exception as e:
        st.error(f"An error occurred: {e}")

# Step 2: Display matched genes and allow docking site input after DEG analysis
deg_results_path = "results/deg_results.csv"

if os.path.exists(deg_results_path):
    st.success("DEG Analysis Completed!")
    st.subheader("Matched Genes for Docking")

    # Load matched genes
    matched_genes = pd.read_csv(deg_results_path)

    # Add columns for user input
    matched_genes["Docking Site Center"] = ""
    matched_genes["Docking Site Size"] = ""

    # Allow user to input docking parameters
    for idx, row in matched_genes.iterrows():
        st.write(f"Gene: {row['Gene']}, Protein: {row['Protein']}")
        use_auto = st.checkbox(f"Auto-detect docking site for {row['Protein']}", value=True, key=f"auto_{idx}")
        if not use_auto:
            center = st.text_input(f"Center (x, y, z) for {row['Protein']}", value="10, 10, 10", key=f"center_{idx}")
            size = st.text_input(f"Size (x, y, z) for {row['Protein']}", value="20, 20, 20", key=f"size_{idx}")
            matched_genes.at[idx, "Docking Site Center"] = center
            matched_genes.at[idx, "Docking Site Size"] = size

    # Save user inputs
    if st.button("Save Docking Parameters"):
        matched_genes.to_csv("results/docking_parameters.csv", index=False)
        st.success("Docking parameters saved. Resume the pipeline.")

    if st.button("Resume Pipeline"):
        st.write("Resuming pipeline for molecular docking...")
    cmd = [
        "nextflow", "run", ".",
        f"--ligand_cid={ligand_cid}",
        f"--db_connection_string={config['database']['connection_string']}",
        f"--docking_params=results/docking_parameters.csv"
    ]
    st.write("Executing the following command:")
    st.code(" ".join(cmd))

    # Run the Nextflow command
    try:
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        for line in process.stdout:
            st.text(line)  # Stream output to the GUI
        process.wait()
        if process.returncode == 0:
            st.success("Molecular docking completed successfully!")
        else:
            st.error(f"Pipeline run failed with return code {process.returncode}")
            st.error(process.stderr.read())
    except Exception as e:
        st.error(f"An error occurred: {e}")

import streamlit as st
import yaml
import subprocess
import os
import pandas as pd
import psycopg2

def fetch_therapeutic_targets(db_connection_string, experiment_id):
    """
    Query the database to fetch therapeutic targets for a given experiment.

    Args:
        db_connection_string (str): Database connection string.
        experiment_id (int): Experiment ID.

    Returns:
        pd.DataFrame: DataFrame containing therapeutic targets.
    """
    try:
        conn = psycopg2.connect(db_connection_string)
        query = """
            SELECT tt.disease_gene_id, dg.gene_name
            FROM therapeutic_targets tt
            INNER JOIN disease_genes dg ON tt.disease_gene_id = dg.disease_gene_id
            WHERE tt.experiment_id = %s;
        """
        targets = pd.read_sql_query(query, conn, params=(experiment_id,))
        conn.close()
        return targets
    except Exception as e:
        st.error(f"Error fetching therapeutic targets: {e}")
        return pd.DataFrame()

def run_nextflow_workflow(workflow, params):
    """
    Run a Nextflow workflow with the given parameters.

    Args:
        workflow (str): Workflow file name.
        params (list): List of workflow parameters.

    Returns:
        bool: True if the workflow completed successfully, False otherwise.
    """
    cmd = ["nextflow", "run", workflow] + params

    st.write("Executing the following command:")
    st.code(" ".join(cmd))

    try:
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        for line in process.stdout:
            st.text(line)
        process.wait()
        if process.returncode == 0:
            st.success(f"{workflow} completed successfully!")
            return True
        else:
            st.error(f"{workflow} failed with return code {process.returncode}")
            st.error(process.stderr.read())
            return False
    except Exception as e:
        st.error(f"An error occurred: {e}")
        return False

def display_docking_input(therapeutic_targets):
    """
    Display input fields for docking site parameters for therapeutic targets.

    Args:
        therapeutic_targets (pd.DataFrame): DataFrame of therapeutic targets.
    """
    therapeutic_targets["Docking Site Center"] = ""
    therapeutic_targets["Docking Site Size"] = ""

    for idx, row in therapeutic_targets.iterrows():
        st.write(f"Gene: {row['gene_name']}, Disease Gene ID: {row['disease_gene_id']}")
        use_auto = st.checkbox(f"Auto-detect docking site for {row['disease_gene_id']}", value=True, key=f"auto_{idx}")
        if not use_auto:
            center = st.text_input(f"Center (x, y, z) for {row['disease_gene_id']}", value="10, 10, 10", key=f"center_{idx}")
            size = st.text_input(f"Size (x, y, z) for {row['disease_gene_id']}", value="20, 20, 20", key=f"size_{idx}")
            therapeutic_targets.at[idx, "Docking Site Center"] = center
            therapeutic_targets.at[idx, "Docking Site Size"] = size

    if st.button("Save Docking Parameters"):
        docking_params_path = "results/docking_parameters.csv"
        therapeutic_targets.to_csv(docking_params_path, index=False)
        st.success("Docking parameters saved. Ready for molecular docking.")

def main():
    with open("config.yaml", "r") as file:
        config = yaml.safe_load(file)

    # GUI defaults
    geo_id = st.text_input("GEO ID:", config["pipeline"]["default_geo_id"])
    samples = st.text_input("Samples:", config["pipeline"]["default_samples"])
    compound_name = st.text_input("Natural Compound (Name):", config["pipeline"]["default_compound"])
    ligand_cid = st.text_input("Ligand PubChem CID:", config["pipeline"]["default_cid"])
    description = st.text_input("Description:", config["pipeline"]["default_description"])

    if st.button("Run Pre-Docking Analysis"):
        params = [
            f"--geo_id={geo_id}",
            f"--samples={samples}",
            f"--compound_name={compound_name}",
            f"--description={description}",
            f"--db_connection_string={config['database']['connection_string']}"
        ]
        if run_nextflow_workflow("workflow1.nf", params):
            st.session_state["pre_docking_completed"] = True

    if st.session_state.get("pre_docking_completed", False):
        experiment_id_path = "results/experiment_id.txt"
        if os.path.exists(experiment_id_path):
            with open(experiment_id_path, "r") as f:
                experiment_id = int(f.read().strip())
            therapeutic_targets = fetch_therapeutic_targets(config["database"]["connection_string"], experiment_id)
            display_docking_input(therapeutic_targets)

    if os.path.exists("results/docking_parameters.csv") and st.button("Run Molecular Docking"):
        params = [
            f"--ligand_cid={ligand_cid}",
            f"--db_connection_string={config['database']['connection_string']}",
            f"--docking_params=results/docking_parameters.csv"
        ]
        run_nextflow_workflow("workflow2.nf", params)

# Entry point
if __name__ == "__main__":
    main()

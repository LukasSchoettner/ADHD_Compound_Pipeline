import streamlit as st
import yaml
import subprocess
import os
import pandas as pd
import psycopg2
from sqlalchemy import create_engine

def fetch_therapeutic_targets(db_connection_string, experiment_id):
    """
    Query the database to fetch therapeutic targets for a given experiment.
    """
    try:
        engine = create_engine(db_connection_string)
        query = """
            SELECT tt.disease_gene_id, dg.gene_name
            FROM therapeutic_targets tt
            INNER JOIN disease_genes dg ON tt.disease_gene_id = dg.disease_gene_id
            WHERE tt.experiment_id = %(exp_id)s;
        """
        params_dict = {"exp_id": experiment_id}
        targets = pd.read_sql_query(query, con=engine, params=params_dict)
        engine.dispose()
        return targets
    except Exception as e:
        st.error(f"Error fetching therapeutic targets: {e}")
        return pd.DataFrame()

def run_nextflow_workflow(workflow_file, params):
    """
    Run a Nextflow workflow with the given parameters.
    """
    cmd = ["nextflow", "run", workflow_file] + params
    st.write("Executing Nextflow command:")
    st.code(" ".join(cmd))
    try:
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        for line in process.stdout:
            st.text(line)
        process.wait()
        if process.returncode == 0:
            st.success(f"{workflow_file} completed successfully!")
            return True
        else:
            st.error(f"{workflow_file} failed with return code {process.returncode}")
            st.error(process.stderr.read())
            return False
    except Exception as e:
        st.error(f"An error occurred: {e}")
        return False

def display_docking_input(therapeutic_targets):
    """
    Display input fields for docking site parameters for each therapeutic target.
    Saves them to a CSV if the user clicks "Save Docking Parameters."
    """
    # Create columns for auto_detect, center, size
    therapeutic_targets["auto_detect"] = True  # default True
    therapeutic_targets["center"] = ""
    therapeutic_targets["size"] = ""

    for idx, row in therapeutic_targets.iterrows():
        st.write(f"Gene: {row['gene_name']}, Disease Gene ID: {row['disease_gene_id']}")

        # Let user pick auto or manual
        use_auto = st.checkbox(f"Auto-detect docking site for {row['disease_gene_id']}",
                               value=True, key=f"auto_{idx}")
        therapeutic_targets.at[idx, "auto_detect"] = use_auto

        if not use_auto:
            center = st.text_input(
                f"Center (x, y, z) for {row['disease_gene_id']}",
                value="10, 10, 10",
                key=f"center_{idx}"
            )
            size = st.text_input(
                f"Size (x, y, z) for {row['disease_gene_id']}",
                value="20, 20, 20",
                key=f"size_{idx}"
            )
            therapeutic_targets.at[idx, "center"] = center
            therapeutic_targets.at[idx, "size"] = size

    if st.button("Save Docking Parameters"):
        docking_params_path = "results/docking_parameters.csv"
        # Only write needed columns
        # or rename them to consistent naming in your CSV
        therapeutic_targets[[
            "disease_gene_id",
            "gene_name",
            "auto_detect",
            "center",
            "size"
        ]].to_csv(docking_params_path, index=False)
        st.success("Docking parameters saved. Ready for molecular docking.")


def main():
    # 1. Load config
    with open("config.yaml", "r") as f:
        config = yaml.safe_load(f)

    st.title("ADHD Compound Pipeline")

    # 2. Basic pipeline inputs
    geo_id = st.text_input("GEO ID:", value=config["pipeline"]["default_geo_id"])
    description = st.text_input("Description:", value=config["pipeline"]["default_description"])
    compound_name = st.text_input("Natural Compound (Name):", value=config["pipeline"]["default_compound"])

    # 3. Fields for Control and Compound samples
    control_samples_str = st.text_input("Control samples (comma-separated):",
                                        value="GSM2286316,GSM2286317")
    compound_samples_str = st.text_input("Compound samples (comma-separated):",
                                         value="GSM2286238,GSM2286239")

    # Merge them for the Nextflow pipeline
    control_list = [s.strip() for s in control_samples_str.split(",") if s.strip()]
    compound_list = [s.strip() for s in compound_samples_str.split(",") if s.strip()]

    samples_merged = ",".join(control_list + compound_list)
    groups_merged = ",".join(["Control"] * len(control_list) +
                             ["Compound"] * len(compound_list))

    # 4. "Run Pre-Docking Analysis" - calls workflow1.nf
    if st.button("Run Pre-Docking Analysis"):
        params = [
            f"--geo_id={geo_id}",
            f"--samples={samples_merged}",
            f"--groups={groups_merged}",
            f"--compound_name={compound_name}",
            f"--description={description}",
            f"--db_connection_string={config['database']['connection_string']}"
        ]
        success = run_nextflow_workflow("workflow1.nf", params)
        if success:
            st.session_state["pre_docking_completed"] = True

    # 5. If pre-docking completed, show docking input
    if st.session_state.get("pre_docking_completed", False):
        experiment_id_path = "results/experiment_id.txt"
        if os.path.exists(experiment_id_path):
            with open(experiment_id_path, "r") as f:
                experiment_id = int(f.read().strip())

            # Fetch the targets from DB, let user specify docking boxes
            therapeutic_targets = fetch_therapeutic_targets(config["database"]["connection_string"],
                                                            experiment_id)
            display_docking_input(therapeutic_targets)

            # If docking params exist, user can run "Run Molecular Docking"
            docking_params_csv = "results/docking_parameters.csv"
            if os.path.exists(docking_params_csv) and st.button("Run Molecular Docking"):
                # 6. Call workflow2.nf, passing the CSV file path
                params = [
                    f"--db_connection_string={config['database']['connection_string']}",
                    f"--ligand_cid={config['pipeline']['default_cid']}",
                    f"--experiment_id={experiment_id}",
                    f"--docking_params={docking_params_csv}",
                    f"--output_dir={config['paths']['results']}"
                ]
                run_nextflow_workflow("workflow2.nf", params)

if __name__ == "__main__":
    main()

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

def fetch_docking_results(db_connection_string, experiment_id, energy_cutoff):
    """
    Query the database to fetch docking results based on energy cutoff.
    Includes uniprot_id so we can locate the correct PDBs.
    """
    try:
        engine = create_engine(db_connection_string)
        query = """
            SELECT dr.docking_id,
                   dr.uniprot_id,    -- Make sure docking_results has this column
                   dr.ligand_cid,
                   dr.deg_id,
                   dr.binding_energy,
                   dg.gene_name
            FROM docking_results dr
            INNER JOIN degs dg ON dr.deg_id = dg.deg_id
            WHERE dr.experiment_id = %(exp_id)s
              AND dr.binding_energy <= %(cutoff)s;
        """
        params_dict = {"exp_id": experiment_id, "cutoff": energy_cutoff}
        results = pd.read_sql_query(query, con=engine, params=params_dict)
        engine.dispose()
        return results
    except Exception as e:
        st.error(f"Error fetching docking results: {e}")
        return pd.DataFrame()


def display_md_input(docking_results, md_param_path):
    """
    Display input fields for MD simulation parameters for each docking result.
    """
    # Set default simulation parameters
    docking_results["generic_sim_time"] = 100
    docking_results["generic_temp"] = 300
    docking_results["generic_pressure"] = 1
    docking_results["generic_solvent"] = "water"

    individual_params = {}

    for idx, row in docking_results.iterrows():
        # Display uniprot_id along with other info
        st.write(
            f"Gene: {row['gene_name']}, "
            f"Docking ID: {row['docking_id']}, "
            f"UniProt ID: {row.get('uniprot_id', 'N/A')}, "
            f"Binding Energy: {row['binding_energy']}"
        )

        individual = st.checkbox(f"Set individual parameters for Docking ID {row['docking_id']}", key=f"indiv_{idx}")
        if individual:
            sim_time = st.number_input(f"Simulation Time (ps) for {row['docking_id']}", value=100, key=f"sim_time_{idx}")
            temp = st.number_input(f"Temperature (K) for {row['docking_id']}", value=300, key=f"temp_{idx}")
            pressure = st.number_input(f"Pressure (atm) for {row['docking_id']}", value=1, key=f"pressure_{idx}")
            solvent = st.text_input(f"Solvent Type for {row['docking_id']}", value="water", key=f"solvent_{idx}")

            individual_params[row['docking_id']] = {
                "simulation_time": sim_time,
                "temperature": temp,
                "pressure": pressure,
                "solvent": solvent
            }

    if st.button("Save MD Parameters"):
        # Save both generic and individual parameters
        docking_results["individual_params"] = docking_results["docking_id"].map(individual_params)

        # This .to_csv() will include the uniprot_id column (and others)
        docking_results.to_csv(md_param_path, index=False)
        st.success("MD parameters saved. Ready for MD simulation.")


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


def display_docking_input(therapeutic_targets, docking_parameter_path):
    """
    Display input fields for docking site parameters for each therapeutic target.
    Saves them to a CSV if the user clicks "Save Docking Parameters."
    """
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
        # Only write needed columns
        therapeutic_targets[[
            "disease_gene_id",
            "gene_name",
            "auto_detect",
            "center",
            "size"
        ]].to_csv(docking_parameter_path, index=False)
        st.success("Docking parameters saved. Ready for molecular docking.")


def main():
    # 1. Load config
    with open("config.yaml", "r") as f:
        config = yaml.safe_load(f)

    results_path = config["paths"]["results"]
    data_path = config["paths"]["data"]
    experiment_id_path = os.path.join(data_path, "experiment_id.txt")
    docking_parameter_path = os.path.join(data_path, "docking_parameters.csv")
    db_connection_string = config["database"]["connection_string"]

    st.title("ADHD Compound Pipeline")

    # 2. Basic pipeline inputs
    description = st.text_input("Description:", value=config["pipeline"]["default_description"])
    geo_id = st.text_input("GEO ID:", value=config["pipeline"]["default_geo_id"])
    control_samples_str = st.text_input("Control samples (comma-separated):", value="GSM2286316,GSM2286317")
    compound_samples_str = st.text_input("Compound samples (comma-separated):", value="GSM2286238,GSM2286239")
    compound_name = st.text_input("Natural Compound (Name):", value=config["pipeline"]["default_compound"])
    ligand_cid = st.text_input("Compound ID (PubChem):", value=config["pipeline"]["default_cid"])
    adj_p = st.text_input("Adjusted P-Value for DEG analysis:", value=config["pipeline"]["adj_p"])
    raw_p = st.text_input("Raw P-Value for DEG analysis:", value=config["pipeline"]["raw_p"])
    log_fc_up = st.text_input("Log_FC upper minimum for DEG analysis (up-regulation, aka log_fc > x)",
                              value=config["pipeline"]["log_fc_up"])
    log_fc_down = st.text_input("Log_FC lower minimum for DEG analysis (down-regulation, aka log_fc < x)",
                                value=config["pipeline"]["log_fc_down"])

    # 3. Combine them for the Nextflow pipeline
    control_list = [s.strip() for s in control_samples_str.split(",") if s.strip()]
    compound_list = [s.strip() for s in compound_samples_str.split(",") if s.strip()]

    samples_merged = ",".join(control_list + compound_list)
    groups_merged = ",".join(["Control"] * len(control_list) + ["Compound"] * len(compound_list))

    # 4. "Run Pre-Docking Analysis"
    if st.button("Run Pre-Docking Analysis"):
        params = [
            f"--geo_id={geo_id}",
            f"--samples={samples_merged}",
            f"--groups={groups_merged}",
            f"--compound_name={compound_name}",
            f"--description={description}",
            f"--db_connection_string={db_connection_string}",
            f"--adj_p={adj_p}",
            f"--raw_p={raw_p}",
            f"--log_fc_up={log_fc_up}",
            f"--log_fc_down={log_fc_down}"
        ]
        success = run_nextflow_workflow("workflow1.nf", params)
        if success:
            st.session_state["pre_docking_completed"] = True

    # 5. If pre-docking completed, show docking input
    if st.session_state.get("pre_docking_completed", False):
        if os.path.exists(experiment_id_path):
            with open(experiment_id_path, "r") as f:
                experiment_id = int(f.read().strip())

            # Fetch the targets from DB, let user specify docking boxes
            therapeutic_targets = fetch_therapeutic_targets(db_connection_string, experiment_id)
            display_docking_input(therapeutic_targets, docking_parameter_path)

            # If docking params exist, user can run "Run Molecular Docking"
            if os.path.exists(docking_parameter_path) and st.button("Run Molecular Docking"):
                params = [
                    f"--db_connection_string={db_connection_string}",
                    f"--ligand_cid={ligand_cid}",
                    f"--experiment_id={experiment_id}",
                    f"--docking_params={docking_parameter_path}",
                    f"--output_dir={results_path}"
                ]
                success2 = run_nextflow_workflow("workflow2.nf", params)
                if success2:
                    st.session_state["docking_completed"] = True

    # 6. MD simulation
    if st.session_state.get("docking_completed", False):
        st.subheader("Molecular Dynamics Simulation")

        binding_energy_cutoff = st.number_input("Enter Binding Energy Cutoff:", value=-7.0, step=0.1)
        if os.path.exists(experiment_id_path):
            with open(experiment_id_path, "r") as f:
                experiment_id = int(f.read().strip())

            # Now fetch results with uniprot_id
            docking_results = fetch_docking_results(db_connection_string, experiment_id, binding_energy_cutoff)
            if not docking_results.empty:
                # This CSV will store uniprot_id, docking_id, gene_name, etc.
                md_param_path = os.path.join(data_path, "md_parameters.csv")

                # Let user specify or override MD parameters
                display_md_input(docking_results, md_param_path)

                # If the MD param CSV exists, user can run "Run Molecular Dynamics"
                if os.path.exists(md_param_path) and st.button("Run Molecular Dynamics"):
                    params = [
                        f"--db_connection_string={db_connection_string}",
                        f"--experiment_id={experiment_id}",
                        f"--md_param_path={md_param_path}",
                        f"--output_dir={results_path}"
                    ]
                    success3 = run_nextflow_workflow("workflow3.nf", params)
                    if success3:
                        st.success("Molecular Dynamics simulation completed!")


if __name__ == "__main__":
    main()

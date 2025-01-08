import streamlit as st
import yaml
import subprocess
import os

# Load configuration
with open("config.yaml", "r") as file:
    config = yaml.safe_load(file)

# GUI defaults from config
geo_id = st.text_input("GEO ID:", config["pipeline"]["default_geo_id"])
samples = st.text_input("Samples:", config["pipeline"]["default_samples"])
compound = st.text_input("Natural Compound:", config["pipeline"]["default_compound"])

# Checkbox for optional analyses
include_ppi = st.checkbox("Include PPI Network Analysis", value=True)
include_pathway = st.checkbox("Perform Pathway Enrichment Analysis", value=False)

# Button to run the pipeline
run_button = st.button("Run Pipeline")

if run_button:
    st.write(f"Running pipeline with GEO ID: {geo_id}, Samples: {samples}, Compound: {compound}")

    # Build the command
    cmd = [
        "nextflow", "run", ".",
        f"--geo_id={geo_id}",
        f"--samples={samples}",
        f"--natural_compound={compound}"
    ]

    # Add optional parameters
    if include_ppi:
        cmd.append("--ppi_analysis")
    if include_pathway:
        cmd.append("--pathway_analysis")

    # Show the command being run
    st.write("Executing the following command:")
    st.code(" ".join(cmd))

    # Run the pipeline and capture output
    try:
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        for line in process.stdout:
            st.text(line)  # Stream output to the GUI
        process.wait()
        if process.returncode == 0:
            st.success("Pipeline run completed successfully!")
        else:
            st.error(f"Pipeline run failed with return code {process.returncode}")
            st.error(process.stderr.read())
    except Exception as e:
        st.error(f"An error occurred: {e}")

# Display the results
if os.path.exists("results/pathway_enrichment_results.csv"):
    st.write("Pathway Enrichment Results:")
    st.download_button(
        label="Download Results",
        data=open("results/pathway_enrichment_results.csv", "rb").read(),
        file_name="pathway_enrichment_results.csv",
        mime="text/csv"
    )

import streamlit as st
import yaml

# Load configuration
with open("config.yaml", "r") as file:
    config = yaml.safe_load(file)

# GUI defaults from config
geo_id = st.text_input("GEO ID:", config["pipeline"]["default_geo_id"])
samples = st.text_input("Samples:", config["pipeline"]["default_samples"])
compound = st.text_input("Natural Compound:", config["pipeline"]["default_compound"])

run_button = st.button("Run Pipeline")

if run_button:
    st.write(f"Running pipeline with GEO ID: {geo_id}, Samples: {samples}, Compound: {compound}")

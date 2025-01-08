import os
import networkx as nx
import pandas as pd
import requests
import matplotlib.pyplot as plt
import py4cytoscape as p4c
import pickle

def ensure_directories_exist():
    """Ensure necessary directories exist."""
    os.makedirs("../data", exist_ok=True)
    os.makedirs("../results", exist_ok=True)


def build_ppi_network(interactions_file):
    """
    Build a protein-protein interaction (PPI) network from a file.

    Parameters:
        interactions_file (str): Path to the file containing interaction data.

    Returns:
        nx.Graph: A NetworkX graph representing the PPI network.
    """
    interactions = pd.read_csv(interactions_file)

    # Check if required columns exist
    required_columns = {'preferredName_A', 'preferredName_B', 'score'}
    if not required_columns.issubset(interactions.columns):
        raise KeyError(f"Missing required columns in the file. Expected: {required_columns}, Found: {interactions.columns}")

    G = nx.Graph()
    for _, row in interactions.iterrows():
        G.add_edge(row['preferredName_A'], row['preferredName_B'], weight=row['score'])
    return G


def fetch_string_interactions(genes, species=9606, score_threshold=400):
    """
    Fetch PPI data from the STRING database.

    Args:
        genes (list): List of target gene symbols.
        species (int): NCBI species taxonomy ID (default: 9606 for humans).
        score_threshold (int): Minimum interaction score to filter (default: 400).

    Returns:
        pd.DataFrame: DataFrame containing interactions.
    """
    url = "https://string-db.org/api/json/network"
    gene_list = "%0d".join(genes)
    params = {
        "identifiers": gene_list,
        "species": species,
        "required_score": score_threshold,
        "caller_identity": "ADHD_Compound_Pipeline"
    }

    response = requests.get(url, params=params)
    if response.status_code != 200:
        raise Exception(f"Failed to fetch STRING data: {response.status_code}, {response.text}")

    data = response.json()
    return pd.DataFrame(data)


def visualize_ppi_network(ppi_data):
    """
    Visualize the PPI network using networkx and matplotlib.

    Args:
        ppi_data (pd.DataFrame): DataFrame with interaction data.
    """
    G = nx.Graph()
    for _, row in ppi_data.iterrows():
        G.add_edge(row['preferredName_A'], row['preferredName_B'], weight=row['score'])

    plt.figure(figsize=(12, 8))
    pos = nx.spring_layout(G, seed=42)
    nx.draw(G, pos, with_labels=True, node_size=700, font_size=10, font_color='black')
    nx.draw_networkx_edge_labels(G, pos, edge_labels=nx.get_edge_attributes(G, 'weight'))
    plt.title("Protein-Protein Interaction (PPI) Network")
    plt.show()


def load_network_into_cytoscape(ppi_data):
    """
    Load the PPI network into Cytoscape for visualization.

    Args:
        ppi_data (pd.DataFrame): DataFrame with interaction data.
    """
    # Prepare edges
    edges = ppi_data[['preferredName_A', 'preferredName_B', 'score']].copy()
    edges.columns = ['source', 'target', 'interaction']

    # Convert columns to strings
    edges['source'] = edges['source'].astype(str)
    edges['target'] = edges['target'].astype(str)
    edges['interaction'] = edges['interaction'].astype(str)

    # Extract unique nodes and create a node DataFrame with 'id' column
    nodes = pd.DataFrame({
        'id': pd.concat([edges['source'], edges['target']]).unique()
    })

    # Convert node IDs to strings
    nodes['id'] = nodes['id'].astype(str)

    # Handle missing values
    edges.fillna("Unknown", inplace=True)
    nodes.fillna("Unknown", inplace=True)

    # Create network in Cytoscape
    try:
        p4c.create_network_from_data_frames(nodes, edges, title="PPI Network", collection="ADHD_Compound_Pipeline")
        print("Network loaded into Cytoscape!")
    except Exception as e:
        print(f"Could not load network into Cytoscape: {e}")




if __name__ == "__main__":
    ensure_directories_exist()

    # Define input and output files
    interactions_file = "../data/ppi_interactions.csv"
    output_pickle = "../results/ppi_network.gpickle"

    # Fetch interactions from STRING
    target_genes = ["DRD2", "DRD4", "CHRNA3", "CYP1A1", "TNF", "IL6", "KCNJ3"]
    ppi_data = fetch_string_interactions(target_genes)
    ppi_data.to_csv(interactions_file, index=False)

    # Build and save the network
    ppi_network = build_ppi_network(interactions_file)
    with open(output_pickle, "wb") as f:
        pickle.dump(ppi_network, f)

    # Visualize the network
    visualize_ppi_network(ppi_data)

    # Load network into Cytoscape
    try:
        load_network_into_cytoscape(ppi_data)
    except Exception as e:
        print(f"Could not load network into Cytoscape: {e}")

import os
import networkx as nx
import pandas as pd
import requests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sqlalchemy import create_engine, text
import py4cytoscape as p4c
import pickle
import argparse

def fetch_therapeutic_targets(DB_CONNECTION_STRING, experiment_id):
    """
    Fetch therapeutic targets for the given experiment ID.

    Args:
        DB_CONNECTION_STRING (str): Database connection string.
        experiment_id (int): ID of the current experiment.

    Returns:
        list: List of gene symbols associated with therapeutic targets.
    """
    try:
        engine = create_engine(DB_CONNECTION_STRING)
        query = text("""
        SELECT DISTINCT tt.deg_name AS current_symbol
        FROM therapeutic_targets tt
        WHERE tt.experiment_id = :experiment_id;
        """)
        with engine.connect() as connection:
            therapeutic_targets = pd.read_sql_query(query, connection, params={"experiment_id": experiment_id})
        return therapeutic_targets['current_symbol'].tolist()
    except Exception as e:
        print(f"Error fetching therapeutic targets for experiment {experiment_id}: {e}")
        return []

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
    if not genes:
        print("No genes provided for PPI network.")
        return pd.DataFrame()

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

def build_ppi_network(interactions_file):
    """
    Build a protein-protein interaction (PPI) network from a file.

    Args:
        interactions_file (str): Path to the file containing interaction data.

    Returns:
        nx.Graph: A NetworkX graph representing the PPI network.
    """
    interactions = pd.read_csv(interactions_file)

    required_columns = {'preferredName_A', 'preferredName_B', 'score'}
    if not required_columns.issubset(interactions.columns):
        raise KeyError(f"Missing required columns in the file. Expected: {required_columns}, Found: {interactions.columns}")

    G = nx.Graph()
    for _, row in interactions.iterrows():
        G.add_edge(row['preferredName_A'], row['preferredName_B'], weight=row['score'])
    return G

def save_ppi_to_db(ppi_data, experiment_id, DB_CONNECTION_STRING):
    """
    Save PPI network data to the database.

    Args:
        ppi_data (pd.DataFrame): DataFrame with interaction data.
        experiment_id (int): ID of the experiment.
        DB_CONNECTION_STRING (str): Database connection string.
    """
    try:
        engine = create_engine(DB_CONNECTION_STRING)
        query = text("""
        INSERT INTO ppi_network (protein_a, protein_b, interaction_score, experiment_id, ppi_id)
        VALUES (:protein_a, :protein_b, :interaction_score, :experiment_id, DEFAULT)
        """)
        with engine.connect() as connection:
            for _, row in ppi_data.iterrows():
                connection.execute(query, {
                    "protein_a": row['preferredName_A'],
                    "protein_b": row['preferredName_B'],
                    "interaction_score": row['score'],
                    "experiment_id": experiment_id
                })
        print("PPI network results saved to the database.")
    except Exception as e:
        print(f"Error saving PPI network results: {e}")


def visualize_ppi_network(ppi_data, output_path=None, title="Protein-Protein Interaction (PPI) Network",
                          with_labels=True):
    """
    Visualize the PPI network using networkx and matplotlib,
    with node sizes/colors based on their connectivity (degree).

    Args:
        ppi_data (pd.DataFrame): DataFrame with interaction data.
        title (str): Title for the plot.
        output_path (str): If provided, saves the plot to a file.
        with_labels (bool): Whether to draw labels (gene names) on the nodes.
    """
    G = nx.Graph()
    for _, row in ppi_data.iterrows():
        G.add_edge(row["preferredName_A"], row["preferredName_B"], weight=row["score"])

    # Node degrees for sizing
    degrees = dict(G.degree())
    node_sizes = [degrees[node] * 300 for node in G.nodes()]  # scale up

    # Layout
    fig, ax = plt.subplots(figsize=(12, 8))
    pos = nx.spring_layout(G, seed=42)

    # Node color by degree
    node_color = [degrees[node] for node in G.nodes()]
    sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, norm=plt.Normalize(vmin=min(node_color), vmax=max(node_color)))

    # Draw nodes
    nx.draw_networkx_nodes(
        G, pos, node_size=node_sizes, node_color=node_color,
        cmap=plt.cm.viridis, alpha=0.9, ax=ax
    )

    # Draw labels
    if with_labels:
        nx.draw_networkx_labels(G, pos, font_size=9, ax=ax)

    # Draw edges with better visibility
    edges = G.edges(data=True)
    scores = [d["weight"] for (_, _, d) in edges]
    edge_colors = "black"  # Set edge color to gray for better contrast
    edge_widths = [s / 10 for s in scores]  # Scale edge widths based on interaction scores
    nx.draw_networkx_edges(
        G, pos, edge_color=edge_colors, width=edge_widths, alpha=0.7, ax=ax
    )

    # Add colorbar
    fig.colorbar(sm, ax=ax, label="Node Degree")
    ax.set_title(title)

    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        plt.savefig(output_path, dpi=150, bbox_inches="tight")
        print(f"PPI network plot saved to {output_path}")
    plt.show()





def load_network_into_cytoscape(ppi_data):
    """
    Load the PPI network into Cytoscape for visualization.

    Args:
        ppi_data (pd.DataFrame): DataFrame with interaction data.
    """
    edges = ppi_data[['preferredName_A', 'preferredName_B', 'score']].copy()
    edges['preferredName_A'] = edges['preferredName_A'].fillna('').astype(str)
    edges['preferredName_B'] = edges['preferredName_B'].fillna('').astype(str)
    edges['score'] = edges['score'].fillna(0).astype(float)  # Ensure scores are floats


    nodes = pd.DataFrame({
        'id': pd.concat([edges['source'], edges['target']]).unique()
    })

    try:
        p4c.create_network_from_data_frames(nodes, edges, title="PPI Network", collection="ADHD_Compound_Pipeline")
        print("Network loaded into Cytoscape!")
    except Exception as e:
        print(f"Could not load network into Cytoscape: {e}")

def get_experiment_id(experiment_id):
    """
    Parse experiment_id from file if it is a file path.

    Args:
        experiment_id (str): Either the experiment ID as a string or a file path.

    Returns:
        int: Parsed experiment ID.
    """
    if os.path.isfile(experiment_id):
        with open(experiment_id, "r") as f:
            return int(f.read().strip())
    return int(experiment_id)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build PPI network")
    parser.add_argument("--experiment_id", type=str, required=True, help="Experiment ID or file containing the ID")
    parser.add_argument("--db_connection_string", type=str, required=True, help="Database connection string")
    parser.add_argument("--base_dir", type=str, required=True, help="Base directory for relative paths")

    args = parser.parse_args()
    base_dir = args.base_dir

    # Resolve paths
    data_dir = os.path.join(base_dir, "data")
    results_dir = os.path.join(base_dir, "results")
    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)

    # Get the actual experiment ID
    experiment_id = get_experiment_id(args.experiment_id)
    DB_CONNECTION_STRING = args.db_connection_string

    # Fetch therapeutic targets
    target_genes = fetch_therapeutic_targets(DB_CONNECTION_STRING, experiment_id)
    if not target_genes:
        print("No therapeutic targets found. Exiting.")
        exit()

    # Fetch STRING interactions
    ppi_data = fetch_string_interactions(target_genes)
    interactions_file = os.path.join(data_dir, "ppi_interactions.csv")
    ppi_data.to_csv(interactions_file, index=False)

    # Save to the database
    save_ppi_to_db(ppi_data, experiment_id, DB_CONNECTION_STRING)

    # Build the network
    output_pickle = os.path.join(results_dir, "ppi_network.gpickle")
    ppi_network = build_ppi_network(interactions_file)
    with open(output_pickle, "wb") as f:
        pickle.dump(ppi_network, f)

    # Visualize the network
    visualize_ppi_network(
        ppi_data=ppi_data,
        title="ADHD Therapeutic Target PPI",
        output_path=os.path.join(results_dir, "PPI_network.png")
    )
    print(f"Saving PPI network plot to {os.path.join(results_dir, "PPI_network.png")}")


    # Load into Cytoscape
    try:
        load_network_into_cytoscape(ppi_data)
    except Exception as e:
        print(f"Could not load network into Cytoscape: {e}")

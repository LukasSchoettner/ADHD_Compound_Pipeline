import os
import networkx as nx
import pandas as pd
import requests
import matplotlib.pyplot as plt
import psycopg2
import py4cytoscape as p4c
import pickle
import argparse

def ensure_directories_exist():
    """Ensure necessary directories exist."""
    os.makedirs("../data", exist_ok=True)
    os.makedirs("../results", exist_ok=True)


def fetch_degs_disease_matches(DB_CONNECTION_STRING, experiment_id):
    """
    Fetch matching DEGs and disease genes for the given experiment ID.

    Args:
        DB_CONNECTION_STRING (str): Database connection string.
        experiment_id (int): ID of the current experiment.

    Returns:
        list: List of gene symbols matched between DEGs and disease genes.
    """
    try:
        conn = psycopg2.connect(DB_CONNECTION_STRING)
        query = """
        SELECT DISTINCT degs.gene_name
        FROM degs
        INNER JOIN disease_genes
        ON degs.gene_name = disease_genes.gene_name
        WHERE degs.experiment_id = %s;
        """
        matched_genes = pd.read_sql_query(query, conn, params=(experiment_id,))
        print(matched_genes)
        conn.close()
        return matched_genes['gene_name'].tolist()
    except Exception as e:
        print(f"Error fetching matches for experiment {experiment_id}: {e}")
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
        conn = psycopg2.connect(DB_CONNECTION_STRING)
        cursor = conn.cursor()

        for _, row in ppi_data.iterrows():
            cursor.execute(
                """
                INSERT INTO ppi_network (protein_a, protein_b, interaction_score, experiment_id, ppi_id)
                VALUES (%s, %s, %s, %s, DEFAULT)
                """,
                (row['preferredName_A'], row['preferredName_B'], row['score'], experiment_id),
            )

        conn.commit()
        conn.close()
        print("PPI network results saved to the database.")
    except Exception as e:
        print(f"Error saving PPI network results: {e}")


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
    edges = ppi_data[['preferredName_A', 'preferredName_B', 'score']].copy()
    edges.columns = ['source', 'target', 'interaction']

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

    args = parser.parse_args()

    # Get the actual experiment ID
    experiment_id = get_experiment_id(args.experiment_id)
    DB_CONNECTION_STRING = args.db_connection_string

    ensure_directories_exist()

    # Fetch matched DEGs and disease genes
    target_genes = fetch_degs_disease_matches(DB_CONNECTION_STRING, experiment_id)
    if not target_genes:
        print("No matched genes found. Exiting.")
        exit()

    # Fetch STRING interactions
    ppi_data = fetch_string_interactions(target_genes)
    interactions_file = "../data/ppi_interactions.csv"
    ppi_data.to_csv(interactions_file, index=False)

    # Save to the database
    save_ppi_to_db(ppi_data, experiment_id, DB_CONNECTION_STRING)

    # Build the network
    output_pickle = "../results/ppi_network.gpickle"
    ppi_network = build_ppi_network(interactions_file)
    with open(output_pickle, "wb") as f:
        pickle.dump(ppi_network, f)

    # Visualize the network
    visualize_ppi_network(ppi_data)

    # Load into Cytoscape
    try:
        load_network_into_cytoscape(ppi_data)
    except Exception as e:
        print(f"Could not load network into Cytoscape: {e}")
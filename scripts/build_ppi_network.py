import os
import networkx as nx
import pandas as pd
import requests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import psycopg2
import pickle
import argparse


def fetch_therapeutic_targets(db_connection_string, experiment_id):
    """
    Fetch therapeutic targets for the given experiment ID.

    Args:
        db_connection_string (str): Database connection string.
        experiment_id (int): ID of the current experiment.

    Returns:
        list: List of gene symbols associated with therapeutic targets.
    """
    query = """
        SELECT dg.gene_name AS current_symbol
        FROM therapeutic_target tt
        INNER JOIN disease_genes dg ON tt.disease_gene_id = dg.disease_gene_id
        WHERE tt.experiment_id = %s;
    """

    try:
        with psycopg2.connect(db_connection_string) as conn:
            with conn.cursor() as cur:
                cur.execute(query, (experiment_id,))
                results = cur.fetchall()

        # Convert to list after fetching all results
        return [row[0] for row in results]

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


def save_ppi_to_db(ppi_data, experiment_id, db_connection_string):
    """
    Save PPI network data to the database.

    Args:
        ppi_data (pd.DataFrame): DataFrame with interaction data.
        experiment_id (int): ID of the experiment.
        db_connection_string (str): Database connection string.
    """
    query = """
        INSERT INTO ppi_network (protein_a, protein_b, interaction_score, experiment_id, ppi_id)
        VALUES (%s, %s, %s, %s, DEFAULT)
    """

    try:
        with psycopg2.connect(db_connection_string) as conn:
            with conn.cursor() as cur:
                for _, row in ppi_data.iterrows():
                    cur.execute(query, (row['preferredName_A'], row['preferredName_B'], row['score'], experiment_id))
        print("PPI network results saved to the database.")

    except Exception as e:
        print(f"Error saving PPI network results: {e}")


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


def visualize_ppi_network(ppi_data, output_path=None, title="Protein-Protein Interaction (PPI) Network",
                          with_labels=True):
    """
    Visualize the PPI network with improved clustering, edge visibility, and node grouping.

    Args:
        ppi_data (pd.DataFrame): DataFrame with interaction data.
        title (str): Title for the plot.
        output_path (str): If provided, saves the plot to a file.
        with_labels (bool): Whether to draw labels (gene names) on the nodes.
    """
    # Build the graph
    G = nx.Graph()
    for _, row in ppi_data.iterrows():
        G.add_edge(row["preferredName_A"], row["preferredName_B"], weight=row["score"])

    # Node degrees for sizing
    degrees = dict(G.degree())
    node_sizes = [degrees[node] * 300 for node in G.nodes()]

    # Node layout
    pos = nx.spring_layout(G, seed=42, k=0.3)

    # Edge widths and colors
    scores = [d["weight"] for (_, _, d) in G.edges(data=True)]
    edge_widths = [(s / max(scores)) * 5 for s in scores]
    edge_colors = [(s / max(scores)) for s in scores]

    # Plot
    fig, ax = plt.subplots(figsize=(14, 10))
    nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color="blue", alpha=0.7, ax=ax)
    nx.draw_networkx_edges(G, pos, edge_color=edge_colors, width=edge_widths, alpha=0.7, ax=ax)
    if with_labels:
        nx.draw_networkx_labels(G, pos, font_size=8, font_color="black", ax=ax)

    # Title
    ax.set_title(title, fontsize=16)
    plt.tight_layout()
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        print(f"PPI network plot saved to {output_path}")
    plt.show()


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
    db_connection_string = args.db_connection_string

    # Fetch therapeutic targets
    target_genes = fetch_therapeutic_targets(db_connection_string, experiment_id)
    if not target_genes:
        print("No therapeutic targets found. Exiting.")
        exit()

    # Fetch STRING interactions
    ppi_data = fetch_string_interactions(target_genes)
    interactions_file = os.path.join(data_dir, "ppi_interactions.csv")
    ppi_data.to_csv(interactions_file, index=False)

    # Save to the database
    save_ppi_to_db(ppi_data, experiment_id, db_connection_string)

    # Build the network
    output_pickle = os.path.join(results_dir, "ppi_network.gpickle")
    ppi_network = build_ppi_network(interactions_file)
    with open(output_pickle, "wb") as f:
        pickle.dump(ppi_network, f)

    # Visualize the network
    visualize_ppi_network(ppi_data, os.path.join(results_dir, "PPI_network.png"))

    print("PPI pipeline completed successfully.")
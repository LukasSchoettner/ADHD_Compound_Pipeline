import networkx as nx
import pandas as pd

def build_ppi_network(interactions_file):
    """
    Build a protein-protein interaction (PPI) network from a file.

    Parameters:
        interactions_file (str): Path to the file containing interaction data.

    Returns:
        nx.Graph: A NetworkX graph representing the PPI network.
    """
    interactions = pd.read_csv(interactions_file)
    G = nx.Graph()
    for _, row in interactions.iterrows():
        G.add_edge(row['Protein1'], row['Protein2'], weight=row['Score'])
    return G

if __name__ == "__main__":
    interactions_file = "data/ppi_interactions.csv"
    ppi_network = build_ppi_network(interactions_file)
    nx.write_gpickle(ppi_network, "results/ppi_network.gpickle")
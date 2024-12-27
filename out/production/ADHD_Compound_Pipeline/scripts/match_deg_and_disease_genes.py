import pandas as pd

def match_genes(deg_file, disease_genes_file):
    """
    Match differentially expressed genes (DEGs) with disease-associated genes.

    Parameters:
        deg_file (str): Path to the file containing DEGs.
        disease_genes_file (str): Path to the file containing disease-associated genes.

    Returns:
        pd.DataFrame: DataFrame containing matched genes.
    """
    degs = pd.read_csv(deg_file)
    disease_genes = pd.read_csv(disease_genes_file)
    matched_genes = pd.merge(degs, disease_genes, on='Gene')
    return matched_genes

if __name__ == "__main__":
    deg_file = "data/degs.csv"
    disease_genes_file = "data/disease_genes.csv"
    matched_genes = match_genes(deg_file, disease_genes_file)
    matched_genes.to_csv("results/matched_genes.csv", index=False)

import psycopg2
import pandas as pd
import gseapy as gp
import os

def fetch_genes_from_db(table_name, DB_CONNECTION_STRING):
    """
    Fetch genes from a specified table in the database using psycopg2.

    Parameters:
        table_name (str): The name of the table to query.
        DB_CONNECTION_STRING (str): Database connection string.

    Returns:
        pd.DataFrame: DataFrame containing genes from the specified table.
    """
    try:
        # Connect to the database
        conn = psycopg2.connect(DB_CONNECTION_STRING)
        query = f"SELECT gene_name FROM {table_name}"

        # Execute query and fetch data
        genes = pd.read_sql_query(query, conn)
        conn.close()
        return genes
    except Exception as e:
        print(f"Error fetching genes from table {table_name}: {e}")
        return pd.DataFrame()


def get_genes_for_enrichment(DB_CONNECTION_STRING):
    """
    Fetch DEGs and disease genes from the database, and compute their intersection.

    Parameters:
        DB_CONNECTION_STRING (str): Database connection string.

    Returns:
        pd.DataFrame: DataFrame containing intersecting genes.
    """
    # Fetch genes from the database
    degs = fetch_genes_from_db("degs", DB_CONNECTION_STRING)
    disease_genes = fetch_genes_from_db("disease_genes", DB_CONNECTION_STRING)

    # Compute the intersection
    intersection = pd.merge(degs, disease_genes, on="gene_name")
    print(f"Number of intersecting genes: {len(intersection)}")
    return intersection


def perform_pathway_enrichment_from_db(DB_CONNECTION_STRING, output_file, gene_set="KEGG_2021_Human"):
    """
    Perform pathway enrichment analysis using genes fetched from the database.

    Parameters:
        DB_CONNECTION_STRING (str): Database connection string.
        output_file (str): Path to save the enriched pathways results.
        gene_set (str): Gene set database to use (default: KEGG_2021_Human).

    Returns:
        pd.DataFrame: DataFrame containing enriched pathways.
    """
    # Get intersecting genes
    intersecting_genes = get_genes_for_enrichment(DB_CONNECTION_STRING)

    # Perform pathway enrichment if there are genes to analyze
    if intersecting_genes.empty:
        print("No intersecting genes found. Pathway enrichment skipped.")
        return pd.DataFrame()

    # Convert to a list of gene symbols
    gene_list = intersecting_genes['gene'].tolist()

    # Perform pathway enrichment analysis
    try:
        enrichment_results = gp.enrichr(
            gene_list=gene_list,
            gene_sets=[gene_set],
            organism="Human",
            outdir=os.path.dirname(output_file),  # Save results in the output directory
        )

        # Extract and save results
        enriched_df = enrichment_results.results
        enriched_df.to_csv(output_file, index=False)
        print(f"Enrichment analysis completed. Results saved to {output_file}")
        return enriched_df

    except Exception as e:
        print(f"Error during pathway enrichment analysis: {e}")
        return pd.DataFrame()


if __name__ == "__main__":
    # Database connection string
    DB_CONNECTION_STRING = "postgresql://postgres:admin@localhost/adhd_research"

    # Path to save results
    output_file = "../results/enriched_pathways.csv"

    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # Perform pathway enrichment
    enriched_pathways = perform_pathway_enrichment_from_db(DB_CONNECTION_STRING, output_file)

    # Print top results for quick inspection
    if not enriched_pathways.empty:
        print(enriched_pathways.head())

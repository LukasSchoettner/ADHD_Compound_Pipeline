import psycopg2
import pandas as pd
import gseapy as gp
import os
import argparse

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


def perform_pathway_enrichment_from_db(DB_CONNECTION_STRING, experiment_id, gene_set="KEGG_2021_Human"):
    """
    Perform pathway enrichment analysis using genes fetched from the database.

    Args:
        DB_CONNECTION_STRING (str): Database connection string.
        experiment_id (str): ID of the experiment.
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

    gene_list = intersecting_genes['gene_name'].tolist()

    try:
        enrichment_results = gp.enrichr(
            gene_list=gene_list,
            gene_sets=[gene_set],
            organism="Human",
            outdir=None,
        )

        enriched_df = enrichment_results.results
        save_pathway_enrichment_to_db(enriched_df, experiment_id, DB_CONNECTION_STRING)
        return enriched_df

    except Exception as e:
        print(f"Error during pathway enrichment analysis: {e}")
        return pd.DataFrame()


def save_pathway_enrichment_to_db(enriched_df, experiment_id, DB_CONNECTION_STRING):
    """
    Save pathway enrichment results to the database.

    Args:
        enriched_df (pd.DataFrame): DataFrame containing pathway enrichment results.
        experiment_id (str): ID of the experiment.
        DB_CONNECTION_STRING (str): Database connection string.
    """
    try:
        conn = psycopg2.connect(DB_CONNECTION_STRING)
        cursor = conn.cursor()

        for _, row in enriched_df.iterrows():
            cursor.execute(
                """
                INSERT INTO pathway_enrichment (experiment_id, p_value, enrichment_id, pathway_name)
                VALUES (%s, %s, DEFAULT, %s)
                """,
                (experiment_id, row['Adjusted P-value'], row['Term']),
            )

        conn.commit()
        conn.close()
        print("Pathway enrichment results saved to the database.")
    except Exception as e:
        print(f"Error saving pathway enrichment results: {e}")

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



import argparse
import os

# Add the helper function here

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform pathway enrichment analysis")
    parser.add_argument("--experiment_id", type=str, required=True, help="Experiment ID or file containing the ID")
    parser.add_argument("--db_connection_string", type=str, required=True, help="Database connection string")

    args = parser.parse_args()

    # Get the actual experiment ID
    experiment_id = get_experiment_id(args.experiment_id)
    DB_CONNECTION_STRING = args.db_connection_string

    enriched_pathways = perform_pathway_enrichment_from_db(DB_CONNECTION_STRING, experiment_id)
    if not enriched_pathways.empty:
        print("Pathway enrichment completed and saved to the database.")


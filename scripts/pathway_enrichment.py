import psycopg2
import pandas as pd
import gseapy as gp
import os
import argparse


def fetch_genes_from_therapeutic_targets(db_connection_string, experiment_id):
    """
    Fetch gene names from the therapeutic_targets table for a specific experiment.

    Parameters:
        db_connection_string (str): Database connection string.
        experiment_id (int): Experiment ID.

    Returns:
        pd.DataFrame: DataFrame containing gene names.
    """
    try:
        # Connect to the database
        conn = psycopg2.connect(db_connection_string)
        query = """
            SELECT DISTINCT dg.gene_name
            FROM therapeutic_targets tt
            INNER JOIN disease_genes dg ON tt.disease_gene_id = dg.disease_gene_id
            WHERE tt.experiment_id = %s;
        """
        genes = pd.read_sql_query(query, conn, params=(experiment_id,))
        conn.close()
        return genes
    except Exception as e:
        print(f"Error fetching genes from therapeutic_targets for experiment {experiment_id}: {e}")
        return pd.DataFrame()


def perform_pathway_enrichment_from_db(db_connection_string, experiment_id, gene_set="KEGG_2021_Human"):
    """
    Perform pathway enrichment analysis using genes fetched from the therapeutic_targets table.

    Args:
        db_connection_string (str): Database connection string.
        experiment_id (int): ID of the experiment.
        gene_set (str): Gene set database to use (default: KEGG_2021_Human).

    Returns:
        pd.DataFrame: DataFrame containing enriched pathways.
    """
    # Fetch genes from the therapeutic_targets table
    target_genes = fetch_genes_from_therapeutic_targets(db_connection_string, experiment_id)

    # Perform pathway enrichment if there are genes to analyze
    if target_genes.empty:
        print("No genes found in therapeutic_targets. Pathway enrichment skipped.")
        return pd.DataFrame()

    gene_list = target_genes['gene_name'].tolist()

    try:
        enrichment_results = gp.enrichr(
            gene_list=gene_list,
            gene_sets=[gene_set],
            organism="Human",
            outdir=None,
        )

        enriched_df = enrichment_results.results
        save_pathway_enrichment_to_db(enriched_df, experiment_id, db_connection_string)
        return enriched_df

    except Exception as e:
        print(f"Error during pathway enrichment analysis: {e}")
        return pd.DataFrame()


def save_pathway_enrichment_to_db(enriched_df, experiment_id, db_connection_string):
    """
    Save pathway enrichment results to the database.

    Args:
        enriched_df (pd.DataFrame): DataFrame containing pathway enrichment results.
        experiment_id (int): ID of the experiment.
        db_connection_string (str): Database connection string.
    """
    try:
        conn = psycopg2.connect(db_connection_string)
        cursor = conn.cursor()

        for _, row in enriched_df.iterrows():
            cursor.execute(
                """
                INSERT INTO pathway_enrichment (experiment_id, pathway_name, p_value, enrichment_id)
                VALUES (%s, %s, %s, DEFAULT)
                """,
                (experiment_id, row['Term'], row['Adjusted P-value']),
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform pathway enrichment analysis")
    parser.add_argument("--experiment_id", type=str, required=True, help="Experiment ID or file containing the ID")
    parser.add_argument("--db_connection_string", type=str, required=True, help="Database connection string")

    args = parser.parse_args()

    # Get the actual experiment ID
    experiment_id = get_experiment_id(args.experiment_id)
    db_connection_string = args.db_connection_string

    enriched_pathways = perform_pathway_enrichment_from_db(db_connection_string, experiment_id)
    if not enriched_pathways.empty:
        print("Pathway enrichment completed and saved to the database.")

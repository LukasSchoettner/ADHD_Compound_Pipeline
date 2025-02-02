import psycopg2
import pandas as pd
import gseapy as gp
import os
import argparse
import matplotlib.pyplot as plt
import numpy as np

def plot_pathway_enrichment(enriched_df, top_n, output_path):
    """
    Visualize pathway enrichment using a simple barplot of -log10(Adjusted P-value).

    Args:
        enriched_df (pd.DataFrame): DataFrame containing at least 'Term' and 'Adjusted P-value' columns.
        top_n (int): Number of top pathways to show (default = 15).
        output_path (str): If given, saves the plot to this path.
    """
    if enriched_df.empty:
        print("No enrichment data to plot.")
        return

    # Sort by p-value ascending
    df_sorted = enriched_df.sort_values("Adjusted P-value", ascending=True).head(top_n)

    # Create a new column for plotting
    df_sorted["-log10(p)"] = -np.log10(df_sorted["Adjusted P-value"])

    # Plot
    plt.figure(figsize=(8, 6))
    plt.barh(
        y=df_sorted["Term"],
        width=df_sorted["-log10(p)"],
        color="skyblue",
        edgecolor="gray"
    )
    plt.xlabel("-log10(Adjusted P-value)")
    plt.title("Pathway Enrichment")
    plt.gca().invert_yaxis()  # so highest bars on top

    # Optionally save
    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches="tight")
        print(f"Pathway enrichment plot saved to {output_path}")

    plt.show()


def fetch_genes_from_therapeutic_targets(db_connection_string, experiment_id):
    """
    Fetch gene names from the therapeutic_targets table for a specific experiment.
    """
    try:
        # Connect to the database
        conn = psycopg2.connect(db_connection_string)
        cursor = conn.cursor()

        # Execute the query
        query = """
            SELECT DISTINCT dg.gene_name
            FROM therapeutic_target tt
            INNER JOIN disease_genes dg ON tt.disease_gene_id = dg.disease_gene_id
            WHERE tt.experiment_id = %s;
        """
        cursor.execute(query, (experiment_id,))
        rows = cursor.fetchall()  # Fetch all results

        # Convert the results into a DataFrame
        genes_df = pd.DataFrame(rows, columns=["gene_name"])

        # Close the cursor and connection
        cursor.close()
        conn.close()

        return genes_df

    except Exception as e:
        print(f"Error fetching genes from therapeutic_target for experiment {experiment_id}: {e}")
        return pd.DataFrame()


def perform_pathway_enrichment_from_db(
        db_connection_string,
        experiment_id,
        gene_set="KEGG_2021_Human",
        exclude_keywords=None,
        include_neuro=False,
        pval_cutoff=0.05
):
    """
    Perform pathway enrichment analysis using genes fetched from the therapeutic_targets table,
    then post-filter the results to remove irrelevant pathways.

    Args:
        db_connection_string (str): Database connection string.
        experiment_id (int): ID of the experiment.
        gene_set (str): Gene set database to use (default: KEGG_2021_Human).
        exclude_keywords (list): List of substrings to exclude from the 'Term' column
            (e.g. ["hepatitis", "influenza", "cancer"]).
        include_neuro (bool): If True, you can optionally *only* keep terms containing
            "neuro", "brain", "dementia", etc.  (Example usage.)
        pval_cutoff (float): P-value cutoff for filtering (Adjusted P-value).

    Returns:
        pd.DataFrame: DataFrame containing enriched pathways (filtered).
    """
    # Fetch genes
    target_genes = fetch_genes_from_therapeutic_targets(db_connection_string, experiment_id)
    if target_genes.empty:
        print("No genes found in therapeutic_target. Pathway enrichment skipped.")
        return pd.DataFrame()

    gene_list = target_genes["gene_name"].tolist()

    # 1) Perform Enrichr
    try:
        enrichment_results = gp.enrichr(
            gene_list=gene_list,
            gene_sets=[gene_set],
            organism="Human",
            outdir=None,
        )
        enriched_df = enrichment_results.results
    except Exception as e:
        print(f"Error during pathway enrichment analysis: {e}")
        return pd.DataFrame()

    # 2) Basic p-value cutoff
    filtered_df = enriched_df[enriched_df["Adjusted P-value"] <= pval_cutoff]

    # 3) Exclude irrelevant pathways by keyword
    if exclude_keywords is None:
        # You can customize these to your domain:
        exclude_keywords = [
            "hepatitis", "influenza", "hiv", "infectious", "ebola",
            "bacterial", "viral", "cancer", "adenoma", "carcinoma"
            # add or remove as needed
        ]

    def is_excluded(term):
        term_lower = term.lower()
        return any(kw in term_lower for kw in exclude_keywords)

    mask_excluded = filtered_df["Term"].apply(is_excluded)
    # Keep rows where is_excluded == False
    filtered_df = filtered_df[~mask_excluded]

    # 4) (Optional) Restrict to neurological terms only
    if include_neuro:
        # define keywords that must appear somewhere in the Term
        neuro_keywords = ["neuro", "brain", "dementia", "synapse", "cns"]
        def is_neuro(term):
            t = term.lower()
            return any(nk in t for nk in neuro_keywords)
        mask_neuro = filtered_df["Term"].apply(is_neuro)
        filtered_df = filtered_df[mask_neuro]

    # Save to DB
    save_pathway_enrichment_to_db(filtered_df, experiment_id, db_connection_string)
    return filtered_df


def save_pathway_enrichment_to_db(enriched_df, experiment_id, db_connection_string):
    """
    Save pathway enrichment results to the pathway_enrichment table.
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
                (experiment_id, row["Term"], row["Adjusted P-value"]),
            )

        conn.commit()
        conn.close()
        print("Pathway enrichment results saved to the database.")
    except Exception as e:
        print(f"Error saving pathway enrichment results: {e}")


def get_experiment_id(experiment_id):
    """
    Parse experiment_id from file if it is a file path.
    """
    if os.path.isfile(experiment_id):
        with open(experiment_id, "r") as f:
            return int(f.read().strip())
    return int(experiment_id)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform pathway enrichment analysis")
    parser.add_argument("--experiment_id", type=str, required=True, help="Experiment ID or file containing the ID")
    parser.add_argument("--db_connection_string", type=str, required=True, help="Database connection string")
    parser.add_argument("--include_neuro", action="store_true",
                        help="Only keep explicitly neurological pathways (brain/neuro/dementia).")
    parser.add_argument("--base_dir", required=True)
    args = parser.parse_args()
    output_path = args.base_dir + "/results/pathway_enrichment_barplot.png"
    top_n = 15
    # Get the actual experiment ID
    exp_id = get_experiment_id(args.experiment_id)

    # Example usage: you could pass an additional argument for pval cutoff if desired
    enriched_pathways = perform_pathway_enrichment_from_db(
        db_connection_string=args.db_connection_string,
        experiment_id=exp_id,
        # optional overrides:
        # pval_cutoff=0.05,
        # exclude_keywords=["hepatitis","cancer",...],
        include_neuro=args.include_neuro
    )
    if not enriched_pathways.empty:
        print("Pathway enrichment completed and saved to the database.")
        print(enriched_pathways[["Term", "Adjusted P-value"]].head())
        plot_pathway_enrichment(enriched_pathways, top_n, output_path)
    else:
        print("No enriched pathways after filtering.")
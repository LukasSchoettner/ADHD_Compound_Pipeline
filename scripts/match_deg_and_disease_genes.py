import argparse
import os

import psycopg2
import pandas as pd


def create_therapeutic_targets_table(cursor):
    """
    Create the therapeutic_targets table if it doesn't exist, now including uniprot_id and deg_id.
    """
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS therapeutic_targets (
            target_id SERIAL PRIMARY KEY,
            experiment_id INT NOT NULL REFERENCES experiment(experiment_id) ON DELETE CASCADE,
            disease_gene_id TEXT NOT NULL REFERENCES disease_genes(disease_gene_id) ON DELETE CASCADE,
            deg_id INT NOT NULL REFERENCES degs(deg_id) ON DELETE CASCADE,  -- Add deg_id as a foreign key
            uniprot_id VARCHAR(20),  -- store UniProt ID here
            deg_name VARCHAR(255),
            UNIQUE (experiment_id, disease_gene_id, deg_id)  -- Prevent duplicate entries for the same deg_id
        );
    """)


def fetch_and_save_therapeutic_targets(db_connection_string, experiment_id):
    """
    Fetch matching DEGs and disease genes (plus aliases) for the given experiment ID,
    and store them in the therapeutic_targets table with deg_name, deg_id, and uniprot_id.

    Args:
        db_connection_string (str): Database connection string.
        experiment_id (int): ID of the current experiment.
    """
    try:
        with psycopg2.connect(db_connection_string) as conn:
            with conn.cursor() as cursor:
                # 1. Create or ensure the therapeutic_targets table exists
                create_therapeutic_targets_table(cursor)

                # 2. Fetch matches between DEGs and disease_genes (or via gene_aliases).
                #    Include deg_id and uniprot_id in the SELECT.
                query = """
                    SELECT DISTINCT
                        dg.disease_gene_id,
                        dg.uniprot_id,
                        degs.gene_name AS deg_name,
                        degs.deg_id
                    FROM degs
                    INNER JOIN gene_aliases ga
                        ON degs.gene_name = ga.alias
                    INNER JOIN disease_genes dg
                        ON ga.disease_gene_id = dg.disease_gene_id
                    WHERE degs.experiment_id = %s
                
                    UNION
                
                    SELECT DISTINCT
                        dg.disease_gene_id,
                        dg.uniprot_id,
                        degs.gene_name AS deg_name,
                        degs.deg_id
                    FROM degs
                    INNER JOIN disease_genes dg
                        ON degs.gene_name = dg.gene_name
                    WHERE degs.experiment_id = %s;
                """

                matched_genes = pd.read_sql_query(query, conn, params=(experiment_id, experiment_id))

                # 3. Insert the matched genes into therapeutic_targets
                for _, row in matched_genes.iterrows():
                    disease_gene_id = row['disease_gene_id']
                    uniprot_id = row['uniprot_id']
                    deg_name = row['deg_name']
                    deg_id = row['deg_id']

                    cursor.execute("""
                        INSERT INTO therapeutic_targets (experiment_id, disease_gene_id, uniprot_id, deg_name, deg_id)
                        VALUES (%s, %s, %s, %s, %s)
                        ON CONFLICT (experiment_id, disease_gene_id, deg_id) DO NOTHING;
                    """, (experiment_id, disease_gene_id, uniprot_id, deg_name, deg_id))

            conn.commit()

        print(f"Therapeutic targets for experiment {experiment_id} saved successfully.")

    except Exception as e:
        print(f"Error while saving therapeutic targets for experiment {experiment_id}: {e}")


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
    parser = argparse.ArgumentParser(description="Match DEGs and disease genes")
    parser.add_argument("--db_connection_string", required=True,
                        help="Database connection string")
    parser.add_argument("--experiment_id", type=str, required=True,
                        help="Experiment ID or file containing the ID")
    args = parser.parse_args()

    EXPERIMENT_ID = get_experiment_id(args.experiment_id)
    DB_CONNECTION_STRING = args.db_connection_string

    # Fetch and save therapeutic targets (now including uniprot_id)
    fetch_and_save_therapeutic_targets(DB_CONNECTION_STRING, EXPERIMENT_ID)

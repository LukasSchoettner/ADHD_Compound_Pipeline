import argparse
import os

import psycopg2


def fetch_and_save_therapeutic_targets(experiment_id):
    """
    Fetch matching DEGs and disease genes (plus aliases) for the given experiment ID,
    and store them in the therapeutic_targets table with deg_name, deg_id, and uniprot_id.

    Args:
        experiment_id (int): ID of the current experiment.
    """
    try:
        # Connect to the PostgreSQL database
        with psycopg2.connect(
                dbname="adhd_research",
                user="postgres",
                password="admin",
                host="db",
                port="5432"
        ) as conn:
            with conn.cursor() as cursor:
                # Query to fetch matches between DEGs and disease_genes
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

                # Execute the query and fetch results
                cursor.execute(query, (experiment_id, experiment_id))
                matched_genes = cursor.fetchall()

                # Check if any results were returned
                if not matched_genes:
                    print(f"No therapeutic targets found for experiment {experiment_id}.")
                    return

                print(f"Found {len(matched_genes)} therapeutic target matches for experiment {experiment_id}.")

                # Insert the matched genes into therapeutic_target
                insert_query = """
                    INSERT INTO therapeutic_target (experiment_id, disease_gene_id, uniprot_id, deg_name, deg_id)
                    VALUES (%s, %s, %s, %s, %s)
                    ON CONFLICT (experiment_id, disease_gene_id, deg_id) DO NOTHING;
                """

                for row in matched_genes:
                    cursor.execute(insert_query, (experiment_id, *row))

            # Commit the transaction
            conn.commit()

        print(f"Therapeutic targets for experiment {experiment_id} saved successfully.")

    except psycopg2.Error as e:
        print(f"Database error while saving therapeutic targets for experiment {experiment_id}: {e}")
    except Exception as e:
        print(f"Unexpected error: {e}")


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
    fetch_and_save_therapeutic_targets(EXPERIMENT_ID)

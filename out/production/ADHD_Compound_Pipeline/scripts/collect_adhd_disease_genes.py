import os
import psycopg2
import requests
import json

# Set file paths
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
print(BASE_DIR)
DATA_DIR = os.path.join(BASE_DIR, "data")
GENE_FILE = os.path.join(DATA_DIR, "adhd_disease_genes.txt")

# PostgreSQL database connection string
DB_CONNECTION_STRING = "postgresql://postgres:admin@localhost/adhd_research"

# Ensembl API endpoint
ENSEMBL_API = "https://rest.ensembl.org"

# Function to read gene names from the file
def read_genes(file_path):
    with open(file_path, "r") as file:
        content = file.read()
    return [gene.strip() for gene in content.split(",")]

# Function to fetch gene information from Ensembl
def fetch_gene_info(gene_name):
    endpoint = f"/lookup/symbol/homo_sapiens/{gene_name}?content-type=application/json"
    response = requests.get(ENSEMBL_API + endpoint)
    if response.ok:
        return response.json()
    else:
        print(f"Error fetching info for {gene_name}: {response.status_code}, {response.text}")
        return None

# Function to save gene information to PostgreSQL database
def save_to_database(connection_string, gene_data):
    # Connect to PostgreSQL database
    conn = psycopg2.connect(connection_string)
    cursor = conn.cursor()

    # Create table if it doesn't exist
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS disease_genes (
            disease_gene_id TEXT PRIMARY KEY,
            gene_name VARCHAR(255),
            description TEXT,
            biotype VARCHAR(255),
            object_type VARCHAR(255),
            source TEXT,
            start INT,
            end INT
        )
    """)

    # Insert gene data
    for gene in gene_data:
        cursor.execute("""
            INSERT INTO disease_genes (disease_gene_id, gene_name, description, biotype, object_type, source, start, end)
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s)
            ON CONFLICT (disease_gene_id) DO NOTHING
        """, (
            gene["id"],
            gene["display_name"],
            gene.get("description", ""),
            gene.get("biotype", ""),
            gene.get("object_type", ""),
            "Ensembl API",  # Fixed value for source
            gene.get("start", None),
            gene.get("end", None)
        ))

    conn.commit()
    conn.close()

# Main workflow
def main():
    # Read genes from file
    gene_names = read_genes(GENE_FILE)

    # Fetch information for each gene
    gene_data = []
    for gene in gene_names:
        info = fetch_gene_info(gene)
        if info:
            gene_data.append({
                "id": info.get("id"),
                "display_name": info.get("display_name", gene),
                "description": info.get("description", ""),
                "biotype": info.get("biotype", ""),
                "object_type": info.get("object_type", ""),
                "start": info.get("start", None),
                "end": info.get("end", None)
            })

    # Save results to database
    save_to_database(DB_CONNECTION_STRING, gene_data)
    print(f"Saved {len(gene_data)} genes to the database.")

if __name__ == "__main__":
    main()

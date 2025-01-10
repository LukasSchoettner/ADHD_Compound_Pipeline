import os
import psycopg2
import requests
import json
import time

# Set file paths
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # Adjust to parent directory
DATA_DIR = os.path.join(BASE_DIR, "data")
GENE_FILE = os.path.join(DATA_DIR, "adhd_disease_genes.txt")
UNRESOLVED_FILE = os.path.join(DATA_DIR, "unresolved_genes.txt")

# PostgreSQL database connection string
DB_CONNECTION_STRING = "postgresql://postgres:admin@localhost/adhd_research"

# Ensembl and HGNC API endpoints
ENSEMBL_API = "https://rest.ensembl.org"
HGNC_API = "https://rest.genenames.org/fetch/symbol/"
HGNC_PREV_SYMBOL_API = "https://rest.genenames.org/search/prev_symbol/"

session = requests.Session()

# Function to read gene names from the file
def read_genes(file_path):
    with open(file_path, "r") as file:
        content = file.read()
    return [gene.strip() for gene in content.split(",")]

# Function to log unresolved genes
def log_unresolved_gene(gene_name):
    with open(UNRESOLVED_FILE, "a") as file:
        file.write(gene_name + "\n")

# Function to resolve gene names using HGNC
def resolve_gene_name_hgnc(gene_name):
    headers = {"Accept": "application/json"}
    aliases = []
    try:
        response = session.get(HGNC_API + gene_name, headers=headers, timeout=10)
        if response.ok:
            data = response.json()
            if "response" in data and "docs" in data["response"] and len(data["response"]["docs"]) > 0:
                result = data["response"]["docs"][0]
                resolved_name = result.get("symbol")
                aliases.extend(result.get("alias_symbol", []))
                return resolved_name, aliases

        # Check for previous symbols if not resolved
        response = session.get(HGNC_PREV_SYMBOL_API + gene_name, headers=headers, timeout=10)
        if response.ok:
            data = response.json()
            if "response" in data and "docs" in data["response"] and len(data["response"]["docs"]) > 0:
                result = data["response"]["docs"][0]
                resolved_name = result.get("symbol")
                aliases.extend(result.get("alias_symbol", []))
                previous_symbols = result.get("prev_symbol", [])
                aliases.extend(previous_symbols)
                return resolved_name, aliases

    except requests.RequestException as e:
        print(f"HGNC request failed for {gene_name}: {e}")
    log_unresolved_gene(gene_name)
    return gene_name, []  # Return original if no match found

# Function to fetch gene information from Ensembl
def fetch_gene_info(gene_name):
    if gene_name.startswith("ENSG"):  # Handle Ensembl IDs directly
        endpoint = f"/lookup/id/{gene_name}?content-type=application/json"
    else:  # Handle gene symbols
        endpoint = f"/lookup/symbol/homo_sapiens/{gene_name}?content-type=application/json"

    try:
        response = session.get(ENSEMBL_API + endpoint, timeout=20)
        if response.ok:
            return response.json()
        elif response.status_code in [429, 500, 502, 503, 504]:
            print(f"Error for {gene_name}")
    except requests.RequestException as e:
        print(f"Request exception for {gene_name}: {e}")

    return None

# Function to save gene information and aliases to PostgreSQL database
def save_to_database(connection_string, gene_data):
    conn = psycopg2.connect(connection_string)
    cursor = conn.cursor()

    # Create tables if they don't exist
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS disease_genes (
            disease_gene_id TEXT PRIMARY KEY,
            gene_name VARCHAR(255),
            description TEXT,
            biotype VARCHAR(255),
            object_type VARCHAR(255),
            source TEXT,
            start_position INT,
            end_position INT
        )
    """)

    cursor.execute("""
        CREATE TABLE IF NOT EXISTS gene_aliases (
            id SERIAL PRIMARY KEY,
            disease_gene_id TEXT REFERENCES disease_genes(disease_gene_id) ON DELETE CASCADE,
            alias VARCHAR(255) NOT NULL
        )
    """)

    # Insert gene data and aliases
    for gene in gene_data:
        cursor.execute("""
            INSERT INTO disease_genes (disease_gene_id, gene_name, description, biotype, object_type, source, start_position, end_position)
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s)
            ON CONFLICT (disease_gene_id) DO NOTHING
        """, (
            gene["id"],
            gene["display_name"],
            gene.get("description", ""),
            gene.get("biotype", ""),
            gene.get("object_type", ""),
            "Ensembl API",
            gene.get("start", None),
            gene.get("end", None)
        ))

        for alias in gene.get("aliases", []):
            cursor.execute("""
                INSERT INTO gene_aliases (disease_gene_id, alias)
                VALUES (%s, %s)
                ON CONFLICT DO NOTHING
            """, (gene["id"], alias))

    conn.commit()
    conn.close()

# Main workflow
def main():
    # Read genes from file
    gene_names = read_genes(GENE_FILE)

    # Fetch information for each gene
    gene_data = []
    for gene in gene_names:
        resolved_gene_name, aliases = resolve_gene_name_hgnc(gene)  # Resolve gene name using HGNC
        potential_names = [resolved_gene_name] + aliases
        fetched = False
        for name in potential_names:
            info = fetch_gene_info(name)
            if info:
                gene_data.append({
                    "id": info.get("id"),
                    "display_name": info.get("display_name", name),
                    "description": info.get("description", ""),
                    "biotype": info.get("biotype", ""),
                    "object_type": info.get("object_type", ""),
                    "start": info.get("start", None),
                    "end": info.get("end", None),
                    "aliases": aliases
                })
                fetched = True
                break

    # Save results to database
    save_to_database(DB_CONNECTION_STRING, gene_data)
    print(f"Saved {len(gene_data)} genes to the database.")

if __name__ == "__main__":
    main()

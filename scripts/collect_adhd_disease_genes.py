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

session = requests.Session()
x = 0
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
    try:
        response = session.get(HGNC_API + gene_name, headers=headers, timeout=10)
        if response.ok:
            data = response.json()
            if "response" in data and "docs" in data["response"] and len(data["response"]["docs"]) > 0:
                result = data["response"]["docs"][0]
                resolved_name = result.get("symbol")
                aliases = result.get("alias_symbol", [])
                #print(f"Resolved {gene_name} to {resolved_name}, aliases: {aliases}")
                return resolved_name, aliases
        else:
            response = session.get("https://rest.genenames.org/search/prev_symbol/" + gene_name, headers=headers, timeout=10)
            if response.ok:
                data = response.json()
                if "response" in data and "docs" in data["response"] and len(data["response"]["docs"]) > 0:
                    result = data["response"]["docs"][0]
                    resolved_name = result.get("symbol")
                    aliases = result.get("alias_symbol", [])
                    #print(f"Resolved {gene_name} to {resolved_name}, aliases: {aliases}")
                    return resolved_name, aliases

    except requests.RequestException as e:
        print(f"HGNC request failed for {gene_name}: {e}")
    log_unresolved_gene(gene_name)
    return gene_name, []  # Return original if no match found

# Function to fetch gene information from Ensembl with retries and synonyms
def fetch_gene_info(gene_name):
    #print(f"trying to fetch gene {gene_name}")

    if gene_name.startswith("ENSG"):  # Handle Ensembl IDs directly
        #print("using ENSG")
        endpoint = f"/lookup/id/{gene_name}?content-type=application/json"
    else:  # Handle gene symbols
        #print("using symbol")
        endpoint = f"/lookup/symbol/homo_sapiens/{gene_name}?content-type=application/json"

    try:
        response = session.get(ENSEMBL_API + endpoint, timeout=20)
        if response.ok:
            return response.json()
        elif response.status_code in [429, 500, 502, 503, 504]:
            print(f"error for {gene_name}")
    except requests.RequestException as e:
        print(f"Request exception for {gene_name}: {e}")

    #print(f"Failed to fetch info for {gene_name}.")
    return None

# Function to fetch synonyms from Ensembl
def fetch_synonyms_ensembl(gene_name):
    endpoint = f"/xrefs/symbol/homo_sapiens/{gene_name}?content-type=application/json"
    try:
        response = session.get(ENSEMBL_API + endpoint, timeout=10)
        if response.ok:
            data = response.json()
            synonyms = [entry["id"] for entry in data if "id" in entry]
            #print(f"Synonyms for {gene_name}: {synonyms}")
            return synonyms
    except requests.RequestException as e:
        print(f"Ensembl synonym lookup failed for {gene_name}: {e}")
    return []

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
            start_position INT,
            end_position INT
        )
    """)

    # Insert gene data
    sorted_display_names = sorted([item["display_name"] for item in gene_data if "display_name" in item])
    print(len(sorted_display_names))
    for gene in gene_data:
        global x
        x +=1
        #print(gene["display_name"] + "saved to db")
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
            "Ensembl API",  # Fixed value for source
            gene.get("start", None),
            gene.get("end", None)
        ))
    print(x)
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
                    "end": info.get("end", None)
                })
                fetched = True
                break
        if not fetched:
            synonyms = fetch_synonyms_ensembl(gene)

            for synonym in synonyms:
                info = fetch_gene_info(synonym)
                if info:
                    gene_data.append({
                        "id": info.get("id"),
                        "display_name": info.get("display_name", synonym),
                        "description": info.get("description", ""),
                        "biotype": info.get("biotype", ""),
                        "object_type": info.get("object_type", ""),
                        "start": info.get("start", None),
                        "end": info.get("end", None)
                    })
                    break

    # Save results to database
    save_to_database(DB_CONNECTION_STRING, gene_data)
    print(f"Saved {len(gene_data)} genes to the database.")

    # Find missing genes
    missing_genes = find_missing_genes_in_database(DB_CONNECTION_STRING, gene_data)
    print(f"Genes in gene_data but not in the database: {missing_genes}")
    if missing_genes:
        with open("missing_genes.txt", "w") as f:
            f.write("\n".join(missing_genes))
        print(f"Missing genes saved to missing_genes.txt.")

def fetch_database_entries(connection_string):
    conn = psycopg2.connect(connection_string)
    cursor = conn.cursor()
    cursor.execute("SELECT gene_name FROM disease_genes;")
    db_ids = {row[0] for row in cursor.fetchall()}
    conn.close()
    return db_ids

def find_missing_genes_in_database(connection_string, gene_data):
    # Fetch gene IDs from the database
    db_ids = fetch_database_entries(connection_string)  # Function to fetch IDs from the database
    gene_data_ids = {gene["display_name"] for gene in gene_data}

    # Find genes in gene_data but not in the database
    missing_in_db = gene_data_ids - db_ids

    return missing_in_db





if __name__ == "__main__":
    main()

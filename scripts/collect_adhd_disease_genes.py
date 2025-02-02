import os
import requests
import psycopg2

###############################################################################
# Configuration
###############################################################################

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE_DIR, "data")

GENE_FILE = os.path.join(DATA_DIR, "adhd_disease_genes.txt")
UNRESOLVED_FILE = os.path.join(DATA_DIR, "unresolved_genes.txt")

DB_CONNECTION_STRING = "postgresql://postgres:admin@db:5432/adhd_research"

ENSEMBL_API = "https://rest.ensembl.org"
HGNC_API = "https://rest.genenames.org/fetch/symbol/"
HGNC_PREV_SYMBOL_API = "https://rest.genenames.org/search/prev_symbol/"
ALIAS_SYMBOL_API = "https://rest.genenames.org/search/alias_symbol/"

# Requests session to reuse connections
session = requests.Session()

###############################################################################
# File I/O
###############################################################################


def read_genes(file_path):
    """
    Read a single line of comma-separated gene names from file_path.
    Returns a list of gene names.
    """
    with open(file_path, "r") as file:
        content = file.read().strip()
    # If your file is guaranteed to have them all on one line, comma-separated
    return [gene.strip() for gene in content.split(",") if gene.strip()]


def log_unresolved_gene(gene_name):
    """
    Append an unresolved gene name to a text file for reference.
    """
    with open(UNRESOLVED_FILE, "a") as file:
        file.write(gene_name + "\n")

###############################################################################
# API Calls
###############################################################################


def fetch_uniprot_id(symbol):
    """
    Try retrieving a UniProt primaryAccession ID for a given symbol.
    Returns the first matching ID if found, or None otherwise.
    """
    uniprot_api = f"https://rest.uniprot.org/uniprotkb/search?query=gene:{symbol}&fields=accession"
    headers = {"Accept": "application/json"}

    try:
        print(f"[UniProt] Querying for symbol: {symbol}")
        response = session.get(uniprot_api, headers=headers, timeout=10)
        if response.ok:
            data = response.json()
            results = data.get("results", [])
            if results:
                uniprot_id = results[0].get("primaryAccession", None)
                print(f"[UniProt] Found ID for {symbol}: {uniprot_id}")
                return uniprot_id
        else:
            print(f"[UniProt] Non-OK response for {symbol} (HTTP {response.status_code})")
    except requests.RequestException as e:
        print(f"[UniProt] Request failed for {symbol}: {e}")
    return None


def resolve_gene_name_hgnc(gene_name):
    """
    Resolve a gene symbol to its current HGNC symbol and attempt to retrieve a UniProt ID.
    Fallback steps:
      1. Look up current symbol via HGNC_API.
      2. If no UniProt ID, check HGNC previous-symbol endpoint.
      3. If still no UniProt ID, try direct UniProt lookup with the resolved symbol.
      4. If still missing, try each alias for a UniProt ID.

    Returns a tuple: (resolved_symbol, [aliases], uniprot_id).
    """
    headers = {"Accept": "application/json"}
    aliases = []
    uniprot_id = None
    resolved_symbol = gene_name

    try:
        # 1) Check HGNC current symbol
        response = session.get(f"{HGNC_API}{gene_name}", headers=headers, timeout=10)
        if response.ok:
            data = response.json()
            docs = data.get("response", {}).get("docs", [])
            if docs:
                result = docs[0]
                resolved_symbol = result.get("symbol", gene_name)
                aliases.extend(result.get("alias_symbol", []))
                uniprot_ids = result.get("uniprot_ids", [])
                if uniprot_ids:
                    uniprot_id = uniprot_ids[0]

        # 2) If no UniProt, check previous-symbol endpoint
        if not uniprot_id:
            response = session.get(f"{HGNC_PREV_SYMBOL_API}{gene_name}", headers=headers, timeout=10)
            if response.ok:
                data = response.json()
                docs = data.get("response", {}).get("docs", [])
                if docs:
                    result = docs[0]
                    new_symbol = result.get("symbol")
                    if new_symbol:
                        resolved_symbol = new_symbol
                    aliases.extend(result.get("alias_symbol", []))
                    prev_symbols = result.get("prev_symbol", [])
                    aliases.extend(prev_symbols)
                    uniprot_ids = result.get("uniprot_ids", [])
                    if uniprot_ids:
                        uniprot_id = uniprot_ids[0]

        # 3) If still no UniProt, check alias_symbol endpoint
        if not uniprot_id:
            print(f"[HGNC] Checking alias symbols for {gene_name}")
            response = session.get(f"{ALIAS_SYMBOL_API}{gene_name}", headers=headers, timeout=10)
            if response.ok:
                data = response.json()
                docs = data.get("response", {}).get("docs", [])
                if docs:
                    result = docs[0]
                    new_symbol = result.get("symbol")
                    if new_symbol:
                        print(f"[HGNC] Resolved alias {gene_name} -> {new_symbol}")
                        resolved_symbol = new_symbol
                    aliases.extend(result.get("alias_symbol", []))
                    # Also might have prev_symbol, unify them if present
                    prev_symbols = result.get("prev_symbol", [])
                    aliases.extend(prev_symbols)
                    uniprot_ids = result.get("uniprot_ids", [])
                    if uniprot_ids:
                        uniprot_id = uniprot_ids[0]

        # 4) If still no UniProt, try direct UniProt lookup on the resolved name
        if not uniprot_id and resolved_symbol != gene_name:
            uniprot_id = fetch_uniprot_id(resolved_symbol)

        # 5) Finally, if still no ID, try each alias in turn
        if not uniprot_id and aliases:
            for alias in aliases:
                uniprot_id = fetch_uniprot_id(alias)
                if uniprot_id:
                    break

    except requests.RequestException as e:
        print(f"[HGNC] Request failed for {gene_name}: {e}")

    # If we never found a UniProt ID, log the gene
    if not uniprot_id:
        log_unresolved_gene(gene_name)

    return resolved_symbol, aliases, uniprot_id


def fetch_gene_info(gene_name):
    """
    Query Ensembl for gene info. If the name starts with ENSG, assume it's an Ensembl ID,
    otherwise treat it as a symbol.

    Returns the JSON response dict or None if not found / error.
    """
    if gene_name.startswith("ENSG"):
        endpoint = f"/lookup/id/{gene_name}?content-type=application/json"
    else:
        endpoint = f"/lookup/symbol/homo_sapiens/{gene_name}?content-type=application/json"

    url = f"{ENSEMBL_API}{endpoint}"
    try:
        print(f"[Ensembl] Querying {url}")
        response = session.get(url, timeout=20)
        if response.ok:
            return response.json()
        else:
            print(f"[Ensembl] Non-OK response for {gene_name} (HTTP {response.status_code})")
    except requests.RequestException as e:
        print(f"[Ensembl] Request exception for {gene_name}: {e}")
    return None

###############################################################################
# Database
###############################################################################

def create_tables(cursor):
    """
    Create disease_genes and gene_aliases tables if they don't exist.
    You can extend or modify the schema as needed.
    """
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS disease_genes (
            disease_gene_id TEXT PRIMARY KEY,
            gene_name VARCHAR(255),
            description TEXT,
            biotype VARCHAR(255),
            object_type VARCHAR(255),
            source TEXT,
            start_position INT,
            end_position INT,
            uniprot_id VARCHAR(20)
        );
    """)

    cursor.execute("""
        CREATE TABLE IF NOT EXISTS gene_aliases (
            id SERIAL PRIMARY KEY,
            disease_gene_id TEXT REFERENCES disease_genes(disease_gene_id) ON DELETE CASCADE,
            alias VARCHAR(255) NOT NULL,
            CONSTRAINT unique_alias_per_gene UNIQUE (disease_gene_id, alias)
        );
    """)


def save_to_database(gene_data):
    """
    Connect to PostgreSQL using psycopg2 and upsert the given gene_data into
    disease_genes and gene_aliases tables.
    """
    if not gene_data:
        print("[WARNING] No gene data to save to the database.")
        return

    try:
        # Connect to the PostgreSQL database
        connection = psycopg2.connect(
            dbname="adhd_research", user="postgres", password="admin", host="db", port="5432"
        )
        connection.autocommit = False
        cursor = connection.cursor()

        for gene in gene_data:
            # Insert or update the disease_genes table
            cursor.execute(
                """
                INSERT INTO disease_genes (
                    disease_gene_id,
                    gene_name,
                    description,
                    biotype,
                    object_type,
                    source,
                    start_position,
                    end_position,
                    uniprot_id
                ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s)
                ON CONFLICT (disease_gene_id)
                DO UPDATE SET
                    gene_name       = EXCLUDED.gene_name,
                    description     = EXCLUDED.description,
                    biotype         = EXCLUDED.biotype,
                    object_type     = EXCLUDED.object_type,
                    source          = EXCLUDED.source,
                    start_position  = EXCLUDED.start_position,
                    end_position    = EXCLUDED.end_position,
                    uniprot_id      = EXCLUDED.uniprot_id
                """,
                (
                    gene["id"],
                    gene["display_name"],
                    gene.get("description", ""),
                    gene.get("biotype", ""),
                    gene.get("object_type", ""),
                    "Ensembl API",
                    gene.get("start", None),
                    gene.get("end", None),
                    gene.get("uniprot_id", None),
                )
            )

            # Insert aliases into the gene_aliases table
            for alias in gene.get("aliases", []):
                cursor.execute(
                    """
                    INSERT INTO gene_aliases (disease_gene_id, alias)
                    VALUES (%s, %s)
                    ON CONFLICT (disease_gene_id, alias) DO NOTHING
                    """,
                    (gene["id"], alias)
                )

        # Commit the transaction
        connection.commit()
        print(f"[INFO] Successfully saved {len(gene_data)} genes to the database.")

    except Exception as e:
        # Rollback the transaction in case of an error
        if connection:
            connection.rollback()
        print(f"[ERROR] Database operation failed: {e}")

    finally:
        # Close the database connection
        if cursor:
            cursor.close()
        if connection:
            connection.close()




###############################################################################
# Main
###############################################################################

def main():
    # 1. Read genes from the comma-separated file
    gene_names = read_genes(GENE_FILE)
    gene_data = []

    # 2. Process each gene
    for gene in gene_names:
        print(f"\n=== Processing gene: {gene} ===")
        # a) Try to resolve the gene name to the current HGNC symbol
        resolved_name, aliases, uniprot_id = resolve_gene_name_hgnc(gene)

        # b) Attempt to fetch Ensembl info by the resolved name + aliases
        potential_names = [resolved_name] + aliases
        info_found = None

        for name in potential_names:
            info = fetch_gene_info(name)
            if info:  # If we get a valid result from Ensembl, use it
                # Possibly do a final UniProt check if not resolved yet
                if not uniprot_id and name != resolved_name:
                    uniprot_id = fetch_uniprot_id(name)

                gene_data.append({
                    "id": info.get("id"),
                    "display_name": resolved_name,  # Store the official HGNC symbol
                    "description": info.get("description", ""),
                    "biotype": info.get("biotype", ""),
                    "object_type": info.get("object_type", ""),
                    "start": info.get("start"),
                    "end": info.get("end"),
                    "aliases": aliases,
                    "uniprot_id": uniprot_id,
                })
                info_found = True
                break  # No need to check other aliases once found

        # If info_found is not set, that means Ensembl had no record;

    # 3. Save all fetched gene data to PostgreSQL
    save_to_database(gene_data)
    print(f"\n[INFO] Saved {len(gene_data)} genes to the database.")


if __name__ == "__main__":
    main()
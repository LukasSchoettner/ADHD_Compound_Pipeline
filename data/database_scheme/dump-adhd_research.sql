--
-- 1) EXPERIMENT
--
CREATE TABLE experiment (
    experiment_id      SERIAL PRIMARY KEY,
    geo_id             VARCHAR(255),
    compound           TEXT,       -- e.g. main compound name or ID
    description        TEXT,
    created_at         TIMESTAMP DEFAULT now(),
    status            TEXT,        -- e.g. 'Running' / 'Completed'
    adj_p             NUMERIC,     -- or FLOAT
    raw_p             NUMERIC,
    log_fc_up         NUMERIC,
    log_fc_down       NUMERIC,
    degs_found        INTEGER      -- number of DEGs found
);

--
-- 2) DEGS
--
CREATE TABLE degs (
    deg_id            SERIAL PRIMARY KEY,
    experiment_id     INT REFERENCES experiment (experiment_id),
    gene_name         VARCHAR(255),
    log_fold_change   NUMERIC,
    p_value           NUMERIC,
    degs              TEXT,
    pubchem_id        VARCHAR
);

--
-- 3) DISEASE_GENES
--
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

--
-- 4) GENE_ALIASES
--
CREATE TABLE IF NOT EXISTS gene_aliases (
            id SERIAL PRIMARY KEY,
            disease_gene_id TEXT REFERENCES disease_genes(disease_gene_id) ON DELETE CASCADE,
            alias VARCHAR(255) NOT NULL,
            CONSTRAINT unique_alias_per_gene UNIQUE (disease_gene_id, alias)
);

--
-- 5) NATURAL_COMPOUND
--
CREATE TABLE natural_compound (
    compound_id         SERIAL PRIMARY KEY,
    compound_name       VARCHAR(255),
    chemical_structure  VARCHAR(255),
    source             VARCHAR(255),
    pubchem_id         VARCHAR(255),
    molecular_weight   NUMERIC
);

--
-- 6) COMPOUND_TARGET
--
CREATE TABLE compound_target (
    target_id     SERIAL PRIMARY KEY,
    compound_id   INT REFERENCES natural_compound(compound_id),
    gene_name     VARCHAR(255)
);

--
-- 7) DOCKING_RESULTS
--
CREATE TABLE docking_results (
    docking_id         SERIAL PRIMARY KEY,
    experiment_id      INT REFERENCES experiment(experiment_id),
    ligand_cid         VARCHAR(255),         -- external ligand ID
    deg_id             INT,          -- Possibly references degs(deg_id)? 
                                     -- If so: REFERENCES degs(deg_id)
    binding_energy     NUMERIC,
    rmsd_lower_bound   NUMERIC,
    rmsd_upper_bound   NUMERIC,
    docking_pose_rank  INT,
    center_x           NUMERIC,
    center_y           NUMERIC,
    center_z           NUMERIC,
    size_x             NUMERIC,
    size_y             NUMERIC,
    size_z             NUMERIC,
    log_file_path      TEXT,
    created_at         TIMESTAMP DEFAULT now(),
    updated_at         TIMESTAMP DEFAULT now(),
    uniprot_id         VARCHAR(255)
);

--
-- 8) MOLECULAR_DYNAMICS
--
CREATE TABLE molecular_dynamics (
    md_id           SERIAL PRIMARY KEY,
    experiment_id   INT REFERENCES experiment(experiment_id),
    observed_effect TEXT,
    compound_id     INT REFERENCES natural_compound(compound_id),
    simulation_time NUMERIC,
    temperature     NUMERIC,
    pressure        NUMERIC,
    solvent_type    VARCHAR(255),
    status          VARCHAR(255)
);

--
-- 9) PATHWAY_ENRICHMENT
--
CREATE TABLE pathway_enrichment (
    enrichment_id   SERIAL PRIMARY KEY,
    experiment_id   INT REFERENCES experiment(experiment_id),
    pathway_name    TEXT,
    p_value         NUMERIC
);

--
-- 10) PPI_NETWORK
--
CREATE TABLE ppi_network (
    ppi_id           SERIAL PRIMARY KEY,
    protein_a        VARCHAR(255),
    protein_b        VARCHAR(255),
    interaction_score NUMERIC,
    experiment_id    INT REFERENCES experiment(experiment_id)
);

--
-- 11) THERAPEUTIC_TARGET
--
CREATE TABLE therapeutic_targets (
    target_id       SERIAL PRIMARY KEY,
    experiment_id   INT REFERENCES experiment(experiment_id),
    disease_gene_id INT REFERENCES disease_genes(disease_gene_id),
    deg_name        VARCHAR(255),
    uniprot_id      VARCHAR(255),
    deg_id          INT REFERENCES degs(deg_id)
);

-- -----------------------------------------------------------------------------
-- Sequences
-- -----------------------------------------------------------------------------

-- DROP SEQUENCE IF EXISTS compound_target_target_id_seq;
CREATE SEQUENCE IF NOT EXISTS compound_target_target_id_seq
    INCREMENT BY 1
    MINVALUE 1
    MAXVALUE 2147483647
    START 1
    CACHE 1
    NO CYCLE;

-- DROP SEQUENCE IF EXISTS degs_deg_id_seq;
CREATE SEQUENCE IF NOT EXISTS degs_deg_id_seq
    INCREMENT BY 1
    MINVALUE 1
    MAXVALUE 2147483647
    START 1
    CACHE 1
    NO CYCLE;

-- DROP SEQUENCE IF EXISTS docking_results_docking_id_seq;
CREATE SEQUENCE IF NOT EXISTS docking_results_docking_id_seq
    INCREMENT BY 1
    MINVALUE 1
    MAXVALUE 2147483647
    START 1
    CACHE 1
    NO CYCLE;

-- DROP SEQUENCE IF EXISTS experiment_experiment_id_seq;
CREATE SEQUENCE IF NOT EXISTS experiment_experiment_id_seq
    INCREMENT BY 1
    MINVALUE 1
    MAXVALUE 2147483647
    START 1
    CACHE 1
    NO CYCLE;

-- DROP SEQUENCE IF EXISTS gene_aliases_id_seq;
CREATE SEQUENCE IF NOT EXISTS gene_aliases_id_seq
    INCREMENT BY 1
    MINVALUE 1
    MAXVALUE 2147483647
    START 1
    CACHE 1
    NO CYCLE;

-- DROP SEQUENCE IF EXISTS molecular_dynamics_md_id_seq;
CREATE SEQUENCE IF NOT EXISTS molecular_dynamics_md_id_seq
    INCREMENT BY 1
    MINVALUE 1
    MAXVALUE 2147483647
    START 1
    CACHE 1
    NO CYCLE;

-- DROP SEQUENCE IF EXISTS natural_compound_compound_id_seq;
CREATE SEQUENCE IF NOT EXISTS natural_compound_compound_id_seq
    INCREMENT BY 1
    MINVALUE 1
    MAXVALUE 2147483647
    START 1
    CACHE 1
    NO CYCLE;

-- DROP SEQUENCE IF EXISTS pathway_enrichment_enrichment_id_seq;
CREATE SEQUENCE IF NOT EXISTS pathway_enrichment_enrichment_id_seq
    INCREMENT BY 1
    MINVALUE 1
    MAXVALUE 2147483647
    START 1
    CACHE 1
    NO CYCLE;

-- DROP SEQUENCE IF EXISTS ppi_network_ppi_id_seq;
CREATE SEQUENCE IF NOT EXISTS ppi_network_ppi_id_seq
    INCREMENT BY 1
    MINVALUE 1
    MAXVALUE 2147483647
    START 1
    CACHE 1
    NO CYCLE;

-- DROP SEQUENCE IF EXISTS therapeutic_target_target_id_seq;
CREATE SEQUENCE IF NOT EXISTS therapeutic_target_target_id_seq
    INCREMENT BY 1
    MINVALUE 1
    MAXVALUE 2147483647
    START 1
    CACHE 1
    NO CYCLE;

-- -----------------------------------------------------------------------------
-- Tables
-- -----------------------------------------------------------------------------

-- compound_target definition
-- DROP TABLE IF EXISTS compound_target;
CREATE TABLE IF NOT EXISTS compound_target (
    target_id SERIAL4 NOT NULL,
    compound_id INT4 NULL,
    disease_gene_id VARCHAR(20) NULL,  -- Changed from INT4 to VARCHAR(20)
    gene_name VARCHAR(255) NULL,
    CONSTRAINT compound_target_pkey PRIMARY KEY (target_id)
);

-- degs definition
-- DROP TABLE IF EXISTS degs;
CREATE TABLE IF NOT EXISTS degs (
    deg_id SERIAL4 NOT NULL,
    experiment_id INT4 NULL,
    gene_name VARCHAR(255) NULL,
    log_fold_change NUMERIC NULL,
    p_value NUMERIC NULL,
    degs TEXT NULL,
    pubchem_id VARCHAR NULL,
    CONSTRAINT degs_pkey PRIMARY KEY (deg_id)
);

-- disease_genes definition
-- DROP TABLE IF EXISTS disease_genes;
CREATE TABLE IF NOT EXISTS disease_genes (
    disease_gene_id VARCHAR(20) NOT NULL,  -- Changed from serial4 to VARCHAR(20)
    gene_name VARCHAR(255) NULL,
    description VARCHAR(255) NULL,
    "source" VARCHAR(255) NULL,
    biotype VARCHAR(255) NULL,
    object_type VARCHAR(255) NULL,
    start_position INT4 NULL,
    end_position INT4 NULL,
    uniprot_id VARCHAR(255) NULL,
    CONSTRAINT disease_genes_pkey PRIMARY KEY (disease_gene_id),
    CONSTRAINT unique_gene_name UNIQUE (gene_name)
);

-- docking_results definition
-- DROP TABLE IF EXISTS docking_results;
CREATE TABLE IF NOT EXISTS docking_results (
    docking_id SERIAL4 NOT NULL,
    experiment_id INT4 NULL,
    ligand_cid VARCHAR(255) NULL,
    deg_id INT4 NULL,
    binding_energy NUMERIC NULL,
    rmsd_lower_bound NUMERIC NULL,
    rmsd_upper_bound NUMERIC NULL,
    docking_pose_rank INT4 NULL,
    center_x NUMERIC NULL,
    center_y NUMERIC NULL,
    center_z NUMERIC NULL,
    size_x NUMERIC NULL,
    size_y NUMERIC NULL,
    size_z NUMERIC NULL,
    log_file_path TEXT NULL,
    created_at TIMESTAMP DEFAULT NOW() NULL,
    updated_at TIMESTAMP DEFAULT NOW() NULL,
    uniprot_id VARCHAR(255) NULL,
    CONSTRAINT docking_results_pkey PRIMARY KEY (docking_id)
);

-- experiment definition
-- DROP TABLE IF EXISTS experiment;
CREATE TABLE IF NOT EXISTS experiment (
    experiment_id SERIAL4 NOT NULL,
    geo_id VARCHAR(255) NULL,
    compound TEXT NULL,
    description TEXT NULL,
    created_at TIMESTAMP DEFAULT NOW() NULL,
    status TEXT NULL,
    adj_p NUMERIC NULL,
    raw_p NUMERIC NULL,
    log_fc_up NUMERIC NULL,
    log_fc_down NUMERIC NULL,
    degs_found INT4 NULL,
    CONSTRAINT experiment_pkey PRIMARY KEY (experiment_id)
);

-- gene_aliases definition
-- DROP TABLE IF EXISTS gene_aliases;
CREATE TABLE IF NOT EXISTS gene_aliases (
    id SERIAL4 NOT NULL,
    disease_gene_id VARCHAR(20) NULL,  -- Changed from INT4 to VARCHAR(20)
    alias VARCHAR(255) NOT NULL,
    CONSTRAINT gene_aliases_pkey PRIMARY KEY (id),
    CONSTRAINT unique_alias_per_gene UNIQUE (disease_gene_id, alias)
);

-- molecular_dynamics definition
-- DROP TABLE IF EXISTS molecular_dynamics;
CREATE TABLE IF NOT EXISTS molecular_dynamics (
    md_id SERIAL4 NOT NULL,
    experiment_id INT4 NULL,
    observed_effect TEXT NULL,
    compound_id INT4 NULL,
    simulation_time NUMERIC NULL,
    temperature NUMERIC NULL,
    pressure NUMERIC NULL,
    solvent_type VARCHAR(255) NULL,
    status VARCHAR(255) NULL,
    CONSTRAINT molecular_dynamics_pkey PRIMARY KEY (md_id)
);

-- natural_compound definition
-- DROP TABLE IF EXISTS natural_compound;
CREATE TABLE IF NOT EXISTS natural_compound (
    compound_id SERIAL4 NOT NULL,
    compound_name VARCHAR(255) NULL,
    chemical_structure VARCHAR(255) NULL,
    "source" VARCHAR(255) NULL,
    pubchem_id VARCHAR(255) NULL,
    molecular_weight NUMERIC NULL,
    CONSTRAINT natural_compound_pkey PRIMARY KEY (compound_id)
);

-- pathway_enrichment definition
-- DROP TABLE IF EXISTS pathway_enrichment;
CREATE TABLE IF NOT EXISTS pathway_enrichment (
    enrichment_id SERIAL4 NOT NULL,
    experiment_id INT4 NULL,
    pathway_name TEXT NULL,
    p_value NUMERIC NULL,
    CONSTRAINT pathway_enrichment_pkey PRIMARY KEY (enrichment_id)
);

-- ppi_network definition
-- DROP TABLE IF EXISTS ppi_network;
CREATE TABLE IF NOT EXISTS ppi_network (
    ppi_id SERIAL4 NOT NULL,
    protein_a VARCHAR(255) NULL,
    protein_b VARCHAR(255) NULL,
    interaction_score NUMERIC NULL,
    experiment_id INT4 NULL,
    CONSTRAINT ppi_network_pkey PRIMARY KEY (ppi_id)
);

-- therapeutic_target definition
-- DROP TABLE IF EXISTS therapeutic_target;
CREATE TABLE IF NOT EXISTS therapeutic_target (
    target_id SERIAL4 NOT NULL,
    experiment_id INT4 NULL,
    disease_gene_id VARCHAR(20) NULL,  -- Changed from INT4 to VARCHAR(20)
    deg_name VARCHAR(255) NULL,
    uniprot_id VARCHAR(255) NULL,
    deg_id INT4 NULL,
    CONSTRAINT therapeutic_target_pkey PRIMARY KEY (target_id)
);

-- -----------------------------------------------------------------------------
-- Foreign Key Constraints
-- -----------------------------------------------------------------------------

-- Add foreign keys for compound_target
ALTER TABLE IF EXISTS compound_target
    ADD CONSTRAINT fk_compound_target_disease_gene
        FOREIGN KEY (disease_gene_id)
        REFERENCES disease_genes(disease_gene_id)
        ON DELETE SET NULL
        ON UPDATE CASCADE;

-- Add foreign keys for gene_aliases
ALTER TABLE IF EXISTS gene_aliases
    ADD CONSTRAINT fk_gene_aliases_disease_gene
        FOREIGN KEY (disease_gene_id)
        REFERENCES disease_genes(disease_gene_id)
        ON DELETE CASCADE
        ON UPDATE CASCADE;

-- Add foreign keys for therapeutic_target
ALTER TABLE IF EXISTS therapeutic_target
    ADD CONSTRAINT fk_therapeutic_target_disease_gene
        FOREIGN KEY (disease_gene_id)
        REFERENCES disease_genes(disease_gene_id)
        ON DELETE SET NULL
        ON UPDATE CASCADE;

ALTER TABLE IF EXISTS degs
    ADD CONSTRAINT fk_degs_experiment
        FOREIGN KEY (experiment_id)
        REFERENCES experiment(experiment_id)
        ON DELETE SET NULL
        ON UPDATE CASCADE;

ALTER TABLE IF EXISTS docking_results
    ADD CONSTRAINT fk_docking_results_experiment
        FOREIGN KEY (experiment_id)
        REFERENCES experiment(experiment_id)
        ON DELETE SET NULL
        ON UPDATE CASCADE;

ALTER TABLE IF EXISTS docking_results
    ADD CONSTRAINT fk_docking_results_degs
        FOREIGN KEY (deg_id)
        REFERENCES degs(deg_id)
        ON DELETE SET NULL
        ON UPDATE CASCADE;

ALTER TABLE IF EXISTS molecular_dynamics
    ADD CONSTRAINT fk_molecular_dynamics_experiment
        FOREIGN KEY (experiment_id)
        REFERENCES experiment(experiment_id)
        ON DELETE SET NULL
        ON UPDATE CASCADE;

ALTER TABLE IF EXISTS molecular_dynamics
    ADD CONSTRAINT fk_molecular_dynamics_compound
        FOREIGN KEY (compound_id)
        REFERENCES natural_compound(compound_id)
        ON DELETE SET NULL
        ON UPDATE CASCADE;

ALTER TABLE IF EXISTS ppi_network
    ADD CONSTRAINT fk_ppi_network_experiment
        FOREIGN KEY (experiment_id)
        REFERENCES experiment(experiment_id)
        ON DELETE SET NULL
        ON UPDATE CASCADE;

ALTER TABLE IF EXISTS pathway_enrichment
    ADD CONSTRAINT fk_pathway_enrichment_experiment
        FOREIGN KEY (experiment_id)
        REFERENCES experiment(experiment_id)
        ON DELETE SET NULL
        ON UPDATE CASCADE;

ALTER TABLE IF EXISTS therapeutic_target
    ADD CONSTRAINT fk_therapeutic_target_degs
        FOREIGN KEY (deg_id)
        REFERENCES degs(deg_id)
        ON DELETE SET NULL
        ON UPDATE CASCADE;

-- Add any additional foreign key constraints here similarly.

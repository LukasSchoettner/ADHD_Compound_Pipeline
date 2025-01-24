PGDMP      6                 }            adhd_research    16.6    16.6 B    �           0    0    ENCODING    ENCODING        SET client_encoding = 'UTF8';
                      false            �           0    0 
   STDSTRINGS 
   STDSTRINGS     (   SET standard_conforming_strings = 'on';
                      false            �           0    0 
   SEARCHPATH 
   SEARCHPATH     8   SELECT pg_catalog.set_config('search_path', '', false);
                      false            �           1262    16484    adhd_research    DATABASE     u   CREATE DATABASE adhd_research WITH TEMPLATE = template0 ENCODING = 'UTF8' LOCALE_PROVIDER = libc LOCALE = 'C.UTF-8';
    DROP DATABASE adhd_research;
                postgres    false                        2615    2200    public    SCHEMA        CREATE SCHEMA public;
    DROP SCHEMA public;
                pg_database_owner    false            �           0    0    SCHEMA public    COMMENT     6   COMMENT ON SCHEMA public IS 'standard public schema';
                   pg_database_owner    false    4            �            1259    16587    degs    TABLE     �   CREATE TABLE public.degs (
    deg_id integer NOT NULL,
    experiment_id integer NOT NULL,
    gene_name character varying(255),
    log_fold_change double precision,
    p_value double precision,
    degs character varying(50)
);
    DROP TABLE public.degs;
       public         heap    postgres    false    4            �            1259    16586    degs_deg_id_seq    SEQUENCE     �   CREATE SEQUENCE public.degs_deg_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;
 &   DROP SEQUENCE public.degs_deg_id_seq;
       public          postgres    false    222    4            �           0    0    degs_deg_id_seq    SEQUENCE OWNED BY     C   ALTER SEQUENCE public.degs_deg_id_seq OWNED BY public.degs.deg_id;
          public          postgres    false    221            �            1259    16486    disease_genes    TABLE     O  CREATE TABLE public.disease_genes (
    disease_gene_id text NOT NULL,
    gene_name character varying(255) NOT NULL,
    description text,
    source character varying(255),
    biotype character varying,
    object_type character varying,
    start_position integer,
    end_position integer,
    uniprot_id character varying(20)
);
 !   DROP TABLE public.disease_genes;
       public         heap    postgres    false    4            �            1259    16485 !   disease_genes_disease_gene_id_seq    SEQUENCE     �   CREATE SEQUENCE public.disease_genes_disease_gene_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;
 8   DROP SEQUENCE public.disease_genes_disease_gene_id_seq;
       public          postgres    false    216    4            �           0    0 !   disease_genes_disease_gene_id_seq    SEQUENCE OWNED BY     g   ALTER SEQUENCE public.disease_genes_disease_gene_id_seq OWNED BY public.disease_genes.disease_gene_id;
          public          postgres    false    215            �            1259    24922    docking_results    TABLE     �  CREATE TABLE public.docking_results (
    docking_id integer NOT NULL,
    experiment_id integer,
    ligand_cid character varying(50),
    deg_id integer,
    binding_energy double precision NOT NULL,
    rmsd_lower_bound double precision,
    rmsd_upper_bound double precision,
    docking_pose_rank integer DEFAULT 1,
    center_x double precision,
    center_y double precision,
    center_z double precision,
    size_x double precision,
    size_y double precision,
    size_z double precision,
    log_file_path character varying(255),
    created_at timestamp without time zone DEFAULT now(),
    updated_at timestamp without time zone DEFAULT now(),
    uniprot_id character varying
);
 #   DROP TABLE public.docking_results;
       public         heap    postgres    false    4            �            1259    24921    docking_results_docking_id_seq    SEQUENCE     �   CREATE SEQUENCE public.docking_results_docking_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;
 5   DROP SEQUENCE public.docking_results_docking_id_seq;
       public          postgres    false    230    4            �           0    0    docking_results_docking_id_seq    SEQUENCE OWNED BY     a   ALTER SEQUENCE public.docking_results_docking_id_seq OWNED BY public.docking_results.docking_id;
          public          postgres    false    229            �            1259    16577 
   experiment    TABLE     �  CREATE TABLE public.experiment (
    experiment_id integer NOT NULL,
    geo_id character varying(50) NOT NULL,
    compound character varying(255) NOT NULL,
    description text,
    created_at timestamp without time zone DEFAULT now(),
    status text,
    adj_p numeric(10,5) DEFAULT 0.05 NOT NULL,
    raw_p numeric(10,5) DEFAULT 0.01 NOT NULL,
    log_fc_up numeric(10,5) DEFAULT 1.5 NOT NULL,
    log_fc_down numeric(10,5) DEFAULT '-1.5'::numeric NOT NULL,
    degs_found integer DEFAULT 0 NOT NULL
);
    DROP TABLE public.experiment;
       public         heap    postgres    false    4            �            1259    16576    experiment_experiment_id_seq    SEQUENCE     �   CREATE SEQUENCE public.experiment_experiment_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;
 3   DROP SEQUENCE public.experiment_experiment_id_seq;
       public          postgres    false    4    220            �           0    0    experiment_experiment_id_seq    SEQUENCE OWNED BY     ]   ALTER SEQUENCE public.experiment_experiment_id_seq OWNED BY public.experiment.experiment_id;
          public          postgres    false    219            �            1259    16728    gene_aliases    TABLE     �   CREATE TABLE public.gene_aliases (
    id integer NOT NULL,
    disease_gene_id text NOT NULL,
    alias character varying(255) NOT NULL
);
     DROP TABLE public.gene_aliases;
       public         heap    postgres    false    4            �            1259    16727    gene_aliases_id_seq    SEQUENCE     �   CREATE SEQUENCE public.gene_aliases_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;
 *   DROP SEQUENCE public.gene_aliases_id_seq;
       public          postgres    false    4    226            �           0    0    gene_aliases_id_seq    SEQUENCE OWNED BY     K   ALTER SEQUENCE public.gene_aliases_id_seq OWNED BY public.gene_aliases.id;
          public          postgres    false    225            �            1259    33122    molecular_dynamics    TABLE     $  CREATE TABLE public.molecular_dynamics (
    md_id integer NOT NULL,
    experiment_id integer,
    pubchem_id character varying,
    observed_effect text,
    simulation_time integer,
    temperature double precision,
    pressure double precision,
    solvent_type text,
    status text
);
 &   DROP TABLE public.molecular_dynamics;
       public         heap    postgres    false    4            �            1259    33121    molecular_dynamics_md_id_seq    SEQUENCE     �   CREATE SEQUENCE public.molecular_dynamics_md_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;
 3   DROP SEQUENCE public.molecular_dynamics_md_id_seq;
       public          postgres    false    4    233            �           0    0    molecular_dynamics_md_id_seq    SEQUENCE OWNED BY     ]   ALTER SEQUENCE public.molecular_dynamics_md_id_seq OWNED BY public.molecular_dynamics.md_id;
          public          postgres    false    232            �            1259    33113    natural_compounds    TABLE     �   CREATE TABLE public.natural_compounds (
    pubchem_id character varying NOT NULL,
    compound_name character varying,
    chemical_structure text,
    source character varying,
    molecular_weight double precision
);
 %   DROP TABLE public.natural_compounds;
       public         heap    postgres    false    4            �            1259    16614    pathway_enrichment    TABLE     �   CREATE TABLE public.pathway_enrichment (
    enrichment_id integer NOT NULL,
    experiment_id integer NOT NULL,
    pathway_name character varying(255),
    p_value double precision
);
 &   DROP TABLE public.pathway_enrichment;
       public         heap    postgres    false    4            �            1259    16613 $   pathway_enrichment_enrichment_id_seq    SEQUENCE     �   CREATE SEQUENCE public.pathway_enrichment_enrichment_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;
 ;   DROP SEQUENCE public.pathway_enrichment_enrichment_id_seq;
       public          postgres    false    224    4            �           0    0 $   pathway_enrichment_enrichment_id_seq    SEQUENCE OWNED BY     m   ALTER SEQUENCE public.pathway_enrichment_enrichment_id_seq OWNED BY public.pathway_enrichment.enrichment_id;
          public          postgres    false    223            �            1259    16523    ppi_network    TABLE     �   CREATE TABLE public.ppi_network (
    ppi_id integer NOT NULL,
    protein_a character varying(255) NOT NULL,
    protein_b character varying(255) NOT NULL,
    interaction_score double precision,
    experiment_id integer
);
    DROP TABLE public.ppi_network;
       public         heap    postgres    false    4            �            1259    16522    ppi_network_ppi_id_seq    SEQUENCE     �   CREATE SEQUENCE public.ppi_network_ppi_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;
 -   DROP SEQUENCE public.ppi_network_ppi_id_seq;
       public          postgres    false    4    218            �           0    0    ppi_network_ppi_id_seq    SEQUENCE OWNED BY     Q   ALTER SEQUENCE public.ppi_network_ppi_id_seq OWNED BY public.ppi_network.ppi_id;
          public          postgres    false    217            �            1259    16805    therapeutic_targets    TABLE     �   CREATE TABLE public.therapeutic_targets (
    target_id integer NOT NULL,
    experiment_id integer NOT NULL,
    disease_gene_id text NOT NULL,
    deg_name character varying(255),
    uniprot_id character varying,
    deg_id integer NOT NULL
);
 '   DROP TABLE public.therapeutic_targets;
       public         heap    postgres    false    4            �            1259    16804 !   therapeutic_targets_target_id_seq    SEQUENCE     �   CREATE SEQUENCE public.therapeutic_targets_target_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;
 8   DROP SEQUENCE public.therapeutic_targets_target_id_seq;
       public          postgres    false    4    228            �           0    0 !   therapeutic_targets_target_id_seq    SEQUENCE OWNED BY     g   ALTER SEQUENCE public.therapeutic_targets_target_id_seq OWNED BY public.therapeutic_targets.target_id;
          public          postgres    false    227            �           2604    16590    degs deg_id    DEFAULT     j   ALTER TABLE ONLY public.degs ALTER COLUMN deg_id SET DEFAULT nextval('public.degs_deg_id_seq'::regclass);
 :   ALTER TABLE public.degs ALTER COLUMN deg_id DROP DEFAULT;
       public          postgres    false    221    222    222            �           2604    16566    disease_genes disease_gene_id    DEFAULT     �   ALTER TABLE ONLY public.disease_genes ALTER COLUMN disease_gene_id SET DEFAULT nextval('public.disease_genes_disease_gene_id_seq'::regclass);
 L   ALTER TABLE public.disease_genes ALTER COLUMN disease_gene_id DROP DEFAULT;
       public          postgres    false    216    215    216            �           2604    24925    docking_results docking_id    DEFAULT     �   ALTER TABLE ONLY public.docking_results ALTER COLUMN docking_id SET DEFAULT nextval('public.docking_results_docking_id_seq'::regclass);
 I   ALTER TABLE public.docking_results ALTER COLUMN docking_id DROP DEFAULT;
       public          postgres    false    230    229    230            �           2604    16580    experiment experiment_id    DEFAULT     �   ALTER TABLE ONLY public.experiment ALTER COLUMN experiment_id SET DEFAULT nextval('public.experiment_experiment_id_seq'::regclass);
 G   ALTER TABLE public.experiment ALTER COLUMN experiment_id DROP DEFAULT;
       public          postgres    false    219    220    220            �           2604    16731    gene_aliases id    DEFAULT     r   ALTER TABLE ONLY public.gene_aliases ALTER COLUMN id SET DEFAULT nextval('public.gene_aliases_id_seq'::regclass);
 >   ALTER TABLE public.gene_aliases ALTER COLUMN id DROP DEFAULT;
       public          postgres    false    225    226    226            �           2604    33125    molecular_dynamics md_id    DEFAULT     �   ALTER TABLE ONLY public.molecular_dynamics ALTER COLUMN md_id SET DEFAULT nextval('public.molecular_dynamics_md_id_seq'::regclass);
 G   ALTER TABLE public.molecular_dynamics ALTER COLUMN md_id DROP DEFAULT;
       public          postgres    false    233    232    233            �           2604    16617     pathway_enrichment enrichment_id    DEFAULT     �   ALTER TABLE ONLY public.pathway_enrichment ALTER COLUMN enrichment_id SET DEFAULT nextval('public.pathway_enrichment_enrichment_id_seq'::regclass);
 O   ALTER TABLE public.pathway_enrichment ALTER COLUMN enrichment_id DROP DEFAULT;
       public          postgres    false    223    224    224            �           2604    16526    ppi_network ppi_id    DEFAULT     x   ALTER TABLE ONLY public.ppi_network ALTER COLUMN ppi_id SET DEFAULT nextval('public.ppi_network_ppi_id_seq'::regclass);
 A   ALTER TABLE public.ppi_network ALTER COLUMN ppi_id DROP DEFAULT;
       public          postgres    false    218    217    218            �           2604    16808    therapeutic_targets target_id    DEFAULT     �   ALTER TABLE ONLY public.therapeutic_targets ALTER COLUMN target_id SET DEFAULT nextval('public.therapeutic_targets_target_id_seq'::regclass);
 L   ALTER TABLE public.therapeutic_targets ALTER COLUMN target_id DROP DEFAULT;
       public          postgres    false    227    228    228            �           2606    16592    degs degs_pkey 
   CONSTRAINT     P   ALTER TABLE ONLY public.degs
    ADD CONSTRAINT degs_pkey PRIMARY KEY (deg_id);
 8   ALTER TABLE ONLY public.degs DROP CONSTRAINT degs_pkey;
       public            postgres    false    222            �           2606    16568     disease_genes disease_genes_pkey 
   CONSTRAINT     k   ALTER TABLE ONLY public.disease_genes
    ADD CONSTRAINT disease_genes_pkey PRIMARY KEY (disease_gene_id);
 J   ALTER TABLE ONLY public.disease_genes DROP CONSTRAINT disease_genes_pkey;
       public            postgres    false    216                       2606    24930 $   docking_results docking_results_pkey 
   CONSTRAINT     j   ALTER TABLE ONLY public.docking_results
    ADD CONSTRAINT docking_results_pkey PRIMARY KEY (docking_id);
 N   ALTER TABLE ONLY public.docking_results DROP CONSTRAINT docking_results_pkey;
       public            postgres    false    230            �           2606    16585    experiment experiment_pkey 
   CONSTRAINT     c   ALTER TABLE ONLY public.experiment
    ADD CONSTRAINT experiment_pkey PRIMARY KEY (experiment_id);
 D   ALTER TABLE ONLY public.experiment DROP CONSTRAINT experiment_pkey;
       public            postgres    false    220            �           2606    16735    gene_aliases gene_aliases_pkey 
   CONSTRAINT     \   ALTER TABLE ONLY public.gene_aliases
    ADD CONSTRAINT gene_aliases_pkey PRIMARY KEY (id);
 H   ALTER TABLE ONLY public.gene_aliases DROP CONSTRAINT gene_aliases_pkey;
       public            postgres    false    226                       2606    33129 *   molecular_dynamics molecular_dynamics_pkey 
   CONSTRAINT     k   ALTER TABLE ONLY public.molecular_dynamics
    ADD CONSTRAINT molecular_dynamics_pkey PRIMARY KEY (md_id);
 T   ALTER TABLE ONLY public.molecular_dynamics DROP CONSTRAINT molecular_dynamics_pkey;
       public            postgres    false    233                       2606    33119 (   natural_compounds natural_compounds_pkey 
   CONSTRAINT     n   ALTER TABLE ONLY public.natural_compounds
    ADD CONSTRAINT natural_compounds_pkey PRIMARY KEY (pubchem_id);
 R   ALTER TABLE ONLY public.natural_compounds DROP CONSTRAINT natural_compounds_pkey;
       public            postgres    false    231            �           2606    16619 *   pathway_enrichment pathway_enrichment_pkey 
   CONSTRAINT     s   ALTER TABLE ONLY public.pathway_enrichment
    ADD CONSTRAINT pathway_enrichment_pkey PRIMARY KEY (enrichment_id);
 T   ALTER TABLE ONLY public.pathway_enrichment DROP CONSTRAINT pathway_enrichment_pkey;
       public            postgres    false    224            �           2606    16530    ppi_network ppi_network_pkey 
   CONSTRAINT     ^   ALTER TABLE ONLY public.ppi_network
    ADD CONSTRAINT ppi_network_pkey PRIMARY KEY (ppi_id);
 F   ALTER TABLE ONLY public.ppi_network DROP CONSTRAINT ppi_network_pkey;
       public            postgres    false    218            �           2606    16814 I   therapeutic_targets therapeutic_targets_experiment_id_disease_gene_id_key 
   CONSTRAINT     �   ALTER TABLE ONLY public.therapeutic_targets
    ADD CONSTRAINT therapeutic_targets_experiment_id_disease_gene_id_key UNIQUE (experiment_id, disease_gene_id);
 s   ALTER TABLE ONLY public.therapeutic_targets DROP CONSTRAINT therapeutic_targets_experiment_id_disease_gene_id_key;
       public            postgres    false    228    228                        2606    16812 ,   therapeutic_targets therapeutic_targets_pkey 
   CONSTRAINT     q   ALTER TABLE ONLY public.therapeutic_targets
    ADD CONSTRAINT therapeutic_targets_pkey PRIMARY KEY (target_id);
 V   ALTER TABLE ONLY public.therapeutic_targets DROP CONSTRAINT therapeutic_targets_pkey;
       public            postgres    false    228                       2606    24980 9   therapeutic_targets therapeutic_targets_unique_constraint 
   CONSTRAINT     �   ALTER TABLE ONLY public.therapeutic_targets
    ADD CONSTRAINT therapeutic_targets_unique_constraint UNIQUE (experiment_id, disease_gene_id, deg_id);
 c   ALTER TABLE ONLY public.therapeutic_targets DROP CONSTRAINT therapeutic_targets_unique_constraint;
       public            postgres    false    228    228    228            
           2606    16593    degs degs_experiment_id_fkey    FK CONSTRAINT     �   ALTER TABLE ONLY public.degs
    ADD CONSTRAINT degs_experiment_id_fkey FOREIGN KEY (experiment_id) REFERENCES public.experiment(experiment_id) ON DELETE CASCADE;
 F   ALTER TABLE ONLY public.degs DROP CONSTRAINT degs_experiment_id_fkey;
       public          postgres    false    220    222    3318                       2606    24936 +   docking_results docking_results_deg_id_fkey    FK CONSTRAINT     �   ALTER TABLE ONLY public.docking_results
    ADD CONSTRAINT docking_results_deg_id_fkey FOREIGN KEY (deg_id) REFERENCES public.degs(deg_id) ON DELETE CASCADE;
 U   ALTER TABLE ONLY public.docking_results DROP CONSTRAINT docking_results_deg_id_fkey;
       public          postgres    false    230    3320    222                       2606    24931 2   docking_results docking_results_experiment_id_fkey    FK CONSTRAINT     �   ALTER TABLE ONLY public.docking_results
    ADD CONSTRAINT docking_results_experiment_id_fkey FOREIGN KEY (experiment_id) REFERENCES public.experiment(experiment_id) ON DELETE CASCADE;
 \   ALTER TABLE ONLY public.docking_results DROP CONSTRAINT docking_results_experiment_id_fkey;
       public          postgres    false    220    3318    230                       2606    16736 .   gene_aliases gene_aliases_disease_gene_id_fkey    FK CONSTRAINT     �   ALTER TABLE ONLY public.gene_aliases
    ADD CONSTRAINT gene_aliases_disease_gene_id_fkey FOREIGN KEY (disease_gene_id) REFERENCES public.disease_genes(disease_gene_id) ON DELETE CASCADE;
 X   ALTER TABLE ONLY public.gene_aliases DROP CONSTRAINT gene_aliases_disease_gene_id_fkey;
       public          postgres    false    216    3314    226                       2606    33130 8   molecular_dynamics molecular_dynamics_experiment_id_fkey    FK CONSTRAINT     �   ALTER TABLE ONLY public.molecular_dynamics
    ADD CONSTRAINT molecular_dynamics_experiment_id_fkey FOREIGN KEY (experiment_id) REFERENCES public.experiment(experiment_id);
 b   ALTER TABLE ONLY public.molecular_dynamics DROP CONSTRAINT molecular_dynamics_experiment_id_fkey;
       public          postgres    false    233    220    3318                       2606    33135 5   molecular_dynamics molecular_dynamics_pubchem_id_fkey    FK CONSTRAINT     �   ALTER TABLE ONLY public.molecular_dynamics
    ADD CONSTRAINT molecular_dynamics_pubchem_id_fkey FOREIGN KEY (pubchem_id) REFERENCES public.natural_compounds(pubchem_id);
 _   ALTER TABLE ONLY public.molecular_dynamics DROP CONSTRAINT molecular_dynamics_pubchem_id_fkey;
       public          postgres    false    231    3334    233                       2606    16620 8   pathway_enrichment pathway_enrichment_experiment_id_fkey    FK CONSTRAINT     �   ALTER TABLE ONLY public.pathway_enrichment
    ADD CONSTRAINT pathway_enrichment_experiment_id_fkey FOREIGN KEY (experiment_id) REFERENCES public.experiment(experiment_id) ON DELETE CASCADE;
 b   ALTER TABLE ONLY public.pathway_enrichment DROP CONSTRAINT pathway_enrichment_experiment_id_fkey;
       public          postgres    false    3318    220    224            	           2606    16682 %   ppi_network ppi_network_experiment_fk    FK CONSTRAINT     �   ALTER TABLE ONLY public.ppi_network
    ADD CONSTRAINT ppi_network_experiment_fk FOREIGN KEY (experiment_id) REFERENCES public.experiment(experiment_id);
 O   ALTER TABLE ONLY public.ppi_network DROP CONSTRAINT ppi_network_experiment_fk;
       public          postgres    false    3318    220    218                       2606    24974 3   therapeutic_targets therapeutic_targets_deg_id_fkey    FK CONSTRAINT     �   ALTER TABLE ONLY public.therapeutic_targets
    ADD CONSTRAINT therapeutic_targets_deg_id_fkey FOREIGN KEY (deg_id) REFERENCES public.degs(deg_id) ON DELETE CASCADE;
 ]   ALTER TABLE ONLY public.therapeutic_targets DROP CONSTRAINT therapeutic_targets_deg_id_fkey;
       public          postgres    false    222    228    3320                       2606    16820 <   therapeutic_targets therapeutic_targets_disease_gene_id_fkey    FK CONSTRAINT     �   ALTER TABLE ONLY public.therapeutic_targets
    ADD CONSTRAINT therapeutic_targets_disease_gene_id_fkey FOREIGN KEY (disease_gene_id) REFERENCES public.disease_genes(disease_gene_id) ON DELETE CASCADE;
 f   ALTER TABLE ONLY public.therapeutic_targets DROP CONSTRAINT therapeutic_targets_disease_gene_id_fkey;
       public          postgres    false    3314    228    216                       2606    16815 :   therapeutic_targets therapeutic_targets_experiment_id_fkey    FK CONSTRAINT     �   ALTER TABLE ONLY public.therapeutic_targets
    ADD CONSTRAINT therapeutic_targets_experiment_id_fkey FOREIGN KEY (experiment_id) REFERENCES public.experiment(experiment_id) ON DELETE CASCADE;
 d   ALTER TABLE ONLY public.therapeutic_targets DROP CONSTRAINT therapeutic_targets_experiment_id_fkey;
       public          postgres    false    3318    220    228           
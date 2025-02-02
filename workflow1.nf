nextflow.enable.dsl = 2

workflow {
    // Create an experiment and capture its ID
    def experiment_id_channel = create_experiment(params.geo_id,
                                                  params.compound_name,
                                                  params.description,
                                                  params.raw_p,
                                                  params.adj_p,
                                                  params.log_fc_up,
                                                  params.log_fc_down)

    // Run analyze_geo and wait for its completion
    def geo_analysis_done = analyze_geo(params.geo_id,
                                        params.samples,
                                        params.groups,
                                        experiment_id_channel,
                                        params.ligand_cid,
                                        params.raw_p,
                                        params.adj_p,
                                        params.log_fc_up,
                                        params.log_fc_down)

    def deg_match_done = match_deg_and_disease_genes(experiment_id_channel, geo_analysis_done)

    // Run only after genes were matched, otherwise tt table wont be filled
    build_ppi_network(experiment_id_channel, deg_match_done)
    perform_pathway_enrichment(params.db_connection_string, experiment_id_channel, deg_match_done)

}


process create_experiment {
    publishDir 'data'
    input:
    val geo_id
    val compound_name
    val description
    val raw_p
    val adj_p
    val log_fc_up
    val log_fc_down

    output:
    path "experiment_id.txt"


    script:
    """
    # Debug: Log current step
    echo "Step: Starting psql command" > create_experiment_debug.log

    # Run the command
    PGPASSWORD=admin psql -h db -p 5432 -U postgres -d adhd_research -t -A -c \\
    "INSERT INTO experiment (geo_id, compound, raw_p, adj_p, log_fc_up, log_fc_down, description, status) VALUES ('$geo_id', '$compound_name', '$raw_p', '$adj_p', '$log_fc_up', '$log_fc_down', '$description', 'Running') RETURNING experiment_id;" | grep -Eo '^[0-9]+' > experiment_id.txt
    """
}


process analyze_geo {
    publishDir 'data'

    input:
    val geo_id
    val samples
    val groups
    val experiment_id
    val ligand_cid
    val raw_p
    val adj_p
    val log_fc_up
    val log_fc_down

    output:
    path "analyze_geo_done.txt"

    script:
    """
    Rscript ${params.scripts_dir}/analyze_geo.R \\
        --geo_id $geo_id --samples "$samples" \\
        --experiment_id $experiment_id \\
        --groups $groups \\
        --db_connection_string "${params.db_connection_string}" \\
        --ligand_cid $ligand_cid \\
        --raw_p $raw_p \\
        --adj_p $adj_p \\
        --log_fc_up $log_fc_up \\
        --log_fc_down $log_fc_down \\
        --results_dir "${params.results_dir}"
    echo "Analyze GEO completed." > analyze_geo_done.txt
    """
}


process match_deg_and_disease_genes {

    input:
    val experiment_id
    val geo_analysis_done

    script:
    """
    python3 ${params.scripts_dir}/match_deg_and_disease_genes.py \\
        --experiment_id $experiment_id \\
        --db_connection_string "${params.db_connection_string}"
    """

}

process build_ppi_network {
    publishDir 'data'

    input:
    val experiment_id
    path geo_analysis_done

    output:
    path "ppi_done.txt"

    script:
    """
    python3 ${params.scripts_dir}/build_ppi_network.py \\
        --experiment_id $experiment_id \\
        --db_connection_string "${params.db_connection_string}" \\
        --base_dir "${params.project_dir}"
    echo "PPI network completed." > ppi_done.txt
    """
}

process perform_pathway_enrichment {
    publishDir 'data'

    input:
    val db_connection_string
    val experiment_id
    path ppi_done

    script:
    """
    python3 ${params.scripts_dir}/pathway_enrichment.py \\
        --experiment_id $experiment_id \\
        --db_connection_string $db_connection_string \\
        --base_dir "${params.project_dir}"
    """
}

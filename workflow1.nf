nextflow.enable.dsl = 2
//stop hereeeeeeeeee
workflow {
    // Create an experiment and capture its ID
    create_experiment(params.geo_id, params.compound_name, params.description)

    // Channel to pass the experiment ID
    def experiment_id_channel = Channel.fromPath("results/experiment_id.txt")

    // Run analyze_geo and wait for its completion
    def geo_analysis_done = analyze_geo(params.geo_id, params.samples, params.groups, experiment_id_channel)

    match_deg_and_disease_genes(experiment_id_channel, geo_analysis_done)

    // Run build_ppi_network only after analyze_geo completes
    def ppi_done = build_ppi_network(experiment_id_channel, geo_analysis_done)

    // Run perform_pathway_enrichment only after build_ppi_network completes
    perform_pathway_enrichment(params.db_connection_string, experiment_id_channel, geo_analysis_done)

    // Run molecular_docking

}


process create_experiment {
    publishDir 'results'
    input:
    val geo_id
    val compound_name
    val description

    output:
    path "experiment_id.txt"


    script:
    """
    # Debug: Log current step
    echo "Step: Starting psql command" > create_experiment_debug.log

    # Run the command
    psql ${params.db_connection_string} -t -A -c \\
    "INSERT INTO experiment (geo_id, compound, description, status) VALUES ('$geo_id', '$compound_name', '$description', 'Running') RETURNING experiment_id;" | grep -Eo '^[0-9]+' > experiment_id.txt
    """

    // Pass the value via the output file
    // experiment_id = file("experiment_id.txt")
}


process analyze_geo {
    publishDir 'results'

    input:
    val geo_id
    val samples
    val groups
    val experiment_id

    output:
    path "analyze_geo_done.txt"

    script:
    """
    Rscript /home/scmbag/Desktop/ADHD_Compound_Pipeline/scripts/analyze_geo.R --geo_id $geo_id --samples "$samples" --experiment_id $experiment_id --groups $groups --db_connection_string "${params.db_connection_string}"
    echo "Analyze GEO completed." > analyze_geo_done.txt
    """
}

process match_deg_and_disease_genes {

    input:
    val experiment_id
    val geo_analysis_done

    script:
    """
    python3 /home/scmbag/Desktop/ADHD_Compound_Pipeline/scripts/match_deg_and_disease_genes.py --experiment_id $experiment_id --db_connection_string "${params.db_connection_string}"
    """

}

process build_ppi_network {
    publishDir 'results'

    input:
    val experiment_id
    path geo_analysis_done

    output:
    path "ppi_done.txt"

    script:
    """
    python3 /home/scmbag/Desktop/ADHD_Compound_Pipeline/scripts/build_ppi_network.py --experiment_id $experiment_id --db_connection_string "${params.db_connection_string}"
    echo "PPI network completed." > ppi_done.txt
    """
}

process perform_pathway_enrichment {
    publishDir 'results'

    input:
    val db_connection_string
    val experiment_id
    path ppi_done

    script:
    """
    python3 /home/scmbag/Desktop/ADHD_Compound_Pipeline/scripts/pathway_enrichment.py --experiment_id $experiment_id --db_connection_string $db_connection_string
    """
}


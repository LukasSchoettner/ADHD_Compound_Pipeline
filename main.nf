nextflow.enable.dsl = 2

workflow {
    def experiment_id_channel = create_experiment(params.geo_id, params.compound, params.description)

    def scripts_channel = Channel.fromPath("scripts")
    analyze_geo(params.geo_id, params.samples, scripts_channel)
}

// Define the individual processes

process create_experiment {
    input:
    val geo_id
    val compound
    val description

    output:
    path "experiment_id.txt"

    script:
    """
    # Debug: Log current step
    echo "Step: Starting psql command" > create_experiment_debug.log

    # Run the command
    psql ${params.db_connection_string} -t -A -c \\
    "INSERT INTO experiment (geo_id, compound, description, status) VALUES ('$geo_id', '$compound', '$description', 'Running') RETURNING experiment_id;" | grep -Eo '^[0-9]+' > experiment_id.txt

    # Debug: Check if file is created
    if [ -f experiment_id.txt ]; then
        echo "Step: File created successfully" >> create_experiment_debug.log
    else
        echo "Step: File creation failed" >> create_experiment_debug.log
        exit 1
    fi
    """
}

process analyze_geo {
    input:
    val geo_id
    val samples
    path scripts

    output:
    path "deg_results_final.csv"
    path "volcano_plot_final.png"

    script:
    """
    Rscript scripts/analyze_geo.R --geo_id $geo_id --samples "$samples" --db_connection_string "${params.db_connection_string}"
    """
}

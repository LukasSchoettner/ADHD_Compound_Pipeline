nextflow.enable.dsl = 2

workflow {
    create_experiment(params.geo_id, params.compound, params.description)

    def scripts_channel = Channel.fromPath("scripts")
    def experiment_id = Channel.fromPath("results/experiment_id.txt")
    analyze_geo(params.geo_id, params.samples, experiment_id)
}


process create_experiment {
    publishDir 'results'
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
    """

    // Pass the value via the output file
    // experiment_id = file("experiment_id.txt")
}


process analyze_geo {
    publishDir 'results'
    input:
    val geo_id
    val samples
    val experiment_id

    script:
    """
    Rscript /home/scmbag/Desktop/ADHD_Compound_Pipeline/scripts/analyze_geo.R --geo_id $geo_id --samples "$samples" --experiment_id $experiment_id --db_connection_string "${params.db_connection_string}"
    """
}

nextflow.enable.dsl = 2

workflow {
    // Load docking parameters
    def docking_params = Channel.fromPath(params.docking_params)

    // Load experiment ID
    def experiment_id_channel = Channel.fromPath("data/experiment_id.txt")

    // Perform molecular docking
    molecular_docking(params.ligand_cid, params.db_connection_string, docking_params, experiment_id_channel)
}

process molecular_docking {
    publishDir 'results'

    input:
    val ligand_cid
    val db_connection_string
    path docking_params
    path experiment_id

    output:
    path "docking_done.txt"

    script:
    """
    python3 ${params.scripts_dir}/molecular_docking.py \\
        --ligand_cid $ligand_cid \\
        --db_connection_string $db_connection_string \\
        --experiment_id \$(cat $experiment_id) \\
        --docking_params $docking_params \\
        --output_dir ${params.docking_results_dir} \\
        --project_dir ${params.project_dir}
    echo "Docking completed." > docking_done.txt
    """
}

process molecular_dynamics {
  input:
    val db_connection
    val experiment_id
    val input_structure

  """
  python ${params.scripts_dir}/molecular_dynamics.py \
    --db_connection_string $db_connection \
    --experiment_id $experiment_id \
    --input_structure $input_structure \
    --simulation_time 100 \
    --temperature 310 \
    --pressure 1.0 \
    --solvent_type water \
    --output_dir md_results
  """
}


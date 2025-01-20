nextflow.enable.dsl = 2

workflow {
    // Load MD parameters
    def md_params = Channel.fromPath(params.md_param_path)

    // Load experiment ID
    def experiment_id_channel = Channel.fromPath("data/experiment_id.txt")

    // Perform molecular docking
    molecular_dynamics(params.db_connection_string, md_params, experiment_id_channel)
}

process molecular_dynamics {
  input:
    val db_connection
    val experiment_id
    val input_structure
    path md_params

  """
  python ${params.scripts_dir}/molecular_dynamics.py \\
    --db_connection_string $db_connection \\
    --experiment_id $experiment_id \\
    --md_param_file $md_params
    --output_dir md_results_path
  """
}


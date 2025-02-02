nextflow.enable.dsl = 2

workflow {
    // Load MD parameters
    def md_param_path = Channel.fromPath(params.md_param_path)

    // Load experiment ID
    def experiment_id_channel = Channel.fromPath("data/experiment_id.txt")

    // Perform molecular docking
    molecular_dynamics(params.db_connection_string, experiment_id_channel, md_param_path, params.output_dir, params.data_dir)
}

process molecular_dynamics {
  input:
    val db_connection
    val experiment_id
    path md_param_path
    path output_dir
    path data_dir

  """
  python ${params.scripts_dir}/molecular_dynamics.py \\
    --db_connection_string $db_connection \\
    --experiment_id $experiment_id \\
    --md_param_file $md_param_path \\
    --output_dir $output_dir \\
    --data_dir $data_dir \\
    --project_dir ${params.project_dir}
  """
}

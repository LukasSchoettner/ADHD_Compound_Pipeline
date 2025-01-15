nextflow.enable.dsl = 2

workflow {
    // Load docking parameters
    def docking_params = Channel.fromPath(params.docking_params)

    // Load experiment ID
    def experiment_id_channel = Channel.fromPath("results/experiment_id.txt")

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
    python3 /home/scmbag/Desktop/ADHD_Compound_Pipeline/scripts/molecular_docking.py \\
        --ligand_cid $ligand_cid \\
        --db_connection_string $db_connection_string \\
        --experiment_id \$(cat $experiment_id) \\
        --docking_params $docking_params \\
        --output_dir "/home/scmbag/Desktop/ADHD_Compound_Pipeline/results/molecular_docking"
    echo "Docking completed." > docking_done.txt
    """
}

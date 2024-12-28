nextflow.enable.dsl = 2

workflow {
    // Step 1: Create the experiment
    create_experiment(params.geo_id, params.compound, params.description)

    // Step 2: Perform GEO analysis
    analyze_geo(params.geo_id, params.samples, "experiment_id.txt")

    // Step 3: Build PPI network
    build_ppi_network(params.geo_id, "experiment_id.txt")

    // Step 4: Collect ADHD-related genes
    collect_adhd_genes(params.api_url, "experiment_id.txt")

    // Step 5: Match DEGs and ADHD genes
    match_deg_and_disease_genes("DEG_results_final.csv", "adhd_genes.json", "experiment_id.txt")

    // Step 6: Molecular docking
    molecular_docking(params.compound, "experiment_id.txt")

    // Step 7: Molecular dynamics
    molecular_dynamics("docking_results.csv", "experiment_id.txt")

    // Step 8: Pathway enrichment
    pathway_enrichment("matched_genes.csv", "experiment_id.txt")

    // Step 9: Synergy analysis
    synergy_analysis("pathway_enrichment_results.csv", "experiment_id.txt")
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
    psql -c "INSERT INTO experiment (geo_id, compound, description, status) VALUES ('$geo_id', '$compound', '$description', 'Running') RETURNING experiment_id" -t > experiment_id.txt
    """
}

process analyze_geo {
    input:
    val geo_id
    val samples
    path "experiment_id.txt"

    output:
    path "DEG_results_final.csv"
    path "volcano_plot_final.png"

    script:
    """
    experiment_id=\$(cat experiment_id.txt)
    Rscript scripts/analyze_geo.R --geo_id $geo_id --samples "$samples" --db_connection_string "${params.db_connection_string}" --experiment_id \$experiment_id
    """
}

process build_ppi_network {
    input:
    val geo_id
    path "experiment_id.txt"

    output:
    path "ppi_network.csv"

    script:
    """
    experiment_id=\$(cat experiment_id.txt)
    python scripts/build_ppi_network.py --geo_id $geo_id --db_connection_string "${params.db_connection_string}" --experiment_id \$experiment_id
    """
}

process collect_adhd_genes {
    input:
    val api_url
    path "experiment_id.txt"

    output:
    path "adhd_genes.json"

    script:
    """
    experiment_id=\$(cat experiment_id.txt)
    python scripts/collect_adhd_genes.py --api_url $api_url --experiment_id \$experiment_id
    """
}

process match_deg_and_disease_genes {
    input:
    path "DEG_results_final.csv"
    path "adhd_genes.json"
    path "experiment_id.txt"

    output:
    path "matched_genes.csv"

    script:
    """
    experiment_id=\$(cat experiment_id.txt)
    python scripts/match_deg_and_disease_genes.py --degs DEG_results_final.csv --adhd_genes adhd_genes.json --db_connection_string "${params.db_connection_string}" --experiment_id \$experiment_id
    """
}

process molecular_docking {
    input:
    val compound
    path "experiment_id.txt"

    output:
    path "docking_results.csv"

    script:
    """
    experiment_id=\$(cat experiment_id.txt)
    bash scripts/molecular_docking.sh --compound $compound --db_connection_string "${params.db_connection_string}" --experiment_id \$experiment_id
    """
}

process molecular_dynamics {
    input:
    path "docking_results.csv"
    path "experiment_id.txt"

    output:
    path "dynamics_simulation.log"

    script:
    """
    experiment_id=\$(cat experiment_id.txt)
    bash scripts/molecular_dynamics.sh --input docking_results.csv --db_connection_string "${params.db_connection_string}" --experiment_id \$experiment_id
    """
}

process pathway_enrichment {
    input:
    path "matched_genes.csv"
    path "experiment_id.txt"

    output:
    path "pathway_enrichment_results.csv"

    script:
    """
    experiment_id=\$(cat experiment_id.txt)
    Rscript scripts/pathway_enrichment.R --input matched_genes.csv --db_connection_string "${params.db_connection_string}" --experiment_id \$experiment_id
    """
}

process synergy_analysis {
    input:
    path "pathway_enrichment_results.csv"
    path "experiment_id.txt"

    output:
    path "synergy_results.csv"

    script:
    """
    experiment_id=\$(cat experiment_id.txt)
    Rscript scripts/synergy_analysis.R --input pathway_enrichment_results.csv --db_connection_string "${params.db_connection_string}" --experiment_id \$experiment_id
    """
}

nextflow.enable.dsl = 2

workflow {
    // Define workflow steps
    analyze_geo(params.geo_id, params.samples)
    build_ppi_network(params.geo_id)
    collect_adhd_genes(params.api_url)
    match_deg_and_disease_genes()
    molecular_docking(params.compound)
    molecular_dynamics()
    pathway_enrichment()
    synergy_analysis()
}

// Define the individual processes

process analyze_geo {
    input:
    val geo_id
    val samples

    output:
    path "DEG_results_final.csv"  // Output from R script
    path "volcano_plot_final.png" // Output plot

    script:
    """
    Rscript scripts/analyze_geo.R --geo_id $geo_id --samples "$samples"
    """
}

process build_ppi_network {
    input:
    val geo_id

    output:
    path "ppi_network.csv"  // Example output file

    script:
    """
    python scripts/build_ppi_network.py --geo_id $geo_id
    """
}

process collect_adhd_genes {
    input:
    val api_url

    output:
    path "adhd_genes.json"

    script:
    """
    python scripts/collect_adhd_genes.py --api_url $api_url
    """
}

process match_deg_and_disease_genes {
    input:
    path "DEG_results_final.csv"  // Input from analyze_geo process
    path "adhd_genes.json"        // Input from collect_adhd_genes process

    output:
    path "matched_genes.csv"

    script:
    """
    python scripts/match_deg_and_disease_genes.py --degs DEG_results_final.csv --adhd_genes adhd_genes.json
    """
}

process molecular_docking {
    input:
    val compound

    output:
    path "docking_results.csv"

    script:
    """
    bash scripts/molecular_docking.sh --compound $compound
    """
}

process molecular_dynamics {
    input:
    path "docking_results.csv"

    output:
    path "dynamics_simulation.log"

    publishDir "${params.results_dir}/molecular_dynamics", mode: 'copy'

    script:
    """
    bash scripts/molecular_dynamics.sh --input docking_results.csv
    """
}


process pathway_enrichment {
    input:
    path "matched_genes.csv"  // Input from match_deg_and_disease_genes process

    output:
    path "pathway_enrichment_results.csv"

    script:
    """
    Rscript scripts/pathway_enrichment.R --input matched_genes.csv
    """
}

process synergy_analysis {
    input:
    path "pathway_enrichment_results.csv"  // Input from pathway_enrichment process

    output:
    path "synergy_results.csv"

    script:
    """
    Rscript scripts/synergy_analysis.R --input pathway_enrichment_results.csv
    """
}

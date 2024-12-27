nextflow.enable.dsl = 2

workflow {
    // Define inputs
    params.geo_id = "GSE85871"
    params.samples = "GSM2286316,GSM2286317,GSM2286238,GSM2286239"
    params.compound = "Gastrodin"

    // Define workflow steps
    analyze_geo(params)
    build_ppi_network(params)
    collect_adhd_genes(params)
    match_deg_and_disease_genes(params)
    molecular_docking(params)
    molecular_dynamics(params)
    pathway_enrichment(params)
    synergy_analysis(params)
}

// Define the individual processes

process analyze_geo {
    input:
    val geo_id from params.geo_id
    val samples from params.samples

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
    val geo_id from params.geo_id

    script:
    """
    python scripts/build_ppi_network.py --geo_id $geo_id
    """
}

process collect_adhd_genes {
    script:
    """
    python scripts/collect_adhd_genes.py
    """
}

process match_deg_and_disease_genes {
    script:
    """
    python scripts/match_deg_and_disease_genes.py
    """
}

process molecular_docking {
    input:
    val compound from params.compound

    script:
    """
    bash scripts/molecular_docking.sh --compound $compound
    """
}

process molecular_dynamics {
    script:
    """
    bash scripts/molecular_dynamics.sh
    """
}

process pathway_enrichment {
    script:
    """
    Rscript scripts/pathway_enrichment.R
    """
}

process synergy_analysis {
    script:
    """
    Rscript scripts/synergy_analysis.R
    """
}

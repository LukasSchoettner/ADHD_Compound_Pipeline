#!/usr/bin/env Rscript

###############################################################################
# Load required libraries
###############################################################################
library(argparse)
library(GEOquery)
library(limma)
library(Biobase)
library(ggplot2)
library(AnnotationDbi)
library(dplyr)
library(DBI)
library(RPostgres)
library(tibble)         # for rownames_to_column

###############################################################################
# Helper Functions
###############################################################################

# 1) Load Annotation Package
load_annotation_package <- function(platform_id) {
  platform_package_map <- list(
    "GPL570" = "hgu133plus2.db",
    "GPL96"  = "hgu133a.db",
    "GPL571" = "hgu133a2.db"
  )
  if (platform_id %in% names(platform_package_map)) {
    annotation_package <- platform_package_map[[platform_id]]
    cat("Using annotation package:", annotation_package, "\n")

    if (!require(annotation_package, character.only = TRUE)) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install(annotation_package)
      library(annotation_package, character.only = TRUE)
    }
    return(annotation_package)
  } else {
    stop("Unsupported platform or no annotation package specified for platform: ", platform_id)
  }
}

# 2) Fetch annotation data from GEO if no local annotation package is available
get_annotation_from_geo <- function(platform_id) {
  cat("Fetching annotation data from GEO for platform:", platform_id, "\n")
  gpl <- getGEO(platform_id, AnnotGPL = TRUE)
  platform_table <- Table(gpl)
  # Check that 'Gene Symbol' column is present
  if (!all(c("ID", "Gene Symbol") %in% colnames(platform_table))) {
    stop("The platform table does not contain 'ID' and 'Gene Symbol' columns.")
  }
  probe_to_gene <- platform_table[, c("ID", "Gene Symbol")]
  colnames(probe_to_gene) <- c("PROBEID", "SYMBOL")
  return(probe_to_gene)
}

# 3) Map probes to gene symbols
map_probes_to_symbols <- function(expr_data, platform_id) {
  # For Affymetrix arrays with known annotation packages
  if (platform_id %in% c("GPL570", "GPL96", "GPL571")) {
    annotation_package <- load_annotation_package(platform_id)
    probe_ids <- rownames(expr_data)
    gene_symbols <- mapIds(
      get(annotation_package),
      keys = probe_ids,
      column = "SYMBOL",
      keytype = "PROBEID",
      multiVals = "first"
    )
    rownames(expr_data) <- gene_symbols

  # For other platforms, try pulling from GEO's own table
  } else {
    probe_to_gene <- get_annotation_from_geo(platform_id)
    expr_data_df <- as.data.frame(expr_data)
    expr_data_df$PROBEID <- rownames(expr_data_df)
    merged_df <- merge(expr_data_df, probe_to_gene, by = "PROBEID")
    gene_symbols <- merged_df$SYMBOL
    # Drop columns we no longer need
    merged_df <- merged_df[, !(names(merged_df) %in% c("PROBEID", "SYMBOL"))]
    expr_data <- as.matrix(merged_df)
    rownames(expr_data) <- gene_symbols
  }

  # Remove any rows with NA gene symbols
  expr_data <- expr_data[!is.na(rownames(expr_data)), ]

  # If multiple probes map to the same gene symbol, we summarize them by mean
  # (You could also use median or max, depending on your preference.)
  expr_data <- expr_data %>%
    as.data.frame() %>%
    rownames_to_column("gene_symbol") %>%
    group_by(gene_symbol) %>%
    summarize(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
    as.data.frame()

  rownames(expr_data) <- expr_data$gene_symbol
  expr_data <- expr_data[, !colnames(expr_data) %in% "gene_symbol"]

  return(as.matrix(expr_data))
}

# 4) Run the limma differential expression pipeline
run_differential_expression <- function(expr_data, pheno_data) {
  # Create design matrix: one column per group factor
  design <- model.matrix(~ 0 + factor(pheno_data$group))
  group_levels <- levels(factor(pheno_data$group))
  colnames(design) <- group_levels

  cat("Design matrix:\n")
  print(design)

  # We'll define a single contrast: "Compound - Control"
  # Adjust as needed if you have more complex designs
  if (!all(c("Control", "Compound") %in% group_levels)) {
    stop("For this example, we expect 'Control' and 'Compound' in the groups, but got:", paste(group_levels, collapse = ", "))
  }

  contrast_matrix <- makeContrasts(
    Compound_vs_Control = Compound - Control,
    levels = design
  )

  cat("Contrast matrix:\n")
  print(contrast_matrix)

  # Fit the model
  fit <- lmFit(expr_data, design)
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)

  # Return the entire fit2 so the caller can do topTable
  return(fit2)
}

###############################################################################
# Main Script
###############################################################################

# Create argument parser
parser <- ArgumentParser(description = "Script for GEO data analysis (limma-based)")
parser$add_argument("--geo_id", type = "character", required = TRUE,
                    help = "GEO ID of the dataset (e.g. GSE85871)")
parser$add_argument("--samples", type = "character", required = TRUE,
                    help = "Comma-separated sample IDs (e.g. GSM1234,GSM1235)")
parser$add_argument("--groups", type = "character", required = TRUE,
                    help = "Comma-separated group labels, matching the number of samples")
parser$add_argument("--db_connection_string", type = "character", required = TRUE,
                    help = "Database connection string")
parser$add_argument("--experiment_id", type = "character", required = TRUE,
                    help = "Experiment ID or path to .txt file containing it")
args <- parser$parse_args()

cat("Parsed Arguments:\n")
cat(" GEO ID:", args$geo_id, "\n")
cat(" Samples:", args$samples, "\n")
cat(" experiment_id:", args$experiment_id, "\n")
cat(" groups:", args$groups, "\n")

# 1) Determine the final experiment_id integer
if (file.exists(args$experiment_id)) {
  experiment_id <- as.integer(readLines(args$experiment_id))
  if (is.na(experiment_id)) {
    stop("Error: experiment_id from file could not be parsed as an integer.")
  }
} else {
  experiment_id <- as.integer(args$experiment_id)
  if (is.na(experiment_id)) {
    stop("Error: experiment_id must be an integer or a valid .txt file.")
  }
}
cat("Final experiment_id:", experiment_id, "\n")

# 2) Download the GEO dataset
gse_list <- getGEO(args$geo_id, GSEMatrix = TRUE)
if (length(gse_list) > 1) {
  gse <- gse_list[[1]]
} else {
  gse <- gse_list[[1]]
}

# 3) Determine platform (GPL)
if ("annotation" %in% slotNames(gse)) {
  platform_id <- annotation(gse)
} else if ("platform_id" %in% names(Meta(gse))) {
  platform_id <- Meta(gse)$platform_id
} else {
  stop("Unable to determine platform info from the GEO dataset.")
}
cat("Platform ID:", platform_id, "\n")

# Extract expression and phenotype data
expr_data <- exprs(gse)
pheno_data <- pData(gse)

# 4) Subset to user-specified samples
samples_vec <- unlist(strsplit(args$samples, ","))
if (!all(samples_vec %in% colnames(expr_data))) {
  missing_samps <- samples_vec[!samples_vec %in% colnames(expr_data)]
  stop("These sample(s) do not exist in the dataset:", paste(missing_samps, collapse = ", "))
}
expr_data <- expr_data[, samples_vec]
pheno_data <- pheno_data[samples_vec, , drop = FALSE]

# 5) Assign user-specified groups
groups_vec <- unlist(strsplit(args$groups, ","))
if (length(groups_vec) != length(samples_vec)) {
  stop("Number of groups must match number of samples! Provided groups: ", length(groups_vec),
       ", samples: ", length(samples_vec))
}
pheno_data$group <- groups_vec

# Check length match
if (length(samples_vec) != length(groups_vec)) {
  stop("Number of samples does not match number of group labels!")
}

# 6) Map probes to gene symbols
expr_data <- map_probes_to_symbols(expr_data, platform_id)

# 7) Run the differential expression pipeline with limma
fit2 <- run_differential_expression(expr_data, pheno_data)
deg_results <- topTable(fit2, adjust.method = "none", number = Inf)
deg_results$adj.P.Value <- p.adjust(deg_results$P.Value, method = "BH")

# Mark significance (simple example)
deg_results$significant <- "Not changed"
deg_results$significant[deg_results$P.Value <= 0.05 & deg_results$logFC > 1]  <- "Up regulated"
deg_results$significant[deg_results$P.Value <= 0.05 & deg_results$logFC < -1] <- "Down regulated"

deg_results_filtered <- deg_results[deg_results$significant %in% c("Up regulated", "Down regulated"), ]

cat("Number of DEGs passing filters:", nrow(deg_results_filtered), "\n")

# 8) Save DEGs to DB (if any)
if (nrow(deg_results_filtered) > 0) {
  degs_to_save <- data.frame(
    gene_name       = rownames(deg_results_filtered),
    log_fold_change = deg_results_filtered$logFC,
    p_value         = deg_results_filtered$P.Value,
    degs            = deg_results_filtered$significant,
    experiment_id   = experiment_id
  )
  rownames(degs_to_save) <- NULL

  # Connect to DB (keeping credentials inline as requested)
  con <- dbConnect(
    RPostgres::Postgres(),
    host = "/run/postgresql",
    port = 5432,
    dbname = "adhd_research",
    user = "postgres",
    password = "admin"
  )

  # Optionally wrap in a transaction
  dbBegin(con)
  tryCatch({
    dbWriteTable(con, name = "degs", value = degs_to_save, row.names = FALSE, append = TRUE)
    dbCommit(con)
  }, error = function(e) {
    dbRollback(con)
    stop("DB insertion failed: ", e$message)
  })

  dbDisconnect(con)
} else {
  cat("No DEGs found to save.\n")
}

# 9) Volcano plot
if (nrow(deg_results) > 0) {
  volcano_plot <- ggplot(deg_results,
                         aes(x = logFC, y = -log10(adj.P.Value), color = significant)) +
    geom_point(alpha = 0.8) +
    scale_color_manual(values = c("Down regulated" = "blue",
                                  "Not changed" = "gray",
                                  "Up regulated" = "red")) +
    labs(title = "Volcano plot", x = "Log2 (Fold change)", y = "-Log10 (Adjusted P-value)") +
    theme_minimal()

  ggsave("volcano_plot_final.png", volcano_plot, width = 8, height = 6)
  cat("Volcano plot saved: volcano_plot_final.png\n")
} else {
  cat("No data for volcano plot.\n")
}

cat("Analysis completed successfully.\n")

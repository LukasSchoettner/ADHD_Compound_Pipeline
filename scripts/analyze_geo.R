#!/usr/bin/env Rscript

source("/app/renv/activate.R")
.libPaths(c("/opt/conda/envs/myenv/lib/R/library", .libPaths()))
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
library(tibble)


# 1) Create an annotation data frame for your probes

get_annotation_df <- function(probe_ids, platform_id) {
  # returns a data.frame with:
  #   probe_id    gene_symbol
  # for each probe in `probe_ids`.

  # If the platform has a local annotation package (e.g. GPL570)
  if (platform_id %in% c("GPL570","GPL96","GPL571")) {
    # Load or install the package
    platform_package <- list(
      "GPL570" = "hgu133plus2.db",
      "GPL96"  = "hgu133a.db",
      "GPL571" = "hgu133a2.db"
    )[[platform_id]]

    if (!require(platform_package, character.only=TRUE)) {
      if (!requireNamespace("BiocManager", quietly=TRUE))
        install.packages("BiocManager")
      BiocManager::install(platform_package)
      library(platform_package, character.only=TRUE)
    }

    # Lookup
    gene_syms <- mapIds(
      get(platform_package),
      keys = probe_ids,
      column = "SYMBOL",
      keytype = "PROBEID",
      multiVals = "first"
    )

    # If we fail to find a gene symbol, it stays NA.
    anno_df <- data.frame(
      probe_id = probe_ids,
      gene_symbol = unname(gene_syms),
      stringsAsFactors = FALSE
    )

  } else {
    # Otherwise fetch annotation from GEO for the entire platform
    gpl <- getGEO(platform_id, AnnotGPL=TRUE)
    gpl_table <- Table(gpl)

    # We assume 'ID' / 'Gene Symbol' are the relevant columns:
    if (!all(c("ID","Gene Symbol") %in% colnames(gpl_table))) {
      stop("Cannot find 'ID' and 'Gene Symbol' in the platform table.")
    }

    # Subset to the probes we actually have
    subset_table <- gpl_table[gpl_table$ID %in% probe_ids, c("ID","Gene Symbol")]
    colnames(subset_table) <- c("probe_id","gene_symbol")

    # Some might be missing -> remain NA
    anno_df <- data.frame(
      probe_id   = probe_ids,
      gene_symbol = NA_character_,
      stringsAsFactors = FALSE
    )
    # Merge in the known symbols from subset_table
    anno_df <- merge(anno_df, subset_table,
                     by="probe_id",
                     all.x=TRUE, sort=FALSE)
    # We end up with columns: probe_id, gene_symbol.x, gene_symbol.y
    # So pick final gene_symbol
    # If gene_symbol.y is not NA, use that; else gene_symbol.x
    anno_df$gene_symbol <- ifelse(
      is.na(anno_df$`gene_symbol.y`) | anno_df$`gene_symbol.y`=="",
      anno_df$`gene_symbol.x`,
      anno_df$`gene_symbol.y`
    )
    # Clean up columns
    anno_df <- anno_df[, c("probe_id","gene_symbol")]
  }

  # If any gene_symbol is empty, keep it NA
  anno_df$gene_symbol[anno_df$gene_symbol==""] <- NA
  return(anno_df)
}


# 2) Run the limma pipeline (Contrast = Compound - Control)

run_limma <- function(expr_data, pheno_data) {
  design <- model.matrix(~ 0 + factor(pheno_data$group))
  group_levels <- levels(factor(pheno_data$group))

  # Reorder columns so "Control" is first, "Compound" is second
  if (all(c("Control","Compound") %in% group_levels)) {
    design <- design[, c("factor(pheno_data$group)Control",
                         "factor(pheno_data$group)Compound")]
    colnames(design) <- c("Control","Compound")
  } else {
    stop("Expected 'Control' and 'Compound' in the sample groups, got: ",
         paste(group_levels,collapse=", "))
  }

  contrast_matrix <- makeContrasts(
    Compound_vs_Control = Compound - Control,
    levels = design
  )

  # Fit
  fit <- lmFit(expr_data, design)
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  return(fit2)
}


###############################################################################
# Main Script
###############################################################################
parser <- ArgumentParser(description="Option C: Best-probe-per-gene after limma")
parser$add_argument("--geo_id", type="character", required=TRUE,
                    help="GEO ID of the dataset (e.g. GSE85871)")
parser$add_argument("--samples", type="character", required=TRUE,
                    help="Comma-separated sample IDs (e.g. GSM1234,GSM1235)")
parser$add_argument("--groups", type="character", required=TRUE,
                    help="Comma-separated group labels, matching #samples (e.g. Compound,Control)")
parser$add_argument("--db_connection_string", type="character", required=TRUE,
                    help="Database connection string")
parser$add_argument("--experiment_id", type="character", required=TRUE,
                    help="Experiment ID (integer) or path to file containing it")
parser$add_argument("--ligand_cid", type="character", required=TRUE,
                    help="PubChem id of ligand")

# Thresholds
parser$add_argument("--raw_p",      type="character", required=TRUE,
                    help="Raw P-Value cut-off for DEG analysis")
parser$add_argument("--adj_p",      type="character", required=TRUE,
                    help="Adjusted P-Value cut-off for DEG analysis")
parser$add_argument("--log_fc_up",  type="character", required=TRUE,
                    help="Minimum logFC for up-regulation")
parser$add_argument("--log_fc_down",type="character", required=TRUE,
                    help="Maximum logFC for down-regulation (should be negative)")
parser$add_argument("--results_dir",type="character", required=TRUE,
                    help="Path to results folder")

args <- parser$parse_args()

cat("\nParsed arguments:\n")
print(args)

# 1) Determine experiment_id
if (file.exists(args$experiment_id)) {
  experiment_id <- as.integer(readLines(args$experiment_id))
  if (is.na(experiment_id)) {
    stop("Error: experiment_id from file could not be parsed as integer.")
  }
} else {
  experiment_id <- as.integer(args$experiment_id)
  if (is.na(experiment_id)) {
    stop("Error: experiment_id must be an integer or valid .txt file.")
  }
}
cat("Final experiment_id:", experiment_id, "\n")

# 2) Download the GEO dataset
gse_list <- getGEO(args$geo_id, GSEMatrix=TRUE)
if (length(gse_list) > 1) {
  gse <- gse_list[[1]]
} else {
  gse <- gse_list[[1]]
}

# 3) Determine platform
if ("annotation" %in% slotNames(gse)) {
  platform_id <- annotation(gse)
} else if ("platform_id" %in% names(Meta(gse))) {
  platform_id <- Meta(gse)$platform_id
} else {
  stop("Cannot find platform info in GEO object.")
}
cat("Platform ID:", platform_id, "\n")

expr_data <- exprs(gse)
pheno_data <- pData(gse)

# 4) Subset to user-specified samples
samples_vec <- unlist(strsplit(args$samples, ","))
if (!all(samples_vec %in% colnames(expr_data))) {
  missing_samps <- samples_vec[!samples_vec %in% colnames(expr_data)]
  stop("These samples do not exist in the dataset:",paste(missing_samps,collapse=", "))
}
expr_data <- expr_data[, samples_vec]
pheno_data <- pheno_data[samples_vec, , drop=FALSE]

# 5) Assign groups
groups_vec <- unlist(strsplit(args$groups, ","))
if (length(groups_vec)!=length(samples_vec)) {
  stop("Number of groups must match number of samples! Provided groups: ",
       length(groups_vec), ", samples: ", length(samples_vec))
}
pheno_data$group <- groups_vec

# (Optional) log2 check like GEO2R if needed:
# qx <- quantile(expr_data, c(0,0.25,0.5,0.75,0.99,1.0), na.rm=TRUE)
# LogC <- (qx[5]>100)||((qx[6]-qx[1]>50)&&(qx[2]>0))
# if (LogC) {
#   expr_data[expr_data<=0] <- NaN
#   expr_data <- log2(expr_data)
# }

# 6) Build an annotation data frame for your probes
probe_ids <- rownames(expr_data)
anno_df <- get_annotation_df(probe_ids, platform_id)

# 7) Run limma with each probe as its own row
fit2 <- run_limma(expr_data, pheno_data)

# 8) Extract all results
deg_results <- topTable(fit2, adjust.method="none", number=Inf)

# add columns for raw + adjusted p
deg_results$adj.P.Value <- p.adjust(deg_results$P.Value, method="BH")

# Convert threshold args
raw_p_val     <- as.numeric(args$raw_p)
adj_p_val     <- as.numeric(args$adj_p)
log_fc_up_val <- as.numeric(args$log_fc_up)
log_fc_down_val <- as.numeric(args$log_fc_down)

# 9) Convert rownames (which are probe IDs) into a column for merging
deg_results <- deg_results %>%
  tibble::rownames_to_column("probe_id")

# Merge in the gene symbols
# Note: some may be NA if the platform annotation didn't have them
deg_results <- left_join(deg_results, anno_df, by="probe_id")

# If gene_symbol is NA, fallback to the probe_id
deg_results$gene_symbol[ is.na(deg_results$gene_symbol) ] <- deg_results$probe_id[ is.na(deg_results$gene_symbol) ]


# 10) For each gene_symbol, keep the single "best" probe by raw p‐value

deg_results <- deg_results %>%
  group_by(gene_symbol) %>%
  slice_min(order_by = P.Value, n=1) %>%
  ungroup()

# Now we have exactly one row per gene


# 11) Mark significance using fallback logic:

deg_results$significant <- "Not changed"

# First see how many pass the adjusted threshold
adj_pass_count <- sum(
  deg_results$adj.P.Value <= adj_p_val &
    (deg_results$logFC >= log_fc_up_val | deg_results$logFC <= log_fc_down_val)
)
cat("Number of genes passing adjusted p-value threshold:", adj_pass_count, "\n")

if (adj_pass_count>0) {
  # Use adjusted p
  deg_results$significant[
    deg_results$adj.P.Value <= adj_p_val & deg_results$logFC >= log_fc_up_val
  ] <- "Up regulated"

  deg_results$significant[
    deg_results$adj.P.Value <= adj_p_val & deg_results$logFC <= log_fc_down_val
  ] <- "Down regulated"

} else {
  # Fallback to raw p
  cat("No genes pass adjusted threshold. Using raw p-value fallback.\n")
  deg_results$significant[
    deg_results$P.Value <= raw_p_val & deg_results$logFC >= log_fc_up_val
  ] <- "Up regulated"

  deg_results$significant[
    deg_results$P.Value <= raw_p_val & deg_results$logFC <= log_fc_down_val
  ] <- "Down regulated"
}

###############################################################################
# 12) Filter up/down for saving
###############################################################################
deg_results_filtered <- deg_results %>%
  filter(significant %in% c("Up regulated","Down regulated"))

cat("Number of final DEGs to save:", nrow(deg_results_filtered), "\n")

###############################################################################
# 13) Save to DB
###############################################################################
if (nrow(deg_results_filtered)>0) {
  degs_to_save <- data.frame(
    gene_name       = deg_results_filtered$gene_symbol,
    log_fold_change = deg_results_filtered$logFC,
    p_value         = deg_results_filtered$P.Value,
    degs            = deg_results_filtered$significant,
    experiment_id   = experiment_id,
    pubchem_id      = args$ligand_cid
  )

  # Connect
  con <- dbConnect(
    RPostgres::Postgres(),
    host="db",
    port=5432,
    dbname="adhd_research",
    user="postgres",
    password="admin"
  )

  dbBegin(con)
  tryCatch({
    dbWriteTable(con, "degs", degs_to_save,
                 row.names=FALSE, append=TRUE)
    dbExecute(con,
      "UPDATE experiment SET degs_found=$1 WHERE experiment_id=$2",
      params=list(nrow(degs_to_save), experiment_id)
    )
    dbCommit(con)
    cat("DEGs saved to DB. degs_found updated.\n")
  }, error=function(e) {
    dbRollback(con)
    stop("DB insertion failed: ", e$message)
  })
  dbDisconnect(con)
} else {
  cat("No DEGs found. Setting degs_found=0 in 'experiment' table.\n")
  con <- dbConnect(
    RPostgres::Postgres(),
    host="db",
    port=5432,
    dbname="adhd_research",
    user="postgres",
    password="admin"
  )
  dbExecute(con,
    "UPDATE experiment SET degs_found=0 WHERE experiment_id=$1",
    params=list(experiment_id)
  )
  dbDisconnect(con)
}

###############################################################################
# 14) Volcano Plot
###############################################################################
if (nrow(deg_results)>0) {
  # deg_results includes the "best probe" for each gene
  volcano_plot <- ggplot(deg_results,
    aes(x=logFC, y=-log10(adj.P.Value), color=significant)) +
    geom_point(alpha=0.8) +
    scale_color_manual(values=c("Up regulated"="red",
                                "Down regulated"="blue",
                                "Not changed"="gray")) +
    labs(title="Volcano plot: best-probe-per-gene",
         x="Log2(Fold Change)",
         y="-Log10(Adjusted P-value)") +
    theme_minimal()

  ggsave(file.path(args$results_dir, "volcano_plot_final.png"), volcano_plot, width=8, height=6)
  cat("Saved volcano_plot_final.png to", args$results_dir, "\n")
} else {
  cat("No data to plot.\n")
}

cat("Analysis completed successfully.\n")
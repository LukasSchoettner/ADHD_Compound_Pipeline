# Load required libraries
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

# Helper Function: Load Annotation Package
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
    stop("Unsupported platform or annotation package not specified for this platform.")
  }
}

# Helper Function: Get Annotation Data from GEO
get_annotation_from_geo <- function(platform_id) {
  cat("Fetching annotation data from GEO for platform:", platform_id, "\n")
  gpl <- getGEO(platform_id, AnnotGPL = TRUE)
  platform_table <- Table(gpl)
  probe_to_gene <- platform_table[, c("ID", "Gene Symbol")]
  colnames(probe_to_gene) <- c("PROBEID", "SYMBOL")
  return(probe_to_gene)
}

# Create argument parser
parser <- ArgumentParser(description = "Script for GEO data analysis")
parser$add_argument("--geo_id", type = "character", required = TRUE, help = "GEO ID of the dataset")
parser$add_argument("--samples", type = "character", required = TRUE, help = "Comma-separated sample IDs")
parser$add_argument("--db_connection_string", type = "character", required = TRUE, help = "Database connection string")
parser$add_argument("--experiment_id", type = "character", required = TRUE, help = "Experiment ID for database integration")
args <- parser$parse_args()

# Step 1: Download the GEO dataset
gse <- getGEO(args$geo_id, GSEMatrix = TRUE)

# Step 2: Handle cases where `gse` is a list
if (length(gse) > 1) {
  gse <- gse[[1]]  # Select the first element
} else {
  gse <- gse[[1]]
}

# Check if experiment_id is a file and read its value
if (file.exists(args$experiment_id)) {
  experiment_id <- as.integer(readLines(args$experiment_id))
  if (is.na(experiment_id)) {
    stop("Error: experiment_id could not be parsed as an integer.")
  }
} else {
  experiment_id <- as.integer(args$experiment_id)
  if (is.na(experiment_id)) {
    stop("Error: experiment_id must be an integer or a valid .txt file.")
  }
}

# Step 2: Extract platform (GPL) information
if ("annotation" %in% slotNames(gse)) {
  platform_id <- annotation(gse)
} else if ("platform_id" %in% names(Meta(gse))) {
  platform_id <- Meta(gse)$platform_id
} else {
  stop("Unable to determine platform information from the GEO dataset.")
}

cat("Platform ID:", platform_id, "\n")
str(gse)

expr_data <- exprs(gse)
pheno_data <- pData(gse)

samples <- unlist(strsplit(args$samples, split = ","))
expr_data <- expr_data[, samples]
pheno_data <- pheno_data[samples, ]

# Step 3: Map probe IDs to gene symbols
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
} else {
  probe_to_gene <- get_annotation_from_geo(platform_id)
  expr_data <- as.data.frame(expr_data)
  expr_data$PROBEID <- rownames(expr_data)
  expr_data <- merge(expr_data, probe_to_gene, by = "PROBEID")
  gene_symbols <- expr_data$SYMBOL
  expr_data <- expr_data[, !(names(expr_data) %in% c("PROBEID", "SYMBOL"))]
}

rownames(expr_data) <- gene_symbols
expr_data <- expr_data[!is.na(rownames(expr_data)), ]

expr_data <- expr_data %>%
  as.data.frame() %>%
  rownames_to_column("gene_symbol") %>%
  group_by(gene_symbol) %>%
  summarize(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
  as.data.frame()

rownames(expr_data) <- expr_data$gene_symbol
expr_data <- expr_data[, -which(names(expr_data) == "gene_symbol")]

# Step 4: Define groups for samples
pheno_data$group <- c("Control", "Control", "Compound", "Compound")

# Create the design matrix
design <- model.matrix(~ 0 + factor(pheno_data$group))
colnames(design) <- c("Control", "Compound")
print(design)  # Debug: Check the design matrix

# Create contrast matrix with a fixed contrast
contrast_matrix <- makeContrasts(Control_vs_Compound = Compound - Control, levels = design)

# Fit the model
fit <- lmFit(expr_data, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Step 5: Get DEGs and prepare for database insertion
deg_results <- topTable(fit2, adjust.method = "none", number = Inf)
deg_results$adj.P.Value <- p.adjust(deg_results$P.Value, method = "BH")

usable_threshold <- 0.05

# Add a column for the p-value type used ("Adjusted" or "Raw")
deg_results$p_value_used <- ifelse(deg_results$adj.P.Value <= usable_threshold,
                                   "Adjusted", "Raw")

# Use adjusted p-values if they are below the threshold, else use raw p-values
deg_results$final_p_value <- ifelse(deg_results$adj.P.Value <= usable_threshold,
                                    deg_results$adj.P.Value, deg_results$P.Value)


# Use raw p-value threshold of 0.01
deg_results$significant <- "Not changed"
deg_results$significant[deg_results$P.Value <= 0.05 & deg_results$logFC > 1] <- "Up regulated"
deg_results$significant[deg_results$P.Value <= 0.05 & deg_results$logFC < -1] <- "Down regulated"

# Filter DEGs for saving
deg_results_filtered <- deg_results[deg_results$significant %in% c("Up regulated", "Down regulated"), ]

if (nrow(deg_results_filtered) > 0) {
  degs_to_save <- data.frame(
    gene_name = rownames(deg_results_filtered),
    log_fold_change = deg_results_filtered$logFC,
    p_value = deg_results_filtered$P.Value,
    degs = deg_results_filtered$significant,
    experiment_id = experiment_id
  )

  con <- dbConnect(
    RPostgres::Postgres(),
    host = "/run/postgresql",
    port = 5432,
    dbname = "adhd_research",
    user = "postgres",
    password = "admin"
  )

  dbWriteTable(con, name = "degs", value = degs_to_save, row.names = FALSE, append = TRUE)
  dbDisconnect(con)
} else {
  cat("No DEGs found to save.\n")
}

# Step 7: Plot a volcano plot
if (nrow(deg_results) > 0) {
  volcano_plot <- ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Value), color = significant)) +
    geom_point(alpha = 0.8) +
    scale_color_manual(values = c("blue", "gray", "red")) +
    labs(title = "Volcano plot", x = "Log2 (Fold change)", y = "-Log10 (Adjusted P-value)") +
    theme_minimal()

  ggsave("volcano_plot_final.png", volcano_plot, width = 8, height = 6)
  print(volcano_plot)
} else {
  cat("No data for volcano plot.\n")
}

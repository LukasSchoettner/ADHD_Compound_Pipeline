# Load required libraries
library(optparse)
library(GEOquery)
library(limma)
library(Biobase)
library(ggplot2)
library(DBI)
library(RPostgres)

# Parse command-line arguments
option_list <- list(
  make_option(c("--geo_id"), type = "character", help = "GEO ID to analyze"),
  make_option(c("--samples"), type = "character", help = "Samples to include, comma-separated"),
  make_option(c("--experiment_id"), type = "character", help = "Experiment ID for database integration")
  make_option(c("--db_connection_string"), type = "character", help = "Database connection string")
)
opt <- parse_args(OptionParser(option_list = option_list))

geo_id <- opt$geo_id
samples <- unlist(strsplit(opt$samples, ","))
experiment_id <- opt$experiment_id
db_connection_string <- opt$db_connection_string

# Step 1: Download the GEO dataset
gse <- getGEO(geo_id, GSEMatrix = TRUE)

# Step 2: Handle cases where `gse` is a list
if (length(gse) > 1) {
  gse <- gse[[1]]  # Select the first element
} else {
  gse <- gse[[1]]
}

# Extract expression data and phenotype data
expr_data <- exprs(gse)
pheno_data <- pData(gse)

# Subset the data to include only the samples of interest
expr_data <- expr_data[, samples]
pheno_data <- pheno_data[samples, ]

# Define the groups for the samples
pheno_data$group <- c("Control", "Control", "Gastrodin", "Gastrodin")

# Create the design matrix using meaningful group names
design <- model.matrix(~ 0 + factor(pheno_data$group))
colnames(design) <- c("Control", "Gastrodin")

# Fit the linear model
fit <- lmFit(expr_data, design)

# Define contrasts (Control vs Gastrodin)
contrast_matrix <- makeContrasts(
  Control_vs_Gastrodin = Gastrodin - Control,
  levels = design
)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Get top table of results
deg_results <- topTable(fit2, adjust.method = "none", number = Inf)

# Adjust p-values using the chosen method (e.g., "BH" for Benjamini-Hochberg)
deg_results$adj.P.Value <- p.adjust(deg_results$P.Value, method = "BH")

# Define a usability threshold for adjusted p-values (e.g., 0.05)
usable_threshold <- 0.05

# Add a column for the p-value type used ("Adjusted" or "Raw")
deg_results$p_value_used <- ifelse(deg_results$adj.P.Value <= usable_threshold,
                                   "Adjusted", "Raw")

# Use adjusted p-values if they are below the threshold, else use raw p-values
deg_results$final_p_value <- ifelse(deg_results$adj.P.Value <= usable_threshold,
                                    deg_results$adj.P.Value, deg_results$P.Value)

# Filter DEGs based on the final p-value and fold change
deg_results$significant <- "Not changed"
deg_results$significant[deg_results$final_p_value < 0.01 & deg_results$logFC > 1] <- "Up regulated"
deg_results$significant[deg_results$final_p_value < 0.01 & deg_results$logFC < -1] <- "Down regulated"

# Prepare the data for insertion into the database
degs_to_save <- data.frame(
  deg_id = seq_len(nrow(deg_results)),                     # Unique identifier
  gene_name = rownames(deg_results),                      # Gene names
  log_fold_change = deg_results$logFC,                    # Log fold change
  p_value = deg_results$final_p_value,                    # Final p-value
  degs = deg_results$significant,                         # DEG classification
  experiment_id = experiment_id                           # Experiment ID
)

con <- dbConnect(RPostgres::Postgres(), db_connection_string)

# Save data to the database
dbWriteTable(
  con,
  name = "degs",       # Table name in the database
  value = degs_to_save, # Data frame to write
  row.names = FALSE,   # Do not save row names
  append = TRUE        # Append to the existing table
)

# Disconnect from the database
dbDisconnect(con)

# Save the results with additional information
write.csv(deg_results, "DEG_results_final.csv", row.names = TRUE)

# Plot a volcano plot with the final p-values
volcano_plot <- ggplot(deg_results, aes(x = logFC, y = -log10(final_p_value), color = significant)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("blue", "gray", "red")) +
  labs(
    title = "Volcano plot",
    x = "Log2 (Fold change)",
    y = "-Log10 (Final P-value)"
  ) +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom"
  )

# Save the updated plot
ggsave("volcano_plot_final.png", volcano_plot, width = 8, height = 6)

# Display the plot
print(volcano_plot)

# Load required libraries
library(GEOquery)
library(limma)
library(Biobase)
library(ggplot2)
library(pheatmap)

# Step 1: Download the GEO dataset
gse <- getGEO("GSE85871", GSEMatrix = TRUE)

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
samples_of_interest <- c("GSM2286316", "GSM2286317", "GSM2286238", "GSM2286239")
expr_data <- expr_data[, samples_of_interest]
pheno_data <- pheno_data[samples_of_interest, ]

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

# Include additional info about the p-value type used to determine significance
deg_results$used_for_DEG <- ifelse(deg_results$significant != "Not changed",
                                   deg_results$p_value_used, NA)

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

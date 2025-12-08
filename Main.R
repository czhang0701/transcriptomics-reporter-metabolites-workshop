# ============================================================================
# MAIN WORKSHOP SCRIPT
# Transcriptomics to Reporter Metabolites - 3-Hour Workshop
# ============================================================================
#
# This script contains all code from the 4 core modules.
# Run line-by-line or section-by-section during the workshop.
#
# BEFORE STARTING:
# 1. Set your working directory to the workshop_materials folder
# 2. Install all required packages (see scripts/00_install_packages.R)
# 3. Verify installation (see scripts/00_check_setup.R)
#
# ============================================================================

# SET WORKING DIRECTORY
# Adjust this path to where you extracted the workshop materials
setwd("C:/path/to/workshop_materials")  # CHANGE THIS PATH!

# Verify you're in the correct directory
getwd()
list.files()  # Should see: data/, scripts/, tutorials/, etc.

# ============================================================================
# INSTALLATION & SETUP (30 minutes)
# ============================================================================

# If packages not installed yet, run:
# source("scripts/00_install_packages.R")

# Verify installation:
source("scripts/00_check_setup.R")

# ============================================================================
# MODULE 1: DIFFERENTIAL EXPRESSION ANALYSIS (20 minutes)
# Using DESeq2 to identify differentially expressed genes
# ============================================================================

cat("\n=== MODULE 1: DIFFERENTIAL EXPRESSION ===\n")

# Load required libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)

## Step 1: Load Data ----
cat("\n--- Step 1: Loading Data ---\n")

# Load count matrix
Counts <- read.table("data/Counts_selected.txt",
                     header = TRUE,
                     row.names = 1,
                     check.names = FALSE)

# Load metadata
Metadata <- read.table("data/patients.txt",
                       header = TRUE,
                       row.names = 1)

# Check dimensions
cat("Count matrix dimensions:", dim(Counts), "\n")
cat("Metadata dimensions:", dim(Metadata), "\n")

# Preview data
head(Counts[, 1:5])
head(Metadata)

## Step 2: Create DESeq2 Object ----
cat("\n--- Step 2: Creating DESeq2 Object ---\n")

# Ensure sample order matches
Metadata <- Metadata[colnames(Counts), , drop = FALSE]
all(rownames(Metadata) == colnames(Counts))  # Should be TRUE

# Create DESeq2 dataset
deSeqData <- DESeqDataSetFromMatrix(
  countData = Counts,
  colData = Metadata,
  design = ~ Type
)

cat("DESeq2 object created successfully\n")
print(deSeqData)

## Step 3: Pre-filtering (Optional) ----
cat("\n--- Step 3: Pre-filtering Low Count Genes ---\n")

# Keep genes with at least 10 reads total
keep <- rowSums(counts(deSeqData)) >= 10
deSeqData <- deSeqData[keep, ]

cat("Genes remaining after filtering:", nrow(deSeqData), "\n")

## Step 4: Run DESeq2 Analysis ----
cat("\n--- Step 4: Running DESeq2 (this may take 1-2 minutes) ---\n")

deSeqAnalysis <- DESeq(deSeqData)

cat("DESeq2 analysis complete!\n")

## Step 5: Extract Results ----
cat("\n--- Step 5: Extracting Results ---\n")

# Get results (Late vs Early)
res <- results(
  deSeqAnalysis,
  contrast = c("Type", "Late", "Early"),
  alpha = 0.05
)

# Summary of results
summary(res)

# Add gene symbols to results
cat("\nAdding gene symbols to results...\n")
ensembl2gene <- read.table("data/Ensembl2gene.tsv",
                          header = TRUE,
                          sep = "\t")
res$gene_symbol <- ensembl2gene$Gene[match(rownames(res), ensembl2gene$Ensembl)]

# Order by adjusted p-value
res_ordered <- res[order(res$padj), ]
cat("\nTop 10 Differentially Expressed Genes:\n")
print(res_ordered[1:10, c("gene_symbol", "baseMean", "log2FoldChange", "padj")])

## Step 6: Count Significant Genes ----
cat("\n--- Step 6: Summary Statistics ---\n")

# Count significant genes
sig_genes <- sum(res$padj < 0.05, na.rm = TRUE)
cat("Significant genes (padj < 0.05):", sig_genes, "\n")

# Up and down regulated
up_genes <- sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE)
down_genes <- sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm = TRUE)

cat("Up-regulated genes:", up_genes, "\n")
cat("Down-regulated genes:", down_genes, "\n")

## Step 7: Volcano Plot ----
cat("\n--- Step 7: Creating Volcano Plot ---\n")

# Prepare data for plotting
res_df <- as.data.frame(res)

# Categorize genes by regulation direction
res_df$regulation <- "Not Significant"
res_df$regulation[res_df$padj < 0.05 & res_df$log2FoldChange > 0] <- "Up-regulated"
res_df$regulation[res_df$padj < 0.05 & res_df$log2FoldChange < 0] <- "Down-regulated"
res_df$regulation <- factor(res_df$regulation,
                            levels = c("Up-regulated", "Down-regulated", "Not Significant"))

# Label top 10 genes
res_df$label <- ""
top_10_indices <- order(res_df$padj)[1:10]
res_df$label[top_10_indices] <- res_df$gene_symbol[top_10_indices]

# Create volcano plot
library(ggrepel)
volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = regulation), alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up-regulated" = "red",
                                 "Down-regulated" = "blue",
                                 "Not Significant" = "gray")) +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 20) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: Late vs Early Stage",
    subtitle = "Top 10 genes labeled",
    x = "Log2 Fold Change (Late vs Early)",
    y = "-Log10 Adjusted P-value",
    color = "Regulation"
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", alpha = 0.5)

print(volcano_plot)

# Save plot
ggsave("figures/Module1_volcano_plot.pdf", volcano_plot, width = 8, height = 6)
cat("Volcano plot saved to figures/Module1_volcano_plot.pdf\n")

## Step 8: Save Results ----
cat("\n--- Step 8: Saving Results ---\n")

# Save all results
write.table(
  as.data.frame(res),
  file = "data/DESeq_output.txt",
  sep = "\t",
  quote = FALSE,
  row.names = TRUE
)

cat("Results saved to data/DESeq_output.txt\n")

cat("\n=== MODULE 1 COMPLETE ===\n")


# ============================================================================
# MODULE 2: GENE SET ENRICHMENT ANALYSIS (15 minutes)
# Using PIANO to identify enriched pathways
# ============================================================================

cat("\n=== MODULE 2: GENE SET ENRICHMENT ANALYSIS ===\n")

# Load required libraries
library(piano)

## Step 1: Load DESeq2 Results ----
cat("\n--- Step 1: Loading DESeq2 Results ---\n")

# Read results if not already in memory
if (!exists("res")) {
  res <- read.table("data/DESeq_output.txt",
                    header = TRUE,
                    row.names = 1)
}

cat("Loaded", nrow(res), "genes\n")

## Step 2: Map Ensembl IDs to Gene Symbols ----
cat("\n--- Step 2: Mapping Gene IDs ---\n")

# Load mapping file
ensembl2gene <- read.table("data/Ensembl2gene.tsv",
                          header = TRUE,
                          sep = "\t")

# Match Ensembl IDs to gene symbols
res$gene_symbol <- ensembl2gene$Gene[match(rownames(res), ensembl2gene$Ensembl)]

# Remove genes without symbols
res_with_symbols <- res[!is.na(res$gene_symbol), ]
cat("Genes with symbols:", nrow(res_with_symbols), "\n")

## Step 3: Prepare Gene Statistics ----
cat("\n--- Step 3: Preparing Statistics for PIANO ---\n")

# Create p-value vector (use gene symbols as names)
pval_vector <- res_with_symbols$pvalue
names(pval_vector) <- res_with_symbols$gene_symbol

# Create fold change direction vector
fc_vector <- res_with_symbols$log2FoldChange
names(fc_vector) <- res_with_symbols$gene_symbol

# Remove NAs
pval_vector <- pval_vector[!is.na(pval_vector)]
fc_vector <- fc_vector[!is.na(fc_vector)]

cat("Genes for GSEA:", length(pval_vector), "\n")

## Step 4: Load Gene Sets ----
cat("\n--- Step 4: Loading Gene Sets ---\n")

# Load MSigDB gene sets (GO Biological Processes)
gsc <- loadGSC("data/c5.bp.v6.2.symbols.gmt")

cat("Number of gene sets loaded:", length(gsc$gsc), "\n")

## Step 5: Run Gene Set Analysis ----
cat("\n--- Step 5: Running PIANO (this may take 2-3 minutes) ---\n")

# Run GSA with reporter genes method
gsaRes <- runGSA(
  geneLevelStats = pval_vector,
  directions = fc_vector,
  gsc = gsc,
  geneSetStat = "reporter",      # Reporter genes method
  signifMethod = "geneSampling",  # Gene sampling for significance
  nPerm = 1000,                   # Number of permutations
  gsSizeLim = c(5, 500),         # Gene set size limits
  adjMethod = "fdr"               # Multiple testing correction
)

cat("PIANO analysis complete!\n")

## Step 6: Extract Top Pathways ----
cat("\n--- Step 6: Viewing Top Enriched Pathways ---\n")

# Get consensus scores (direction-aware)
gsaSummary <- GSAsummaryTable(gsaRes)

# View top 20 pathways
cat("\nTop 20 Enriched Pathways:\n")
print(head(gsaSummary, 20))

## Step 7: Identify Metabolic Pathways ----
cat("\n--- Step 7: Identifying Metabolic Pathways ---\n")

# Filter for metabolic pathways
metabolic_keywords <- c("METABOL", "GLYCOLYSIS", "OXIDATIVE_PHOSPHORYLATION",
                       "TCA", "CITRIC_ACID", "FATTY_ACID", "LIPID")

metabolic_pathways <- gsaSummary[
  grepl(paste(metabolic_keywords, collapse = "|"), rownames(gsaSummary)),
]

cat("\nTop Metabolic Pathways:\n")
print(head(metabolic_pathways, 15))

## Step 8: Visualize Top Pathways with Heatmap ----
cat("\n--- Step 8: Creating Pathway Enrichment Heatmap ---\n")

# Get top 20 pathways by significance
top_pathways <- head(gsaSummary, 20)

# Create a matrix for heatmap (p-values and direction)
# Using distinct directional p-values
heatmap_data <- data.frame(
  Up = -log10(top_pathways$`p adj (dist.dir.up)`),
  Down = -log10(top_pathways$`p adj (dist.dir.dn)`)
)

# Use pathway names from the Name column (not rownames/IDs)
if ("Name" %in% colnames(top_pathways)) {
  pathway_names <- top_pathways$Name
} else {
  pathway_names <- rownames(top_pathways)
}

# Clean up pathway names (remove "GO_" or "GOBP_" prefix and make readable)
pathway_names <- gsub("^GO_", "", pathway_names)
pathway_names <- gsub("^GOBP_", "", pathway_names)
pathway_names <- gsub("_", " ", pathway_names)
pathway_names <- substr(pathway_names, 1, 60)  # Truncate long names

# Set cleaned names as rownames
rownames(heatmap_data) <- pathway_names

# Create heatmap
library(pheatmap)

# Save heatmap to PDF
pdf("figures/Module2_pathway_heatmap.pdf", width = 10, height = 12)
pathway_heatmap <- pheatmap(
  as.matrix(heatmap_data),
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  color = colorRampPalette(c("white", "orange", "red"))(50),
  main = "Top 20 Enriched Pathways (GSEA)",
  fontsize_row = 10,
  fontsize_col = 12,
  angle_col = 0,
  cellwidth = 80,
  cellheight = 20,
  border_color = "gray90",
  na_col = "white",
  legend = TRUE
)
dev.off()

cat("Pathway heatmap saved to figures/Module2_pathway_heatmap.pdf\n")

# Also display in R (optional)
print(pathway_heatmap)

## Step 9: Save Results ----
cat("\n--- Step 9: Saving Results ---\n")

# Save full GSEA results
write.table(
  gsaSummary,
  file = "data/Piano_output.txt",
  sep = "\t",
  quote = FALSE,
  row.names = TRUE
)

cat("GSEA results saved to data/Piano_output.txt\n")

cat("\n=== MODULE 2 COMPLETE ===\n")


# ============================================================================
# MODULE 3: CO-EXPRESSION NETWORK ANALYSIS (15 minutes)
# Building gene correlation networks
# ============================================================================

cat("\n=== MODULE 3: CO-EXPRESSION NETWORK ANALYSIS ===\n")

# Load required libraries
library(Hmisc)
library(corrplot)
library(reshape2)

## Step 1: Load FPKM Data ----
cat("\n--- Step 1: Loading Expression Data ---\n")

# Load FPKM (normalized expression)
FPKM <- read.table("data/FPKM_selected.txt",
                   header = TRUE,
                   row.names = 1,
                   check.names = FALSE)

cat("FPKM matrix dimensions:", dim(FPKM), "\n")

## Step 2: Filter for Significant Genes ----
cat("\n--- Step 2: Filtering Significant Genes ---\n")

# Get significant genes from Module 1
sig_gene_ids <- rownames(res)[res$padj < 0.05 & !is.na(res$padj)]

# Subset FPKM to significant genes only
FPKM_sig <- FPKM[rownames(FPKM) %in% sig_gene_ids, ]

cat("Significant genes for network:", nrow(FPKM_sig), "\n")

# For faster computation, limit to top genes by variance
gene_vars <- apply(FPKM_sig, 1, var)
top_genes <- names(sort(gene_vars, decreasing = TRUE)[1:min(200, nrow(FPKM_sig))])
FPKM_network <- FPKM_sig[top_genes, ]

cat("Genes for network analysis:", nrow(FPKM_network), "\n")

## Step 3: Calculate Correlations ----
cat("\n--- Step 3: Calculating Spearman Correlations ---\n")

# Transpose (genes as columns)
FPKM_transposed <- t(FPKM_network)

# Calculate Spearman correlation with p-values
cor_result <- rcorr(FPKM_transposed, type = "spearman")

cat("Correlation matrix calculated\n")
cat("Correlation matrix dimensions:", dim(cor_result$r), "\n")

## Step 4: Convert to Long Format ----
cat("\n--- Step 4: Preparing Network Data ---\n")

# Get correlation matrix and p-values
cor_matrix <- cor_result$r
pval_matrix <- cor_result$P

# Convert to long format
cor_long <- melt(cor_matrix, varnames = c("Gene1", "Gene2"), value.name = "Correlation")
pval_long <- melt(pval_matrix, varnames = c("Gene1", "Gene2"), value.name = "Pvalue")

# Combine
network_data <- merge(cor_long, pval_long, by = c("Gene1", "Gene2"))

# Remove self-correlations
network_data <- network_data[network_data$Gene1 != network_data$Gene2, ]

# Remove duplicate pairs (keep only upper triangle)
network_data <- network_data[as.character(network_data$Gene1) < as.character(network_data$Gene2), ]

cat("Total gene pairs:", nrow(network_data), "\n")

## Step 5: FDR Correction ----
cat("\n--- Step 5: Multiple Testing Correction ---\n")

# FDR correction
network_data$FDR <- p.adjust(network_data$Pvalue, method = "fdr")

cat("FDR correction applied\n")

## Step 6: Filter Significant Correlations ----
cat("\n--- Step 6: Filtering Significant Correlations ---\n")

# Filter for significant and strong correlations
significant_network <- network_data[
  network_data$FDR < 0.01 &
  abs(network_data$Correlation) > 0.7,
]

cat("Significant correlations (FDR < 0.01, |r| > 0.7):", nrow(significant_network), "\n")

# Positive and negative correlations
pos_cor <- sum(significant_network$Correlation > 0)
neg_cor <- sum(significant_network$Correlation < 0)

cat("Positive correlations:", pos_cor, "\n")
cat("Negative correlations:", neg_cor, "\n")

## Step 7: Visualize Top Correlations ----
cat("\n--- Step 7: Visualizing Network ---\n")

# Get top genes by number of connections
gene_connections <- table(c(as.character(significant_network$Gene1),
                           as.character(significant_network$Gene2)))
top_connected_genes <- names(sort(gene_connections, decreasing = TRUE)[1:min(50, length(gene_connections))])

# Subset correlation matrix for visualization
cor_subset <- cor_matrix[top_connected_genes, top_connected_genes]

# Create correlation plot
pdf("figures/Module3_correlation_network.pdf", width = 10, height = 10)
corrplot(cor_subset,
         method = "color",
         type = "upper",
         tl.col = "black",
         tl.cex = 0.6,
         title = "Co-expression Network (Top 50 Connected Genes)")
dev.off()

cat("Correlation plot saved to figures/Module3_correlation_network.pdf\n")

## Step 8: Save Network Data ----
cat("\n--- Step 8: Saving Network Data ---\n")

# Save significant correlations
write.table(
  significant_network,
  file = "data/Coexpression_network.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Network data saved to data/Coexpression_network.txt\n")
cat("Can be imported into Cytoscape for visualization\n")

cat("\n=== MODULE 3 COMPLETE ===\n")


# ============================================================================
# MODULE 4: REPORTER METABOLITE ANALYSIS (60 minutes) ⭐ MAIN FOCUS
# Identifying key metabolites driving metabolic changes
# ============================================================================

cat("\n=== MODULE 4: REPORTER METABOLITE ANALYSIS ===\n")

# Load required libraries
library(XML)

# Load custom reporter metabolites functions
source("scripts/parse_sbml_model.R")
source("scripts/reporter_metabolites.R")

## Step 1: Understand Reporter Metabolites ----
cat("\n--- Step 1: What Are Reporter Metabolites? ---\n")
cat("
Reporter metabolites are metabolites around which
transcriptional regulation is significantly enriched.

Algorithm (Patil & Nielsen, 2005):
1. Convert gene P-values to Z-scores
2. For each metabolite, aggregate Z-scores of neighboring genes
3. Correct for background using random sampling
4. Calculate metabolite-level P-values

They connect:
  Genes → Enzymes → Reactions → Metabolites → Pathways

This helps identify KEY metabolic hubs in your biological system.
\n")

## Step 2: Load or Parse Metabolic Model ----
cat("\n--- Step 2: Loading Metabolic Model ---\n")

# Check if pre-parsed model exists
model_file <- "data/parsed_model.rds"

if (file.exists(model_file)) {
  cat("Loading pre-parsed model (fast)...\n")
  model <- loadModel(model_file)
} else {
  cat("Parsing SBML model for the first time...\n")
  cat("This will take 3-5 minutes but only needs to be done once!\n\n")

  model <- parseSBMLModel("data/Reference_model.xml")

  # Save for future use
  saveModel(model, model_file)
  cat("\nModel saved! Next time this will load instantly.\n")
}

## Step 3: Prepare Gene Statistics ----
cat("\n--- Step 3: Preparing Gene Statistics ---\n")

# Get gene symbols and statistics from DESeq2 results
genes_for_reporter <- res$gene_symbol[!is.na(res$gene_symbol)]
pvals_for_reporter <- res$pvalue[!is.na(res$gene_symbol)]
fc_for_reporter <- res$log2FoldChange[!is.na(res$gene_symbol)]

# Remove any remaining NAs
valid <- !is.na(pvals_for_reporter) & !is.na(fc_for_reporter)
genes_for_reporter <- genes_for_reporter[valid]
pvals_for_reporter <- pvals_for_reporter[valid]
fc_for_reporter <- fc_for_reporter[valid]

cat("Genes for reporter analysis:", length(genes_for_reporter), "\n")
cat("Genes in model:", sum(genes_for_reporter %in% model$genes), "\n")

## Step 4: Run Reporter Metabolite Analysis ----
cat("\n--- Step 4: Running Reporter Metabolite Analysis ---\n")
cat("This will take 5-10 minutes for full analysis...\n\n")

# Run reporter metabolites algorithm
reporter_results <- reporterMetabolites(
  model = model,
  genes = genes_for_reporter,
  genePValues = pvals_for_reporter,
  geneFoldChanges = fc_for_reporter,
  printResults = TRUE  # Print top 20 to console
)

cat("\nReporter metabolite analysis complete!\n")

## Step 5: Interpret Results ----
cat("\n--- Step 5: Biological Interpretation ---\n")

# Extract results for "all genes" test
reporter_all <- reporter_results[[1]]

# Show top significant metabolites
cat("\nTop 10 Significant Reporter Metabolites:\n")
cat(sprintf("%-20s %-40s %12s %10s\n", "ID", "Name", "P-value", "N Genes"))
cat(strrep("-", 85), "\n")
for (i in 1:min(10, length(reporter_all$mets))) {
  cat(sprintf("%-20s %-40s %.2e %10d\n",
              substr(reporter_all$mets[i], 1, 20),
              substr(reporter_all$metNames[i], 1, 40),
              reporter_all$metPValues[i],
              reporter_all$metNGenes[i]))
}

# Look for specific metabolic signatures
cat("\nKey Metabolic Changes Detected:\n\n")

metabolite_names <- tolower(paste(reporter_all$metNames, collapse = " "))

if (grepl("glucose|glycolysis|lactate|pyruvate", metabolite_names)) {
  cat("✓ GLYCOLYSIS/GLUCOSE metabolism\n")
  cat("  → Changes in glucose utilization\n")
  cat("  → Potential Warburg effect\n\n")
}

if (grepl("citrate|succinate|malate|fumarate|tca|krebs", metabolite_names)) {
  cat("✓ TCA CYCLE metabolites\n")
  cat("  → Altered energy metabolism\n")
  cat("  → Mitochondrial reprogramming\n\n")
}

if (grepl("fatty|lipid|acyl|palmit", metabolite_names)) {
  cat("✓ FATTY ACID metabolism\n")
  cat("  → Lipid metabolism changes\n")
  cat("  → Energy source shift\n\n")
}

if (grepl("amino|glutam|leucine|alanine", metabolite_names)) {
  cat("✓ AMINO ACID metabolism\n")
  cat("  → Protein/nitrogen metabolism\n")
  cat("  → Potential nutrient stress response\n\n")
}

## Step 6: Connect to Earlier Results ----
cat("\n--- Step 6: Connecting Results Across Modules ---\n")

cat("
INTEGRATION ACROSS MODULES:

Module 1 (DESeq2):
  → Found", sig_genes, "differentially expressed genes

Module 2 (GSEA):
  → Identified enriched biological pathways

Module 3 (Networks):
  → Found", nrow(significant_network), "significant gene correlations

Module 4 (Reporter Metabolites): ⭐
  → Identified", length(reporter_all$mets), "metabolites with scores
  → Found", sum(reporter_all$metPValues < 0.05), "significant reporter metabolites (P < 0.05)
  → KEY METABOLIC HUBS identified!

BIOLOGICAL INSIGHT:
These reporter metabolites represent the KEY NODES in metabolic
reprogramming. They are the metabolites around which most
transcriptional changes occur.

Algorithm: For each metabolite, we:
1. Identified all genes encoding enzymes that produce/consume it
2. Aggregated their differential expression Z-scores
3. Corrected for gene set size using random sampling
4. Calculated metabolite-level significance

NEXT STEPS:
1. Validate with metabolomics data (measure actual metabolite levels)
2. Target high-scoring metabolites experimentally
3. Consider as therapeutic targets
")

## Step 7: Save Results ----
cat("\n--- Step 7: Saving Reporter Metabolite Results ---\n")

# Save detailed results using custom function
saveReporterMetabolites(reporter_results, "data/Reporter_Metabolites_output.txt")

# Also create a simple summary table for the "all genes" test
reporter_df <- data.frame(
  Metabolite_ID = reporter_all$mets,
  Metabolite_Name = reporter_all$metNames,
  Z_Score = reporter_all$metZScores,
  P_Value = reporter_all$metPValues,
  N_Genes = reporter_all$metNGenes,
  Mean_Z = reporter_all$meanZ,
  Std_Z = reporter_all$stdZ,
  stringsAsFactors = FALSE
)

write.table(
  reporter_df,
  file = "data/Reporter_Metabolites_summary.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Reporter metabolite summary saved to data/Reporter_Metabolites_summary.txt\n")

## Step 8: Visualization ----
cat("\n--- Step 8: Visualizing Reporter Metabolites ---\n")

# Get top 20 metabolites by Z-score
top_n <- min(20, length(reporter_all$mets))
top_indices <- 1:top_n

# Create data for visualization
viz_data <- data.frame(
  Metabolite = reporter_all$metNames[top_indices],
  Z_Score = reporter_all$metZScores[top_indices],
  Neg_Log_P = -log10(reporter_all$metPValues[top_indices]),
  N_Genes = reporter_all$metNGenes[top_indices]
)

# Clean up metabolite names for display
viz_data$Metabolite <- gsub("\\[.*\\]", "", viz_data$Metabolite)  # Remove compartment
viz_data$Metabolite <- substr(viz_data$Metabolite, 1, 40)  # Truncate long names

# Create bar plot of Z-scores
pdf("figures/Module4_reporter_metabolites.pdf", width = 12, height = 10)
par(mar = c(5, 15, 4, 2))
barplot(
  viz_data$Z_Score,
  names.arg = viz_data$Metabolite,
  horiz = TRUE,
  las = 1,
  col = ifelse(viz_data$Z_Score > 0, "firebrick", "steelblue"),
  main = "Top 20 Reporter Metabolites (Ordered by Z-score)",
  xlab = "Aggregated Z-score",
  cex.names = 0.8
)
abline(v = 0, lty = 2, col = "gray50")
dev.off()

cat("Reporter metabolite plot saved to figures/Module4_reporter_metabolites.pdf\n")

# If up/down results exist, create comparison heatmap
if (length(reporter_results) == 3) {
  cat("\nCreating up/down regulation heatmap...\n")

  # Get top metabolites from each test
  top_all <- reporter_results[[1]]$mets[1:min(15, length(reporter_results[[1]]$mets))]
  top_up <- reporter_results[[2]]$mets[1:min(15, length(reporter_results[[2]]$mets))]
  top_down <- reporter_results[[3]]$mets[1:min(15, length(reporter_results[[3]]$mets))]

  # Union of top metabolites
  all_top_mets <- unique(c(top_all, top_up, top_down))

  # Create matrix
  heatmap_data <- matrix(0, nrow = length(all_top_mets), ncol = 3)
  colnames(heatmap_data) <- c("All Genes", "Up-regulated", "Down-regulated")
  rownames(heatmap_data) <- all_top_mets

  for (i in 1:length(all_top_mets)) {
    met <- all_top_mets[i]

    # All genes
    idx <- which(reporter_results[[1]]$mets == met)
    if (length(idx) > 0) {
      heatmap_data[i, 1] <- reporter_results[[1]]$metZScores[idx]
    }

    # Up-regulated
    idx <- which(reporter_results[[2]]$mets == met)
    if (length(idx) > 0) {
      heatmap_data[i, 2] <- reporter_results[[2]]$metZScores[idx]
    }

    # Down-regulated
    idx <- which(reporter_results[[3]]$mets == met)
    if (length(idx) > 0) {
      heatmap_data[i, 3] <- reporter_results[[3]]$metZScores[idx]
    }
  }

  # Get metabolite names
  met_names <- sapply(all_top_mets, function(met) {
    idx <- which(reporter_all$mets == met)
    if (length(idx) > 0) {
      name <- reporter_all$metNames[idx]
      name <- gsub("\\[.*\\]", "", name)
      return(substr(name, 1, 40))
    }
    return(met)
  })
  rownames(heatmap_data) <- met_names

  pdf("figures/Module4_reporter_metabolites_heatmap.pdf", width = 10, height = 12)
  pheatmap(
    heatmap_data,
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    color = colorRampPalette(c("blue", "white", "red"))(50),
    main = "Reporter Metabolites: Up vs Down Regulation",
    fontsize_row = 9,
    fontsize_col = 11,
    display_numbers = FALSE,
    breaks = seq(-max(abs(heatmap_data)), max(abs(heatmap_data)), length.out = 51)
  )
  dev.off()

  cat("Heatmap saved to figures/Module4_reporter_metabolites_heatmap.pdf\n")
}

cat("\n=== MODULE 4 COMPLETE ===\n")


# ============================================================================
# WORKSHOP COMPLETE! 🎉
# ============================================================================

cat("\n")
cat("============================================================================\n")
cat("                    WORKSHOP COMPLETE! 🎉                                   \n")
cat("============================================================================\n")
cat("\n")
cat("You have successfully:\n")
cat("  ✓ Analyzed differential gene expression (DESeq2)\n")
cat("  ✓ Identified enriched pathways (PIANO GSEA)\n")
cat("  ✓ Built co-expression networks\n")
cat("  ✓ Identified REPORTER METABOLITES ⭐\n")
cat("  ✓ Connected genes to metabolic reprogramming\n")
cat("\n")
cat("Output files generated:\n")
cat("  • data/DESeq_output.txt\n")
cat("  • data/Piano_output.txt\n")
cat("  • data/Coexpression_network.txt\n")
cat("  • data/Reporter_Metabolites_output.txt ⭐\n")
cat("  • figures/ (all plots)\n")
cat("\n")
cat("Next steps:\n")
cat("  1. Explore the detailed tutorials in tutorials/ folder\n")
cat("  2. Try the practice dataset (scripts/create_practice_dataset.R)\n")
cat("  3. Apply to your own data!\n")
cat("  4. (Optional) Explore MOFA2 bonus content\n")
cat("\n")
cat("Questions? Check the tutorials or README.md\n")
cat("\n")
cat("============================================================================\n")

# ============================================================================
# BONUS: Quick Session Info
# ============================================================================

cat("\nSession Information:\n")
sessionInfo()

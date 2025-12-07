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

# Order by adjusted p-value
res_ordered <- res[order(res$padj), ]
head(res_ordered, 10)

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
res_df$significant <- ifelse(res_df$padj < 0.05, "Significant", "Not Significant")

# Create volcano plot
volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = significant), alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Not Significant" = "gray", "Significant" = "red")) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: Late vs Early Stage",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue")

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

## Step 8: Save Results ----
cat("\n--- Step 8: Saving Results ---\n")

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
library(piano)

## Step 1: Understand Reporter Metabolites ----
cat("\n--- Step 1: What Are Reporter Metabolites? ---\n")
cat("
Reporter metabolites are metabolites around which
transcriptional regulation is significantly enriched.

They connect:
  Genes → Enzymes → Reactions → Metabolites → Pathways

This helps identify KEY metabolic hubs in your biological system.
\n")

## Step 2: Load Metabolic Model ----
cat("\n--- Step 2: Loading Genome-Scale Metabolic Model ---\n")

# Load metabolic model (SBML format)
# This contains metabolite-gene associations
library(XML)

model <- xmlParse("data/Reference_model.xml")
cat("Metabolic model loaded\n")

# In practice, we'll use pre-processed metabolite-gene associations
# For this workshop, we'll use example gene sets

## Step 3: Create Metabolite-Gene Associations ----
cat("\n--- Step 3: Creating Metabolite-Gene Gene Sets ---\n")

# For demonstration, we'll use GO Biological Process as proxy
# In real analysis, you would extract from metabolic model

# We already have gene sets from Module 2
# Let's focus on metabolic pathways

# Filter for metabolic gene sets
metabolic_gene_sets <- gsc$gsc[grepl("METABOL|GLYCOLYSIS|TCA|FATTY|LIPID|OXIDATIVE",
                                     names(gsc$gsc))]

cat("Number of metabolic gene sets:", length(metabolic_gene_sets), "\n")

# Create new gene set collection
gsc_metabolites <- list(gsc = metabolic_gene_sets)
class(gsc_metabolites) <- "GSC"

## Step 4: Prepare Gene Statistics ----
cat("\n--- Step 4: Preparing Gene Statistics ---\n")

# Use the same statistics from Module 2
# (already have pval_vector and fc_vector)

cat("Genes for reporter analysis:", length(pval_vector), "\n")

## Step 5: Run Reporter Metabolite Analysis ----
cat("\n--- Step 5: Running Reporter Analysis (2-3 minutes) ---\n")

# Run GSA with reporter method (same as Module 2, but focused on metabolites)
reporter_results <- runGSA(
  geneLevelStats = pval_vector,
  directions = fc_vector,
  gsc = gsc_metabolites,
  geneSetStat = "reporter",      # Reporter genes method
  signifMethod = "geneSampling",
  nPerm = 1000,
  gsSizeLim = c(3, 100),         # Metabolite gene sets are usually smaller
  adjMethod = "fdr"
)

cat("Reporter metabolite analysis complete!\n")

## Step 6: Extract Reporter Metabolites ----
cat("\n--- Step 6: Identifying Top Reporter Metabolites ---\n")

# Get results summary
reporter_summary <- GSAsummaryTable(reporter_results)

# View top reporter metabolites
cat("\nTop 20 Reporter Metabolites:\n")
print(head(reporter_summary, 20))

## Step 7: Biological Interpretation ----
cat("\n--- Step 7: Biological Interpretation ---\n")

# Look for specific metabolic signatures
cat("\nKey Metabolic Changes Detected:\n\n")

# Glycolysis
if (any(grepl("GLYCOLYSIS", rownames(reporter_summary)))) {
  cat("✓ GLYCOLYSIS pathways detected\n")
  cat("  → Warburg effect: increased lactate production\n")
  cat("  → Common in cancer cells\n\n")
}

# TCA cycle
if (any(grepl("TCA|CITRIC_ACID", rownames(reporter_summary)))) {
  cat("✓ TCA CYCLE changes detected\n")
  cat("  → Altered energy metabolism\n")
  cat("  → Potential metabolic reprogramming\n\n")
}

# Oxidative phosphorylation
if (any(grepl("OXIDATIVE_PHOSPHORYLATION", rownames(reporter_summary)))) {
  cat("✓ OXIDATIVE PHOSPHORYLATION changes\n")
  cat("  → Mitochondrial function altered\n")
  cat("  → May indicate metabolic shift\n\n")
}

# Fatty acid metabolism
if (any(grepl("FATTY_ACID", rownames(reporter_summary)))) {
  cat("✓ FATTY ACID metabolism changes\n")
  cat("  → Lipid metabolism reprogramming\n")
  cat("  → Energy source utilization changes\n\n")
}

## Step 8: Connect to Earlier Results ----
cat("\n--- Step 8: Connecting Results Across Modules ---\n")

cat("
INTEGRATION ACROSS MODULES:

Module 1 (DESeq2):
  → Found", sig_genes, "differentially expressed genes

Module 2 (GSEA):
  → Identified", nrow(metabolic_pathways), "enriched metabolic pathways

Module 3 (Networks):
  → Found", nrow(significant_network), "significant gene correlations

Module 4 (Reporter Metabolites): ⭐
  → Identified", nrow(reporter_summary), "metabolite-associated processes
  → KEY METABOLIC HUBS identified!

BIOLOGICAL INSIGHT:
These reporter metabolites represent the KEY NODES in metabolic
reprogramming. They are the metabolites around which most
transcriptional changes occur.

NEXT STEPS:
1. Validate with metabolomics data
2. Test hypotheses experimentally
3. Consider as therapeutic targets
")

## Step 9: Save Results ----
cat("\n--- Step 9: Saving Reporter Metabolite Results ---\n")

# Save reporter results
write.table(
  reporter_summary,
  file = "data/Reporter_Metabolites_output.txt",
  sep = "\t",
  quote = FALSE,
  row.names = TRUE
)

cat("Reporter metabolite results saved to data/Reporter_Metabolites_output.txt\n")

## Step 10: Visualization ----
cat("\n--- Step 10: Visualizing Reporter Metabolites ---\n")

# Create bar plot of top reporters
top_reporters <- head(reporter_summary, 15)
top_reporters$pathway <- rownames(top_reporters)

# Plot p-values (as -log10)
reporter_plot <- ggplot(top_reporters, aes(x = reorder(pathway, -`p adj (dist.dir.up)`),
                                           y = -log10(`p adj (dist.dir.up)`))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Top 15 Reporter Metabolites",
    x = "Metabolic Process",
    y = "-Log10 Adjusted P-value"
  ) +
  theme(axis.text.y = element_text(size = 8))

print(reporter_plot)

# Save plot
ggsave("figures/Module4_reporter_metabolites.pdf", reporter_plot, width = 10, height = 8)
cat("Reporter metabolite plot saved to figures/Module4_reporter_metabolites.pdf\n")

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

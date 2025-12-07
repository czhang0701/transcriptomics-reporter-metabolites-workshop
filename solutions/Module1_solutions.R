################################################################################
# Module 1: Differential Expression Analysis - Solutions
# Integrative Multi-Omics Analysis Workshop
################################################################################

# Load required libraries
library(DESeq2)
library(ggplot2)

# Set working directory (adjust to your path)
setwd("path/to/workshop_materials")

################################################################################
# Load and prepare data (from tutorial)
################################################################################

# Import counts data
Counts <- as.matrix(read.csv(file = "../Counts_selected.txt",
                             sep = "\t",
                             row.names = 1))

# Import metadata
Metadata <- read.csv(file = "../patients.txt",
                    sep = "\t",
                    row.names = 1)

# Create DESeqDataSet
deSeqData <- DESeqDataSetFromMatrix(
  countData = Counts,
  colData = Metadata,
  design = ~ Type
)

# Pre-filter low counts
keep <- rowSums(counts(deSeqData) >= 10) >= 5
deSeqData <- deSeqData[keep,]

# Run DESeq2
deSeqAnalysis <- DESeq(deSeqData)

# Extract results
res <- results(
  deSeqAnalysis,
  contrast = c("Type", "Late", "Early"),
  lfcThreshold = 0,
  altHypothesis = "greaterAbs",
  alpha = 0.05
)

################################################################################
# Exercise Solutions
################################################################################

### Question 1: How many genes have padj < 0.05?

# Solution:
sig_genes <- sum(res$padj < 0.05, na.rm = TRUE)
cat("Number of significant genes (padj < 0.05):", sig_genes, "\n")

# Expected: ~2000-3000 genes depending on filtering


### Question 2: How many genes have padj < 0.05 AND |log2FoldChange| > 1?

# Solution:
sig_lfc_genes <- sum(res$padj < 0.05 & abs(res$log2FoldChange) > 1, na.rm = TRUE)
cat("Significant genes with |LFC| > 1:", sig_lfc_genes, "\n")

# Alternative with explicit conditions:
sig_up <- sum(res$padj < 0.05 & res$log2FoldChange > 1, na.rm = TRUE)
sig_down <- sum(res$padj < 0.05 & res$log2FoldChange < -1, na.rm = TRUE)
cat("  Up-regulated:", sig_up, "\n")
cat("  Down-regulated:", sig_down, "\n")
cat("  Total:", sig_up + sig_down, "\n")


### Question 3: What is the most up-regulated gene in Late stage?

# Solution:
# Order by log2FoldChange (descending)
res_ordered <- res[order(-res$log2FoldChange), ]
most_upregulated <- res_ordered[1, ]

cat("\nMost up-regulated gene in Late stage:\n")
cat("Gene ID:", rownames(res_ordered)[1], "\n")
cat("Log2 Fold Change:", most_upregulated$log2FoldChange, "\n")
cat("Adjusted p-value:", most_upregulated$padj, "\n")

# Optional: Add gene symbol if mapping available
# (This would require loading Ensembl2gene.tsv)


### Question 4: What is the most down-regulated gene in Late stage?

# Solution:
# Order by log2FoldChange (ascending)
res_ordered_down <- res[order(res$log2FoldChange), ]
most_downregulated <- res_ordered_down[1, ]

cat("\nMost down-regulated gene in Late stage:\n")
cat("Gene ID:", rownames(res_ordered_down)[1], "\n")
cat("Log2 Fold Change:", most_downregulated$log2FoldChange, "\n")
cat("Adjusted p-value:", most_downregulated$padj, "\n")


### Question 5: Find a gene with high significance (padj < 0.001) but small fold change (|LFC| < 0.5)

# Solution:
# Filter for highly significant genes with small effect
small_effect_sig <- res[
  !is.na(res$padj) &
  res$padj < 0.001 &
  abs(res$log2FoldChange) < 0.5,
]

cat("\nHighly significant genes with small fold changes:\n")
print(head(small_effect_sig))

# Interpretation: These genes are consistently different between groups
# but the magnitude of change is small. Could be:
# 1. Highly expressed housekeeping genes with tight regulation
# 2. Genes with low biological variability
# 3. Important for subtle regulatory differences

cat("\nNumber of such genes:", nrow(small_effect_sig), "\n")


################################################################################
# Bonus Challenge: Create filtered subset
################################################################################

# Create subset with:
# - Significant genes (padj < 0.05)
# - Strong effect (|log2FoldChange| > 1.5)
# - Well-expressed genes (baseMean > 100)

bonus_subset <- res[
  !is.na(res$padj) &
  res$padj < 0.05 &
  abs(res$log2FoldChange) > 1.5 &
  res$baseMean > 100,
]

# Sort by adjusted p-value
bonus_subset <- bonus_subset[order(bonus_subset$padj), ]

cat("\n=== Bonus Challenge Results ===\n")
cat("Number of high-confidence differentially expressed genes:", nrow(bonus_subset), "\n")
cat("\nTop 10 genes:\n")
print(head(bonus_subset, 10))

# Save to file
write.table(
  as.data.frame(bonus_subset),
  file = "high_confidence_DEGs.txt",
  sep = "\t",
  row.names = TRUE,
  col.names = NA
)

cat("\nSubset saved to: high_confidence_DEGs.txt\n")


################################################################################
# Additional Explorations
################################################################################

### Explore distribution of fold changes

# Create histogram of log2 fold changes
pdf("LFC_distribution.pdf", width = 8, height = 6)
hist(res$log2FoldChange,
     breaks = 50,
     main = "Distribution of Log2 Fold Changes",
     xlab = "Log2 Fold Change (Late vs Early)",
     col = "lightblue",
     border = "darkblue")
abline(v = 0, col = "red", lwd = 2, lty = 2)
dev.off()


### Compare mean expression between groups

# Get normalized counts
norm_counts <- counts(deSeqAnalysis, normalized = TRUE)

# Calculate mean expression per group
early_samples <- which(Metadata$Type == "Early")
late_samples <- which(Metadata$Type == "Late")

mean_early <- rowMeans(norm_counts[, early_samples])
mean_late <- rowMeans(norm_counts[, late_samples])

# Create comparison data frame
comparison_df <- data.frame(
  Gene = rownames(res),
  Mean_Early = mean_early,
  Mean_Late = mean_late,
  Log2FC = res$log2FoldChange,
  padj = res$padj,
  Significant = ifelse(is.na(res$padj), "No",
                      ifelse(res$padj < 0.05, "Yes", "No"))
)

# Plot mean expression comparison
pdf("Mean_expression_comparison.pdf", width = 8, height = 8)
ggplot(comparison_df, aes(x = Mean_Early, y = Mean_Late, color = Significant)) +
  geom_point(alpha = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(values = c("Yes" = "red", "No" = "gray")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "Mean Expression: Early vs Late",
    x = "Mean Expression (Early)",
    y = "Mean Expression (Late)"
  )
dev.off()


################################################################################
# Summary Statistics
################################################################################

cat("\n=== Summary Statistics ===\n")
cat("Total genes tested:", nrow(res), "\n")
cat("Significant genes (FDR < 0.05):", sum(res$padj < 0.05, na.rm = TRUE), "\n")
cat("  Up in Late:", sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE), "\n")
cat("  Down in Late:", sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm = TRUE), "\n")
cat("\nHighly significant (FDR < 0.001):", sum(res$padj < 0.001, na.rm = TRUE), "\n")
cat("Large fold change (|LFC| > 2):", sum(abs(res$log2FoldChange) > 2, na.rm = TRUE), "\n")
cat("\nMedian log2 fold change:", median(res$log2FoldChange, na.rm = TRUE), "\n")
cat("Mean absolute log2 fold change:", mean(abs(res$log2FoldChange), na.rm = TRUE), "\n")


################################################################################
# Interpretation Notes
################################################################################

cat("\n=== Interpretation Guidelines ===\n")
cat("1. Log2 Fold Change Interpretation:\n")
cat("   LFC = 1  → 2-fold increase\n")
cat("   LFC = 2  → 4-fold increase\n")
cat("   LFC = -1 → 2-fold decrease\n")
cat("   LFC = -2 → 4-fold decrease\n\n")

cat("2. Statistical Significance:\n")
cat("   Use padj (FDR-adjusted), NOT pvalue\n")
cat("   Typical threshold: padj < 0.05\n")
cat("   Stringent threshold: padj < 0.01\n\n")

cat("3. Biological Significance:\n")
cat("   Consider both statistical significance AND effect size\n")
cat("   Small LFC with low padj → consistent but subtle change\n")
cat("   Large LFC with high padj → variable, may be unreliable\n")
cat("   Ideal: Large LFC with low padj → strong candidate\n\n")

cat("Analysis complete!\n")

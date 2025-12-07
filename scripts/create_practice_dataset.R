################################################################################
# Create Small Practice Dataset
# Integrative Multi-Omics Analysis Workshop
################################################################################

cat("========================================\n")
cat("Creating Practice Dataset\n")
cat("========================================\n\n")

# Set seed for reproducibility
set.seed(12345)

# Set working directory (adjust as needed)
# setwd("path/to/workshop_materials")

################################################################################
# 1. Create Small Count Matrix (10 samples, 500 genes)
################################################################################

cat("Step 1: Creating count matrix...\n")

# Sample information
n_samples <- 10
n_genes <- 500

sample_names <- paste0("Sample_", sprintf("%02d", 1:n_samples))
gene_names <- paste0("GENE_", sprintf("%04d", 1:n_genes))

# Assign groups (5 Control, 5 Treatment)
groups <- factor(rep(c("Control", "Treatment"), each = 5))

# Generate realistic count data
# Control samples: baseline expression
# Treatment samples: some genes up, some down

# Baseline mean expression (varies by gene)
base_mean <- rnbinom(n_genes, mu = 100, size = 10)

# Create count matrix
counts_practice <- matrix(0, nrow = n_genes, ncol = n_samples)
rownames(counts_practice) <- gene_names
colnames(counts_practice) <- sample_names

for (i in 1:n_genes) {
  for (j in 1:n_samples) {
    # Add some noise
    counts_practice[i, j] <- rnbinom(1,
                                     mu = base_mean[i],
                                     size = 10)
  }
}

# Add differential expression to some genes
# 50 up-regulated genes
up_genes <- 1:50
for (i in up_genes) {
  # 2-fold increase in Treatment
  counts_practice[i, 6:10] <- rnbinom(5,
                                      mu = base_mean[i] * 2,
                                      size = 10)
}

# 50 down-regulated genes
down_genes <- 51:100
for (i in down_genes) {
  # 2-fold decrease in Treatment
  counts_practice[i, 6:10] <- rnbinom(5,
                                      mu = base_mean[i] * 0.5,
                                      size = 10)
}

cat("✓ Created count matrix:", nrow(counts_practice), "genes ×",
    ncol(counts_practice), "samples\n\n")

################################################################################
# 2. Create Metadata
################################################################################

cat("Step 2: Creating metadata...\n")

metadata_practice <- data.frame(
  Sample = sample_names,
  Group = groups,
  Batch = factor(rep(c("A", "B"), times = 5)),
  row.names = sample_names
)

cat("✓ Created metadata\n")
print(metadata_practice)
cat("\n")

################################################################################
# 3. Create FPKM-like Normalized Data
################################################################################

cat("Step 3: Creating FPKM normalized data...\n")

# Simulate gene lengths
gene_lengths <- sample(500:5000, n_genes, replace = TRUE)

# Simple FPKM calculation
# FPKM = (counts * 1e9) / (gene_length * total_counts)

fpkm_practice <- matrix(0, nrow = n_genes, ncol = n_samples)
rownames(fpkm_practice) <- gene_names
colnames(fpkm_practice) <- sample_names

for (j in 1:n_samples) {
  total_counts <- sum(counts_practice[, j])
  for (i in 1:n_genes) {
    fpkm_practice[i, j] <- (counts_practice[i, j] * 1e9) /
                           (gene_lengths[i] * total_counts)
  }
}

cat("✓ Created FPKM matrix\n\n")

################################################################################
# 4. Create Gene Annotation
################################################################################

cat("Step 4: Creating gene annotation...\n")

# Create mapping to gene symbols
gene_symbols <- c(
  # Metabolic genes (first 100)
  paste0("MET", 1:100),
  # Other genes
  paste0("GENE", 101:500)
)

gene_annotation <- data.frame(
  GeneID = gene_names,
  Symbol = gene_symbols,
  Length = gene_lengths,
  row.names = gene_names
)

cat("✓ Created gene annotation\n\n")

################################################################################
# 5. Create Gene Set for GSEA (Metabolic Pathway)
################################################################################

cat("Step 5: Creating gene sets...\n")

# Create a simple GMT file with metabolic pathways
gene_sets <- list(
  "GLYCOLYSIS" = c("MET1", "MET2", "MET3", "MET4", "MET5",
                  "MET6", "MET7", "MET8", "MET9", "MET10"),
  "TCA_CYCLE" = c("MET11", "MET12", "MET13", "MET14", "MET15",
                 "MET16", "MET17", "MET18"),
  "OXIDATIVE_PHOSPHORYLATION" = c("MET19", "MET20", "MET21", "MET22",
                                 "MET23", "MET24", "MET25"),
  "FATTY_ACID_METABOLISM" = c("MET26", "MET27", "MET28", "MET29",
                             "MET30", "MET31", "MET32"),
  "AMINO_ACID_METABOLISM" = c("MET33", "MET34", "MET35", "MET36",
                             "MET37", "MET38", "MET39", "MET40"),
  "NUCLEOTIDE_METABOLISM" = c("MET41", "MET42", "MET43", "MET44",
                             "MET45", "MET46"),
  "GLUTAMINE_METABOLISM" = c("MET47", "MET48", "MET49", "MET50"),
  "LACTATE_METABOLISM" = c("MET1", "MET2", "MET3"),  # Overlap with glycolysis
  "PENTOSE_PHOSPHATE" = c("MET51", "MET52", "MET53", "MET54"),
  "PYRUVATE_METABOLISM" = c("MET1", "MET11", "MET12")  # Connects pathways
)

# Write GMT file
gmt_file <- "data/practice_pathways.gmt"
gmt_conn <- file(gmt_file, "w")

for (pathway_name in names(gene_sets)) {
  genes <- gene_sets[[pathway_name]]
  line <- paste(c(pathway_name, "Practice_pathway", genes), collapse = "\t")
  writeLines(line, gmt_conn)
}

close(gmt_conn)

cat("✓ Created", length(gene_sets), "gene sets in GMT format\n\n")

################################################################################
# 6. Create Reporter Metabolite Gene Sets
################################################################################

cat("Step 6: Creating metabolite-gene associations...\n")

metabolite_genes <- list(
  "Glucose" = c("MET1", "MET2", "MET3"),
  "Pyruvate" = c("MET1", "MET11", "MET12"),
  "Lactate" = c("MET1", "MET2"),
  "Acetyl_CoA" = c("MET11", "MET26", "MET27"),
  "Citrate" = c("MET11", "MET12", "MET13"),
  "Succinate" = c("MET14", "MET15"),
  "ATP" = c("MET19", "MET20", "MET21", "MET22"),
  "NADH" = c("MET19", "MET20", "MET11", "MET1"),
  "Glutamine" = c("MET47", "MET48", "MET49", "MET50"),
  "Fatty_Acids" = c("MET26", "MET27", "MET28", "MET29", "MET30")
)

# Save as RDS for easy loading
saveRDS(metabolite_genes, "data/practice_metabolite_genes.rds")

cat("✓ Created", length(metabolite_genes), "metabolite gene sets\n\n")

################################################################################
# 7. Save All Practice Files
################################################################################

cat("Step 7: Saving practice dataset files...\n")

# Create data directory if needed
if (!dir.exists("data")) {
  dir.create("data")
}

# Save count matrix
write.table(counts_practice,
           "data/practice_counts.txt",
           sep = "\t",
           quote = FALSE,
           col.names = NA)
cat("✓ Saved: practice_counts.txt\n")

# Save metadata
write.table(metadata_practice,
           "data/practice_metadata.txt",
           sep = "\t",
           quote = FALSE,
           col.names = NA)
cat("✓ Saved: practice_metadata.txt\n")

# Save FPKM
write.table(fpkm_practice,
           "data/practice_fpkm.txt",
           sep = "\t",
           quote = FALSE,
           col.names = NA)
cat("✓ Saved: practice_fpkm.txt\n")

# Save annotation
write.table(gene_annotation,
           "data/practice_gene_annotation.txt",
           sep = "\t",
           quote = FALSE,
           col.names = NA)
cat("✓ Saved: practice_gene_annotation.txt\n")

################################################################################
# 8. Create Quick Analysis Script
################################################################################

cat("\nStep 8: Creating quick analysis script...\n")

quick_analysis <- '
################################################################################
# Quick Practice Analysis
################################################################################

library(DESeq2)

# Load practice data
counts <- as.matrix(read.table("data/practice_counts.txt", sep = "\\t", row.names = 1))
metadata <- read.table("data/practice_metadata.txt", sep = "\\t", row.names = 1)

# Run DESeq2
dds <- DESeqDataSetFromMatrix(counts, metadata, design = ~ Group)
dds <- DESeq(dds)
res <- results(dds, contrast = c("Group", "Treatment", "Control"))

# View results
summary(res)
head(res[order(res$padj), ])

# How many significant genes?
cat("Significant genes (padj < 0.05):", sum(res$padj < 0.05, na.rm = TRUE), "\\n")
cat("Expected: ~100 (50 up, 50 down)\\n")
'

writeLines(quick_analysis, "scripts/quick_practice_analysis.R")
cat("✓ Saved: quick_practice_analysis.R\n\n")

################################################################################
# 9. Summary
################################################################################

cat("========================================\n")
cat("Practice Dataset Created!\n")
cat("========================================\n\n")

cat("Dataset Summary:\n")
cat("- Samples: 10 (5 Control, 5 Treatment)\n")
cat("- Genes: 500\n")
cat("- Differentially expressed: ~100 genes\n")
cat("  - Up-regulated: 50 genes\n")
cat("  - Down-regulated: 50 genes\n")
cat("- Pathways: 10 metabolic pathways\n")
cat("- Metabolites: 10 with gene associations\n\n")

cat("Files created:\n")
cat("✓ data/practice_counts.txt\n")
cat("✓ data/practice_metadata.txt\n")
cat("✓ data/practice_fpkm.txt\n")
cat("✓ data/practice_gene_annotation.txt\n")
cat("✓ data/practice_pathways.gmt\n")
cat("✓ data/practice_metabolite_genes.rds\n")
cat("✓ scripts/quick_practice_analysis.R\n\n")

cat("To use:\n")
cat("1. Run: source(\"scripts/quick_practice_analysis.R\")\n")
cat("2. Should find ~100 significant genes\n")
cat("3. Use for practice before working with full TCGA dataset\n\n")

cat("Practice dataset creation complete!\n")
'
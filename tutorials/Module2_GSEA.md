# Module 2: Gene Set Enrichment Analysis with PIANO

**Duration**: 45 minutes
**Difficulty**: Intermediate

---

## 📖 Learning Objectives

By the end of this module, you will be able to:
- Understand gene set enrichment analysis (GSEA) concepts
- Convert Ensembl IDs to gene symbols for pathway analysis
- Use PIANO to identify enriched metabolic pathways
- Interpret directional enrichment results
- Understand distinct-directional vs mixed-directional enrichment
- Connect pathway results to metabolic reprogramming

---

## 🧬 Biological Context

### What is Gene Set Enrichment Analysis?

**Gene Set Enrichment Analysis (GSEA)** asks: *Are genes in a specific pathway or biological process coordinately up- or down-regulated?*

**Example**:
- Individual gene analysis: Gene A (p=0.04), Gene B (p=0.06), Gene C (p=0.08)
- Each gene alone is marginally significant
- But all three are in "Glycolysis" pathway
- **GSEA conclusion**: Glycolysis pathway is significantly enriched!

### Why Use GSEA?

1. **Increased power**: Detects coordinated small changes
2. **Biological interpretation**: Identifies affected processes
3. **Mechanistic insights**: Suggests functional consequences
4. **Hypothesis generation**: Points to follow-up experiments

### GSEA vs Over-Representation Analysis (ORA)

| Feature | ORA | GSEA (PIANO) |
|---------|-----|--------------|
| **Input** | List of significant genes | All genes with statistics |
| **Threshold** | Requires p-value cutoff | Uses full distribution |
| **Direction** | Up/down separated | Directional enrichment |
| **Power** | Lower | Higher |
| **Example** | Fisher's exact test | Reporter genes method |

---

## 🔬 PIANO: Platform for Integrated Analysis of Omics data

### What is PIANO?

**PIANO** is a Bioconductor package for gene set analysis with several advantages:

- **Multiple methods**: Reporter genes, Fisher's exact, GSEA, etc.
- **Directional analysis**: Distinguishes up vs down vs mixed regulation
- **Permutation testing**: Robust significance assessment
- **Consensus scoring**: Combines multiple statistical approaches

### Key Concepts

**Gene Set Statistics**:
- **Distinct-directional up**: Genes mostly up-regulated
- **Distinct-directional down**: Genes mostly down-regulated
- **Non-directional**: Genes significantly changed (either direction)
- **Mixed-directional**: Both up and down in same pathway

**Example Interpretation**:
- Pathway with 10 genes, 9 up-regulated, 1 down → **Distinct-directional up**
- Pathway with 10 genes, 5 up, 5 down → **Mixed-directional**

---

## 💻 Step-by-Step Analysis

### Step 1: Load Required Libraries

```r
library(piano)
library(DESeq2)

# Set working directory (adjust to your path)
setwd("path/to/workshop_materials")
```

---

### Step 2: Load DESeq2 Results from Module 1

```r
# Load differential expression results
DESeq_results <- read.delim("data/DESeq output.txt",
                            row.names = 1,
                            stringsAsFactors = FALSE)

# Check the data
head(DESeq_results)
dim(DESeq_results)

# Summary
summary(DESeq_results)
```

**Expected output**:
```
                baseMean log2FoldChange     lfcSE      stat    pvalue      padj
ENSG00000000003  4367.04        -0.3824  0.20839 -1.835141 0.0664848 0.3105794

[1] 15234     6
```

---

### Step 3: Map Ensembl IDs to Gene Symbols

**Why?** MSigDB gene sets use gene symbols (like TP53), not Ensembl IDs (like ENSG00000141510)

```r
# Load gene mapping
ensembl2gene <- read.delim("data/Ensembl2gene.tsv",
                          row.names = 2,
                          stringsAsFactors = FALSE)

# Preview mapping
head(ensembl2gene)

# Add gene symbols to DESeq results
DESeq_results$Gene <- ensembl2gene[rownames(DESeq_results), "Gene"]

# Check how many genes have symbols
cat("Genes with symbols:", sum(!is.na(DESeq_results$Gene)), "\n")
cat("Genes without symbols:", sum(is.na(DESeq_results$Gene)), "\n")

# Remove genes without symbols
DESeq_with_genes <- DESeq_results[!is.na(DESeq_results$Gene), ]

cat("\nAfter filtering:", nrow(DESeq_with_genes), "genes\n")
```

**Expected output**:
```
                Gene
ENSG00000000003 TSPAN6
ENSG00000000005   TNMD

Genes with symbols: 14523
Genes without symbols: 711

After filtering: 14523 genes
```

---

### Step 4: Prepare Data for PIANO

PIANO requires:
1. **Gene statistics** (p-values and fold changes)
2. **Gene set collection** (GMT file)

```r
# Extract gene symbols as row names
rownames(DESeq_with_genes) <- DESeq_with_genes$Gene

# Keep only required columns
DESeq_for_piano <- DESeq_with_genes[, c('log2FoldChange', 'pvalue')]

# Handle NA values
# Set NA p-values to 1 (non-significant)
# Set NA fold changes to 0 (no change)
DESeq_for_piano$pvalue[is.na(DESeq_for_piano$pvalue)] <- 1
DESeq_for_piano$log2FoldChange[is.na(DESeq_for_piano$log2FoldChange)] <- 0

# Check for duplicates (some Ensembl IDs map to same symbol)
duplicated_genes <- duplicated(rownames(DESeq_for_piano))
cat("Duplicated gene symbols:", sum(duplicated_genes), "\n")

# Remove duplicates (keep first occurrence)
DESeq_for_piano <- DESeq_for_piano[!duplicated_genes, ]

cat("Final genes for PIANO:", nrow(DESeq_for_piano), "\n")
```

---

### Step 5: Create Matrices for PIANO

```r
# PIANO expects matrices with gene names as row names

# P-values matrix
pval_matrix <- as.matrix(DESeq_for_piano[, 'pvalue', drop = FALSE])
rownames(pval_matrix) <- rownames(DESeq_for_piano)
colnames(pval_matrix) <- "Late_vs_Early"

# Fold changes matrix (for directionality)
fc_matrix <- as.matrix(DESeq_for_piano[, 'log2FoldChange', drop = FALSE])
rownames(fc_matrix) <- rownames(DESeq_for_piano)
colnames(fc_matrix) <- "Late_vs_Early"

# Preview
head(pval_matrix)
head(fc_matrix)

cat("\nData prepared for PIANO!\n")
cat("Genes:", nrow(pval_matrix), "\n")
```

---

### Step 6: Load Gene Set Collection

```r
# Load MSigDB GO Biological Processes gene sets
gsc <- loadGSC("data/c5.bp.v6.2.symbols.gmt")

# Check gene set collection
cat("Number of gene sets:", length(gsc$gsc), "\n")

# Preview first few gene sets
cat("\nFirst 5 gene sets:\n")
names(gsc$gsc)[1:5]

# Check size distribution
gene_set_sizes <- sapply(gsc$gsc, length)
summary(gene_set_sizes)

# Plot gene set sizes
hist(gene_set_sizes,
     breaks = 50,
     main = "Gene Set Size Distribution",
     xlab = "Number of genes per set",
     col = "skyblue")
```

**Expected output**:
```
Number of gene sets: 5192

First 5 gene sets:
[1] "GO_POSITIVE_REGULATION_OF_VIRAL_TRANSCRIPTION"
[2] "GO_CARDIAC_CHAMBER_DEVELOPMENT"
[3] "GO_DNA_DEPENDENT_DNA_REPLICATION_MAINTENANCE_OF_FIDELITY"
[4] "GO_CIRCADIAN_RHYTHM"
[5] "GO_PHOSPHATIDYLSERINE_ACYL_CHAIN_REMODELING"

   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
   5.00   14.00   29.00   56.82   67.00  995.00
```

---

### Step 7: Run PIANO Gene Set Analysis

```r
# Run GSA with reporter genes method
cat("Running PIANO gene set analysis...\n")
cat("This may take a few minutes...\n\n")

gsaRes <- runGSA(
  geneLevelStats = pval_matrix,
  directions = fc_matrix,
  gsc = gsc,
  geneSetStat = "reporter",      # Reporter genes method (Patil & Nielsen)
  signifMethod = "geneSampling",  # Gene sampling for p-values
  nPerm = 1000,                   # Number of permutations
  gsSizeLim = c(5, 500),         # Gene set size limits (5-500 genes)
  adjMethod = "fdr"               # FDR correction
)

cat("\n✓ PIANO analysis complete!\n")
```

**Parameters explained**:

- `geneSetStat = "reporter"`: Uses reporter algorithm (aggregates p-values)
- `signifMethod = "geneSampling"`: Permutes genes to get null distribution
- `nPerm = 1000`: 1000 permutations (more = slower but more accurate)
- `gsSizeLim = c(5, 500)`: Only pathways with 5-500 genes
- `adjMethod = "fdr"`: False Discovery Rate correction

**This will take 2-5 minutes...**

---

### Step 8: View Results Summary

```r
# Get summary table
gsa_summary <- GSAsummaryTable(gsaRes, save = FALSE)

# View top pathways
head(gsa_summary, 20)

# Save full results
GSAsummaryTable(gsaRes,
               save = TRUE,
               file = "data/Piano_output.txt")

cat("\n✓ Results saved to: Piano_output.txt\n")
```

**Key columns in output**:

- `Name`: Gene set/pathway name
- `Genes (tot)`: Total genes in pathway
- `Stat (dist.dir.up)`: Statistic for distinct-directional up
- `Stat (dist.dir.dn)`: Statistic for distinct-directional down
- `p (dist.dir.up)`: P-value for up-regulation
- `p adj (dist.dir.up)`: FDR-adjusted p-value for up
- `Genes (up)`: Number of up-regulated genes
- `Genes (down)`: Number of down-regulated genes

---

### Step 9: Extract Significant Pathways

```r
# Extract pathways with FDR < 0.05

# Up-regulated pathways
sig_up <- gsa_summary[gsa_summary$`p adj (dist.dir.up)` < 0.05, ]
sig_up <- sig_up[order(sig_up$`p adj (dist.dir.up)`), ]

cat("Significantly UP-regulated pathways:", nrow(sig_up), "\n")
if (nrow(sig_up) > 0) {
  cat("\nTop 10 up-regulated:\n")
  print(head(sig_up[, c("Name", "Genes (tot)", "Genes (up)",
                        "p adj (dist.dir.up)")], 10))
}

# Down-regulated pathways
sig_down <- gsa_summary[gsa_summary$`p adj (dist.dir.dn)` < 0.05, ]
sig_down <- sig_down[order(sig_down$`p adj (dist.dir.dn)`), ]

cat("\n\nSignificantly DOWN-regulated pathways:", nrow(sig_down), "\n")
if (nrow(sig_down) > 0) {
  cat("\nTop 10 down-regulated:\n")
  print(head(sig_down[, c("Name", "Genes (tot)", "Genes (down)",
                          "p adj (dist.dir.dn)")], 10))
}
```

---

### Step 10: Visualize Top Pathways

```r
library(ggplot2)

# Prepare data for plotting
# Select top 10 up and top 10 down pathways

n_top <- 10

if (nrow(sig_up) > 0) {
  top_up <- head(sig_up, n_top)
  top_up$Direction <- "Up in Late"
  top_up$NegLogFDR <- -log10(top_up$`p adj (dist.dir.up)`)
}

if (nrow(sig_down) > 0) {
  top_down <- head(sig_down, n_top)
  top_down$Direction <- "Down in Late"
  top_down$NegLogFDR <- -log10(top_down$`p adj (dist.dir.dn)`)
  # Make negative for plotting
  top_down$NegLogFDR <- -top_down$NegLogFDR
}

# Combine for plotting
if (exists("top_up") && exists("top_down")) {
  plot_data <- rbind(
    data.frame(
      Pathway = top_up$Name,
      Score = top_up$NegLogFDR,
      Direction = top_up$Direction
    ),
    data.frame(
      Pathway = top_down$Name,
      Score = top_down$NegLogFDR,
      Direction = top_down$Direction
    )
  )

  # Simplify pathway names (remove GO_ prefix)
  plot_data$Pathway <- gsub("GO_", "", plot_data$Pathway)
  plot_data$Pathway <- gsub("_", " ", plot_data$Pathway)

  # Create barplot
  ggplot(plot_data, aes(x = reorder(Pathway, Score), y = Score, fill = Direction)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = c("Up in Late" = "red", "Down in Late" = "blue")) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8)) +
    labs(
      title = "Top Enriched Pathways (Late vs Early)",
      x = "Pathway",
      y = "-log10(FDR) × Direction",
      fill = ""
    ) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black")

  ggsave("figures/enriched_pathways.pdf", width = 10, height = 8)
  cat("\n✓ Plot saved to: figures/enriched_pathways.pdf\n")
}
```

---

## 📊 Interpreting Results

### Understanding Directional Enrichment

**Distinct-directional UP** (in Late stage):
- Most genes in pathway are up-regulated
- Suggests **increased pathway activity**
- Example: Cell cycle → increased proliferation

**Distinct-directional DOWN** (in Late stage):
- Most genes in pathway are down-regulated
- Suggests **decreased pathway activity**
- Example: Differentiation → loss of specialized functions

**Mixed-directional**:
- Both up and down genes in same pathway
- Suggests **complex regulation** or **feedback loops**
- Example: Apoptosis (pro- and anti-apoptotic genes)

---

### Biological Interpretation Example

**Finding**: "GLYCOLYSIS" pathway significantly up-regulated (FDR < 0.001)

**Interpretation**:
1. **Observation**: 8/10 glycolysis genes up-regulated in Late stage
2. **Mechanism**: Late stage tumors have enhanced glucose metabolism
3. **Context**: Warburg effect - cancer cells prefer glycolysis
4. **Consequence**: Produces lactate (see Module 4 - Reporter Metabolites!)
5. **Clinical**: Could target glycolysis therapeutically

---

### Connecting to Metabolic Context

Look for metabolic pathway patterns:

**Energy metabolism**:
- Glycolysis ↑ + Oxidative phosphorylation ↓ = Warburg effect
- TCA cycle ↑ = Oxidative metabolism
- Fatty acid oxidation ↑ = Alternative fuel source

**Biosynthesis**:
- Nucleotide synthesis ↑ = DNA replication (proliferation)
- Amino acid synthesis ↑ = Protein production
- Lipid synthesis ↑ = Membrane biosynthesis

**Redox balance**:
- Pentose phosphate pathway ↑ = NADPH production
- Glutathione metabolism ↑ = Antioxidant defense

---

## ✏️ Hands-On Exercises

### Exercise 1: Pathway Statistics

Using the PIANO results, answer:

1. How many pathways are significantly up-regulated (FDR < 0.05)?
2. How many are significantly down-regulated?
3. What is the most significantly enriched pathway overall?
4. How many genes are in this pathway?

**Code to help**:
```r
# Question 1
n_sig_up <- sum(gsa_summary$`p adj (dist.dir.up)` < 0.05)

# Question 2
n_sig_down <- sum(gsa_summary$`p adj (dist.dir.dn)` < 0.05)

# Question 3
min_pval_up <- min(gsa_summary$`p adj (dist.dir.up)`, na.rm = TRUE)
min_pval_down <- min(gsa_summary$`p adj (dist.dir.dn)`, na.rm = TRUE)
# Compare which is smaller

# Question 4
# Get pathway name with smallest p-value and check "Genes (tot)"
```

---

### Exercise 2: Focus on Metabolism

Filter results to find metabolic pathways:

1. Find all pathways containing "METABOL" in name
2. Which metabolic pathways are up-regulated?
3. Which are down-regulated?
4. What biological story do these tell?

**Code to help**:
```r
# Find metabolic pathways
metabolic <- gsa_summary[grep("METABOL", gsa_summary$Name), ]

# Separate by direction
metabolic_up <- metabolic[metabolic$`p adj (dist.dir.up)` < 0.05, ]
metabolic_down <- metabolic[metabolic$`p adj (dist.dir.dn)` < 0.05, ]

# View
print(metabolic_up[, c("Name", "Genes (tot)", "p adj (dist.dir.up)")])
print(metabolic_down[, c("Name", "Genes (tot)", "p adj (dist.dir.dn)")])
```

---

### Exercise 3: Gene-Level Investigation

Pick one significantly enriched pathway and investigate:

1. **Which specific genes drive the enrichment?**
2. **What are their fold changes?**
3. **Are they all changing in same direction?**

**Code to help**:
```r
# Pick a pathway (example: first significant up pathway)
pathway_name <- sig_up$Name[1]

# Get genes in this pathway
pathway_genes <- gsc$gsc[[pathway_name]]

# Get their statistics
pathway_stats <- DESeq_for_piano[pathway_genes, ]
pathway_stats <- pathway_stats[order(pathway_stats$log2FoldChange), ]

# View
print(pathway_stats)

# Visualize
barplot(pathway_stats$log2FoldChange,
       names.arg = rownames(pathway_stats),
       las = 2,
       col = ifelse(pathway_stats$pvalue < 0.05, "red", "gray"),
       main = paste("Gene Expression in", pathway_name),
       ylab = "Log2 Fold Change")
abline(h = 0)
```

---

## 🎯 Advanced Topics

### 1. Consensus Scoring

PIANO can combine multiple methods:

```r
# Run with multiple gene set statistics
gsaRes_consensus <- runGSA(
  geneLevelStats = pval_matrix,
  directions = fc_matrix,
  gsc = gsc,
  geneSetStat = c("mean", "median", "sum", "maxmean", "reporter"),
  signifMethod = "geneSampling",
  nPerm = 1000,
  gsSizeLim = c(5, 500)
)

# View consensus scores
consensus <- consensusScores(gsaRes_consensus, method = "mean")
```

---

### 2. Network View of Pathways

Create pathway-pathway similarity network:

```r
# Calculate overlap between significant pathways
sig_pathway_names <- c(sig_up$Name[1:10], sig_down$Name[1:10])

# Extract genes for each pathway
pathway_gene_lists <- gsc$gsc[sig_pathway_names]

# Calculate Jaccard similarity
# (For demonstration - full implementation would be more complex)
```

---

### 3. Compare with Practice Dataset

Run same analysis on practice data:

```r
# Load practice DESeq results (after running Module 1 on practice data)
# practice_results <- read.delim("practice_DESeq_output.txt")

# Load practice gene sets
# practice_gsc <- loadGSC("data/practice_pathways.gmt")

# Run PIANO on practice data
# Should find enriched metabolic pathways as designed
```

---

## 🎯 Key Takeaways

1. ✅ **GSEA** detects coordinated changes in gene sets
2. ✅ **More powerful** than individual gene analysis
3. ✅ **Directional enrichment** distinguishes up vs down regulation
4. ✅ **Reporter genes method** aggregates gene-level evidence
5. ✅ **FDR correction** essential for multiple testing
6. ✅ **Gene symbols required** for MSigDB gene sets
7. ✅ **Metabolic pathways** provide context for Module 4

---

## 🔗 Connecting to Other Modules

**From Module 1** (Differential Expression):
- Individual gene statistics → Input for GSEA
- P-values and fold changes used

**To Module 3** (Co-expression):
- Pathways suggest which genes to focus on
- Network analysis can reveal pathway structure

**To Module 4** (Reporter Metabolites):
- Enriched metabolic pathways point to key metabolites
- Glycolysis enrichment → Expect lactate/pyruvate as reporters
- Provides biological context for metabolic interpretation

---

## 📚 Additional Resources

### Papers
- **PIANO**: Väremo et al. (2013) Nucleic Acids Res - "Enriching the gene set analysis"
- **GSEA**: Subramanian et al. (2005) PNAS - "Gene set enrichment analysis"
- **Reporter genes**: Patil & Nielsen (2005) PNAS - "Transcriptional regulation of metabolism"

### Databases
- **MSigDB**: http://www.gsea-msigdb.org/gsea/msigdb/
- **GO**: http://geneontology.org/
- **KEGG**: https://www.genome.jp/kegg/

### Software
- **PIANO**: https://bioconductor.org/packages/release/bioc/html/piano.html
- **clusterProfiler**: Alternative R package for enrichment
- **Enrichr**: Web-based enrichment tool

---

## ⏭️ Next Module

Continue to **Module 3: Co-expression Network Analysis** to explore gene-gene relationships and identify metabolic gene modules.

**Solutions**: Available in `solutions/Module2_solutions.R`

---

**Great progress!** You now have:
- ✅ Differential expression (Module 1)
- ✅ Pathway enrichment (Module 2)
- Next: Network structure (Module 3)
- Then: Reporter metabolites (Module 4) ⭐

**Almost ready for the main focus on metabolic reprogramming!** 🔬

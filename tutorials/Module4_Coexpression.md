# Module 4: Advanced Co-expression Network Analysis

**Duration**: 30 minutes
**Difficulty**: Intermediate

---

## 📖 Learning Objectives

By the end of this module, you will be able to:
- Understand co-expression network concepts
- Calculate gene-gene correlation matrices
- Apply multiple hypothesis testing correction
- Filter networks for significant correlations
- Visualize correlation networks
- Export networks to Cytoscape
- Identify metabolic gene modules
- Interpret biological meaning of co-expression

---

## 🧬 Biological Context

### What is Co-expression?

**Co-expression** = Genes with similar expression patterns across samples

**Example**:
```
        Sample1  Sample2  Sample3  Sample4  Sample5
GeneA      10       15       20       12       18
GeneB       9       16       19       11       17
GeneC     100       50       30       90       40
```
- **GeneA and GeneB**: Highly correlated (co-expressed)
- **GeneC**: Different pattern (not co-expressed with A/B)

### Why Study Co-expression?

1. **Functional relationships**: Co-expressed genes often work together
2. **Regulatory networks**: Shared transcription factors
3. **Pathway discovery**: Identify genes in same pathway
4. **Novel gene function**: Predict function from network neighbors
5. **Module detection**: Find groups of coordinated genes

### Co-expression vs Physical Interaction

| Co-expression | Protein-Protein Interaction |
|---------------|----------------------------|
| mRNA level | Protein level |
| Temporal association | Physical binding |
| Indirect relationship | Direct interaction |
| High-throughput | Targeted |
| Hypothesis generating | Hypothesis testing |

---

## 🔬 Network Analysis Concepts

### Types of Correlations

**Pearson correlation**:
- Linear relationships
- Assumes normal distribution
- Sensitive to outliers
- Range: -1 to +1

**Spearman correlation**:
- Rank-based (monotonic relationships)
- Non-parametric
- Robust to outliers
- **Recommended for gene expression** ✓

**Correlation coefficient interpretation**:
- |r| > 0.7: Strong correlation
- |r| 0.4-0.7: Moderate correlation
- |r| < 0.4: Weak correlation

### Network Terminology

- **Node**: Gene
- **Edge**: Correlation between genes
- **Weight**: Correlation strength
- **Module/Cluster**: Groups of highly connected genes
- **Hub gene**: Highly connected node (many edges)

---

## 💻 Step-by-Step Analysis

### Step 1: Load Required Libraries

```r
library(corrplot)      # Correlation visualization
library(Hmisc)         # Correlation with p-values
library(reshape2)      # Data reshaping
library(ggplot2)       # Advanced plotting

# Set working directory
setwd("path/to/workshop_materials")
```

---

### Step 2: Load and Prepare Expression Data

We'll use FPKM data (not counts) because:
- Already normalized
- Better for correlation analysis
- Accounts for gene length and library size

```r
# Load FPKM expression data
FPKM_data <- read.csv("data/FPKM_selected.txt",
                      sep = "\t",
                      row.names = 1)

# Check data
dim(FPKM_data)
head(FPKM_data[, 1:5])

cat("Genes:", nrow(FPKM_data), "\n")
cat("Samples:", ncol(FPKM_data), "\n")
```

---

### Step 3: Select Genes for Analysis

**Strategy**: Focus on subset for demonstration
- Full dataset: ~20,000 genes → ~200 million pairwise correlations!
- We'll use 113 pre-selected genes for speed

```r
# Load selected genes
target_genes <- read.delim("data/Genes_selected.txt",
                          sep = "\t",
                          row.names = 1)

cat("Target genes:", nrow(target_genes), "\n")
head(target_genes)

# Subset FPKM data to target genes
FPKM_subset <- FPKM_data[rownames(target_genes), ]

dim(FPKM_subset)
```

---

### Step 4: Map to Gene Symbols

```r
# Load gene mapping
ensembl2gene <- read.delim("data/Ensembl2gene.tsv",
                          row.names = 2,
                          stringsAsFactors = FALSE)

# Map Ensembl IDs to gene symbols
gene_symbols <- ensembl2gene[rownames(FPKM_subset), "Gene"]

# Update row names
rownames(FPKM_subset) <- gene_symbols

# Check
head(rownames(FPKM_subset))
```

---

### Step 5: Filter Lowly Expressed Genes

**Why?** Lowly expressed genes have unreliable correlations

```r
# Calculate mean expression per gene
mean_expression <- rowMeans(FPKM_subset)

# Plot distribution
hist(log2(mean_expression + 1),
     breaks = 30,
     main = "Distribution of Mean Expression",
     xlab = "log2(Mean FPKM + 1)",
     col = "lightblue")
abline(v = log2(5 + 1), col = "red", lwd = 2, lty = 2)
legend("topright", "Cutoff: FPKM = 5", col = "red", lty = 2)

# Filter: Keep genes with mean FPKM > 5
FPKM_filtered <- FPKM_subset[mean_expression > 5, ]

cat("\nBefore filtering:", nrow(FPKM_subset), "genes\n")
cat("After filtering:", nrow(FPKM_filtered), "genes\n")
cat("Removed:", nrow(FPKM_subset) - nrow(FPKM_filtered), "lowly expressed genes\n")
```

---

### Step 6: Transpose Data for Correlation

**Important**: Correlation functions expect samples in rows, genes in columns

```r
# Transpose: genes in columns, samples in rows
FPKM_transposed <- t(as.matrix(FPKM_filtered))

dim(FPKM_transposed)
# Should show: [samples x genes]

# Check
head(FPKM_transposed[, 1:5])
```

---

### Step 7: Calculate Correlation Matrix with P-values

```r
cat("Calculating correlations...\n")
cat("This may take a minute...\n\n")

# Calculate Spearman correlation with p-values
cor_result <- rcorr(FPKM_transposed, type = "spearman")

# Extract correlation coefficients and p-values
cor_matrix <- cor_result$r
pval_matrix <- cor_result$P

cat("✓ Correlation matrix calculated\n")
cat("Dimensions:", dim(cor_matrix), "\n")

# Preview
cor_matrix[1:5, 1:5]
```

**What we get**:
- `cor_matrix`: Correlation coefficients (r values)
- `pval_matrix`: P-values for each correlation

---

### Step 8: Visualize Full Correlation Matrix

```r
# Create correlation plot
# Note: May be crowded with many genes

pdf("figures/correlation_matrix.pdf", width = 12, height = 12)

corrplot(cor_matrix,
         method = "square",           # Square cells
         type = "full",              # Show full matrix
         order = "hclust",           # Hierarchical clustering
         tl.cex = 0.5,               # Text label size
         tl.col = "black",           # Text color
         p.mat = pval_matrix,        # P-values for significance
         sig.level = 0.01,           # Significance threshold
         insig = "blank",            # Hide non-significant
         col = colorRampPalette(c("blue", "white", "red"))(200),
         addgrid.col = NA            # No grid
)

dev.off()

cat("✓ Correlation plot saved to: figures/correlation_matrix.pdf\n")
```

**Interpretation**:
- **Red**: Positive correlation (co-expressed)
- **Blue**: Negative correlation (oppositely expressed)
- **White/blank**: Not significant
- **Clusters**: Groups of co-expressed genes

---

### Step 9: Multiple Testing Correction

**Problem**: With N genes, we have N*(N-1)/2 pairwise tests!

Example: 100 genes = 4,950 tests → Many false positives

**Solution**: FDR correction

```r
# Extract upper triangle (avoid counting each pair twice)
upper_triangle_cor <- cor_matrix
upper_triangle_cor[lower.tri(upper_triangle_cor, diag = TRUE)] <- NA

upper_triangle_pval <- pval_matrix
upper_triangle_pval[lower.tri(upper_triangle_pval, diag = TRUE)] <- NA

# Reshape to long format
cor_long <- melt(upper_triangle_cor)
pval_long <- melt(upper_triangle_pval)

# Combine
correlation_data <- data.frame(
  Gene1 = cor_long$Var1,
  Gene2 = cor_long$Var2,
  Correlation = cor_long$value,
  Pvalue = pval_long$value
)

# Remove NA (diagonal and lower triangle)
correlation_data <- correlation_data[!is.na(correlation_data$Correlation), ]

# Remove self-correlations if any remain
correlation_data <- correlation_data[
  as.character(correlation_data$Gene1) != as.character(correlation_data$Gene2),
]

cat("Total pairwise correlations:", nrow(correlation_data), "\n")

# Apply FDR correction
correlation_data$FDR <- p.adjust(correlation_data$Pvalue, method = "fdr")

# Check how many significant
sig_05 <- sum(correlation_data$FDR < 0.05, na.rm = TRUE)
sig_01 <- sum(correlation_data$FDR < 0.01, na.rm = TRUE)

cat("Significant correlations (FDR < 0.05):", sig_05, "\n")
cat("Significant correlations (FDR < 0.01):", sig_01, "\n")
```

---

### Step 10: Filter for Significant Correlations

```r
# Keep only significant correlations
significant_cors <- correlation_data[
  !is.na(correlation_data$FDR) & correlation_data$FDR < 0.05,
]

# Order by absolute correlation strength
significant_cors <- significant_cors[
  order(-abs(significant_cors$Correlation)),
]

cat("\nSignificant correlations:\n")
cat("Total:", nrow(significant_cors), "\n")
cat("Positive:", sum(significant_cors$Correlation > 0), "\n")
cat("Negative:", sum(significant_cors$Correlation < 0), "\n\n")

# View top correlations
cat("Top 10 strongest correlations:\n")
print(head(significant_cors, 10))
```

---

### Step 11: Save Network for Cytoscape

```r
# Prepare edge list for Cytoscape
# Format: Source, Target, Correlation, P-value, FDR, CorrelationType

cytoscape_network <- data.frame(
  Source = as.character(significant_cors$Gene1),
  Target = as.character(significant_cors$Gene2),
  Correlation = significant_cors$Correlation,
  Pvalue = significant_cors$Pvalue,
  FDR = significant_cors$FDR,
  Type = ifelse(significant_cors$Correlation > 0, "Positive", "Negative"),
  stringsAsFactors = FALSE
)

# Save
write.table(cytoscape_network,
           "data/Coexpression_network.txt",
           sep = "\t",
           row.names = FALSE,
           quote = FALSE)

cat("✓ Network saved for Cytoscape: Coexpression_network.txt\n")
cat("  Nodes:", length(unique(c(cytoscape_network$Source,
                               cytoscape_network$Target))), "\n")
cat("  Edges:", nrow(cytoscape_network), "\n")
```

---

### Step 12: Create Simple Network Visualization in R

```r
library(igraph)  # If available

# Create network visualization
# Note: Install igraph if not available: install.packages("igraph")

if (requireNamespace("igraph", quietly = TRUE)) {
  library(igraph)

  # Create graph from edge list
  g <- graph_from_data_frame(
    cytoscape_network[, c("Source", "Target", "Correlation")],
    directed = FALSE
  )

  # Set edge colors
  E(g)$color <- ifelse(E(g)$Correlation > 0, "red", "blue")

  # Set edge width by correlation strength
  E(g)$width <- abs(E(g)$Correlation) * 5

  # Plot
  pdf("figures/coexpression_network.pdf", width = 12, height = 12)

  plot(g,
       vertex.size = 5,
       vertex.label.cex = 0.6,
       vertex.label.color = "black",
       vertex.color = "lightblue",
       edge.arrow.size = 0,
       layout = layout_with_fr(g),  # Fruchterman-Reingold layout
       main = "Co-expression Network (FDR < 0.05)")

  legend("bottomright",
         legend = c("Positive correlation", "Negative correlation"),
         col = c("red", "blue"),
         lty = 1,
         lwd = 2)

  dev.off()

  cat("✓ Network plot saved: figures/coexpression_network.pdf\n")

} else {
  cat("Note: igraph not available. Install for network plots.\n")
  cat("      install.packages('igraph')\n")
}
```

---

### Step 13: Identify Network Modules

```r
# Detect communities/modules in network

if (requireNamespace("igraph", quietly = TRUE)) {

  # Community detection using Louvain algorithm
  communities <- cluster_louvain(g)

  cat("\nNetwork modules detected:", length(communities), "\n")

  # Get module membership
  module_membership <- data.frame(
    Gene = V(g)$name,
    Module = membership(communities)
  )

  # How many genes per module?
  module_sizes <- table(module_membership$Module)
  print(module_sizes)

  # Save module information
  write.table(module_membership,
             "data/network_modules.txt",
             sep = "\t",
             row.names = FALSE,
             quote = FALSE)

  # Visualize with module colors
  pdf("figures/coexpression_modules.pdf", width = 12, height = 12)

  plot(communities, g,
       vertex.size = 5,
       vertex.label.cex = 0.6,
       edge.arrow.size = 0,
       main = "Co-expression Network Modules")

  dev.off()

  cat("✓ Modules saved: data/network_modules.txt\n")
  cat("✓ Module plot saved: figures/coexpression_modules.pdf\n")
}
```

---

## 📊 Interpreting Results

### Understanding Correlations

**High positive correlation (r > 0.7)**:
- Genes likely in same pathway
- May share transcription factors
- Could be functionally related
- Example: Glycolysis enzymes

**High negative correlation (r < -0.7)**:
- Genes in opposing processes
- Mutual exclusion
- Example: Proliferation vs differentiation genes

**Moderate correlation (0.4 < |r| < 0.7)**:
- Indirect relationship
- Different but related pathways
- Shared upstream regulators

---

### Biological Interpretation

**Example Module**:
Genes: LDHA, LDHB, PKM, HK2, SLC2A1

**Analysis**:
1. **Function**: All involved in glycolysis
2. **Co-expression**: r > 0.8 for all pairs
3. **Biological meaning**: Coordinated glycolytic program
4. **Clinical relevance**: Warburg effect in cancer
5. **Follow-up**: Check if lactate is reporter metabolite (Module 3!)

---

### Hub Genes

**Hub gene** = Highly connected gene (many significant correlations)

```r
# Calculate node degree (number of connections)
if (exists("g")) {
  node_degrees <- degree(g)

  # Top hub genes
  top_hubs <- sort(node_degrees, decreasing = TRUE)[1:10]

  cat("\nTop 10 hub genes:\n")
  print(top_hubs)

  # Hub genes often:
  # - Central to biological process
  # - Potential drug targets
  # - Master regulators
}
```

---

## ✏️ Hands-On Exercises

### Exercise 1: Correlation Statistics

Using the correlation results:

1. **What is the strongest positive correlation?**
2. **What is the strongest negative correlation?**
3. **What percentage of all correlations are significant (FDR < 0.05)?**
4. **What is the median absolute correlation for significant pairs?**

**Code to help**:
```r
# Question 1
max_pos <- significant_cors[which.max(significant_cors$Correlation), ]
print(max_pos)

# Question 2
min_neg <- significant_cors[which.min(significant_cors$Correlation), ]
print(min_neg)

# Question 3
percent_sig <- (nrow(significant_cors) / nrow(correlation_data)) * 100
cat("Percentage significant:", percent_sig, "%\n")

# Question 4
median_abs_cor <- median(abs(significant_cors$Correlation))
print(median_abs_cor)
```

---

### Exercise 2: Metabolic Gene Focus

Focus on metabolic genes from enriched pathways (Module 2):

1. **Find correlations involving glycolysis genes**
2. **Are they positively or negatively correlated?**
3. **Do they form a tight module?**

**Code to help**:
```r
# Define glycolysis genes (examples)
glycolysis_genes <- c("HK1", "HK2", "GPI", "PFKM", "ALDOA",
                     "TPI1", "GAPDH", "PGK1", "PKM", "LDHA", "LDHB")

# Find in our network
glycolysis_in_network <- glycolysis_genes[glycolysis_genes %in%
                                         rownames(FPKM_filtered)]

cat("Glycolysis genes in network:", glycolysis_in_network, "\n")

# Extract correlations between glycolysis genes
glycolysis_cors <- significant_cors[
  significant_cors$Gene1 %in% glycolysis_in_network &
  significant_cors$Gene2 %in% glycolysis_in_network,
]

print(glycolysis_cors)

# Are they mostly positive?
table(glycolysis_cors$Type)
```

---

### Exercise 3: Module Enrichment

For each module, check if enriched for specific functions:

1. **Extract genes in largest module**
2. **Look up their functions** (use Google/GeneCards)
3. **Identify common theme**

**Code to help**:
```r
if (exists("module_membership")) {
  # Find largest module
  largest_module <- names(which.max(table(module_membership$Module)))

  # Get genes
  module_genes <- module_membership[
    module_membership$Module == largest_module,
    "Gene"
  ]

  cat("Largest module (", largest_module, "):\n", sep = "")
  print(module_genes)

  # You can then search these genes online
}
```

---

## 🎯 Advanced Topics

### 1. Partial Correlations

Account for confounding:

```r
# Partial correlation removes effect of confounders
# Requires ppcor package
# install.packages("ppcor")

if (requireNamespace("ppcor", quietly = TRUE)) {
  library(ppcor)

  # Calculate partial correlations
  # (accounts for all other genes)
  # WARNING: Computationally intensive
}
```

---

### 2. Differential Co-expression

Compare networks between conditions:

```r
# Build separate networks for Early vs Late
# Compare edge differences
# Genes that change correlation = rewired regulation
```

---

### 3. Weighted Gene Co-expression Network Analysis (WGCNA)

More sophisticated approach:

```r
# WGCNA package
# Soft thresholding instead of hard cutoff
# Automatic module detection
# install.packages("WGCNA")
```

---

## 🎯 Key Takeaways

1. ✅ **Co-expression** = Similar expression patterns
2. ✅ **Spearman correlation** best for gene expression
3. ✅ **FDR correction** essential for multiple testing
4. ✅ **Network modules** reveal functional groups
5. ✅ **Hub genes** are central regulators
6. ✅ **Export to Cytoscape** for advanced visualization
7. ✅ **Metabolic modules** connect to reporter metabolites (Module 4)

---

## 🔗 Connecting to Other Modules

**From Module 2** (GSEA):
- Enriched pathways suggest gene groups
- Co-expression confirms pathway coherence
- Example: Glycolysis genes co-expressed → validates enrichment

**To Module 3** (Reporter Metabolites):
- Co-expressed metabolic genes → shared metabolites
- Network modules → metabolic subsystems
- Hub genes → key metabolic enzymes
- Provides gene-level context for metabolite-level findings

**Network + Pathways + Metabolites** = Complete picture!

---

## 📚 Additional Resources

### Papers
- **Co-expression**: Stuart et al. (2003) Science - "A gene-coexpression network"
- **WGCNA**: Langfelder & Horvath (2008) BMC Bioinformatics
- **Network biology**: Barabási & Oltvai (2004) Nature Reviews Genetics

### Software
- **Cytoscape**: https://cytoscape.org/
- **WGCNA**: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/
- **igraph**: https://igraph.org/

### Databases
- **STRING**: Protein interaction/co-expression database
- **GeneMANIA**: Gene function prediction from networks

---

## 🔄 Cytoscape Import Guide

**Optional**: Visualize in Cytoscape for publication-quality networks

1. **Open Cytoscape**
2. **Import Network**:
   - File → Import → Network from File
   - Select: `Coexpression_network.txt`
   - Source = Source, Target = Target

3. **Import Edge Attributes**:
   - Correlation, FDR already included

4. **Style Network**:
   - Edge color by Type (red=positive, blue=negative)
   - Edge width by absolute Correlation
   - Node size by degree

5. **Layout**:
   - Layout → Prefuse Force Directed Layout

6. **Export**:
   - File → Export → Network as Image

---

## ⏭️ Next Module

**Note**: In the 3-hour workshop, Module 3 (Reporter Metabolites) comes BEFORE Module 4 (this module).

This tutorial shows the advanced co-expression analysis that builds upon:
- Differential expression (Module 1)
- Enriched pathways (Module 2)
- Reporter metabolites (Module 3) ⭐

Together, these modules provide a **complete picture of metabolic reprogramming**!

**Solutions**: Available in `solutions/Module3_solutions.R`

---

**Excellent progress!** You now understand:
- ✅ Individual gene changes (Module 1)
- ✅ Pathway-level changes (Module 2)
- ✅ Gene-gene relationships (Module 3)
- Ready for: **Metabolic network integration (Module 4)** 🎯

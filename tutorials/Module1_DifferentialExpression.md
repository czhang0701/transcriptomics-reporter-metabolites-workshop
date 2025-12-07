# Module 1: Differential Expression Analysis with DESeq2

**Duration**: 45 minutes
**Difficulty**: Intermediate

---

## 📖 Learning Objectives

By the end of this module, you will be able to:
- Load and explore RNA-seq count data and metadata
- Create a DESeq2 dataset object
- Perform differential expression analysis
- Interpret results and identify significant genes
- Generate quality control visualizations

---

## 🧬 Biological Context

We are analyzing RNA-seq data from **TCGA tumor samples** comparing:
- **Early stage** tumors (better prognosis)
- **Late stage** tumors (worse prognosis)

**Research Question**: Which genes are differentially expressed between early and late stage tumors?

---

## 📊 Understanding the Data

### Input Files

1. **Counts_selected.txt**
   - Gene expression counts (raw, unnormalized)
   - Rows: Genes (Ensembl IDs)
   - Columns: Samples (TCGA barcodes)
   - Values: Integer counts from RNA-seq

2. **patients.txt**
   - Sample metadata
   - Column "Type": Early or Late stage classification

### Why Use Counts?

DESeq2 requires **raw counts** (not normalized values like FPKM or TPM) because:
- It models count data using negative binomial distribution
- It performs its own normalization based on library size
- Raw counts preserve the statistical properties needed for hypothesis testing

---

## 💻 Step-by-Step Analysis

### Step 1: Load Required Libraries

```r
# Load DESeq2 for differential expression analysis
library(DESeq2)

# Set working directory (adjust to your path)
setwd("C:/work/For course/MSc Microbiome KCL 2025/Integerative_Rui/Materials")
```

**What does this do?**
- Loads the DESeq2 package with all necessary functions
- Sets the working directory to access data files

---

### Step 2: Import Count Data and Metadata

```r
# Import counts data
Counts <- as.matrix(read.csv(file = "Counts_selected.txt",
                             sep = "\t",
                             row.names = 1))

# Import metadata
Metadata <- read.csv(file = "patients.txt",
                    sep = "\t",
                    row.names = 1)

# Preview the data
head(Counts[, 1:5])  # First 5 samples
head(Metadata)
```

**Expected Output:**
```
                    TCGA.MR.A8JO.01A TCGA.NI.A8LF.01A TCGA.RC.A6M6.01A
ENSG00000000003              3699              5548              3413
ENSG00000000005                 2                31                 0
ENSG00000000419               734               813              2200

                     Type
TCGA.MR.A8JO.01A     Late
TCGA.NI.A8LF.01A     Late
TCGA.RC.A6M6.01A     Late
```

**Key Points:**
- Counts are integers (not decimals)
- Sample IDs match between counts and metadata
- Metadata has sample IDs as row names

---

### Step 3: Explore the Data Structure

```r
# Check dimensions
dim(Counts)      # Number of genes × samples
dim(Metadata)    # Number of samples × variables

# Check sample names match
all(colnames(Counts) == rownames(Metadata))

# Summary of metadata
table(Metadata$Type)
```

**Expected Output:**
```
[1] 20000   50    # 20,000 genes, 50 samples
[1] 50   1        # 50 samples, 1 variable
[1] TRUE          # Sample names match!

Early  Late
  25    25       # Balanced design
```

**Why is this important?**
- Ensures samples are correctly matched between counts and metadata
- Confirms balanced experimental design
- Identifies the comparison groups

---

### Step 4: Create DESeq2 Dataset Object

```r
# Create DESeqDataSet object
deSeqData <- DESeqDataSetFromMatrix(
  countData = Counts,
  colData = Metadata,
  design = ~ Type
)

# View the object
deSeqData
```

**Parameters Explained:**
- `countData`: Matrix of gene counts
- `colData`: Metadata with sample information
- `design = ~ Type`: Statistical model formula
  - `~` means "modeled by"
  - `Type` is the variable of interest (Early vs Late)

**Expected Output:**
```
class: DESeqDataSet
dim: 20000 50
metadata(1): version
assays(1): counts
rownames(20000): ENSG00000000003 ENSG00000000005 ...
rowData names(0):
colnames(50): TCGA.MR.A8JO.01A TCGA.NI.A8LF.01A ...
colData names(1): Type
```

---

### Step 5: Pre-filtering (Optional but Recommended)

```r
# Remove genes with very low counts
# Keep genes with at least 10 counts in at least 5 samples
keep <- rowSums(counts(deSeqData) >= 10) >= 5
deSeqData <- deSeqData[keep,]

# Check how many genes remain
nrow(deSeqData)
```

**Why filter?**
- Reduces multiple testing burden
- Improves statistical power
- Removes noise from lowly expressed genes
- Speeds up analysis

---

### Step 6: Run Differential Expression Analysis

```r
# Run the DESeq2 pipeline
deSeqAnalysis <- DESeq(deSeqData)
```

**What happens internally?**
1. **Normalization**: Calculates size factors to account for library size
2. **Dispersion estimation**: Models gene-wise variability
3. **Statistical testing**: Fits negative binomial model for each gene

**This may take a few minutes...**

---

### Step 7: Extract Results

```r
# Extract results for Late vs Early comparison
res <- results(
  deSeqAnalysis,
  contrast = c("Type", "Late", "Early"),  # Compare Late to Early
  lfcThreshold = 0,                       # No fold-change threshold
  altHypothesis = "greaterAbs",          # Two-sided test
  alpha = 0.05                            # FDR threshold
)

# View summary
summary(res)

# Preview results
head(res)
```

**Parameters Explained:**
- `contrast = c("Type", "Late", "Early")`:
  - Variable name: "Type"
  - Numerator: "Late" (test condition)
  - Denominator: "Early" (reference condition)
  - **Positive log2FC** = higher in Late
  - **Negative log2FC** = higher in Early

- `lfcThreshold = 0`: Test against LFC = 0 (any change)
- `altHypothesis = "greaterAbs"`: Two-sided test (up OR down)
- `alpha = 0.05`: FDR cutoff for significance

**Expected Output:**
```
out of 15000 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 1200, 8%
LFC < 0 (down)     : 1100, 7.3%
outliers [1]       : 0, 0%
low counts [2]     : 3000, 20%

                baseMean log2FoldChange     lfcSE      stat    pvalue      padj
ENSG00000000003  4367.04        -0.3824  0.20839 -1.835141 0.0664848 0.3105794
```

**Column Definitions:**
- `baseMean`: Average normalized counts across all samples
- `log2FoldChange`: Log2(Late/Early) - **key result!**
- `lfcSE`: Standard error of log2FC
- `stat`: Test statistic
- `pvalue`: Unadjusted p-value
- `padj`: **Adjusted p-value (FDR)** - **use this for significance!**

---

### Step 8: Save Results

```r
# Save complete results to file
write.table(
  as.data.frame(res),
  file = "DESeq_output.txt",
  sep = "\t",
  row.names = TRUE,
  col.names = NA
)
```

---

## 📈 Quality Control Visualizations

### 1. MA Plot (Mean vs Log Fold Change)

```r
# Create MA plot
plotMA(res, ylim = c(-5, 5), alpha = 0.05)
```

**Interpretation:**
- **X-axis**: Mean expression level
- **Y-axis**: Log2 fold change
- **Blue points**: Significant genes (padj < 0.05)
- **Gray points**: Non-significant genes
- **Red triangles**: Points outside y-axis limits

**What to look for:**
- Even distribution above/below zero
- More variability at low expression (left side) - normal!
- No systematic bias

---

### 2. Volcano Plot (Fold Change vs Significance)

```r
# Create volcano plot
library(ggplot2)

# Prepare data
res_df <- as.data.frame(res)
res_df$significant <- res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1

ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = significant), alpha = 0.5) +
  scale_color_manual(values = c("gray", "red")) +
  theme_minimal() +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title = "Volcano Plot: Late vs Early",
       x = "Log2 Fold Change",
       y = "-Log10 P-value")
```

**Interpretation:**
- **X-axis**: Magnitude of change
- **Y-axis**: Statistical significance
- **Top corners**: Most significant differential genes
- **Red points**: Significant with |LFC| > 1

---

### 3. Sample Distance Heatmap

```r
# Transform data for visualization
vsd <- vst(deSeqData, blind = FALSE)

# Calculate sample distances
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)

# Plot heatmap
library(pheatmap)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         annotation_col = Metadata)
```

**Interpretation:**
- Darker colors = more similar samples
- Samples should cluster by Type (Early/Late)
- Identifies outlier samples

---

## 🔍 Interpreting Results

### How many genes are differentially expressed?

```r
# Count significant genes
sum(res$padj < 0.05, na.rm = TRUE)

# Up-regulated in Late stage
sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE)

# Down-regulated in Late stage
sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm = TRUE)
```

### Find the most significant genes

```r
# Order by adjusted p-value
res_ordered <- res[order(res$padj), ]
head(res_ordered, 10)
```

### Find the largest fold changes

```r
# Order by absolute fold change
res_ordered_fc <- res[order(-abs(res$log2FoldChange)), ]
head(res_ordered_fc, 10)
```

---

## ✏️ Hands-On Exercise

**Task**: Answer the following questions using the DESeq2 results:

1. How many genes have padj < 0.05?
2. How many genes have padj < 0.05 AND |log2FoldChange| > 1?
3. What is the most up-regulated gene in Late stage?
4. What is the most down-regulated gene in Late stage?
5. Find a gene with high significance (padj < 0.001) but small fold change (|LFC| < 0.5)

**Bonus Challenge**:
Create a subset of results containing only:
- Significant genes (padj < 0.05)
- Strong effect (|log2FoldChange| > 1.5)
- Well-expressed genes (baseMean > 100)

Save this subset to a new file.

---

## 🎯 Key Takeaways

1. ✅ DESeq2 requires **raw counts**, not normalized data
2. ✅ Sample IDs must **match exactly** between counts and metadata
3. ✅ Use **padj**, not pvalue, for determining significance
4. ✅ **Log2FoldChange** direction:
   - Positive = higher in Late (numerator)
   - Negative = higher in Early (denominator)
5. ✅ Pre-filtering low-count genes improves power
6. ✅ Quality control plots are essential before interpretation

---

## 📚 Additional Resources

- [DESeq2 Paper (Love et al., 2014)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)
- [DESeq2 Vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
- [Understanding padj vs pvalue](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2907892/)

---

## ⏭️ Next Module

Continue to **Module 2: Gene Set Enrichment Analysis** to understand which biological pathways are affected by these gene expression changes.

Solutions to exercises are available in `solutions/Module1_solutions.R`

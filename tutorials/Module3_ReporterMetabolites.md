# Module 3: Reporter Metabolite Analysis ⭐

**Duration**: 60 minutes
**Difficulty**: Intermediate-Advanced

---

## 📖 Learning Objectives

By the end of this module, you will be able to:
- Understand the concept of reporter metabolites
- Link gene expression changes to metabolic networks
- Use genome-scale metabolic models (GEMs) to interpret transcriptomics
- Identify key metabolic reprogramming events
- Interpret reporter metabolite results in biological context

---

## 🧬 Biological Context

### What are Reporter Metabolites?

**Reporter metabolites** are metabolites whose surrounding genes in a metabolic network show coordinated differential expression.

**Concept**:
- Metabolic reactions are catalyzed by enzymes (proteins)
- Enzymes are encoded by genes
- If many genes associated with a metabolite are differentially expressed...
- → That metabolite is likely a "reporter" of metabolic reprogramming

### Why is This Useful?

1. **Mechanistic insight**: Links gene expression to metabolism
2. **Upstream regulators**: Identifies metabolites controlling pathways
3. **Drug targets**: Highlights metabolic vulnerabilities
4. **Biomarkers**: Suggests metabolites to measure experimentally

---

## 🔬 The Reporter Metabolite Method

### Step-by-Step Logic

1. **Input**: Differential expression results (gene-level p-values)
2. **Network**: Genome-scale metabolic model (gene-metabolite associations)
3. **Aggregation**: For each metabolite, collect p-values of associated genes
4. **Statistics**: Combine gene p-values into metabolite-level score
5. **Output**: Ranked list of reporter metabolites

### Statistical Method

**Reporter Score** (Patil & Nielsen, 2005):

For each metabolite *m*:
- Collect p-values of genes encoding enzymes that produce/consume *m*
- Convert p-values to Z-scores: `Z = Φ⁻¹(1 - p)`
- Aggregate Z-scores: `Z_agg = mean(Z) / sqrt(k)` where k = number of genes
- Correct for background using random metabolites

**Interpretation**:
- **Positive reporter score** → Metabolite's genes are significantly changed
- **High absolute score** → Strong evidence for metabolic reprogramming
- **Sign of fold changes** → Direction of regulation

---

## 💻 Analysis Workflow

### Approach 1: Using PIANO (R-based)

PIANO can perform reporter metabolite analysis using SBML models.

#### Step 1: Load Required Libraries

```r
library(piano)
library(DESeq2)

# Optional: for SBML parsing
# library(rsbml)  # Can be difficult to install

# Set working directory
setwd("path/to/workshop_materials")
```

---

#### Step 2: Prepare Differential Expression Results

```r
# Load DESeq2 results from Module 1
DESeq_results <- read.delim("data/DESeq output.txt",
                            row.names = 1,
                            stringsAsFactors = FALSE)

# Check data
head(DESeq_results)
dim(DESeq_results)

# Remove genes with NA values
DESeq_results_clean <- DESeq_results[!is.na(DESeq_results$log2FoldChange), ]

# Keep only pvalue and log2FoldChange columns
DESeq_for_reporter <- DESeq_results_clean[, c('pvalue', 'log2FoldChange')]

dim(DESeq_for_reporter)
```

**Expected output**:
```
                baseMean log2FoldChange     lfcSE      stat    pvalue      padj
ENSG00000000003  4367.04        -0.3824  0.20839 -1.835141 0.0664848 0.3105794

[1] 15234     6  # Original dimensions
[1] 15234     2  # After filtering to pvalue and LFC
```

---

#### Step 3: Load Gene-Metabolite Associations from SBML Model

```r
# Option 1: Load from SBML file (requires rsbml - can be tricky)
# This step is demonstrated but may not work on all systems

# Check if SBML model exists
if (file.exists("data/Reference_model.xml")) {
  cat("✓ SBML model found\n")

  # Try to load gene sets from SBML
  # NOTE: This may fail if rsbml is not properly installed
  tryCatch({
    # Load gene set collection from SBML
    metabolite_gsc <- loadGSC("data/Reference_model.xml")

    cat("✓ Loaded", length(metabolite_gsc$gsc), "metabolite gene sets\n")

    # Preview first gene set
    names(metabolite_gsc$gsc)[1:5]

  }, error = function(e) {
    cat("✗ Could not load SBML model\n")
    cat("Error:", e$message, "\n")
    cat("Will use alternative approach...\n")
  })
} else {
  cat("✗ SBML model not found\n")
}
```

---

#### Step 4: Alternative Approach - Pre-computed Gene Sets

Since SBML parsing can be challenging, we'll use a pre-processed gene set format:

```r
# For this workshop, we'll create a simplified example
# In practice, you would extract this from the SBML model

# Example: Create gene sets for a few metabolites
# Format: List where each element is a metabolite with associated genes

# This would normally come from parsing Reference_model.xml
# Here we demonstrate the concept with example data

example_metabolite_genes <- list(
  "Pyruvate" = c("LDHA", "LDHB", "PKM", "PDHA1", "PDHB"),
  "Lactate" = c("LDHA", "LDHB", "SLC16A1", "SLC16A3"),
  "Glucose" = c("SLC2A1", "SLC2A3", "HK1", "HK2", "GCK"),
  "Citrate" = c("CS", "ACLY", "SLC25A1"),
  "Glutamine" = c("GLUL", "GLS", "GLS2", "SLC1A5", "SLC38A1"),
  "Fatty_Acids" = c("FASN", "ACACA", "SCD", "ACSL1", "ACSL4"),
  "ATP" = c("ATP5A1", "ATP5B", "ATP5C1", "ATP5D", "ATP5E"),
  "NAD" = c("NAMPT", "NMNAT1", "NMNAT2", "NMNAT3"),
  "Succinate" = c("SDHA", "SDHB", "SDHC", "SDHD", "SUCLA2"),
  "Acetyl_CoA" = c("ACLY", "ACSS2", "PDHA1", "PDHB")
)

cat("Created", length(example_metabolite_genes), "example metabolite gene sets\n")
```

---

#### Step 5: Map Ensembl IDs to Gene Symbols

```r
# Load gene ID mapping
ensembl2gene <- read.delim("data/Ensembl2gene.tsv",
                          row.names = 2,
                          stringsAsFactors = FALSE)

# Map DESeq results to gene symbols
DESeq_for_reporter$Gene <- ensembl2gene[rownames(DESeq_for_reporter), "Gene"]

# Remove genes without symbols
DESeq_for_reporter <- DESeq_for_reporter[!is.na(DESeq_for_reporter$Gene), ]

# Set gene symbols as row names
rownames(DESeq_for_reporter) <- DESeq_for_reporter$Gene
DESeq_for_reporter <- DESeq_for_reporter[, c('pvalue', 'log2FoldChange')]

cat("Mapped", nrow(DESeq_for_reporter), "genes to symbols\n")
head(DESeq_for_reporter)
```

---

#### Step 6: Run Reporter Metabolite Analysis

```r
# Prepare data for PIANO
pval_vector <- as.matrix(DESeq_for_reporter[, 'pvalue'])
fc_vector <- as.matrix(DESeq_for_reporter[, 'log2FoldChange'])

rownames(pval_vector) <- rownames(DESeq_for_reporter)
rownames(fc_vector) <- rownames(DESeq_for_reporter)

# Replace any NA values
pval_vector[is.na(pval_vector)] <- 1
fc_vector[is.na(fc_vector)] <- 0

# Convert example gene sets to PIANO format
# Create a GSC object manually
metabolite_gsc_list <- list()
for (met_name in names(example_metabolite_genes)) {
  metabolite_gsc_list[[met_name]] <- example_metabolite_genes[[met_name]]
}

# Prepare GSC object
gsc_metabolites <- list(gsc = metabolite_gsc_list)
class(gsc_metabolites) <- "GSC"

# Run reporter analysis with PIANO
reporter_results <- runGSA(
  geneLevelStats = pval_vector,
  directions = fc_vector,
  gsc = gsc_metabolites,
  geneSetStat = "reporter",      # Use reporter metabolite method
  signifMethod = "geneSampling",  # Gene sampling for significance
  nPerm = 1000,                   # Number of permutations
  gsSizeLim = c(3, 100),         # Gene set size limits
  adjMethod = "fdr"               # FDR correction
)

cat("\n✓ Reporter metabolite analysis complete!\n")
```

**What this does**:
- `geneSetStat = "reporter"`: Uses reporter algorithm (Patil & Nielsen)
- `nPerm = 1000`: Generates null distribution via permutation
- `gsSizeLim = c(3, 100)`: Only considers metabolites with 3-100 genes
- `adjMethod = "fdr"`: FDR correction for multiple testing

---

#### Step 7: Extract and Interpret Results

```r
# Get summary table
reporter_summary <- GSAsummaryTable(reporter_results,
                                   save = FALSE)

# View top reporter metabolites
head(reporter_summary, 10)

# Save results to file
GSAsummaryTable(reporter_results,
               save = TRUE,
               file = "Reporter_Metabolites_output.txt")

cat("\n✓ Results saved to: Reporter_Metabolites_output.txt\n")
```

**Expected output columns**:
- `Name`: Metabolite name
- `Genes (tot)`: Total genes associated with metabolite
- `Stat (dist.dir.up)`: Distinct directional up-regulation
- `Stat (dist.dir.dn)`: Distinct directional down-regulation
- `p (dist.dir.up)`: P-value for up-regulation
- `p adj (dist.dir.up)`: FDR-adjusted p-value
- `Genes (up)`: Number of up-regulated genes
- `Genes (down)`: Number of down-regulated genes

---

#### Step 8: Visualize Top Reporter Metabolites

```r
library(ggplot2)

# Extract reporter scores and p-values
# Create a data frame for plotting
reporter_df <- data.frame(
  Metabolite = rownames(reporter_summary),
  DistDirUp = reporter_summary$`Stat (dist.dir.up)`,
  DistDirDn = reporter_summary$`Stat (dist.dir.dn)`,
  PvalUp = reporter_summary$`p adj (dist.dir.up)`,
  PvalDn = reporter_summary$`p adj (dist.dir.dn)`,
  GenesUp = reporter_summary$`Genes (up)`,
  GenesDn = reporter_summary$`Genes (down)`
)

# Calculate overall reporter score
reporter_df$ReporterScore <- with(reporter_df,
  ifelse(PvalUp < PvalDn, DistDirUp, -DistDirDn))

reporter_df$Significant <- with(reporter_df,
  ifelse(PvalUp < 0.05 | PvalDn < 0.05, "Yes", "No"))

# Order by absolute reporter score
reporter_df <- reporter_df[order(-abs(reporter_df$ReporterScore)), ]

# Plot top 10 reporter metabolites
top_n <- 10
plot_data <- head(reporter_df, top_n)

ggplot(plot_data, aes(x = reorder(Metabolite, ReporterScore),
                      y = ReporterScore,
                      fill = Significant)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("No" = "gray70", "Yes" = "red")) +
  theme_minimal() +
  labs(
    title = "Top Reporter Metabolites (Late vs Early)",
    x = "Metabolite",
    y = "Reporter Score",
    fill = "FDR < 0.05"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed")

ggsave("figures/reporter_metabolites_barplot.pdf", width = 8, height = 6)
```

---

### Approach 2: Using RAVEN Toolbox (MATLAB/Octave)

For more comprehensive analysis with the full SBML model:

#### Step 1: Prepare Input File

We already created this in Module 1:

```r
# This was created earlier
# DESeq output for RAVEN.txt contains:
# - Row names: Ensembl IDs
# - Column 1: pvalue
# - Column 2: log2FoldChange
# - No column names
```

#### Step 2: Run RAVEN (Outside of R)

**In MATLAB/Octave** (not covered in detail in this workshop):

```matlab
% Load model
model = importModel('Reference_model.xml');

% Load gene data
geneData = readtable('DESeq output for RAVEN.txt');

% Run reporter metabolite analysis
[reporterMetabolites, scores, pvals] = reporterMetabolites(...
    model, geneData(:,1), geneData(:,2));

% Export results
writetable(reporterMetabolites, 'RAVEN_reporter_output.txt');
```

**Note**: RAVEN setup is beyond this workshop scope, but this shows the workflow.

---

## 📊 Interpreting Results

### Understanding Reporter Scores

**High positive score**:
- Genes producing/consuming this metabolite are UP-regulated
- Metabolite production likely increased
- Suggests enhanced pathway flux

**High negative score**:
- Genes are DOWN-regulated
- Metabolite production likely decreased
- Suggests suppressed pathway

**Low score**:
- Genes show mixed or no coordinated regulation
- Not a reporter metabolite

---

### Biological Interpretation Example

**Scenario**: Lactate is top reporter metabolite (score > 5, FDR < 0.001)

**Associated genes**: LDHA, LDHB (up-regulated), SLC16A1 (up-regulated)

**Biological interpretation**:
1. **Lactate dehydrogenase** (LDHA/B) enzymes increased
2. **Lactate transporter** (SLC16A1) increased
3. → **Warburg effect**: Cancer cells produce lactate even with oxygen
4. → Suggests glycolytic metabolism in Late stage tumors
5. → Potential therapeutic target with LDHA inhibitors

---

### Common Reporter Metabolites in Cancer

| Metabolite | Direction | Biological Meaning |
|------------|-----------|-------------------|
| **Lactate** | ↑ | Glycolytic metabolism (Warburg effect) |
| **Glutamine** | ↑ | Glutamine addiction, anaplerosis |
| **Fatty Acids** | ↑ | Lipid synthesis for membranes |
| **ATP** | ↑ | Energy production |
| **NAD/NADH** | ↑ | Redox balance |
| **Succinate** | ↓ | TCA cycle disruption |
| **Citrate** | ↓ | Reduced oxidative metabolism |

---

## 🔍 Hands-On Exercise

### Exercise 1: Identify Top Reporter Metabolites

Using the results from the analysis above:

1. **How many reporter metabolites have FDR < 0.05?**
2. **What is the top up-regulated reporter metabolite?**
3. **What is the top down-regulated reporter metabolite?**
4. **How many genes are associated with the top reporter?**

**Code to answer**:
```r
# Question 1
sig_metabolites <- sum(reporter_df$PvalUp < 0.05 | reporter_df$PvalDn < 0.05)

# Question 2
top_up <- reporter_df[which.max(reporter_df$DistDirUp), ]

# Question 3
top_down <- reporter_df[which.max(reporter_df$DistDirDn), ]

# Question 4
total_genes <- with(top_up, GenesUp + GenesDn)
```

---

### Exercise 2: Gene-Level Investigation

For the top reporter metabolite:

1. **List all genes associated with this metabolite**
2. **What percentage are up-regulated? Down-regulated?**
3. **Look up the biological function** (use Google/Wikipedia)
4. **Propose a biological hypothesis** for why this metabolite is important

**Code to help**:
```r
# Get top metabolite name
top_metabolite <- reporter_df$Metabolite[1]

# Get associated genes
associated_genes <- example_metabolite_genes[[top_metabolite]]

# Get their fold changes
gene_fc <- DESeq_for_reporter[associated_genes, ]
print(gene_fc)

# Count up vs down
n_up <- sum(gene_fc$log2FoldChange > 0, na.rm = TRUE)
n_down <- sum(gene_fc$log2FoldChange < 0, na.rm = TRUE)

cat("Up-regulated:", n_up, "/", nrow(gene_fc), "\n")
cat("Down-regulated:", n_down, "/", nrow(gene_fc), "\n")
```

---

### Exercise 3: Pathway Context

**Research question**: How do reporter metabolites connect to enriched pathways from Module 2?

1. Compare reporter metabolites to GSEA results
2. Find overlap between metabolic genes and enriched pathways
3. Create integrated interpretation

**Example**:
- GSEA shows "Glycolysis" pathway enriched
- Reporter analysis shows "Lactate" and "Pyruvate" as reporters
- → **Integrated interpretation**: Late stage tumors have enhanced glycolysis, producing lactate from pyruvate

---

## 🎯 Advanced Topics

### 1. Network Visualization

Visualize metabolite-gene networks:

```r
# Create edge list for Cytoscape
edges <- data.frame()

for (met in names(example_metabolite_genes)) {
  genes <- example_metabolite_genes[[met]]
  for (gene in genes) {
    edges <- rbind(edges, data.frame(
      Source = met,
      Target = gene,
      Type = "metabolite_gene"
    ))
  }
}

# Add gene fold changes
edges$FoldChange <- DESeq_for_reporter[edges$Target, "log2FoldChange"]

# Save for Cytoscape
write.table(edges,
           "metabolite_gene_network.txt",
           sep = "\t",
           row.names = FALSE,
           quote = FALSE)

cat("✓ Network file created for Cytoscape\n")
```

---

### 2. Subsystem Analysis

Group metabolites by pathway:

```r
# Define subsystems
subsystems <- list(
  "Glycolysis" = c("Glucose", "Pyruvate", "Lactate"),
  "TCA_Cycle" = c("Citrate", "Succinate", "ATP"),
  "Amino_Acid" = c("Glutamine"),
  "Lipid" = c("Fatty_Acids", "Acetyl_CoA"),
  "Cofactors" = c("NAD", "ATP")
)

# Aggregate scores by subsystem
subsystem_scores <- sapply(subsystems, function(mets) {
  scores <- reporter_df[reporter_df$Metabolite %in% mets, "ReporterScore"]
  mean(scores, na.rm = TRUE)
})

print(subsystem_scores)
```

---

### 3. Compare with Measured Metabolomics

If you have metabolomics data:

```r
# Example: Compare predicted vs measured changes
# predicted_changes <- reporter_df$ReporterScore
# measured_changes <- metabolomics_fc

# correlation <- cor.test(predicted, measured)
#
# ggplot(comparison, aes(x = Predicted, y = Measured)) +
#   geom_point() +
#   geom_smooth(method = "lm") +
#   labs(title = "Predicted vs Measured Metabolite Changes")
```

---

## 🎯 Key Takeaways

1. ✅ **Reporter metabolites** = metabolites with coordinated gene regulation
2. ✅ Links **gene expression** to **metabolic network**
3. ✅ Provides **mechanistic** interpretation of transcriptomics
4. ✅ Requires **genome-scale model** (gene-metabolite associations)
5. ✅ Complements **pathway analysis** with network-level view
6. ✅ **Direction matters**: up/down-regulation suggests flux changes
7. ✅ **Validation**: Compare with measured metabolomics when available

---

## 📚 Additional Resources

### Papers
- **Original method**: Patil & Nielsen (2005) PNAS - "Uncovering transcriptional regulation of metabolism"
- **Application**: Bordbar et al. (2010) Mol Syst Biol - "Model-driven multi-omic data analysis"

### Software
- **PIANO**: https://bioconductor.org/packages/release/bioc/html/piano.html
- **RAVEN Toolbox**: https://github.com/SysBioChalmers/RAVEN
- **iMAT**: Alternative method for metabolic integration

### Models
- **Human-GEM**: https://github.com/SysBioChalmers/Human-GEM
- **Recon3D**: http://vmh.life
- **BiGG Models**: http://bigg.ucsd.edu

---

## ⏭️ Next Steps

**Main Workshop Path**:
- Workshop complete! Review and practice

**Bonus Advanced Path**:
- Continue to **Module 5: MOFA2 Introduction** for multi-omics integration
- Learn how to integrate transcriptomics with proteomics, metabolomics, etc.

**Solutions**: Available in `solutions/Module4_solutions.R`

---

**Congratulations on completing the core workshop!** 🎉

You now have skills in:
- ✅ Differential expression analysis
- ✅ Gene set enrichment analysis
- ✅ Co-expression networks
- ✅ Reporter metabolite analysis

These form a complete transcriptomics-to-metabolism analysis pipeline!

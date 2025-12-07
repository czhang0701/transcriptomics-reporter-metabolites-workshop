# Complete Workshop Guide
## Transcriptomics to Reporter Metabolites - 3 Hour Workshop

**Duration**: 3 hours (including installation)
**Main Goal**: Complete reporter metabolite analysis
**Bonus**: MOFA2 for fast learners (optional, can do after workshop)

---

## 📅 Workshop Timeline

### **Total: 3 Hours**

```
0:00-0:30 (30m) Installation & Setup (Windows & Mac)
0:30-0:50 (20m) Module 1: Differential Expression
0:50-1:05 (15m) Module 2: Gene Set Enrichment
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
1:05-1:15 (10m) BREAK
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
1:15-1:30 (15m) Module 3: Co-expression Networks
1:30-2:30 (60m) ⭐ Module 4: Reporter Metabolites (MAIN)
2:30-2:50 (20m) Interpretation & Discussion
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
2:50-3:00 (10m) Wrap-up & Next Steps
```

**BONUS** (Optional - for fast students or after workshop):
- MOFA2 materials available for self-study
- Can work on during workshop if finished early
- Full tutorials provided for after workshop

---

## 🎯 What You'll Learn

**Guaranteed (Everyone)**:
- ✅ Install R, RStudio, and packages
- ✅ Run differential expression analysis
- ✅ Perform pathway enrichment
- ✅ Build co-expression networks
- ✅ **Identify reporter metabolites** ⭐

**Bonus (Fast Learners/After Workshop)**:
- Multi-omics integration with MOFA2
- Advanced factor analysis
- Self-paced learning materials

---

# PART 1: Installation & Setup (30 minutes)

## 💻 Step 1: Install R (10 minutes)

### **Windows Users**

1. **Download R**:
   - Go to: https://cran.r-project.org/bin/windows/base/
   - Click **"Download R-4.x.x for Windows"**
   - Save the `.exe` file

2. **Install R**:
   - Double-click the downloaded `.exe` file
   - Click **"Next"** through all prompts (accept all defaults)
   - Click **"Finish"**

3. **Verify**:
   - Press `Windows + R` key
   - Type `cmd` and press Enter
   - Type: `R --version`
   - Should show: `R version 4.x.x`

---

### **Mac Users**

1. **Check Your Mac Type**:
   - Click Apple menu () → **About This Mac**
   - Note if you have:
     - **Apple Silicon** (M1, M2, M3) OR
     - **Intel processor**

2. **Download R**:
   - Go to: https://cran.r-project.org/bin/macosx/
   - **Apple Silicon (M1/M2/M3)**: Download `R-4.x.x-arm64.pkg`
   - **Intel Mac**: Download `R-4.x.x-x86_64.pkg`

3. **Install R**:
   - Double-click the downloaded `.pkg` file
   - Click **"Continue"** through installer
   - Enter your Mac password when prompted
   - Click **"Close"**

4. **Verify**:
   - Open **Terminal** (Cmd+Space, type "Terminal")
   - Type: `R --version`
   - Should show: `R version 4.x.x`

**Mac Issue**: If "command not found":
```bash
echo 'export PATH="/Library/Frameworks/R.framework/Resources/bin:$PATH"' >> ~/.zshrc
source ~/.zshrc
```

---

## 🖥️ Step 2: Install RStudio (5 minutes)

### **Both Windows and Mac**

1. **Download RStudio**:
   - Go to: https://posit.co/download/rstudio-desktop/
   - Scroll to **"All Installers"**
   - **Windows**: Download `.exe` file
   - **Mac**: Download `.dmg` file

2. **Install**:
   - **Windows**: Double-click `.exe` → Next → Next → Finish
   - **Mac**: Double-click `.dmg` → Drag RStudio to Applications

3. **Launch RStudio**:
   - **Windows**: Start Menu → RStudio
   - **Mac**: Applications → RStudio

4. **Check**:
   - Bottom-left panel should show: `R version 4.x.x`

---

## 📦 Step 3: Install R Packages (15 minutes)

**CRITICAL**: This takes 10-15 minutes. Start immediately!

### **Copy and Paste into RStudio Console**

```r
# Install BiocManager (package manager)
install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install(c("DESeq2", "piano"))

# Install CRAN packages
install.packages(c("corrplot", "Hmisc", "reshape2",
                   "ggplot2", "pheatmap", "tidyverse"))
```

**Press Enter and wait...**

### **During Installation** (you'll see a lot of text):

**When asked "Update all/some/none?"**
- Type `a` (for "all") and press Enter

**When asked "Do you want to install from sources?"**
- Type `n` (for "no") and press Enter

**Red text is OK!** - Not all red text means error

**This will take 10-15 minutes** - be patient!

---

### **Platform-Specific Issues During Installation**

#### **Windows Issues**

**Issue**: "WARNING: Rtools is required"
- Download: https://cran.r-project.org/bin/windows/Rtools/
- Install Rtools43
- Restart RStudio

**Issue**: "Permission denied"
- Close RStudio
- Right-click RStudio icon → "Run as administrator"
- Try again

---

#### **Mac Issues**

**Issue**: "clang: error" or "compilation failed"
- Open **Terminal**
- Type: `xcode-select --install`
- Click "Install" in popup
- Wait 10-15 minutes
- Return to RStudio and try again

**Issue**: "gfortran is required"
- Go to: https://mac.r-project.org/tools/
- **Intel Mac**: Download `gfortran-12.2-universal.pkg`
- **Apple Silicon**: Download `gfortran-ARM-12.2-Monterey.dmg`
- Install and restart RStudio

---

### **Verify Installation**

After packages finish installing, run this:

```r
# Check if packages installed
packages <- c("DESeq2", "piano", "corrplot", "Hmisc", "reshape2", "ggplot2")

for (pkg in packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat("✓", pkg, "\n")
  } else {
    cat("✗", pkg, "NOT INSTALLED\n")
  }
}
```

**Expected Output**:
```
✓ DESeq2
✓ piano
✓ corrplot
✓ Hmisc
✓ reshape2
✓ ggplot2
```

**If you see any ✗**: Ask instructor for help during break!

---

## 📥 Step 4: Download Workshop Materials

### **Option 1: Download ZIP (Easiest)**

1. Go to: [GitHub URL will be provided by instructor]
2. Click green **"Code"** button
3. Click **"Download ZIP"**
4. Extract ZIP:
   - **Windows**: Right-click → "Extract All" → Choose location
   - **Mac**: Double-click ZIP file

### **Remember the location!**

**Windows example**: `C:\Users\YourName\Downloads\workshop_materials`
**Mac example**: `/Users/YourName/Downloads/workshop_materials`

---

### **Set Working Directory in RStudio**

```r
# ADJUST THIS PATH TO YOUR LOCATION!

# Windows example:
setwd("C:/Users/YourName/Downloads/workshop_materials")

# Mac example:
setwd("/Users/YourName/Downloads/workshop_materials")

# Check it worked:
getwd()

# List files:
list.files()
```

**You should see folders**: data, scripts, tutorials, etc.

---

# PART 2: Core Analysis (1 hour 55 minutes)

## Module 1: Differential Expression (20 minutes)

### **Load Libraries**

```r
library(DESeq2)
library(ggplot2)
```

---

### **Load Data**

```r
# Import counts data
Counts <- as.matrix(read.csv(file = "data/Counts_selected.txt",
                             sep = "\t",
                             row.names = 1))

# Import metadata
Metadata <- read.csv(file = "data/patients.txt",
                    sep = "\t",
                    row.names = 1)

# Check data
head(Counts[, 1:5])
head(Metadata)
```

---

### **Create DESeq Object and Run Analysis**

```r
# Create DESeq object
deSeqData <- DESeqDataSetFromMatrix(
  countData = Counts,
  colData = Metadata,
  design = ~ Type
)

# Run DESeq2 (takes 1-2 minutes)
deSeqAnalysis <- DESeq(deSeqData)

# Extract results
res <- results(
  deSeqAnalysis,
  contrast = c("Type", "Late", "Early"),
  alpha = 0.05
)

# View summary
summary(res)
```

---

### **Quick Volcano Plot**

```r
# Prepare data
res_df <- as.data.frame(res)
res_df$significant <- res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1

# Plot
ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = significant), alpha = 0.5) +
  scale_color_manual(values = c("gray", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot: Late vs Early",
       x = "Log2 Fold Change", y = "-Log10 P-value")
```

---

### **Save Results**

```r
# Save for next modules
write.table(as.data.frame(res),
           file = "data/DESeq_output.txt",
           sep = "\t",
           row.names = TRUE,
           col.names = NA)

cat("✓ DESeq2 complete! Significant genes:",
    sum(res$padj < 0.05, na.rm = TRUE), "\n")
```

---

## Module 2: Gene Set Enrichment (15 minutes)

### **Load and Prepare Data**

```r
library(piano)

# Load DESeq results
DESeq_results <- read.delim("data/DESeq_output.txt",
                            row.names = 1,
                            stringsAsFactors = FALSE)

# Map to gene symbols
ensembl2gene <- read.delim("data/Ensembl2gene.tsv",
                          row.names = 2,
                          stringsAsFactors = FALSE)

DESeq_results$Gene <- ensembl2gene[rownames(DESeq_results), "Gene"]
DESeq_results <- DESeq_results[!is.na(DESeq_results$Gene), ]
rownames(DESeq_results) <- DESeq_results$Gene
```

---

### **Prepare for PIANO**

```r
# Keep needed columns
DESeq_for_piano <- DESeq_results[, c('log2FoldChange', 'pvalue')]

# Handle NA values
DESeq_for_piano$pvalue[is.na(DESeq_for_piano$pvalue)] <- 1
DESeq_for_piano$log2FoldChange[is.na(DESeq_for_piano$log2FoldChange)] <- 0

# Remove duplicates
DESeq_for_piano <- DESeq_for_piano[!duplicated(rownames(DESeq_for_piano)), ]

# Create matrices
pval_matrix <- as.matrix(DESeq_for_piano[, 'pvalue', drop = FALSE])
fc_matrix <- as.matrix(DESeq_for_piano[, 'log2FoldChange', drop = FALSE])
```

---

### **Run PIANO (takes 3-5 minutes)**

```r
# Load gene sets
gsc <- loadGSC("data/c5.bp.v6.2.symbols.gmt")

# Run enrichment
cat("Running PIANO (takes 3-5 min)...\n")
gsaRes <- runGSA(
  geneLevelStats = pval_matrix,
  directions = fc_matrix,
  gsc = gsc,
  geneSetStat = "reporter",
  nPerm = 1000,
  gsSizeLim = c(5, 500),
  adjMethod = "fdr"
)

# Save results
GSAsummaryTable(gsaRes, save = TRUE, file = "data/Piano_output.txt")
cat("✓ GSEA complete!\n")
```

---

## 🕐 BREAK (10 minutes)

Stretch, coffee, ask questions!

---

## Module 3: Co-expression Networks (15 minutes)

### **Load and Prepare Data**

```r
library(corrplot)
library(Hmisc)
library(reshape2)

# Load FPKM data
FPKM_data <- read.csv("data/FPKM_selected.txt",
                      sep = "\t",
                      row.names = 1)

# Load target genes
target_genes <- read.delim("data/Genes_selected.txt",
                          sep = "\t",
                          row.names = 1)

# Subset to target genes
FPKM_subset <- FPKM_data[rownames(target_genes), ]

# Map to gene symbols
gene_symbols <- ensembl2gene[rownames(FPKM_subset), "Gene"]
rownames(FPKM_subset) <- gene_symbols
```

---

### **Calculate Correlations**

```r
# Filter low expression
mean_expression <- rowMeans(FPKM_subset)
FPKM_filtered <- FPKM_subset[mean_expression > 5, ]

# Transpose (samples in rows, genes in columns)
FPKM_transposed <- t(as.matrix(FPKM_filtered))

# Calculate correlations
cor_result <- rcorr(FPKM_transposed, type = "spearman")
cor_matrix <- cor_result$r
pval_matrix <- cor_result$P
```

---

### **Extract Significant Correlations**

```r
# Get upper triangle
upper_tri_cor <- cor_matrix
upper_tri_cor[lower.tri(upper_tri_cor, diag = TRUE)] <- NA

upper_tri_pval <- pval_matrix
upper_tri_pval[lower.tri(upper_tri_pval, diag = TRUE)] <- NA

# Reshape
cor_long <- melt(upper_tri_cor)
pval_long <- melt(upper_tri_pval)

# Combine
correlation_data <- data.frame(
  Gene1 = cor_long$Var1,
  Gene2 = cor_long$Var2,
  Correlation = cor_long$value,
  Pvalue = pval_long$value
)

# Remove NA
correlation_data <- correlation_data[!is.na(correlation_data$Correlation), ]
correlation_data <- correlation_data[
  as.character(correlation_data$Gene1) != as.character(correlation_data$Gene2),
]

# FDR correction
correlation_data$FDR <- p.adjust(correlation_data$Pvalue, method = "fdr")

# Keep significant
significant_cors <- correlation_data[correlation_data$FDR < 0.05, ]
```

---

### **Save Network**

```r
# Save for Cytoscape
write.table(significant_cors,
           "data/Coexpression_network.txt",
           sep = "\t",
           row.names = FALSE,
           quote = FALSE)

cat("✓ Network complete! Significant correlations:",
    nrow(significant_cors), "\n")
```

---

## ⭐ Module 4: Reporter Metabolites (60 minutes)

### **MAIN FOCUS - This is what we came for!**

---

### **Part 1: Concept (10 minutes)**

**Instructor Explains**:

**What are Reporter Metabolites?**
- Metabolites whose surrounding genes show coordinated differential expression
- Links transcriptomics to metabolic networks
- Identifies key metabolic reprogramming events

**How does it work?**
1. Take differential expression results (gene-level p-values)
2. Use genome-scale metabolic model (gene-metabolite associations)
3. For each metabolite, aggregate p-values of associated genes
4. Statistical test: Are this metabolite's genes coordinately changed?

**Why is this useful?**
- Mechanistic interpretation of gene expression
- Identifies metabolic targets
- Connects to Warburg effect, glutamine addiction, etc.
- Suggests what to measure experimentally

---

### **Part 2: Analysis (30 minutes)**

#### **Load DESeq Results**

```r
library(piano)

# Load DESeq results
DESeq_results <- read.delim("data/DESeq_output.txt",
                            row.names = 1,
                            stringsAsFactors = FALSE)

# Map to gene symbols
DESeq_results$Gene <- ensembl2gene[rownames(DESeq_results), "Gene"]
DESeq_results <- DESeq_results[!is.na(DESeq_results$Gene), ]
rownames(DESeq_results) <- DESeq_results$Gene

# Prepare for analysis
DESeq_for_reporter <- DESeq_results[, c('pvalue', 'log2FoldChange')]
DESeq_for_reporter$pvalue[is.na(DESeq_for_reporter$pvalue)] <- 1
DESeq_for_reporter$log2FoldChange[is.na(DESeq_for_reporter$log2FoldChange)] <- 0
```

---

#### **Create Metabolite-Gene Associations**

For this workshop, we'll use example associations (faster than full model):

```r
# Example metabolite-gene associations
metabolite_genes <- list(
  "Pyruvate" = c("LDHA", "LDHB", "PKM", "PDHA1", "PDHB", "PC"),
  "Lactate" = c("LDHA", "LDHB", "SLC16A1", "SLC16A3", "SLC16A7"),
  "Glucose" = c("SLC2A1", "SLC2A3", "HK1", "HK2", "GCK"),
  "Citrate" = c("CS", "ACLY", "SLC25A1", "IDH1", "IDH2"),
  "Glutamine" = c("GLUL", "GLS", "GLS2", "SLC1A5", "SLC38A1"),
  "Glutamate" = c("GLUD1", "GLUD2", "GOT1", "GOT2", "GLS"),
  "Fatty_Acids" = c("FASN", "ACACA", "SCD", "ACSL1", "ACSL4", "CPT1A"),
  "ATP" = c("ATP5A1", "ATP5B", "ATP5C1", "ATP5D", "ATP5E", "ATP5F1"),
  "NAD" = c("NAMPT", "NMNAT1", "NMNAT2", "NMNAT3", "NAPRT"),
  "Succinate" = c("SDHA", "SDHB", "SDHC", "SDHD", "SUCLA2", "SUCLG1"),
  "Acetyl_CoA" = c("ACLY", "ACSS2", "PDHA1", "PDHB", "ACSS1"),
  "Alpha_Ketoglutarate" = c("IDH1", "IDH2", "IDH3A", "OGDH", "GOT1", "GOT2"),
  "Fumarate" = c("FH", "SDHA", "SDHB", "SDHC", "SDHD"),
  "Malate" = c("MDH1", "MDH2", "ME1", "ME2", "ME3"),
  "Oxaloacetate" = c("MDH1", "MDH2", "GOT1", "GOT2", "PC")
)

cat("Created", length(metabolite_genes), "metabolite gene sets\n")
```

---

#### **Prepare for PIANO**

```r
# Convert to GSC format
gsc_metabolites <- list()
for (met_name in names(metabolite_genes)) {
  gsc_metabolites[[met_name]] <- metabolite_genes[[met_name]]
}

# Create GSC object
class(gsc_metabolites) <- "GSC"

# Prepare matrices
pval_vector <- as.matrix(DESeq_for_reporter[, 'pvalue'])
fc_vector <- as.matrix(DESeq_for_reporter[, 'log2FoldChange'])
rownames(pval_vector) <- rownames(DESeq_for_reporter)
rownames(fc_vector) <- rownames(DESeq_for_reporter)
```

---

#### **Run Reporter Metabolite Analysis**

```r
cat("Running reporter metabolite analysis...\n")
cat("This will take 3-5 minutes...\n\n")

reporter_results <- runGSA(
  geneLevelStats = pval_vector,
  directions = fc_vector,
  gsc = gsc_metabolites,
  geneSetStat = "reporter",
  signifMethod = "geneSampling",
  nPerm = 1000,
  gsSizeLim = c(3, 100),
  adjMethod = "fdr"
)

cat("✓ Reporter analysis complete!\n")
```

---

#### **Extract Results**

```r
# Get summary
reporter_summary <- GSAsummaryTable(reporter_results, save = FALSE)

# Save to file
GSAsummaryTable(reporter_results,
               save = TRUE,
               file = "data/Reporter_Metabolites_output.txt")

# View top results
head(reporter_summary, 15)
```

---

### **Part 3: Interpretation (20 minutes)**

#### **Identify Top Reporter Metabolites**

```r
# Prepare for visualization
reporter_df <- data.frame(
  Metabolite = rownames(reporter_summary),
  DistDirUp = reporter_summary$`Stat (dist.dir.up)`,
  DistDirDn = reporter_summary$`Stat (dist.dir.dn)`,
  PvalUp = reporter_summary$`p adj (dist.dir.up)`,
  PvalDn = reporter_summary$`p adj (dist.dir.dn)`,
  GenesUp = reporter_summary$`Genes (up)`,
  GenesDn = reporter_summary$`Genes (down)`
)

# Calculate overall score
reporter_df$ReporterScore <- with(reporter_df,
  ifelse(PvalUp < PvalDn, DistDirUp, -DistDirDn))

reporter_df$Significant <- with(reporter_df,
  ifelse(PvalUp < 0.05 | PvalDn < 0.05, "Yes", "No"))

# Order by score
reporter_df <- reporter_df[order(-abs(reporter_df$ReporterScore)), ]

# View top 10
print(reporter_df[1:10, ])
```

---

#### **Visualize Top Reporter Metabolites**

```r
library(ggplot2)

# Top 10 reporters
plot_data <- head(reporter_df, 10)

ggplot(plot_data, aes(x = reorder(Metabolite, ReporterScore),
                      y = ReporterScore,
                      fill = Significant)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("No" = "gray70", "Yes" = "red")) +
  theme_minimal() +
  labs(
    title = "Top Reporter Metabolites (Late vs Early Stage)",
    x = "Metabolite",
    y = "Reporter Score",
    fill = "FDR < 0.05"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed")

ggsave("figures/reporter_metabolites.pdf", width = 8, height = 6)
```

---

#### **Biological Interpretation**

**Group Discussion** (Instructor leads):

**Example: If Lactate is top reporter metabolite**

1. **Observation**:
   - Lactate has high positive reporter score (e.g., > 5)
   - FDR < 0.001
   - Genes: LDHA, LDHB both up-regulated

2. **Mechanism**:
   - Lactate dehydrogenase (LDHA/LDHB) converts pyruvate → lactate
   - Lactate transporter (SLC16A1) also up-regulated
   - Coordinated upregulation of lactate production and export

3. **Biological Context**:
   - **Warburg effect**: Cancer cells produce lactate even with oxygen
   - Late-stage tumors have enhanced glycolysis
   - Lactate acidifies tumor microenvironment
   - Supports tumor growth and immune evasion

4. **Connection to Earlier Results**:
   - Module 2 showed "Glycolysis" pathway enriched → produces pyruvate
   - Module 4 shows "Lactate" as reporter → pyruvate converted to lactate
   - **Complete story**: Enhanced glycolysis → pyruvate → lactate production

5. **Clinical Implications**:
   - Could target LDHA pharmacologically
   - Measure lactate levels in patient samples
   - Prognostic marker for aggressive disease

---

#### **Investigate Specific Metabolite**

```r
# Pick top metabolite
top_metabolite <- reporter_df$Metabolite[1]
cat("\nInvestigating:", top_metabolite, "\n\n")

# Get associated genes
associated_genes <- metabolite_genes[[top_metabolite]]
cat("Associated genes:", paste(associated_genes, collapse = ", "), "\n\n")

# Get their fold changes
gene_stats <- DESeq_for_reporter[associated_genes, ]
gene_stats <- gene_stats[order(-gene_stats$log2FoldChange), ]

cat("Gene expression changes:\n")
print(gene_stats)

# How many up vs down?
n_up <- sum(gene_stats$log2FoldChange > 0, na.rm = TRUE)
n_down <- sum(gene_stats$log2FoldChange < 0, na.rm = TRUE)

cat("\nUp-regulated:", n_up, "/", nrow(gene_stats), "\n")
cat("Down-regulated:", n_down, "/", nrow(gene_stats), "\n")
```

---

## 💬 Discussion & Wrap-up (20 minutes)

### **What We Learned**

**The Complete Pipeline**:
```
1. DESeq2 → 2,000 significant genes
   ↓
2. GSEA → 50 enriched pathways (e.g., Glycolysis)
   ↓
3. Networks → Gene modules (e.g., glycolysis genes correlated)
   ↓
4. Reporter Metabolites → Key metabolites (e.g., Lactate)
   ↓
5. Biological Interpretation → Warburg effect in cancer
```

---

### **Key Takeaways**

1. **Reporter metabolites** = metabolites with coordinated gene regulation
2. **Connects** transcriptomics to metabolism
3. **Mechanistic** interpretation beyond statistics
4. **Identifies** therapeutic targets
5. **Validates** with metabolomics when available

---

### **Your Results**

You now have:
- ✅ `DESeq_output.txt` - differential expression
- ✅ `Piano_output.txt` - enriched pathways
- ✅ `Coexpression_network.txt` - gene networks
- ✅ `Reporter_Metabolites_output.txt` - **key metabolic changes**

**These are publication-ready results!**

---

## 🕐 FINAL BREAK (10 minutes)

---

## 🎓 Next Steps & Bonus Material (10 minutes)

### **Continue Learning**

**Full Tutorials** (in `tutorials/` folder):
- Detailed versions of Modules 1-4
- Additional exercises
- Advanced topics
- More interpretation

**Practice Dataset**:
```r
# Create practice dataset
source("scripts/create_practice_dataset.R")

# Practice complete workflow
# (smaller, faster for learning)
```

---

### **Apply to Your Own Data**

**Requirements**:
1. Gene expression counts matrix
2. Sample metadata (groups to compare)
3. Gene ID mapping (if using Ensembl IDs)

**Steps**:
1. Replace file paths in code
2. Adjust group names in design formula
3. Run same workflow!

---

### **BONUS: MOFA2 Multi-Omics Integration**

**For Fast Learners** (optional, self-paced):

If you finished early or want to continue learning after workshop:

**Module 5: MOFA2 Introduction** (30 min)
- Located in `tutorials/Module5_MOFA2_Introduction.md`
- Concept of multi-omics factor analysis
- When to use MOFA2
- How it extends reporter metabolites

**What is MOFA2?**
- Integrates MULTIPLE omics types (RNA-seq + proteomics + metabolomics)
- Identifies latent factors explaining variation
- Determines which omics layers drive each factor
- Complements reporter metabolite analysis

**Self-Study Materials**:
- Module 5 tutorial (complete)
- Modules 6-8 tutorials (coming soon)
- Practice with simulated multi-omics data
- Apply to own multi-omics datasets

**Note**: MOFA2 requires:
- Multi-omics data (not just RNA-seq)
- Larger sample sizes (n > 20 recommended)
- More advanced statistical knowledge

---

### **Resources**

**Papers**:
- Reporter metabolites: Patil & Nielsen (2005) PNAS
- DESeq2: Love et al. (2014) Genome Biology
- PIANO: Väremo et al. (2013) Nucleic Acids Research
- MOFA2: Argelaguet et al. (2020) Genome Biology

**Databases**:
- MSigDB: http://www.gsea-msigdb.org/
- Human-GEM: https://github.com/SysBioChalmers/Human-GEM
- MOFA2: https://biofam.github.io/MOFA2/

**Software**:
- Cytoscape: https://cytoscape.org/
- RAVEN Toolbox: https://github.com/SysBioChalmers/RAVEN

---

### **Get Help**

**Questions?**
- Email instructor: [email]
- Check tutorial PDFs in `tutorials/`
- Solutions in `solutions/` folder

**Feedback**:
- What worked well?
- What was confusing?
- What would you like to see added?

---

## ✅ Final Checklist

**By end of workshop, you should have**:

**Files Generated**:
- [ ] `DESeq_output.txt`
- [ ] `Piano_output.txt`
- [ ] `Coexpression_network.txt`
- [ ] `Reporter_Metabolites_output.txt` ⭐

**Skills Learned**:
- [ ] Install R, RStudio, packages
- [ ] Run differential expression
- [ ] Perform pathway enrichment
- [ ] Build co-expression networks
- [ ] **Identify reporter metabolites** ⭐
- [ ] **Interpret metabolic reprogramming** ⭐

**Understanding**:
- [ ] What are reporter metabolites?
- [ ] How to connect genes → pathways → metabolites
- [ ] How to interpret results biologically
- [ ] How to apply to own data

---

## 🎉 Congratulations!

You've completed the workshop and learned:
- ✅ Complete transcriptomics analysis workflow
- ✅ Metabolic network integration
- ✅ Reporter metabolite analysis
- ✅ Biological interpretation of metabolic reprogramming

**Continue exploring with**:
- Full detailed tutorials
- Practice dataset
- MOFA2 bonus materials (optional)
- Your own data!

---

## 📚 Appendix: Quick Reference

### **Common Commands**

```r
# Set working directory
setwd("path/to/workshop_materials")

# Check where you are
getwd()

# List files
list.files()

# Install a package
install.packages("package_name")

# Load a library
library(package_name)

# Get help
?function_name
help(package_name)
```

---

### **File Paths**

**Windows**: Use forward slashes or double backslashes
```r
setwd("C:/Users/Name/Documents/workshop")  # Good
setwd("C:\\Users\\Name\\Documents\\workshop")  # Also good
```

**Mac**: Use forward slashes
```r
setwd("/Users/Name/Documents/workshop")
```

---

### **Troubleshooting**

**Error: "cannot open file"**
- Check working directory: `getwd()`
- Check file exists: `file.exists("data/filename.txt")`
- Fix path or set working directory

**Error: "object not found"**
- Did you run previous code?
- Did you load the library?
- Check spelling

**Error: "package not installed"**
- Install it: `install.packages("package")`
- Or for Bioconductor: `BiocManager::install("package")`

---

**Workshop Complete! 🎊**

**All materials available at**: `workshop_materials/`

**Questions?** Email instructor or check tutorials!

# Workshop Dataset Documentation

This document describes the datasets used in the Integrative Multi-Omics Analysis Workshop.

---

## 📊 Dataset Overview

### Source
**The Cancer Genome Atlas (TCGA)** - Pan-cancer cohort of tumor samples

### Sample Information
- **Total samples**: 50
- **Groups**:
  - Early stage: 25 samples
  - Late stage: 25 samples
- **Data type**: RNA-sequencing (RNA-seq)
- **Sequencing platform**: Illumina

### Biological Question
*Which genes and pathways are differentially expressed between early and late stage tumors?*

---

## 📁 File Descriptions

### 1. Counts_selected.txt

**Description**: Gene expression counts matrix

**Format**:
- Tab-delimited text file
- Rows: Genes (Ensembl gene IDs)
- Columns: Samples (TCGA barcodes)
- Values: Integer counts (unnormalized)

**Dimensions**: ~20,000 genes × 50 samples

**Example**:
```
                    TCGA.MR.A8JO.01A  TCGA.NI.A8LF.01A  ...
ENSG00000000003              3699              5548     ...
ENSG00000000005                 2                31     ...
```

**Usage**: Input for DESeq2 differential expression analysis

**Important Notes**:
- Raw counts (not FPKM, TPM, or other normalized values)
- Contains all genes passing initial quality filters
- Some genes may have zero counts in some samples

---

### 2. FPKM_selected.txt

**Description**: Normalized gene expression values

**Format**:
- Tab-delimited text file
- Rows: Genes (Ensembl gene IDs)
- Columns: Samples (TCGA barcodes)
- Values: FPKM (Fragments Per Kilobase Million) - continuous values

**Dimensions**: ~20,000 genes × 50 samples

**Normalization**: FPKM accounts for:
- Gene length
- Sequencing depth (library size)

**Usage**:
- Co-expression network analysis (Module 3)
- Visualization of expression patterns

**Important Notes**:
- Pre-normalized, do NOT use for DESeq2
- Better for comparing expression patterns across samples
- May contain decimal values and zeros

---

### 3. patients.txt

**Description**: Sample metadata

**Format**:
- Tab-delimited text file
- Rows: Sample IDs (matching Counts and FPKM files)
- Columns: Clinical/experimental variables

**Columns**:
- `Type`: Early or Late stage classification

**Example**:
```
                     Type
TCGA.MR.A8JO.01A     Late
TCGA.NI.A8LF.01A     Late
TCGA.RC.A6M6.01A     Late
```

**Usage**: Design matrix for differential expression analysis

**Important Notes**:
- Sample IDs must exactly match column names in counts matrix
- Row names (sample IDs) are critical for matching

---

### 4. Ensembl2gene.tsv

**Description**: Gene ID to gene symbol mapping

**Format**:
- Tab-delimited text file
- Rows indexed by Ensembl IDs (row 2)
- Column "Gene": Gene symbols (HGNC nomenclature)

**Example**:
```
                    Gene
ENSG00000000003     TSPAN6
ENSG00000000005     TNMD
ENSG00000000419     DPM1
```

**Usage**:
- Convert Ensembl IDs to gene symbols for GSEA
- Add gene symbols to result tables for interpretation

**Important Notes**:
- Not all Ensembl IDs have gene symbols
- Some IDs may map to multiple symbols (use primary)
- Essential for PIANO analysis (requires gene symbols)

---

### 5. c5.bp.v6.2.symbols.gmt

**Description**: Gene Ontology (GO) Biological Processes gene sets

**Source**: [MSigDB](http://www.gsea-msigdb.org/gsea/msigdb/) (Molecular Signatures Database)

**Version**: v6.2

**Format**: GMT (Gene Matrix Transposed)
- Each line = one gene set
- Format: `NAME<tab>DESCRIPTION<tab>GENE1<tab>GENE2<tab>...`

**Example**:
```
GO_APOPTOTIC_PROCESS  http://...  TP53  BCL2  BAX  CASP3  ...
GO_CELL_CYCLE        http://...  CDK1  CCNB1  CDC20  ...
```

**Number of gene sets**: ~5,000 GO Biological Process terms

**Usage**: Gene Set Enrichment Analysis (GSEA) with PIANO

**Important Notes**:
- Uses gene symbols (not Ensembl IDs)
- Each gene set = list of genes involved in a biological process
- Size of gene sets varies (5-500 genes typically)

**Citation**:
Liberzon, A. et al. (2015) Molecular signatures database (MSigDB) 3.0. *Bioinformatics*

---

### 6. Genes_selected.txt

**Description**: Subset of 113 target genes for focused analysis

**Format**:
- Tab-delimited text file
- Two columns: Ensembl ID, Gene symbol

**Example**:
```
                    Gene
ENSG00000026508     CD44
ENSG00000151012     SLC7A11
ENSG00000131069     ACSS2
```

**Usage**: Co-expression network analysis (Module 3)

**Purpose**:
- Reduces computational time for demonstration
- Focuses on biologically relevant genes
- Allows detailed network visualization

**Selection criteria**:
- Highly variable across samples
- Known cancer-related genes
- Diverse functional categories

---

### 7. Reference_model.xml

**Description**: Genome-scale metabolic model (GEM)

**Format**: SBML (Systems Biology Markup Language) XML

**Content**:
- Metabolic reactions
- Gene-protein-reaction (GPR) associations
- Metabolite identifiers
- Stoichiometry

**Usage**: Reporter Metabolite Analysis (Module 4)

**Purpose**:
- Links gene expression to metabolic pathways
- Identifies metabolites with coordinated gene regulation
- Provides mechanistic interpretation

**Important Notes**:
- Large file (~20 MB)
- Requires rsbml or similar parser
- Analysis is computationally intensive (demonstration only)

---

### 8. Pre-computed Output Files

These files are provided as examples of expected outputs:

#### DESeq output.txt
- Complete differential expression results
- Used as input for GSEA and reporter metabolite analysis

#### Piano output.txt
- Gene set enrichment results
- Pathway-level statistics

#### Coexpression_113genes.txt
- Significant gene-gene correlations
- Edge list for network visualization

#### DESeq output for RAVEN.txt
- Simplified format for metabolic analysis
- pvalue and log2FoldChange only

#### Networks.cys
- Cytoscape session file
- Pre-built network visualization

---

## 🔍 Data Quality and Preprocessing

### RNA-seq Processing Pipeline

1. **Raw sequencing reads** (FASTQ)
   ↓
2. **Quality control** (FastQC, trimming)
   ↓
3. **Alignment** to human genome (STAR or HISAT2)
   ↓
4. **Quantification** (featureCounts or HTSeq)
   ↓
5. **Gene-level counts** → `Counts_selected.txt`
   ↓
6. **FPKM normalization** → `FPKM_selected.txt`

### Quality Filters Applied

- **Gene filter**: Expressed in at least 10% of samples
- **Sample filter**: Library size > 5M reads
- **Outlier removal**: Based on PCA
- **Annotation**: Human genome GRCh38/hg38

---

## 📈 Expected Results Summary

### Differential Expression (Module 1)
- **Significant genes** (padj < 0.05): ~2,000-3,000
- **Up-regulated in Late**: ~1,200 genes
- **Down-regulated in Late**: ~1,100 genes
- **Median |log2FC|**: 0.5-1.0

### GSEA (Module 2)
- **Enriched pathways** (FDR < 0.05): ~50-100
- **Top processes**:
  - Cell cycle (up in Late)
  - Immune response (up in Late)
  - Differentiation (down in Late)

### Co-expression (Module 3)
- **Significant correlations** (FDR < 0.05): ~200-400
- **Positive correlations**: ~60%
- **Negative correlations**: ~40%
- **Network modules**: 3-5 distinct clusters

---

## ⚠️ Important Considerations

### Sample Size
- 50 samples is modest for RNA-seq
- Sufficient for demonstration purposes
- Real studies often use 100+ samples

### Batch Effects
- Samples may come from different batches
- In real analysis, check and correct for batch effects
- This dataset is pre-processed to minimize batch effects

### Confounding Variables
- Only "Type" (Early/Late) is provided
- Real data has age, sex, treatment, etc.
- These can confound results if not accounted for

### Multiple Testing
- With ~20,000 genes, FDR correction is essential
- padj (adjusted p-value) accounts for this
- Always use padj, not raw p-value

---

## 📚 Additional Resources

### TCGA Data Access
- [TCGA Portal](https://portal.gdc.cancer.gov/)
- [TCGA Documentation](https://docs.gdc.cancer.gov/)

### File Format Specifications
- [GMT Format](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29)
- [SBML](http://sbml.org/Documents/Specifications)

### Data Processing
- [RNA-seq Best Practices](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4728800/)
- [TCGA mRNA Pipeline](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/)

---

## 📧 Questions?

If you have questions about the data:
- Check CLAUDE.md for analysis workflow
- Review tutorial modules for usage examples
- Contact instructor: [email]

---

**Data Version**: 1.0
**Last Updated**: December 2025

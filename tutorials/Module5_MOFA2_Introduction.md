# Module 5: Introduction to MOFA2

**Duration**: 30 minutes
**Difficulty**: Intermediate

---

## 📖 Learning Objectives

By the end of this module, you will be able to:
- Understand the concept of multi-omics factor analysis
- Explain when and why to use MOFA2
- Distinguish MOFA2 from other dimensionality reduction methods
- Understand the MOFA2 workflow and key concepts

---

## 🧬 What is Multi-Omics Integration?

### The Challenge

Modern biology generates multiple data types ("omics") from the same samples:
- **Transcriptomics**: RNA expression levels
- **Proteomics**: Protein abundance
- **Metabolomics**: Metabolite concentrations
- **Epigenomics**: DNA methylation, histone modifications
- **Genomics**: Mutations, copy number variations

**Problem**: Each omics layer provides a partial view of the biological system. How do we:
1. Integrate these different data types?
2. Identify shared vs. unique biological signals?
3. Understand which omics layers drive phenotypic variation?

---

## 🎯 What is MOFA2?

**MOFA2** (Multi-Omics Factor Analysis version 2) is an **unsupervised machine learning method** that:

1. **Integrates** multiple omics datasets measured on the same samples
2. **Identifies** latent factors that capture variation across and within omics layers
3. **Decomposes** variation into shared (multi-omics) and specific (single-omics) components
4. **Prioritizes** which features (genes, proteins, etc.) drive each factor

### Key Paper
Argelaguet et al., "MOFA+: a statistical framework for comprehensive integration of multi-modal single-cell data" *Genome Biology* (2020)

---

## 🔄 MOFA2 vs Other Methods

### Comparison with Common Approaches

| Method | Multi-Omics? | Unsupervised? | Interpretable? | Variance Explained? |
|--------|--------------|---------------|----------------|---------------------|
| **PCA** | ❌ Single omics | ✅ Yes | ⚠️ Moderate | ✅ Yes |
| **NMF** | ❌ Single omics | ✅ Yes | ✅ Yes | ⚠️ Indirect |
| **CCA** | ✅ Two omics | ✅ Yes | ⚠️ Moderate | ❌ No |
| **PLS** | ✅ Two omics | ❌ Supervised | ⚠️ Moderate | ⚠️ Limited |
| **MOFA2** | ✅ Multi-omics | ✅ Yes | ✅ Yes | ✅ Yes |

### What Makes MOFA2 Unique?

1. **Multi-view learning**: Handles >2 omics layers simultaneously
2. **Variance decomposition**: Quantifies how much each omics contributes to each factor
3. **Sparsity**: Automatically identifies relevant features (genes/proteins)
4. **Missing data**: Handles incomplete measurements across omics layers
5. **Groups**: Can analyze multiple sample groups (e.g., tissues, timepoints)

---

## 🧮 The MOFA2 Model

### Conceptual Framework

```
Multi-Omics Data → Latent Factors → Biological Interpretation

View 1 (mRNA)     ┐
View 2 (Protein)  ├─→ Factors 1, 2, 3, ... K → Biological processes
View 3 (Metabol)  ┘                           Cell types
                                               Disease states
```

### Mathematical Representation

For each omics view *m* and sample *n*:

**Data** ≈ **Weights** × **Factors**

Where:
- **Data** (Y): Original measurements (genes × samples)
- **Factors** (Z): Latent factors (K factors × samples)
- **Weights** (W): Feature loadings (genes × K factors)

### What Are Latent Factors?

**Latent factors** are hidden variables that explain observed variation:
- Each factor captures a **biological process** or **source of variation**
- Factors are **shared** across omics or **specific** to one omics layer
- Similar to principal components but **multi-view aware**

**Example**: Factor 1 might represent:
- Cell cycle progression (visible in mRNA, protein, metabolomics)
- High in proliferating cells, low in quiescent cells

---

## 📊 MOFA2 Workflow Overview

### Step 1: Data Preparation
- Prepare multiple omics matrices (features × samples)
- Select highly variable features per omics
- Normalize data appropriately

### Step 2: Model Training
- Create MOFA object with all omics layers
- Set model parameters (number of factors, sparsity)
- Train model to convergence

### Step 3: Factor Interpretation
- Identify significant factors
- Examine variance explained per omics
- Correlate factors with phenotypes

### Step 4: Feature Analysis
- Identify top-weighted features per factor
- Perform enrichment analysis on feature sets
- Validate key findings

---

## 🎨 Types of Factors

MOFA2 identifies different factor patterns:

### 1. Shared Factors
**Active across multiple omics layers**

```
mRNA:      ████████████ (high variance explained)
Protein:   ████████████ (high variance explained)
Metabol:   ████████     (moderate variance)

→ Represents coordinated multi-omics process
```

**Example**: Cell differentiation affects transcription, translation, and metabolism

### 2. Omics-Specific Factors
**Active in only one omics layer**

```
mRNA:      ████████████ (high variance explained)
Protein:   ░░░░░░░░░░░░ (no variance)
Metabol:   ░░░░░░░░░░░░ (no variance)

→ Represents post-transcriptional regulation
```

**Example**: microRNA regulation affects mRNA but not protein levels

### 3. Partially Shared Factors
**Active in subset of omics**

```
mRNA:      ████████████ (high variance)
Protein:   ████████     (moderate variance)
Metabol:   ░░░░░░░░░░░░ (no variance)

→ Represents processes with translation lag
```

---

## 💡 When to Use MOFA2

### ✅ Ideal Use Cases

1. **Exploratory analysis** of multi-omics data
   - You have ≥2 omics layers on same samples
   - You want to discover major sources of variation
   - No pre-defined hypotheses

2. **Understanding omics relationships**
   - Which omics layers are correlated?
   - Which biological processes are multi-omics vs. omics-specific?
   - What drives sample heterogeneity?

3. **Biomarker discovery**
   - Identify multi-omics signatures of disease
   - Prioritize features for validation
   - Understand which omics is most informative

4. **Quality control**
   - Identify technical artifacts
   - Detect batch effects across omics
   - Find outlier samples

### ❌ When NOT to Use MOFA2

1. **Single omics dataset** → Use PCA or NMF instead
2. **Supervised classification** → Use machine learning classifiers
3. **Causal inference** → MOFA2 finds associations, not causation
4. **Small sample size** (n < 10-15) → Insufficient statistical power
5. **No biological replication** → Cannot distinguish signal from noise

---

## 🔬 Example Application: Cancer Multi-Omics

### Research Question
*What multi-omics factors distinguish cancer subtypes?*

### Input Data (50 tumor samples)
1. **mRNA expression**: 2,000 highly variable genes
2. **Protein abundance**: 500 proteins
3. **Metabolomics**: 200 metabolites

### Expected Output
- **Factor 1**: Proliferation signature (mRNA + protein)
- **Factor 2**: Metabolic reprogramming (all omics)
- **Factor 3**: Immune infiltration (mRNA + protein)
- **Factor 4**: Post-transcriptional regulation (mRNA only)

### Interpretation
- Factors 1-3 are **shared** → coordinate regulation across omics
- Factor 4 is **omics-specific** → regulatory layer uncoupled from protein/metabolite

---

## 📐 Key MOFA2 Parameters

You'll encounter these in the next modules:

### Data Options
- **scale_views**: Scale each omics to unit variance?
- **scale_groups**: Scale sample groups separately?

### Model Options
- **num_factors**: How many factors to infer? (default: 15, auto-prunes)
- **likelihoods**: Data type per omics (gaussian, bernoulli, poisson)
- **sparsity_prior**: Enforce feature sparsity? (default: yes)

### Training Options
- **max_iterations**: Maximum training iterations
- **convergence_mode**: When to stop (fast vs slow)
- **gpu_mode**: Use GPU acceleration?

**Don't worry!** We'll explore these in detail in Module 7.

---

## 🧪 Hands-On Concept Check

### Question 1
You have RNA-seq and proteomics data from 30 samples. After running MOFA2, you find:

```
Factor 1 variance explained:
  - RNA: 25%
  - Protein: 23%

Factor 2 variance explained:
  - RNA: 15%
  - Protein: 2%
```

**Questions**:
1. Is Factor 1 shared or omics-specific?
2. Is Factor 2 shared or omics-specific?
3. What biological scenario could explain Factor 2?

---

### Question 2
You want to analyze:
- mRNA expression (20,000 genes)
- DNA methylation (450,000 CpG sites)
- Proteomics (5,000 proteins)

**Questions**:
1. Should you include all features? Why or why not?
2. How would you select features for each omics?
3. Should the number of features be similar across omics?

---

### Question 3
Rank these scenarios from most to least suitable for MOFA2:

A. 5 samples with RNA-seq and proteomics
B. 100 samples with RNA-seq only
C. 50 samples with RNA-seq, proteomics, and metabolomics
D. 30 samples with RNA-seq and clinical variables

---

## 🎯 Key Takeaways

1. ✅ MOFA2 is an **unsupervised** multi-omics integration method
2. ✅ It identifies **latent factors** explaining variation across omics
3. ✅ Factors can be **shared** (multi-omics) or **omics-specific**
4. ✅ MOFA2 provides **variance decomposition** per omics layer
5. ✅ Requires **≥2 omics** on **same samples** with **n ≥ 15-20**
6. ✅ Best for **exploratory analysis**, not supervised prediction

---

## 📚 Additional Resources

### Papers
- [MOFA+ Paper (Genome Biology, 2020)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02015-1)
- [Original MOFA Paper (Mol Sys Biol, 2018)](https://www.embopress.org/doi/full/10.15252/msb.20178124)

### Documentation
- [MOFA2 Website](https://biofam.github.io/MOFA2/)
- [MOFA2 Tutorials](https://biofam.github.io/MOFA2/tutorials.html)
- [MOFA2 FAQ](https://biofam.github.io/MOFA2/faq.html)

### Videos
- [MOFA2 Webinar (YouTube)](https://www.youtube.com/watch?v=lPQSV-f2YH8)

---

## ⏭️ Next Module

Continue to **Module 6: MOFA2 Data Preparation** where we'll prepare our RNA-seq and other omics data for MOFA2 analysis.

**Solutions to concept checks** are available in `solutions/Module5_solutions.md`

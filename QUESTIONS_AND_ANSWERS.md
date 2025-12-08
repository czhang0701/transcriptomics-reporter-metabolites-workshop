# Workshop Questions and Answers

Quick reference guide for discussion questions throughout the workshop.

---

## Module 1: Differential Expression (DESeq2)

### Q1: Why use negative binomial distribution instead of normal distribution?

**Answer**: RNA-seq count data has three key properties:
- **Discrete values** (counts: 0, 1, 2, 3... not continuous)
- **Non-negative** (can't have negative counts)
- **Overdispersed** (variance > mean)

The normal distribution assumes continuous data and variance = mean. The negative binomial distribution is designed for count data with overdispersion (common in biological data due to technical and biological variability).

---

### Q2: What does a positive log2FoldChange mean in Late vs Early comparison?

**Answer**:
- **Positive log2FC** → Gene is **higher** in Late samples (numerator)
- **Negative log2FC** → Gene is **higher** in Early samples (denominator)

Example: log2FC = +2 means the gene is 2² = 4× higher in Late vs Early.

---

### Q3: Why use adjusted p-values (FDR) instead of raw p-values?

**Answer**: We test ~20,000 genes simultaneously. With p-value < 0.05:
- Expected false positives = 20,000 × 0.05 = **1,000 false discoveries!**
- FDR correction (Benjamini-Hochberg) controls the expected proportion of false discoveries
- FDR < 0.05 means we expect ≤ 5% of significant results to be false positives

---

### Q4: Why do some genes have high fold-change but are NOT significant?

**Answer**: Three main reasons:

1. **High variance**: Gene expression varies a lot between samples → Less confidence
2. **Small sample size**: Few replicates → Can't distinguish signal from noise
3. **Low counts**: Lowly expressed genes have unreliable fold changes

**Key insight**: Significance considers BOTH effect size (fold change) AND confidence (variance, sample size). A large fold change with high uncertainty is not statistically significant.

---

## Module 2: Gene Set Enrichment (GSEA)

### Q1: Why analyze pathways instead of individual genes?

**Answer**:

**Advantages**:
- **Better signal-to-noise**: Many weakly changing genes → Strong pathway signal
- **Biological interpretability**: "Glycolysis is upregulated" vs "PFKFB3 is up"
- **Statistical power**: Aggregating genes increases detection sensitivity
- **Functional context**: Understand coordinated biological processes

**Example**: 10 glycolysis genes each with log2FC = 0.5 (not individually significant) → Pathway is highly significant!

---

### Q2: What does "enriched in up-regulated genes" mean?

**Answer**: Many genes in this pathway are **coordinately upregulated**, suggesting:
- The pathway is **activated** in Late vs Early samples
- The biological process is **more active**
- Transcriptional program is driving increased function

**Example**: "GLYCOLYSIS" enriched in up-regulated genes → Glycolysis is activated in Late-stage cancer

---

### Q3: What is the Warburg effect?

**Answer**:

**Warburg Effect** = Cancer cells preferentially use glycolysis even when oxygen is available (normally cells use oxidative phosphorylation when oxygen is present).

**Evidence in your data**:
- **GLYCOLYSIS** pathways: Enriched in **up-regulated** genes
- **OXIDATIVE_PHOSPHORYLATION**: Enriched in **down-regulated** genes

**Why cancer cells do this**:
- Produces ATP faster (though less efficiently)
- Generates biosynthetic precursors (nucleotides, lipids, amino acids)
- Acidifies tumor microenvironment (lactate export)

---

### Q4: Why do different databases (GO vs KEGG) give different results?

**Answer**:

| Database | Focus | Granularity | Curation |
|----------|-------|-------------|----------|
| **GO Biological Process** | Broad biological concepts | Very specific (1000s of terms) | Community-driven |
| **KEGG** | Metabolic/signaling pathways | Coarse (hundreds of pathways) | Expert-curated |
| **Reactome** | Biochemical reactions | Intermediate | Expert-curated |

**Result**: Same biological process may be represented differently:
- GO: "GLYCOLYTIC_PROCESS" (specific)
- KEGG: "Glycolysis / Gluconeogenesis" (broader, includes reverse pathway)

**Recommendation**: Use multiple databases for comprehensive analysis.

---

## Module 3: Reporter Metabolites ⭐

### Q1: What is the difference between Reporter Metabolites and GSEA?

**Answer**:

| Method | Input | Output | Connectivity |
|--------|-------|--------|--------------|
| **GSEA** | Genes → Predefined gene sets | Enriched pathways | No network structure |
| **Reporter Metabolites** | Genes → Metabolic model (network) | Key metabolites | Uses reaction network topology |

**Key difference**: Reporter metabolites use the **metabolic network structure** to identify metabolites where neighboring genes (enzymes) are coordinately regulated.

**Analogy**:
- GSEA: "Which neighborhoods (pathways) have many changed houses (genes)?"
- Reporter: "Which intersections (metabolites) have many changed roads (reactions) nearby?"

---

### Q2: Why aggregate Z-scores instead of counting significant genes?

**Answer**:

**Z-score aggregation** (what we do):
- Uses **continuous information** (P-values → Z-scores)
- Captures **magnitude** of changes (strong vs weak signals)
- More **statistically powerful** (detects subtle coordinated changes)

**Counting significant genes** (simpler approach):
- Binary classification (significant vs not)
- **Loses information** about effect size
- Sensitive to arbitrary thresholds (P < 0.05 vs P = 0.051)

**Example**:
- Metabolite A: 20 genes, all P = 0.04 (just significant)
- Metabolite B: 10 genes, all P = 1e-10 (extremely significant)
- Counting: Metabolite A "wins" (20 > 10)
- Z-score: Metabolite B wins (stronger signal)

---

### Q3: ⭐ Why might a metabolite appear in BOTH up AND down analyses?

**Answer**:

A metabolite can be significant in both directions when genes producing vs consuming it are regulated oppositely.

**Mechanisms**:

1. **Metabolic flux reversal**:
   - Forward reaction genes UP, reverse reaction genes DOWN

2. **Compartmentalization**:
   - Cytosolic enzymes UP, mitochondrial enzymes DOWN

3. **Competing pathways**:
   - Production pathway UP, consumption pathway DOWN
   - Result: Metabolite accumulates!

4. **Anabolic vs Catabolic**:
   - Biosynthesis genes UP, degradation genes DOWN

5. **Feedback regulation**:
   - High metabolite → Compensatory downregulation of production

**Hypothetical Example**: Lactate
- IF production enzymes (LDHA) are upregulated
- AND export transporters (SLC16A3) are upregulated
- BUT consumption enzymes are downregulated
- THEN: Lactate appears in both up/down analyses!

**Key insight**: This is **biologically meaningful** - it suggests metabolite accumulation or flux reprogramming!

---

### Q4: If a metabolite has high Z-score but only 3 neighboring genes, should you trust it?

**Answer**: **Be cautious!**

**Concerns**:
- **Small sample size** → Low statistical power
- **Sensitive to outliers** → One highly significant gene dominates
- **Less confidence** in biological interpretation

**When to trust**:
- All 3 genes have **very strong** signals (P < 1e-6)
- Genes are **well-known** key enzymes for that metabolite
- Result is **consistent** with literature

**When to doubt**:
- Mixed signals (1 very significant, 2 not)
- Genes have unclear metabolic roles
- No literature support

**Recommendation**: Cross-validate with:
- Other enrichment methods (GSEA)
- Literature review
- Experimental validation

---

### Q5: How would you validate a reporter metabolite experimentally?

**Answer**:

**1. Targeted Metabolomics**:
- LC-MS or GC-MS to measure metabolite concentration
- Compare Early vs Late samples
- Expect: Significant metabolite should show concentration differences

**2. Isotope Tracing**:
- Feed cells ¹³C-glucose or ¹³C-glutamine
- Track labeled metabolite production
- Measures actual metabolic flux

**3. Enzyme Activity Assays**:
- Measure activity of key enzymes producing/consuming metabolite
- Validate that mRNA changes → protein activity changes

**4. Perturbation Experiments**:
- Knockdown/overexpress key enzyme
- Measure effect on metabolite level and cell phenotype

**5. Imaging Mass Spectrometry**:
- Spatial distribution of metabolite in tumor tissue
- See if it correlates with gene expression patterns

---

## Module 4: Co-expression Networks

### Q1: What does positive correlation between two genes mean biologically?

**Answer**:

**Positive correlation** → Genes go up/down together → Suggests:

1. **Co-regulation**: Same transcription factor controls both
2. **Same pathway**: Genes work together in a process
3. **Protein complex**: Gene products physically interact
4. **Similar function**: Genes respond to same signals
5. **Same cell type**: Both expressed in same cells

**Example**: HK2 and PFKP (glycolysis enzymes) positively correlated → Both upregulated together in cancer

---

### Q2: Why filter to top 10% of correlations?

**Answer**:

**Problem with keeping all significant correlations**:
- Network is too **dense** (too many edges)
- Difficult to detect **modules** (communities)
- Weak correlations add **noise**

**Solution - Top 10% filtering**:
- Keep only **strongest** relationships
- Makes modules more **distinct**
- Focuses on **core** regulatory relationships
- Improves **visualization** and interpretation

**Analogy**: In a social network, focus on close friends (strong ties) rather than all acquaintances.

---

### Q3: What is a "module" and why is it biologically meaningful?

**Answer**:

**Module** = Group of genes more connected to **each other** than to rest of network

**Biological meaning**:
- Genes in same **functional pathway**
- **Co-regulated** by same mechanism
- Form **protein complexes**
- Respond to same **biological signal**

**Example**:
- Module 1: Cell cycle genes (all peak during mitosis)
- Module 2: Immune response genes (all induced by interferon)
- Module 3: Metabolic genes (all regulated by nutrient availability)

**Why important**: Modules represent **functional units** - can predict gene function, identify drug targets, understand disease mechanisms.

---

### Q4: What does modularity measure? Is 0.3-0.5 good?

**Answer**:

**Modularity** = How well network divides into distinct modules

- **Range**: -0.5 to +1.0
- **Interpretation**:
  - < 0.3: Poor module structure (random)
  - 0.3-0.5: **Good** module structure
  - 0.5-0.7: **Very good** structure
  - > 0.7: Excellent (rare in biological networks)

**Your result (~0.3-0.5)** → **Good module structure**, suggests biologically meaningful communities.

---

### Q5: Why analyze positive and negative correlations separately?

**Answer**:

**Positive correlations**:
- Genes **co-activated** together
- Same biological process
- Co-regulated by same mechanism

**Negative correlations**:
- Genes show **antagonistic** regulation
- Different cell states (proliferation vs differentiation)
- Metabolic **trade-offs** (glycolysis vs OXPHOS)
- **Feedback inhibition** loops

**Example**:
- Positive module: All cell cycle genes upregulated together
- Negative module: Proliferation genes (UP) vs differentiation genes (DOWN)

**Different biology** → Analyze separately!

---

### Q6: Module with 100 genes but no pathway enrichment - what does this mean?

**Answer**:

**Possible explanations**:

1. **Novel functional module**: Not yet annotated in databases
   - **Action**: Publish! You may have found new biology

2. **Technical artifact**:
   - Batch effects
   - Chromosomal location (genes on same chromosome)
   - **Action**: Check for technical confounders

3. **Database incomplete**:
   - Pathway exists but not in your gene set database
   - Try different databases (GO, KEGG, Reactome)

4. **Weak enrichment**:
   - Module has biological function but not strongly enriched
   - **Action**: Try GO Molecular Function, Cellular Component

5. **Non-functional criteria**:
   - Genes cluster by expression level, not function
   - **Action**: Check if module is just "highly expressed genes"

**Recommendation**: Don't discard! Investigate with literature review, check hub genes.

---

### Q7: Why are hub genes good drug targets?

**Answer**:

**Hub gene** = Highly connected gene in network (many edges)

**Why good drug targets**:

1. **Central to pathway**: Affects many downstream processes
2. **High impact**: Perturbing hub → Large network effect
3. **Essential function**: Often required for cell survival
4. **Biomarkers**: Hub expression predicts pathway activity

**Caution**:
- **Toxicity risk**: Hub in normal cells too → Side effects
- **Redundancy**: Other genes may compensate

**Examples**:
- HIF1A: Hub in hypoxia response → Target for anti-angiogenesis
- MYC: Hub in proliferation → Target for cancer therapy

**Best targets**: Hubs **specific to disease** (cancer) but not essential in normal tissue.

---

### Q8: ⭐ How to integrate network modules with reporter metabolites?

**Answer**:

**Integration strategies**:

**1. Module → Metabolite enrichment**:
- For each module, check if genes are enriched for specific metabolite associations
- Example: Module 1 genes highly enriched for "glucose" neighborhood → Module is glycolysis

**2. Hub genes → Metabolic enzymes**:
- Check if hub genes encode enzymes for reporter metabolites
- Example: LDHA is hub in Module 2 + lactate is reporter metabolite → Lactate production is key

**3. Compare enrichments**:
- Module enriched for "GLYCOLYSIS" pathway
- Reporter metabolites: pyruvate, lactate
- **Convergent evidence** → High confidence in glycolysis reprogramming

**4. Mechanistic stories**:
- Module 1: Upregulated glycolysis genes
- Reporter: Lactate (end product)
- Integration: "Coordinated upregulation of glycolysis genes → lactate accumulation"

**Analogy**:
- Modules = "Factories" (groups of workers)
- Reporter metabolites = "Products" (what factories produce)
- Integration = Connect factories to products!

---

## Integration Across All Modules

### Q: How do the 4 modules complement each other?

**Answer**:

| Module | Question | Output | Strength |
|--------|----------|--------|----------|
| **1. DESeq2** | Which genes change? | 2,000 significant genes | Identifies **WHAT** changes |
| **2. GSEA** | Which pathways change? | 50 enriched pathways | Groups genes into **PROCESSES** |
| **3. Reporter Metabolites** | Which metabolites are key? | 10-20 key metabolites | Finds **METABOLIC HUBS** |
| **4. Networks** | How are genes connected? | 5-10 modules | Reveals **COORDINATED REGULATION** |

**Integration story**:
1. DESeq2: "LDHA, PFKP, HK2 are upregulated"
2. GSEA: "GLYCOLYSIS pathway is enriched"
3. Reporter Metabolites: "Lactate is a key metabolite"
4. Networks: "Glycolysis genes form a co-expression module"

**Conclusion**: "Cancer cells coordinately upregulate glycolysis genes → increased lactate production (Warburg effect)"

**Each module adds unique information** → Together they tell the complete biological story!

---

## Tips for Discussing Results

### Do:
✅ Compare YOUR results with these questions
✅ Discuss unexpected findings with your neighbor
✅ Think about biological mechanisms
✅ Connect results across modules
✅ Question assumptions

### Don't:
❌ Memorize "right answers"
❌ Assume your data will exactly match examples
❌ Accept results without critical thinking
❌ Forget to consider technical artifacts

### Key principle:
**"Trust but verify"** → Use multiple lines of evidence (DESeq2 + GSEA + Reporter + Networks) for strong conclusions.

---

## For After the Workshop

### Next Steps:
1. **Apply to your data**: Use these questions to interrogate your own results
2. **Read original papers**: Patil & Nielsen (2005) for reporter metabolites
3. **Try variations**: Different gene sets, parameters, databases
4. **Validate findings**: Plan experiments to test key predictions

### Resources:
- Detailed explanations: `REFLECTION_QUESTIONS.md`
- Step-by-step tutorials: `tutorials/` folder
- Technical details: `CLAUDE.md`

---

**Remember**: These questions are designed to develop **critical thinking** about your data. There's rarely one "right" answer - biology is complex! The goal is to build intuition and ask better questions.

# Workshop Reflection Questions

## Purpose
These questions are designed to help you think critically about the analyses and develop deeper understanding. Discuss with your neighbors or think through them individually.

---

## Module 1: Differential Expression Analysis (DESeq2)

### Understanding the Results

**Q1.1**: Why do we use negative binomial distribution for RNA-seq count data instead of normal distribution?
<details>
<summary>Hint</summary>
Think about the properties of count data (discrete, non-negative) and overdispersion.
</details>

**Q1.2**: What does a positive log2FoldChange mean in our Late vs Early comparison?
<details>
<summary>Answer</summary>
Positive log2FC means higher expression in Late samples (numerator in Late/Early).
</details>

**Q1.3**: Why do we use adjusted p-values (FDR) instead of raw p-values?
<details>
<summary>Hint</summary>
We're testing 20,000 genes simultaneously. What happens to false positives?
</details>

**Q1.4**: Look at your volcano plot. Why do some genes have high fold-change but are not significant?
<details>
<summary>Think about</summary>
What role does variance play? What about sample size?
</details>

### Critical Thinking

**Q1.5**: If you found 5,000 differentially expressed genes, is that a lot or a little? How would you decide?

**Q1.6**: What biological processes might cause genes to be downregulated in late-stage cancer?

---

## Module 2: Gene Set Enrichment Analysis (GSEA/PIANO)

### Understanding Enrichment

**Q2.1**: Why do we analyze pathways instead of just looking at individual genes?
<details>
<summary>Think about</summary>
Signal-to-noise ratio, biological interpretability, statistical power.
</details>

**Q2.2**: What does it mean when a pathway is enriched in "up-regulated" genes?
<details>
<summary>Answer</summary>
Many genes in this pathway are coordinately upregulated, suggesting the pathway is activated.
</details>

**Q2.3**: Why might different pathway databases (GO vs KEGG) give different results?
<details>
<summary>Consider</summary>
Pathway definitions, gene coverage, granularity, curation methods.
</details>

### Critical Analysis

**Q2.4**: You see "OXIDATIVE_PHOSPHORYLATION" enriched in down-regulated genes. What does this tell you about cancer metabolism?
<details>
<summary>Biological context</summary>
Warburg effect - cancer cells often shift from oxidative phosphorylation to glycolysis.
</details>

**Q2.5**: Why might some highly significant individual genes NOT appear in any enriched pathway?

**Q2.6**: If a pathway has 500 genes but only 10 are in your dataset, should you trust the enrichment?

---

## Module 3: Reporter Metabolites Analysis ⭐

### Understanding the Algorithm

**Q3.1**: What is the difference between Reporter Metabolites and GSEA?
<details>
<summary>Key difference</summary>
GSEA: Gene sets → Pathways
Reporter: Gene sets → Metabolites (network topology-based)
</details>

**Q3.2**: Why do we aggregate Z-scores instead of just counting significant genes?
<details>
<summary>Think about</summary>
Continuous vs binary information, statistical power, magnitude of effects.
</details>

**Q3.3**: What does "background correction" mean and why is it necessary?
<details>
<summary>Hint</summary>
Different metabolites have different numbers of associated genes. How does this bias the analysis?
</details>

### Critical Interpretation

**Q3.4**: ⭐ **Why might a metabolite show up as significant in BOTH up-regulated AND down-regulated analyses?**
<details>
<summary>Biological explanations</summary>

1. **Metabolic flux reversal**: Genes catalyzing forward/reverse reactions regulated differently
2. **Compartmentalization**: Same metabolite in different compartments (cytosol vs mitochondria)
3. **Competing pathways**: Multiple pathways producing/consuming the metabolite
4. **Regulatory complexity**: Anabolic genes up, catabolic genes down (or vice versa)
5. **Feedback loops**: Metabolite accumulation triggers compensatory responses

**Example**: Lactate
- Lactate PRODUCTION genes (LDHA) upregulated → Warburg effect
- Lactate IMPORT genes (MCT1) downregulated → Cancer cells export lactate
- Result: Lactate appears in both up AND down analyses!
</details>

**Q3.5**: If a metabolite has a very high Z-score but only 3 neighboring genes, should you trust it?

**Q3.6**: How would you validate a reporter metabolite finding experimentally?
<details>
<summary>Experimental approaches</summary>

- Targeted metabolomics (LC-MS/GC-MS)
- Isotope tracing experiments
- Enzyme activity assays
- Measure metabolite concentrations in samples
</details>

### Integration Questions

**Q3.7**: Compare your top reporter metabolites with enriched pathways from Module 2. Do they match?

**Q3.8**: If "glucose" is a top reporter metabolite, what would you expect to see in Module 1 (DESeq2)?
<details>
<summary>Expected genes</summary>
GLUT transporters, HK (hexokinase), PFKFB, etc.
</details>

**Q3.9**: Why do we use Ensembl IDs for matching genes to the model instead of gene symbols?

---

## Module 4: Advanced Co-expression Network Analysis

### Understanding Networks

**Q4.1**: What does a positive correlation between two genes tell you biologically?
<details>
<summary>Biological interpretation</summary>

- Co-regulation by same transcription factor
- Part of same protein complex or pathway
- Responding to same environmental stimulus
- Similar degradation rates
</details>

**Q4.2**: Why do we filter to the top 10% of correlations instead of keeping all significant ones?
<details>
<summary>Network analysis rationale</summary>
Too many edges → Dense network → Hard to detect modules. Strong correlations = stronger biological signal.
</details>

**Q4.3**: What is a "module" and why is it biologically meaningful?
<details>
<summary>Module concept</summary>
A module is a group of genes that are more connected to each other than to the rest of the network. Suggests functional unit or co-regulated pathway.
</details>

### Module Detection

**Q4.4**: What does "modularity" measure? Is 0.3 a good modularity score?
<details>
<summary>Modularity explanation</summary>
Modularity measures how well-separated the modules are. Values > 0.3 suggest good module structure.
</details>

**Q4.5**: Why might we analyze positive and negative correlations separately?
<details>
<summary>Biological distinction</summary>

- **Positive**: Co-activation, same pathway
- **Negative**: Antagonistic regulation, different cell states, metabolic trade-offs
</details>

**Q4.6**: If a module is enriched for "CELL_CYCLE" pathways, what does this tell you?

### Critical Analysis

**Q4.7**: You find a module with 100 genes but no significant pathway enrichment. What could this mean?
<details>
<summary>Possibilities</summary>

- Novel functional module not yet annotated
- Technical artifact (batch effect, outlier samples)
- Module defined by non-functional criteria (chromosomal location)
- Pathway database incomplete
</details>

**Q4.8**: Why might hub genes (highly connected) be good drug targets?

**Q4.9**: How would you integrate network modules with reporter metabolites?
<details>
<summary>Integration approach</summary>

1. Check if genes in a module are enriched for specific metabolite associations
2. See if hub genes in modules encode metabolic enzymes
3. Compare module enrichment with reporter metabolite pathways
</details>

---

## Integration Across All Modules

### Big Picture Questions

**I1**: Draw the flow of information through all 4 modules. How does each module build on the previous?
<details>
<summary>Flow</summary>
DESeq2 (genes) → GSEA (pathways) → Reporter Metabolites (metabolic hubs) → Networks (gene relationships)
</details>

**I2**: If you could only do 2 of the 4 modules, which would you choose and why?

**I3**: How do the different modules complement each other? What unique information does each provide?

**I4**: Imagine you're presenting to a clinician. How would you explain the reporter metabolite results in plain language?

### Experimental Design

**I5**: What additional data would help validate your findings from this workshop?
<details>
<summary>Data types</summary>

- Metabolomics data (validate reporter metabolites)
- Proteomics data (check if mRNA changes match protein)
- Additional timepoints (track temporal dynamics)
- Perturbation experiments (knockdown/overexpression)
- Independent validation cohort
</details>

**I6**: If you were to design a follow-up experiment, what would you target?

**I7**: How would this analysis change if you were comparing:
- Different tissues instead of disease stages?
- Treatment vs control instead of Early vs Late?
- Mouse data instead of human data?

---

## Advanced Discussion Topics

### Statistical Considerations

**A1**: We used FDR < 0.05 as a threshold. How would results change with FDR < 0.01? Which is "better"?

**A2**: Why do we use Spearman correlation instead of Pearson for co-expression networks?

**A3**: In reporter metabolites, why is random sampling necessary for background correction?

### Biological Depth

**A4**: The Warburg effect is famous in cancer. Can you identify evidence for it in your results?

**A5**: How might tumor heterogeneity affect your differential expression results?

**A6**: Why might some metabolic changes be adaptive (helping cancer cells) vs. consequential (side effects)?

### Methodological Critique

**A7**: What are the limitations of using only transcriptomics data to infer metabolic changes?
<details>
<summary>Limitations</summary>

- mRNA ≠ protein ≠ enzyme activity
- Post-translational modifications
- Allosteric regulation
- Metabolite concentrations affect flux independent of enzymes
</details>

**A8**: How could batch effects in RNA-seq data affect your conclusions?

**A9**: What assumptions does the reporter metabolites algorithm make? Are they reasonable?

---

## Preparing for Your Own Analysis

**P1**: What would you do differently if you had:
- 1000 samples instead of 50?
- 10 samples instead of 50?
- Paired samples (same patient, before/after)?

**P2**: How would you adapt this workflow for:
- Microbiome data?
- Proteomics data?
- Single-cell RNA-seq?

**P3**: What quality control steps should you add before running this analysis on your own data?

---

## Discussion Prompts for Instructors

**During Module 1**: "Look at your neighbor's volcano plot. Do you see the same patterns?"

**During Module 2**: "Which pathway surprises you the most? Why?"

**During Module 3**: "Compare your top 5 reporter metabolites with your neighbor. Discuss why they might be important in cancer."

**During Module 4**: "If you could validate one module experimentally, which would it be and how?"

**Final Discussion**: "What is the single most interesting biological insight you gained today?"

---

## Take-Home Reflection

After the workshop, write a brief summary (1 paragraph each):

1. **What I learned**: Key technical skills and biological insights
2. **What surprised me**: Unexpected results or concepts
3. **What I'll use**: How this applies to my research
4. **What I still don't understand**: Topics for further study
5. **Next steps**: How I would extend this analysis

---

## Additional Resources for Deep Dives

**Reporter Metabolites**:
- Original paper: Patil & Nielsen (2005) PNAS
- RAVEN toolbox documentation
- Genome-scale metabolic models (HMR database)

**Network Analysis**:
- Louvain/Leiden algorithm papers
- WGCNA package for weighted correlation networks
- Cytoscape tutorials

**Integration**:
- Multi-omics integration methods
- MOFA2 for unsupervised integration
- Network-based stratification

---

**Remember**: The goal is not to get "right answers" but to think critically about your data and develop biological intuition!

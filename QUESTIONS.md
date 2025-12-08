# Workshop Discussion Questions

Questions for reflection and discussion throughout the workshop.

**Purpose**: These questions are designed to encourage critical thinking and deeper understanding. Discuss with your neighbors or think through them individually.

---

## Module 1: Differential Expression (DESeq2)

### Q1: Why use negative binomial distribution instead of normal distribution?

### Q2: What does a positive log2FoldChange mean in Late vs Early comparison?

### Q3: Why use adjusted p-values (FDR) instead of raw p-values?

### Q4: Why do some genes have high fold-change but are NOT significant?

### Q5 (Bonus): If you found 2,000 differentially expressed genes, is that a lot or a little? How would you decide?

---

## Module 2: Gene Set Enrichment (GSEA)

### Q1: Why analyze pathways instead of individual genes?

### Q2: What does "enriched in up-regulated genes" mean?

### Q3: What is the Warburg effect? How would you identify evidence for it in your results?

### Q4: Why might different pathway databases (GO vs KEGG) give different results?

### Q5 (Bonus): Why might some highly significant individual genes NOT appear in any enriched pathway?

---

## Module 3: Reporter Metabolites ⭐

### Q1: What is the difference between Reporter Metabolites and GSEA?

### Q2: Why do we aggregate Z-scores instead of just counting significant genes?

### Q3: ⭐ IMPORTANT: Why might a metabolite show up as significant in BOTH up-regulated AND down-regulated analyses?

**Think about**:
- Metabolic flux reversal (forward vs reverse reactions)
- Compartmentalization (cytosol vs mitochondria)
- Competing pathways (multiple ways to produce/consume)
- Regulatory complexity (anabolic vs catabolic genes)
- Feedback loops (metabolite accumulation triggers responses)

### Q4: If a metabolite has a very high Z-score but only 3 neighboring genes, should you trust it? Why or why not?

### Q5: How would you validate a reporter metabolite finding experimentally?

### Q6: Compare your top reporter metabolites with enriched pathways from Module 2. Do they match? What does this tell you?

---

## Module 4: Co-expression Networks

### Q1: What does a positive correlation between two genes tell you biologically?

### Q2: Why do we filter to the top 10% of correlations instead of keeping all significant ones?

### Q3: What is a "module" and why is it biologically meaningful?

### Q4: What does "modularity" measure? What range of values indicates good module structure?

### Q5: Why might we analyze positive and negative correlations separately?

### Q6: You found a module with 100 genes but no significant pathway enrichment. What could this mean?

### Q7: Why might hub genes (highly connected) be good drug targets?

### Q8: ⭐ INTEGRATION: How would you integrate network modules with reporter metabolites from Module 3?

---

## Integration Across All Modules

### Q1: Draw the flow of information through all 4 modules. How does each module build on the previous?

### Q2: If you could only do 2 of the 4 modules, which would you choose and why?

### Q3: How do the different modules complement each other? What unique information does each provide?

### Q4: Imagine you're presenting to a clinician. How would you explain the reporter metabolite results in plain language?

---

## Experimental Design

### Q1: What additional data would help validate your findings from this workshop?

**Consider**:
- Metabolomics data
- Proteomics data
- Additional timepoints
- Perturbation experiments
- Independent validation cohort

### Q2: If you were to design a follow-up experiment, what would you target?

### Q3: How would this analysis change if you were comparing:
- Different tissues instead of disease stages?
- Treatment vs control instead of Early vs Late?
- Mouse data instead of human data?

---

## Advanced Discussion Topics

### Statistical Considerations

**A1**: We used FDR < 0.05 as a threshold. How would results change with FDR < 0.01? Which is "better"?

**A2**: Why do we use Spearman correlation instead of Pearson for co-expression networks?

**A3**: In reporter metabolites, why is random sampling necessary for background correction?

---

### Biological Depth

**A4**: The Warburg effect is famous in cancer. Can you identify evidence for it in your results?

**A5**: How might tumor heterogeneity affect your differential expression results?

**A6**: Why might some metabolic changes be adaptive (helping cancer cells) vs. consequential (side effects)?

---

### Methodological Critique

**A7**: What are the limitations of using only transcriptomics data to infer metabolic changes?

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

## Take-Home Reflection

After the workshop, write a brief summary (1 paragraph each):

1. **What I learned**: Key technical skills and biological insights
2. **What surprised me**: Unexpected results or concepts
3. **What I'll use**: How this applies to my research
4. **What I still don't understand**: Topics for further study
5. **Next steps**: How I would extend this analysis

---

## Tips for Discussion

### Do:
✅ Compare YOUR results with these questions
✅ Discuss unexpected findings with your neighbor
✅ Think about biological mechanisms
✅ Connect results across modules
✅ Question assumptions

### Don't:
❌ Expect one "right answer"
❌ Assume your data will exactly match examples
❌ Accept results without critical thinking
❌ Forget to consider technical artifacts

---

## Additional Resources

For detailed answers and explanations, see:
- **REFLECTION_QUESTIONS.md** - Hints and expandable answers
- **Tutorial files** in `tutorials/` folder
- Original papers cited in workshop materials

---

**Remember**: These questions are designed to develop **critical thinking** about your data. Biology is complex - the goal is to build intuition and ask better questions!

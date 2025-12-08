# Workshop Questions

Quick reference guide for discussion questions throughout the workshop.

---

## Module 1: Differential Expression (DESeq2)

### Q1: Why use negative binomial distribution instead of normal distribution?


---

### Q2: What does a positive log2FoldChange mean in Late vs Early comparison?


---

### Q3: Why use adjusted p-values (FDR) instead of raw p-values?


---

### Q4: Why do some genes have high fold-change but are NOT significant?


**Key insight**: Significance considers BOTH effect size (fold change) AND confidence (variance, sample size). A large fold change with high uncertainty is not statistically significant.

---

## Module 2: Gene Set Enrichment (GSEA)

### Q1: Why analyze pathways instead of individual genes?


---

### Q2: What does "enriched in up-regulated genes" mean?


---

### Q3: What is the Warburg effect?

---

### Q4: Why do different databases (GO vs KEGG) give different results?


---

## Module 3: Reporter Metabolites ⭐

### Q1: What is the difference between Reporter Metabolites and GSEA?

---

### Q2: Why aggregate Z-scores instead of counting significant genes?

---

### Q3: ⭐ Why might a metabolite appear in BOTH up AND down analyses?


---

### Q4: If a metabolite has high Z-score but only 3 neighboring genes, should you trust it?

---

### Q5: How would you validate a reporter metabolite experimentally?


---

## Module 4: Co-expression Networks

### Q1: What does positive correlation between two genes mean biologically?


---

### Q3: What is a "module" and why is it biologically meaningful?


---

### Q4: What does modularity measure? Is 0.3-0.5 good?


---

### Q5: Why analyze positive and negative correlations separately?


---

### Q6: Module with 100 genes but no pathway enrichment - what does this mean?


---

### Q7: Why are hub genes good drug targets?


---

### Q8: ⭐ How to integrate network modules with reporter metabolites?


---

## Integration Across All Modules

### Q: How do the 4 modules complement each other?



# Instructor Guide
## Transcriptomics to Reporter Metabolites Workshop

This guide is for instructors delivering the 3-hour workshop.

---

## 📅 Workshop Schedule

### 3-Hour Workshop (Current Format)

**Main Focus**: Reporter Metabolite Analysis

| Time | Module | Duration | Type |
|------|--------|----------|------|
| 0:00-0:30 | Installation & Setup | 30 min | Hands-on + Support |
| 0:30-0:50 | Module 1: Differential Expression (DESeq2) | 20 min | Hands-on |
| 0:50-1:05 | Module 2: Gene Set Enrichment (PIANO) | 15 min | Hands-on |
| 1:05-1:15 | **Break** | 10 min | - |
| 1:15-2:15 | ⭐ Module 3: Reporter Metabolites | 60 min | Hands-on (MAIN) |
| 2:15-2:45 | Module 4: Advanced Co-expression Networks | 30 min | Hands-on |
| 2:45-2:55 | Discussion & Interpretation | 10 min | Discussion |
| 2:55-3:00 | Wrap-up & Next Steps | 5 min | Lecture |

**BONUS**: Module 5 (MOFA2) available for self-study or fast learners

---

## 🎯 Learning Outcomes Assessment

### Module Completion Criteria

Participants should be able to:

**Module 1: Differential Expression (DESeq2)**:
- ✅ Load count data and create DESeq2 object
- ✅ Interpret log2FoldChange and padj values
- ✅ Generate volcano plots
- ✅ Identify differentially expressed genes

**Module 2: Gene Set Enrichment (PIANO)**:
- ✅ Map Ensembl IDs to gene symbols
- ✅ Run PIANO enrichment analysis
- ✅ Interpret pathway results
- ✅ Identify enriched biological processes

**Module 3: Reporter Metabolites ⭐ (MAIN)**:
- ✅ Understand reporter metabolite algorithm
- ✅ Load and parse metabolic model (SBML)
- ✅ Calculate metabolite-level Z-scores
- ✅ Perform background correction
- ✅ Identify key reporter metabolites
- ✅ Interpret metabolic reprogramming

**Module 4: Advanced Co-expression Networks**:
- ✅ Calculate gene co-expression correlations
- ✅ Apply FDR correction
- ✅ Detect network modules (Louvain algorithm)
- ✅ Perform functional enrichment on modules
- ✅ Export networks for Cytoscape

**BONUS - Module 5: MOFA2** (Optional):
- ✅ Prepare multi-omics data for MOFA2
- ✅ Train and evaluate MOFA2 models
- ✅ Interpret variance decomposition
- ✅ Identify factor-associated features

---

## 📚 Pre-Workshop Preparation

### 1 Week Before

- [ ] Send welcome email with installation instructions
- [ ] Share INSTALLATION.md
- [ ] Request participants to run `00_check_setup.R`
- [ ] Create Slack/Teams channel for questions
- [ ] Test all code on multiple platforms (Windows/Mac/Linux)

### 3 Days Before

- [ ] Send reminder email
- [ ] Share schedule and learning objectives
- [ ] Provide Zoom/Teams link (if virtual)
- [ ] Share pre-reading materials (optional)

### 1 Day Before

- [ ] Send final reminder
- [ ] Ensure all materials are on GitHub/shared drive
- [ ] Prepare presentation slides
- [ ] Test screen sharing and recording setup
- [ ] Prepare backup plans for common errors

### Day Of

- [ ] Arrive 30 minutes early (in-person) or test setup (virtual)
- [ ] Have backup internet connection ready
- [ ] Prepare printed materials (optional)
- [ ] Set up recording if desired

---

## 🎓 Teaching Tips

### General Approach

1. **Live Coding**:
   - Type code in real-time (don't copy-paste)
   - Explain each line as you type
   - Make intentional mistakes to show debugging
   - Use RStudio shortcuts (Ctrl+Enter, etc.)

2. **Interactive Questions**:
   - Pause every 10-15 minutes for questions
   - Use "raise hand" feature (virtual) or ask verbally
   - Cold call occasionally to check understanding
   - Encourage peer-to-peer explanations

3. **Pacing**:
   - Check completion status regularly
   - Have "stretch goals" for fast learners
   - Prepare abbreviated versions if running behind
   - Build in buffer time

4. **Engagement**:
   - Use polls/quizzes between modules
   - Encourage hypothesis generation before revealing results
   - Ask participants to predict outputs
   - Share real-world examples

---

## 🔧 Common Technical Issues

### Issue 1: Packages Won't Install

**Cause**: Outdated R version or missing system libraries

**Solution**:
```r
# Check R version
R.version.string

# Update Bioconductor
BiocManager::install(version = "3.18")

# Install system deps (Linux)
# Show in Terminal, not R:
sudo apt install libcurl4-openssl-dev libxml2-dev
```

**Prevention**: Emphasize minimum R version in pre-workshop email

---

### Issue 2: Data Files Not Found

**Cause**: Working directory not set correctly

**Solution**:
```r
# Show current directory
getwd()

# Set directory (adjust path)
setwd("path/to/workshop_materials")

# Check files are there
list.files()
```

**Prevention**: Include directory setup in Module 1 introduction

---

### Issue 3: MOFA2 Python Issues

**Cause**: Python not installed or not detected

**Solution**:
```r
# Use R-only mode
library(MOFA2)
# When running model:
model <- run_mofa(model, use_basilisk = TRUE)
```

**Prevention**: Mention Python is optional in installation guide

---

### Issue 4: Out of Memory Errors

**Cause**: Large dataset + limited RAM

**Solution**:
```r
# Use fewer features
top_var <- head(order(-apply(data, 1, var)), 1000)
data_subset <- data[top_var, ]

# Increase memory (Windows)
memory.limit(16000)
```

**Prevention**: Pre-filter datasets to manageable sizes

---

### Issue 5: Slow MOFA2 Training

**Cause**: Too many iterations or factors

**Solution**:
```r
# Reduce iterations for demo
training_opts <- get_default_training_options(model)
training_opts$maxiter <- 500  # Instead of 1000

# Reduce initial factors
model_opts <- get_default_model_options(model)
model_opts$num_factors <- 10  # Instead of 15
```

**Prevention**: Use pre-trained models for demonstration

---

## 📊 Key Discussion Points

### Module 1: Differential Expression

**Questions to Ask**:
- "Why do we use raw counts instead of FPKM for DESeq2?"
- "What does a negative log2FoldChange mean in our contrast?"
- "Why is padj more important than pvalue?"

**Common Confusion**:
- Direction of fold change (Late vs Early)
- Interpretation of NA values in results
- When to use lfcThreshold

---

### Module 5: MOFA2 Introduction

**Questions to Ask**:
- "How is MOFA2 different from running PCA on each omics separately?"
- "What would a factor active in mRNA but not protein suggest biologically?"
- "When would you NOT want to use MOFA2?"

**Common Confusion**:
- Difference between factors and principal components
- Why factors can be omics-specific
- How to determine optimal number of factors

---

## 🎨 Presentation Slides Outline

### Slide Deck 1: Introduction (Day 1)
1. Title slide with workshop info
2. Instructor introduction
3. Learning objectives
4. Schedule overview
5. TCGA dataset overview
6. Biological question
7. Analysis workflow diagram
8. Questions?

### Slide Deck 2: MOFA2 Introduction (Day 2)
1. Recap of Day 1
2. Multi-omics integration challenges
3. MOFA2 conceptual framework
4. Factor types (shared vs specific)
5. Example applications
6. MOFA2 workflow
7. Questions?

---

## ✏️ Exercise Solutions Strategy

### During Workshop

**Option 1**: Give 5-10 minutes for exercises, then solve together
- Pro: Everyone learns the solution
- Con: Slower learners may not attempt

**Option 2**: Provide solutions file, work independently
- Pro: Self-paced learning
- Con: Less interactive

**Option 3**: Break into groups, share solutions
- Pro: Peer learning
- Con: Requires more time

**Recommended**: Combination approach
- Quick exercises: Solve together
- Longer exercises: Independent + group discussion

---

## 📝 Assessment Options

### Formative Assessment (During Workshop)

1. **Quick Polls**:
   - "Which genes are up-regulated in Late stage? A) LFC > 0, B) LFC < 0"
   - "How many omics layers does MOFA2 require? A) 1, B) 2+, C) 3+"

2. **Code Prediction**:
   - "What will this code output?"
   - "Will this plot show up-regulated or down-regulated genes?"

3. **Debugging Exercises**:
   - Provide code with intentional errors
   - Ask participants to find and fix

### Summative Assessment (Post-Workshop)

1. **Take-Home Assignment**:
   - Analyze new dataset with same workflow
   - Submit R script + interpretation
   - Due 1 week after workshop

2. **Short Quiz**:
   - 10 multiple choice questions
   - Covers key concepts from all modules

3. **Project Proposal**:
   - Design multi-omics experiment
   - Describe analysis plan using workshop methods

---

## 🔄 Iterative Improvement

### After Each Workshop

1. **Collect Feedback**:
   - Send post-workshop survey
   - Ask about pace, clarity, usefulness
   - Request suggestions for improvement

2. **Review Common Issues**:
   - Note which concepts were confusing
   - Identify technical problems
   - Update materials accordingly

3. **Update Materials**:
   - Fix errors in tutorials
   - Add clarifications based on questions
   - Update package versions if needed

4. **Share Improvements**:
   - Commit changes to GitHub
   - Update version number
   - Notify previous participants of major updates

---

## 📧 Email Templates

### Pre-Workshop Welcome Email

```
Subject: Welcome to Integrative Multi-Omics Workshop - Setup Instructions

Dear Participants,

Welcome to the Integrative Multi-Omics Analysis Workshop! We're excited to have you join us.

**Workshop Details**:
- Date: [DATE]
- Time: [TIME]
- Location/Link: [LOCATION/ZOOM]

**Before the Workshop**:
1. Install R (≥4.0) and RStudio
2. Download workshop materials from: [GITHUB LINK]
3. Run setup script: source("scripts/00_install_packages.R")
4. Verify setup: source("scripts/00_check_setup.R")

Full instructions: [LINK TO INSTALLATION.MD]

**Need Help?**
- Join our Slack channel: [LINK]
- Email questions to: [EMAIL]

Looking forward to seeing you!

Best regards,
[Your Name]
```

### Post-Workshop Follow-Up

```
Subject: Workshop Materials & Next Steps

Dear Participants,

Thank you for attending the Integrative Multi-Omics Workshop!

**Workshop Materials**:
- Slides: [LINK]
- Recording: [LINK]
- Solutions: Available in solutions/ folder

**Next Steps**:
1. Complete practice exercises
2. Apply methods to your own data
3. Share your feedback: [SURVEY LINK]

**Additional Resources**:
- MOFA2 documentation: https://biofam.github.io/MOFA2/
- DESeq2 vignette: [LINK]
- Office hours: [TIME/LINK]

Please don't hesitate to reach out with questions!

Best regards,
[Your Name]
```

---

## 📖 Recommended Reading

### For Instructors

- Love et al. (2014) "Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2"
- Argelaguet et al. (2020) "MOFA+: a statistical framework for comprehensive integration"
- Ritchie et al. (2015) "limma powers differential expression analyses"

### Optional for Participants

- Conesa et al. (2016) "A survey of best practices for RNA-seq data analysis"
- Hasin et al. (2017) "Multi-omics approaches to disease"

---

## ⏱️ Time Management Tips

1. **Start on time** - don't wait for latecomers
2. **Keep breaks short but firm** (15 min maximum)
3. **Use timers** for exercises
4. **Have "bonus" content** that can be skipped if running late
5. **End on time** - respect participants' schedules

---

## 🎯 Success Metrics

Workshop is successful if:

- [ ] >80% of participants complete all modules
- [ ] >70% rate workshop as "Very Useful"
- [ ] >60% successfully apply methods to their own data
- [ ] <20% report major technical issues
- [ ] Positive written feedback

---

**Good luck with your workshop!**

Questions? Contact: [your.email@institution.edu]

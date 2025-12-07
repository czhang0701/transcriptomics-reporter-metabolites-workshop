# Quick Start Guide
## Get Started in 15 Minutes

This guide gets you up and running quickly. For detailed instructions, see [INSTALLATION.md](INSTALLATION.md).

---

## ⚡ Fast Track Setup

### Step 1: Install R and RStudio (5 min)

**Already have R ≥ 4.0 and RStudio?** Skip to Step 2.

1. Download R: https://cran.r-project.org/
2. Download RStudio: https://posit.co/download/rstudio-desktop/
3. Install both with default settings

---

### Step 2: Install Packages (5-10 min)

Open RStudio and run:

```r
# Install packages
source("https://raw.githubusercontent.com/YOUR_USERNAME/YOUR_REPO/main/scripts/00_install_packages.R")
```

Or manually:

```r
# BiocManager
install.packages("BiocManager")

# Core packages
BiocManager::install(c("DESeq2", "piano", "MOFA2"))
install.packages(c("tidyverse", "pheatmap", "corrplot", "Hmisc", "reshape2"))
```

**Wait for installation to complete** (5-10 minutes)

---

### Step 3: Verify Setup (2 min)

```r
# Check everything works
source("scripts/00_check_setup.R")
```

**Expected**: `✓✓✓ ALL REQUIRED PACKAGES INSTALLED AND WORKING! ✓✓✓`

**If errors**: See [INSTALLATION.md](INSTALLATION.md#troubleshooting)

---

### Step 4: Start Learning! (Now!)

```r
# Set working directory
setwd("path/to/workshop_materials")

# Open first tutorial
file.edit("tutorials/Module1_DifferentialExpression.md")
```

**OR** start with the R script:

```r
# Run Module 1 script (when created)
source("scripts/01_differential_expression.R")
```

---

## 📖 Tutorial Order

Follow in sequence:

1. **Module 1**: Differential Expression (45 min) ✅
2. **Module 2**: GSEA (45 min) ⏳
3. **Module 3**: Co-expression (45 min) ⏳
4. **Module 4**: Reporter Metabolites (30 min) ⏳
5. **Module 5**: MOFA2 Introduction (30 min) ✅
6. **Module 6**: MOFA2 Data Prep (45 min) ⏳
7. **Module 7**: MOFA2 Training (45 min) ⏳
8. **Module 8**: MOFA2 Analysis (45 min) ⏳

**Total time**: 5-6 hours (can be split across days)

---

## 🎯 Learning by Doing

### Option A: Guided Tutorials
- Read tutorials in `tutorials/` folder
- Type code manually
- Complete exercises
- Check solutions in `solutions/`

### Option B: Interactive Scripts
- Run scripts in `scripts/` folder
- Modify parameters
- Experiment with code
- Generate your own plots

### Option C: Hybrid
- Read tutorial for concepts
- Run script for implementation
- Try exercises for practice

---

## 🆘 Quick Troubleshooting

### Package won't install?
```r
# Update R
# Then retry:
BiocManager::install("PackageName", force = TRUE)
```

### Can't find data files?
```r
# Check directory
getwd()

# Set correct directory
setwd("path/to/workshop_materials")

# Verify files
list.files("data/")
```

### Code throws error?
1. Check you've run all previous code
2. Verify data loaded correctly
3. Check package loaded: `library(PackageName)`
4. See full error message
5. Search error online or ask instructor

---

## 📊 What You'll Create

By the end of this workshop:

✅ Differential expression results
✅ Enriched pathway tables
✅ Co-expression networks
✅ Reporter metabolite rankings
✅ MOFA2 trained models
✅ Publication-ready plots
✅ Comprehensive analysis reports

---

## 💡 Pro Tips

1. **Type code, don't copy-paste** - you'll learn better
2. **Read error messages carefully** - they often explain the problem
3. **Save your work frequently** - use R scripts
4. **Experiment!** - change parameters to see what happens
5. **Ask questions** - no question is too basic

---

## 🎓 Prerequisites

### Required Knowledge
- ✅ Basic R (variables, functions, data frames)
- ✅ Basic statistics (p-values, hypothesis testing)
- ✅ Basic genomics (genes, RNA-seq concepts)

### Helpful But Not Required
- Data visualization with ggplot2
- Bioconductor packages
- Multi-omics concepts

### Don't Know R?
**Start here first**:
- [R for Data Science](https://r4ds.had.co.nz/) (Chapters 1-5)
- [Intro to R](https://www.datacamp.com/courses/free-introduction-to-r) (DataCamp, free)

---

## 📁 Quick File Reference

| File/Folder | Description |
|-------------|-------------|
| `README.md` | Full workshop overview |
| `INSTALLATION.md` | Detailed setup guide |
| `tutorials/` | Step-by-step learning modules |
| `scripts/` | Runnable R scripts |
| `data/` | Workshop datasets |
| `solutions/` | Exercise answers |
| `figures/` | Generated plots |

---

## ⏱️ Time Estimates

**Quick learner** (programming background): 4 hours
**Average** (some R experience): 5-6 hours
**Thorough** (want to understand everything): 8-10 hours

**Tip**: Take breaks! Better to spread over 2-3 days.

---

## 🚀 After the Workshop

### Apply to Your Data
1. Prepare your data in same format
2. Modify scripts with your filenames
3. Adjust parameters for your experiment
4. Interpret results in your biological context

### Keep Learning
- Read cited papers
- Explore MOFA2 advanced features
- Try other Bioconductor packages
- Join online communities

### Share Your Work
- Publish analysis code on GitHub
- Write methods section using workshop approaches
- Cite tools appropriately
- Help others learn

---

## 📚 Resources

**Documentation**: See each file's README
**Help**: GitHub Issues or instructor email
**Community**: Bioconductor support forum
**Updates**: Check GitHub repo regularly

---

## ✅ Pre-Workshop Checklist

Before starting, make sure:

- [ ] R (≥4.0) installed
- [ ] RStudio installed
- [ ] Packages installed (`00_check_setup.R` passes)
- [ ] Data files downloaded
- [ ] Working directory set correctly
- [ ] Can create plots in R (`plot(1:10)` works)

---

**Ready?** Open `tutorials/Module1_DifferentialExpression.md` and let's begin!

**Questions?** See [README.md](README.md) for full details or [INSTALLATION.md](INSTALLATION.md) for setup help.

**Good luck and enjoy the workshop! 🎉**

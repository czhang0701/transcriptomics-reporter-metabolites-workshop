# Transcriptomics to Reporter Metabolites Workshop
## 3-Hour Hands-On Workshop

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-%3E%3D4.0-blue)](https://www.r-project.org/)

---

## 📋 Workshop Overview

Learn to analyze RNA-seq data and identify key metabolic reprogramming events through reporter metabolite analysis. This hands-on workshop takes you from raw counts to biological interpretation in just 3 hours!

**Duration**: 3 hours (including installation)
**Platform**: Windows & Mac (detailed instructions for both)
**Level**: Intermediate (basic R knowledge helpful)
**Main Goal**: Complete reporter metabolite analysis

---

## 🎯 What You'll Learn

### During the Workshop (All Students):
- ✅ Install R, RStudio, and packages (Windows & Mac)
- ✅ Run differential expression with DESeq2
- ✅ Perform pathway enrichment with PIANO
- ✅ **Identify reporter metabolites** ⭐ MAIN GOAL
- ✅ **Interpret metabolic reprogramming** ⭐
- ✅ Build advanced co-expression networks with module detection

### Bonus (Fast Learners or After Workshop):
- 📚 MOFA2 multi-omics integration (optional)
- 📚 Advanced tutorials with more details
- 📚 Practice dataset for additional learning

---

## ⏰ Workshop Timeline

```
┌─────────────────────────────────────────────────┐
│ 3-HOUR WORKSHOP                                  │
├─────────────────────────────────────────────────┤
│ 0:00-0:30 (30m) Installation (Windows & Mac)    │
│ 0:30-0:50 (20m) Module 1: Differential Expression│
│ 0:50-1:05 (15m) Module 2: Pathway Enrichment    │
│ ─────────────────────────────────────────────── │
│ 1:05-1:15 (10m) ☕ BREAK                        │
│ ─────────────────────────────────────────────── │
│ 1:15-2:15 (60m) ⭐ Module 3: Reporter Metabolites│
│ 2:15-2:45 (30m) Module 4: Advanced Co-expression│
│ 2:45-2:55 (10m) Discussion & Interpretation     │
│ ─────────────────────────────────────────────── │
│ 2:55-3:00 (5m)  Wrap-up & Next Steps           │
└─────────────────────────────────────────────────┘

BONUS (Optional - during workshop if fast, or after):
├─ MOFA2 tutorials (self-paced)
├─ Advanced topics
└─ Apply to own data
```

---

## 🚀 Getting Started

### **ONE COMPLETE GUIDE**

👉 **[WORKSHOP_GUIDE_COMPLETE.md](WORKSHOP_GUIDE_COMPLETE.md)** 👈

**This single document contains EVERYTHING**:
- ✅ Installation instructions (Windows & Mac)
- ✅ All 4 modules with code
- ✅ Step-by-step analysis
- ✅ Interpretation guidance
- ✅ Troubleshooting
- ✅ Bonus MOFA2 information

**Just follow this ONE guide during the workshop!**

---

## 📁 What's Included

### **Main Workshop Guide**
- `WORKSHOP_GUIDE_COMPLETE.md` - **START HERE!** (Complete 3-hour workshop)

### **Data**
- `data/` - TCGA RNA-seq dataset (50 samples, 20,000 genes)
- `data/README_data.md` - Data documentation

### **Full Tutorials** (Reference/Self-Study)
- `tutorials/Module1_DifferentialExpression.md` - Detailed DESeq2 guide
- `tutorials/Module2_GSEA.md` - Detailed GSEA guide
- `tutorials/Module3_ReporterMetabolites.md` - Detailed reporter metabolite guide ⭐
- `tutorials/Module4_Coexpression.md` - Advanced network analysis guide
- `tutorials/Module5_MOFA2_Introduction.md` - MOFA2 bonus (optional)

### **Scripts**
- `scripts/00_install_packages.R` - Auto-install script
- `scripts/00_check_setup.R` - Verify installation
- `scripts/create_practice_dataset.R` - Practice data generator

### **Solutions**
- `solutions/Module1_solutions.R` - Exercise answers

---

## 💻 System Requirements

### **Software**
- R (≥ 4.0.0)
- RStudio (latest version)

### **Platform Support**
- ✅ Windows 10/11
- ✅ Mac (Intel processors)
- ✅ Mac (Apple Silicon - M1/M2/M3)
- ✅ Linux (should work, less tested)

### **Disk Space**
- 2-3 GB free (for R + packages + data)

### **Internet**
- Required for package installation
- Needed during first 30 minutes

---

## 📚 Workshop Structure

### **Core Modules** (Required - Everyone Does These)

| Time | Module | Focus | Duration |
|------|--------|-------|----------|
| 0:00 | Setup | Install R/RStudio/Packages | 30 min |
| 0:30 | **1** | Differential Expression (DESeq2) | 20 min |
| 0:50 | **2** | Gene Set Enrichment (PIANO) | 15 min |
| 1:15 | **3** | ⭐ Reporter Metabolites (Patil & Nielsen 2005) | 60 min |
| 2:15 | **4** | Advanced Co-expression Networks + Module Detection | 30 min |

**Total**: 2 hours 20 minutes + 20 min breaks + 20 min discussion = **3 hours**

---

### **Bonus Module** (Optional - Fast Students or After Workshop)

| Module | Topic | When | Duration |
|--------|-------|------|----------|
| **5** | MOFA2 Multi-Omics | Self-paced | 30 min+ |

**MOFA2 is completely optional!**
- For students who finish early
- For after-workshop self-study
- Requires multi-omics data (not just RNA-seq)
- Full tutorial provided in `tutorials/Module5_MOFA2_Introduction.md`

---

## 🎓 Learning Outcomes

By the end of this workshop, you will:

### **Conceptual Understanding**
- Understand differential expression analysis
- Know when to use gene set enrichment
- Grasp co-expression network concepts
- **Understand reporter metabolite analysis** ⭐
- **Connect genes to metabolic networks** ⭐

### **Practical Skills**
- Install and configure R environment (both platforms)
- Run DESeq2 for differential expression
- Perform PIANO pathway enrichment
- Calculate gene correlations
- **Identify reporter metabolites** ⭐
- **Interpret metabolic reprogramming** ⭐

### **Deliverables**
- Differential expression results file
- Enriched pathway tables
- Co-expression network file
- **Reporter metabolite rankings** ⭐
- **Biological interpretation** ⭐

---

## 🔬 Scientific Context

### **Biological Question**
*What metabolic changes occur in cancer progression?*

### **Analysis Approach**
```
RNA-seq Data (50 samples: Early vs Late stage)
    ↓
Differential Expression (DESeq2)
    ↓ 2,000 significant genes
Pathway Enrichment (PIANO)
    ↓ 50 enriched pathways (glycolysis, TCA, etc.)
Co-expression Networks
    ↓ Gene modules
Reporter Metabolites ⭐
    ↓
Key Metabolic Reprogramming Identified!
(e.g., Warburg effect: lactate production)
```

---

## 💡 What Makes This Workshop Unique

1. ⭐ **Only workshop** focused on reporter metabolites
2. ⏰ **Realistic 3-hour** format (including installation!)
3. 🖥️ **True cross-platform** support (Windows & Mac detailed instructions)
4. 📄 **Single comprehensive guide** (no jumping between documents)
5. 🎯 **Clear goal** (complete reporter metabolite analysis)
6. 📚 **Flexible learning** (fast students get MOFA2 bonus)
7. 🔬 **Real biological insights** (Warburg effect, cancer metabolism)

---

## 🤝 For Instructors

### **Teaching Materials**
- `WORKSHOP_GUIDE_COMPLETE.md` - Complete workshop script
- `INSTRUCTOR_GUIDE.md` - Teaching tips and strategies
- Email templates and schedules in documentation

### **Pre-Workshop**
- No student setup required!
- Send workshop location/time
- Ensure good internet connection

### **During Workshop**
- Follow WORKSHOP_GUIDE_COMPLETE.md step-by-step
- First 30 min: Help with installation (both platforms)
- Modules 1-3: Quick context for Module 4
- Module 4: Full depth - this is the goal
- Bonus: Direct fast students to MOFA2 materials

---

## 📖 Additional Resources

### **Reporter Metabolites**
- [Patil & Nielsen (2005) PNAS](https://www.pnas.org/doi/10.1073/pnas.0406811102) - Original method
- [PIANO Package](https://bioconductor.org/packages/release/bioc/html/piano.html)

### **RNA-seq Analysis**
- [DESeq2 Vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
- [RNA-seq Best Practices](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4728800/)

### **Metabolic Models**
- [Human-GEM](https://github.com/SysBioChalmers/Human-GEM)
- [VMH Database](http://vmh.life)

### **MOFA2 (Bonus)**
- [MOFA2 Website](https://biofam.github.io/MOFA2/)
- [MOFA2 Paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02015-1)

---

## 📧 Contact

**Questions?** Contact workshop instructor: cheng.5.zhang@kcl.ac.uk

**Issues?** Please report any errors or suggestions

---

## 📄 License

MIT License - Free to use and modify with attribution

Data from TCGA - Subject to TCGA data access policies

---

## 🙏 Acknowledgments

- TCGA Research Network for cancer genomics data
- Patil & Nielsen for reporter metabolite methodology
- PIANO team (Väremo et al.) for excellent R package
- MOFA2 team (Argelaguet et al.) for multi-omics tools
- Bioconductor community

---

## ✅ Quick Start Summary

**For Students**:
1. Come to workshop (bring laptop)
2. Open `WORKSHOP_GUIDE_COMPLETE.md`
3. Follow along for 3 hours
4. Learn reporter metabolite analysis!
5. (Optional) Continue with MOFA2 bonus materials

**For Instructors**:
1. Review `WORKSHOP_GUIDE_COMPLETE.md`
2. Test on your platform
3. Prepare to help with installation (30 min)
4. Follow guide during workshop
5. Direct fast students to MOFA2 materials

---

## 🎯 Workshop Goal

**By the end of 3 hours, ALL students will:**
- ✅ Have R and RStudio working
- ✅ Complete differential expression analysis
- ✅ Understand pathway enrichment
- ✅ Build gene networks
- ✅ **Identify reporter metabolites** ⭐
- ✅ **Interpret metabolic reprogramming in cancer** ⭐

**Fast students will additionally:**
- ✅ Start exploring MOFA2 multi-omics integration
- ✅ Work through advanced tutorials
- ✅ Apply to practice dataset

---

**Ready to uncover metabolic reprogramming? Let's get started!** 🧬→🔬

👉 **[WORKSHOP_GUIDE_COMPLETE.md](WORKSHOP_GUIDE_COMPLETE.md)** 👈

**Everything you need in ONE place!**

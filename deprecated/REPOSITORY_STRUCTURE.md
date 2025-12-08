# Repository Structure
## Transcriptomics to Reporter Metabolites Workshop

**Version**: 3.0.0 Final
**Status**: ✅ Production Ready
**Last Updated**: December 2025

---

## 📁 Clean Repository Structure

```
workshop_materials/
│
├── 📄 README.md                       ⭐ START HERE - Repository overview
├── 📄 WORKSHOP_GUIDE_COMPLETE.md      ⭐ MAIN GUIDE - Complete 3-hour workshop
├── 📄 INSTRUCTOR_GUIDE.md             👨‍🏫 Teaching tips and strategies
├── 📄 LICENSE                         📜 MIT License
├── 📄 .gitignore                      🔧 Git ignore file
├── 📄 REPOSITORY_STRUCTURE.md         📋 This file
│
├── 📂 tutorials/                      📚 Detailed reference materials
│   ├── Module1_DifferentialExpression.md   (12 KB)
│   ├── Module2_GSEA.md                      (19 KB)
│   ├── Module3_Coexpression.md              (20 KB)
│   ├── Module4_ReporterMetabolites.md       (19 KB) ⭐ Main focus
│   └── Module5_MOFA2_Introduction.md        (11 KB) 💎 Bonus
│
├── 📂 scripts/                        💻 R scripts
│   ├── 00_install_packages.R          Auto-install all packages
│   ├── 00_check_setup.R               Verify installation
│   └── create_practice_dataset.R      Generate practice data
│
├── 📂 solutions/                      ✅ Exercise answers
│   └── Module1_solutions.R
│
├── 📂 data/                           📊 RNA-seq dataset
│   ├── README_data.md                 Data documentation
│   └── [TCGA data files...]           Real cancer RNA-seq data
│
├── 📂 figures/                        📈 Output plots (generated during workshop)
│
└── 📂 deprecated/                     ⚠️ Archived old files (DO NOT USE)
    ├── README.md                      Explains why files are deprecated
    └── [old documentation files...]
```

---

## 🎯 Primary User Journey

### **For Students**

1. **Start**: Open `README.md`
2. **Main Workshop**: Follow `WORKSHOP_GUIDE_COMPLETE.md` for 3 hours
3. **Deep Dive**: Explore `tutorials/` folder for detailed explanations
4. **Practice**: Use scripts to generate practice datasets
5. **Check**: Review `solutions/` for exercise answers

### **For Instructors**

1. **Overview**: Read `README.md`
2. **Preparation**: Review `INSTRUCTOR_GUIDE.md`
3. **Teaching**: Follow `WORKSHOP_GUIDE_COMPLETE.md` step-by-step
4. **Support**: Use `scripts/00_check_setup.R` to help students verify installation

---

## 📄 File Descriptions

### **Root Level Files**

| File | Purpose | Size | Essential? |
|------|---------|------|------------|
| `README.md` | Repository landing page, overview | 12 KB | ✅ YES |
| `WORKSHOP_GUIDE_COMPLETE.md` | Complete 3-hour workshop guide | 26 KB | ✅ YES |
| `INSTRUCTOR_GUIDE.md` | Teaching tips, strategies | 12 KB | ✅ YES |
| `LICENSE` | MIT License with attributions | 2 KB | ✅ YES |
| `.gitignore` | Git ignore rules | <1 KB | ✅ YES |
| `REPOSITORY_STRUCTURE.md` | This file | 8 KB | 📚 Reference |

### **Tutorials Folder** (📂 tutorials/)

| File | Content | Size | Workshop Time |
|------|---------|------|---------------|
| `Module1_DifferentialExpression.md` | DESeq2 detailed guide | 12 KB | 20 min |
| `Module2_GSEA.md` | PIANO pathway enrichment | 19 KB | 15 min |
| `Module3_Coexpression.md` | Network analysis | 20 KB | 15 min |
| `Module4_ReporterMetabolites.md` | Reporter metabolites ⭐ | 19 KB | 60 min |
| `Module5_MOFA2_Introduction.md` | MOFA2 bonus 💎 | 11 KB | Optional |

**Total**: 81 KB of comprehensive tutorial content

### **Scripts Folder** (📂 scripts/)

| File | Purpose | When to Use |
|------|---------|-------------|
| `00_install_packages.R` | Auto-install all R packages | Start of workshop |
| `00_check_setup.R` | Verify installation success | After installation |
| `create_practice_dataset.R` | Generate small practice dataset | For extra practice |

### **Solutions Folder** (📂 solutions/)

| File | Content | When to Share |
|------|---------|---------------|
| `Module1_solutions.R` | Exercise answers for Module 1 | After students attempt exercises |

### **Data Folder** (📂 data/)

- `README_data.md` - Documentation of TCGA dataset
- TCGA RNA-seq data files (50 samples, ~20,000 genes)

### **Deprecated Folder** (📂 deprecated/)

Contains **outdated** documentation files from development:
- ❌ DO NOT use for workshop
- 📚 Kept for historical reference only
- See `deprecated/README.md` for details

---

## 🚀 Quick Start Guide

### **Option 1: GitHub Web Interface**

1. Visit repository on GitHub
2. Click **"Code"** → **"Download ZIP"**
3. Extract ZIP file
4. Open `README.md` to get started

### **Option 2: Git Clone**

```bash
git clone [repository-url]
cd workshop_materials
```

### **Then Follow**

1. Open `WORKSHOP_GUIDE_COMPLETE.md`
2. Follow installation instructions (first 30 minutes)
3. Complete all 4 modules (2 hours 30 minutes)
4. Finish with discussion and wrap-up

---

## 📊 Workshop Timeline

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
│ 1:15-1:30 (15m) Module 3: Co-expression         │
│ 1:30-2:30 (60m) ⭐ Module 4: Reporter Metabolites│
│ 2:30-2:50 (20m) Discussion & Interpretation     │
│ ─────────────────────────────────────────────── │
│ 2:50-3:00 (10m) Wrap-up & Next Steps           │
└─────────────────────────────────────────────────┘

BONUS (Optional - during workshop if fast, or after):
├─ MOFA2 tutorials (self-paced)
├─ Advanced topics
└─ Apply to own data
```

---

## 🎓 Learning Outcomes

**By the end of 3 hours, ALL students will:**
- ✅ Have R and RStudio installed and working
- ✅ Complete differential expression analysis (DESeq2)
- ✅ Perform pathway enrichment (PIANO)
- ✅ Build co-expression networks
- ✅ **Identify reporter metabolites** ⭐
- ✅ **Interpret metabolic reprogramming** ⭐

**Fast students will additionally:**
- ✅ Explore MOFA2 multi-omics integration
- ✅ Work through detailed tutorials
- ✅ Apply to practice dataset

---

## 💻 System Requirements

### **Software**
- R (≥ 4.0.0)
- RStudio (latest version)

### **Platforms Supported**
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

## 📝 What Students Need to Bring

- [ ] Laptop (fully charged or with charger)
- [ ] Workshop materials downloaded
- [ ] Internet connection available
- [ ] Enthusiasm to learn! 🎉

**Note**: Students do NOT need to install R/RStudio before workshop. Installation happens during first 30 minutes with instructor help.

---

## 👨‍🏫 What Instructors Need

### **Before Workshop**
- [ ] Review `WORKSHOP_GUIDE_COMPLETE.md`
- [ ] Read `INSTRUCTOR_GUIDE.md`
- [ ] Test all code on your platform
- [ ] Prepare backup pre-run results (optional)
- [ ] Ensure stable internet connection

### **During Workshop**
- [ ] Projector/screen sharing working
- [ ] USB backup of materials
- [ ] `scripts/00_check_setup.R` for troubleshooting
- [ ] Timer for keeping pace

---

## 🔧 Customization

### **To Adapt for Your Institution**

1. **Update contact info** in:
   - `README.md` (line 291)
   - `WORKSHOP_GUIDE_COMPLETE.md`
   - `INSTRUCTOR_GUIDE.md`

2. **Adjust timeline** if needed:
   - All timing in `WORKSHOP_GUIDE_COMPLETE.md`
   - Can extend to 3.5-4 hours if desired

3. **Add your data**:
   - Replace data in `data/` folder
   - Update `data/README_data.md`
   - Adjust gene IDs in Module 1 code

---

## 📈 Repository Statistics

**Content Size**:
- Documentation: ~170 KB
- Scripts: 3 files
- Tutorials: 5 modules (81 KB)
- Solutions: 1 file
- Total: Clean, focused, production-ready

**Lines of Code**:
- R scripts: ~500 lines
- Tutorial code: ~1000 lines
- Total teaching material: Comprehensive

---

## 🌟 What Makes This Repository Special

1. ⭐ **Only workshop** focused on reporter metabolites
2. ⏰ **Realistic 3-hour** format (including installation)
3. 🖥️ **True cross-platform** (Windows & Mac detailed)
4. 📄 **Single comprehensive guide** (no document jumping)
5. 🎯 **Clear goal** (complete Module 4 guaranteed)
6. 📚 **Flexible learning** (bonus MOFA2 content)
7. 🔬 **Real insights** (Warburg effect, cancer metabolism)
8. 🧹 **Clean structure** (deprecated files separated)

---

## 📧 Support and Feedback

### **For Students**
- Questions during workshop: Ask instructor
- Issues after workshop: Open GitHub Issue
- General questions: Email instructor

### **For Instructors**
- Teaching tips: See `INSTRUCTOR_GUIDE.md`
- Technical issues: Check `scripts/00_check_setup.R`
- Contributions: Submit Pull Request

---

## 📖 Citation

If you use this workshop material, please cite:

**Workshop Materials**:
```
Transcriptomics to Reporter Metabolites Workshop (2025)
GitHub: [repository-url]
```

**Key Methods**:
- Patil, K. R., & Nielsen, J. (2005). Reporter metabolites. PNAS, 102(8), 2685-2689.
- Väremo, L., et al. (2013). PIANO. Nucleic Acids Research, 41(8), 4378-4391.

---

## 📜 License

MIT License - Free to use and modify with attribution

See `LICENSE` file for full details.

---

## 🙏 Acknowledgments

- TCGA Research Network for cancer genomics data
- Patil & Nielsen for reporter metabolite methodology
- PIANO team (Väremo et al.) for excellent R package
- MOFA2 team (Argelaguet et al.) for multi-omics tools
- Bioconductor community

---

## ✅ Repository Checklist

**For GitHub Upload**:
- [x] Clean structure (deprecated files moved)
- [x] README.md complete
- [x] WORKSHOP_GUIDE_COMPLETE.md finalized
- [x] All tutorials complete
- [x] Scripts functional
- [x] Solutions provided
- [x] LICENSE included
- [x] .gitignore configured

**Ready to Upload**: ✅ YES

---

**Version**: 3.0.0 Final
**Status**: Production Ready for GitHub
**Folder Structure**: Clean and Organized
**User Path**: Clear and Simple

**START HERE**: README.md → WORKSHOP_GUIDE_COMPLETE.md → Success! 🎉

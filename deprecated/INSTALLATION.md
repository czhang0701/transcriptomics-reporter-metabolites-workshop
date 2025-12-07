# Installation Guide
## Integrative Multi-Omics Analysis Workshop

This guide will help you set up your environment for the workshop.

---

## 📋 System Requirements

### Minimum Requirements
- **RAM**: 8 GB (16 GB recommended)
- **Disk Space**: 5 GB free space
- **Internet**: Required for package installation

### Operating Systems
- ✅ Windows 10/11
- ✅ macOS 10.14+
- ✅ Ubuntu 18.04+ / Debian / Fedora

---

## 🔧 Software Installation

### Step 1: Install R (≥ 4.0.0)

#### Windows
1. Download R from [CRAN](https://cran.r-project.org/bin/windows/base/)
2. Run the installer
3. Accept default settings
4. Verify installation: Open R and check version

#### macOS
1. Download R from [CRAN](https://cran.r-project.org/bin/macosx/)
2. Download the `.pkg` file for your macOS version
3. Run the installer
4. Verify installation in Terminal: `R --version`

#### Linux (Ubuntu/Debian)
```bash
sudo apt update
sudo apt install r-base r-base-dev
R --version
```

---

### Step 2: Install RStudio (Recommended)

1. Download [RStudio Desktop](https://posit.co/download/rstudio-desktop/)
2. Install for your operating system
3. Launch RStudio
4. Verify R is detected (bottom left panel should show R version)

---

### Step 3: Install R Packages

#### Option A: Automated Installation (Recommended)

1. Open RStudio
2. Set working directory to the workshop materials folder:
```r
setwd("path/to/workshop_materials")  # Adjust path
```

3. Run the installation script:
```r
source("scripts/00_install_packages.R")
```

4. Wait for installation to complete (10-30 minutes depending on your system)

#### Option B: Manual Installation

If automated installation fails, install packages manually:

```r
# Install BiocManager
install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install(c("DESeq2", "piano", "MOFA2"))

# Install CRAN packages
install.packages(c(
  "tidyverse", "ggplot2", "pheatmap", "corrplot",
  "Hmisc", "reshape2", "data.table"
))
```

---

### Step 4: Verify Installation

Run the verification script:

```r
source("scripts/00_check_setup.R")
```

**Expected output:**
```
✓ R version 4.x.x
✓ DESeq2 (version x.x.x)
✓ piano (version x.x.x)
✓ MOFA2 (version x.x.x)
...
✓✓✓ ALL REQUIRED PACKAGES INSTALLED AND WORKING! ✓✓✓
```

If you see errors, see the **Troubleshooting** section below.

---

## 📁 Download Workshop Materials

### Option 1: Clone with Git
```bash
git clone https://github.com/yourusername/integrative-omics-workshop.git
cd integrative-omics-workshop
```

### Option 2: Download ZIP
1. Go to the GitHub repository
2. Click "Code" → "Download ZIP"
3. Extract the ZIP file
4. Navigate to the extracted folder

---

## 🗂️ Directory Structure

After setup, your directory should look like:

```
workshop_materials/
├── README.md                    # Workshop overview
├── INSTALLATION.md             # This file
├── data/                       # Workshop datasets
│   ├── Counts_selected.txt
│   ├── FPKM_selected.txt
│   ├── patients.txt
│   └── ...
├── scripts/                    # R scripts
│   ├── 00_install_packages.R
│   ├── 00_check_setup.R
│   └── ...
├── tutorials/                  # Step-by-step guides
│   ├── Module1_DifferentialExpression.md
│   ├── Module2_GSEA.md
│   └── ...
├── solutions/                  # Exercise solutions
└── figures/                    # Generated plots
```

---

## 🐍 Optional: Python Setup for MOFA2

MOFA2 can run in R-only mode, but Python backend provides better performance.

### Install Python (≥ 3.8)

#### Windows
1. Download Python from [python.org](https://www.python.org/downloads/)
2. **Important**: Check "Add Python to PATH" during installation

#### macOS
```bash
brew install python3
```

#### Linux
```bash
sudo apt install python3 python3-pip
```

### Install MOFA2 Python Package

```bash
pip install mofapy2
```

### Configure R to Use Python

In R:
```r
install.packages("reticulate")
library(reticulate)

# Check Python detection
py_config()

# If needed, specify Python path
use_python("/path/to/python3")
```

---

## 🔍 Troubleshooting

### Issue 1: Package Installation Fails

**Error**: `installation of package 'X' had non-zero exit status`

**Solution**:
1. Update R to the latest version
2. Update BiocManager: `BiocManager::install(version = "3.18")`
3. Install system dependencies (Linux):
   ```bash
   sudo apt install libcurl4-openssl-dev libssl-dev libxml2-dev
   ```
4. Try installing packages one by one

---

### Issue 2: DESeq2 Installation Fails

**Solution** (macOS):
```bash
# Install Xcode Command Line Tools
xcode-select --install
```

**Solution** (Linux):
```bash
sudo apt install build-essential gfortran
```

---

### Issue 3: MOFA2 Installation Fails

**Error**: `Python dependencies not found`

**Solution 1** - Use R-only mode:
```r
# MOFA2 can run without Python
library(MOFA2)
# Set use_basilisk = TRUE in run_mofa()
```

**Solution 2** - Install Python dependencies:
```bash
pip install mofapy2 numpy pandas scipy h5py
```

---

### Issue 4: Cannot Load Packages

**Error**: `there is no package called 'X'`

**Solution**:
```r
# Check installed packages
installed.packages()[, "Package"]

# Reinstall specific package
BiocManager::install("PackageName", force = TRUE)
```

---

### Issue 5: Memory Errors

**Error**: `cannot allocate vector of size X Gb`

**Solution**:
1. Close other applications
2. Increase R memory limit (Windows):
   ```r
   memory.limit(size = 16000)  # 16 GB
   ```
3. Use fewer features in analysis
4. Consider using a machine with more RAM

---

### Issue 6: Data Files Not Found

**Error**: `cannot open file 'Counts_selected.txt'`

**Solution**:
```r
# Check current directory
getwd()

# Set correct directory
setwd("path/to/workshop_materials/data")

# Or use full path
Counts <- read.table("path/to/Counts_selected.txt", sep="\t")
```

---

## 💻 Platform-Specific Notes

### Windows
- Use RStudio instead of R GUI for better experience
- File paths use backslashes: `C:\Users\...`
- In R, use forward slashes: `"C:/Users/..."`

### macOS
- May need to install Xcode Command Line Tools
- Some packages require compilation (takes longer)
- Use Terminal for command-line operations

### Linux
- Best performance for bioinformatics tools
- May need to install development libraries
- Check distribution-specific package names

---

## ✅ Pre-Workshop Checklist

Before the workshop starts, ensure:

- [ ] R (≥ 4.0) installed
- [ ] RStudio installed
- [ ] All required packages installed (run `00_check_setup.R`)
- [ ] Workshop materials downloaded
- [ ] Data files accessible
- [ ] Can create plots (test: `plot(1:10)`)
- [ ] Internet connection available (for package updates)

---

## 📧 Getting Help

### During Installation

1. **Check error messages carefully** - they often indicate the solution
2. **Google the error** - many installation issues have documented solutions
3. **Check package documentation** - `?packageName` or visit Bioconductor
4. **Ask for help**:
   - Email instructor: your.email@institution.edu
   - Open GitHub issue with:
     - Your OS and R version
     - Error message
     - Output of `sessionInfo()`

### During the Workshop

- Raise hand for immediate help
- Post in workshop Slack/Teams channel
- Check `solutions/` folder for exercise answers

---

## 🎓 Additional Resources

### R and RStudio
- [RStudio Cheat Sheets](https://posit.co/resources/cheatsheets/)
- [R for Data Science](https://r4ds.had.co.nz/)

### Bioconductor
- [Bioconductor Installation Guide](https://bioconductor.org/install/)
- [BiocManager Documentation](https://cran.r-project.org/web/packages/BiocManager/)

### Package Documentation
- [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
- [MOFA2](https://bioconductor.org/packages/release/bioc/html/MOFA2.html)
- [piano](https://bioconductor.org/packages/release/bioc/html/piano.html)

---

## 🔄 Updates

**Last Updated**: December 2025

Check the GitHub repository for:
- Package version updates
- Bug fixes
- Additional materials

---

**Ready to proceed?** Return to the main [README.md](README.md) to start the workshop!

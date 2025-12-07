# Pre-Workshop Setup Guide
## **COMPLETE THIS BEFORE THE WORKSHOP** ⚠️

**Time Required**: 30-45 minutes
**Must Complete**: At least 24 hours before workshop

---

## ⏰ Why Complete This Early?

Package installation can take 20-40 minutes depending on your internet speed and computer. **DO NOT wait until the workshop starts!**

The workshop is only 3 hours, and we need that time for analysis, not installation.

---

## 📋 Pre-Workshop Checklist

Complete these steps IN ORDER:

- [ ] **Step 1**: Install R (15 min)
- [ ] **Step 2**: Install RStudio (5 min)
- [ ] **Step 3**: Install R packages (20-30 min)
- [ ] **Step 4**: Verify installation (5 min)
- [ ] **Step 5**: Download workshop materials (5 min)

**Total time**: 45-60 minutes

---

## 💻 Step 1: Install R (≥ 4.0.0)

### Windows Users

1. **Download R**:
   - Go to: https://cran.r-project.org/bin/windows/base/
   - Click **"Download R-4.x.x for Windows"** (latest version)
   - Save the `.exe` file

2. **Install R**:
   - Double-click the downloaded `.exe` file
   - Click **"Next"** through all prompts (accept defaults)
   - Click **"Finish"**

3. **Verify Installation**:
   - Press `Windows + R` key
   - Type `cmd` and press Enter
   - In the black window, type: `R --version`
   - Should show: `R version 4.x.x`

**Screenshot**: You should see something like `R version 4.3.2 (2023-10-31)`

---

### Mac Users

1. **Download R**:
   - Go to: https://cran.r-project.org/bin/macosx/
   - Download the `.pkg` file for your macOS version:
     - **macOS 11+** (Big Sur or newer): Download latest R-4.x.x-arm64.pkg OR R-4.x.x-x86_64.pkg
     - Check your Mac: Apple menu → About This Mac
     - **Apple Silicon (M1/M2/M3)**: Use arm64 version
     - **Intel processor**: Use x86_64 version

2. **Install R**:
   - Double-click the downloaded `.pkg` file
   - Click **"Continue"** through the installer
   - Enter your Mac password when prompted
   - Click **"Close"** when done

3. **Verify Installation**:
   - Open **Terminal** (Applications → Utilities → Terminal)
   - Type: `R --version`
   - Should show: `R version 4.x.x`

**Common Mac Issue**: If you get "command not found":
```bash
# Add R to your path (in Terminal):
echo 'export PATH="/Library/Frameworks/R.framework/Resources/bin:$PATH"' >> ~/.zshrc
source ~/.zshrc
```

---

## 🖥️ Step 2: Install RStudio

### Both Windows and Mac

1. **Download RStudio**:
   - Go to: https://posit.co/download/rstudio-desktop/
   - Scroll down to **"All Installers"**
   - **Windows**: Download `RStudio-202x.xx.x-xxx.exe`
   - **Mac**: Download `RStudio-202x.xx.x-xxx.dmg`

2. **Install RStudio**:

   **Windows**:
   - Double-click the `.exe` file
   - Click **"Next"** through all prompts
   - Click **"Finish"**

   **Mac**:
   - Double-click the `.dmg` file
   - Drag **RStudio** icon to **Applications** folder
   - Eject the disk image

3. **Launch RStudio**:
   - **Windows**: Start Menu → RStudio
   - **Mac**: Applications folder → RStudio
   - Or double-click RStudio icon

4. **Verify R is Detected**:
   - Look at bottom-left panel in RStudio
   - Should say: `R version 4.x.x`
   - **If you see an error**: R is not installed correctly, go back to Step 1

**Screenshot**: RStudio has 4 panels. Bottom-left should show R version.

---

## 📦 Step 3: Install R Packages

**THIS IS THE MOST IMPORTANT STEP!** ⚠️

### Option A: Automatic Installation (Recommended)

1. **Download Installation Script**:
   - Download workshop materials first (see Step 5)
   - OR copy the script from below

2. **In RStudio**:
   ```r
   # Set working directory to workshop folder
   # Windows example:
   setwd("C:/Users/YourName/Downloads/workshop_materials")

   # Mac example:
   setwd("/Users/YourName/Downloads/workshop_materials")

   # Run installation script
   source("scripts/00_install_packages.R")
   ```

3. **Wait for Installation**:
   - **THIS WILL TAKE 20-30 MINUTES**
   - You'll see lots of text scrolling
   - Red text is OK - not all are errors!
   - **DO NOT CLOSE RStudio** during installation

4. **When Prompted**:
   - "Update all/some/none?" → Type `a` (all) and press Enter
   - "Do you want to install from sources?" → Type `n` (no) and press Enter

---

### Option B: Manual Installation (If Script Fails)

Copy and paste this into RStudio console:

```r
# Install BiocManager (package manager)
install.packages("BiocManager")

# Install Bioconductor packages (10-15 minutes)
BiocManager::install(c("DESeq2", "piano"), update = FALSE, ask = FALSE)

# Install CRAN packages (10-15 minutes)
install.packages(c(
  "corrplot", "Hmisc", "reshape2",
  "ggplot2", "pheatmap", "tidyverse"
), dependencies = TRUE)
```

**Press Enter and wait...**

---

### Platform-Specific Issues

#### Windows Issues

**Issue 1**: "WARNING: Rtools is required"
- **Solution**:
  - Download Rtools: https://cran.r-project.org/bin/windows/Rtools/
  - Install Rtools43 (or latest version)
  - Restart RStudio
  - Try package installation again

**Issue 2**: "Permission denied"
- **Solution**: Run RStudio as Administrator
  - Right-click RStudio icon → "Run as administrator"

**Issue 3**: Antivirus blocking installation
- **Solution**: Temporarily disable antivirus during installation
  - Re-enable after packages are installed

---

#### Mac Issues

**Issue 1**: "clang: error: unsupported option '-fopenmp'"
- **Solution**: Install Xcode Command Line Tools
  ```bash
  # In Terminal:
  xcode-select --install
  ```
  - Click "Install" in popup window
  - Wait 10-15 minutes for installation
  - Try R package installation again

**Issue 2**: "gfortran is required but missing"
- **Solution**: Install gfortran
  - **Intel Mac**:
    - Download from: https://mac.r-project.org/tools/
    - Install gfortran-12.2-universal.pkg
  - **Apple Silicon (M1/M2)**:
    - Download: gfortran-ARM-12.2-Monterey.dmg
  - Restart RStudio

**Issue 3**: "Cannot install source packages"
- **Solution**: Install binary versions only
  ```r
  # In RStudio:
  options(pkgType = "binary")
  # Then retry installation
  ```

---

## ✅ Step 4: Verify Installation

**Copy and paste this into RStudio console**:

```r
# Create verification script
cat('
cat("\\n=== Checking Installation ===\\n")

# Check R version
r_version <- R.version.string
cat("R version:", r_version, "\\n")

if (getRversion() < "4.0") {
  stop("ERROR: R version too old. Need >= 4.0")
}

# Check required packages
required_packages <- c("DESeq2", "piano", "corrplot",
                      "Hmisc", "reshape2", "ggplot2")

cat("\\nChecking packages:\\n")
all_installed <- TRUE

for (pkg in required_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat("  ✓", pkg, "\\n")
  } else {
    cat("  ✗", pkg, "- NOT INSTALLED\\n")
    all_installed <- FALSE
  }
}

if (all_installed) {
  cat("\\n✓✓✓ ALL PACKAGES INSTALLED SUCCESSFULLY! ✓✓✓\\n")
  cat("You are ready for the workshop!\\n")
} else {
  cat("\\n⚠️  SOME PACKAGES MISSING  ⚠️\\n")
  cat("Please install missing packages before workshop.\\n")
}
', file = "check_setup.R")

# Run verification
source("check_setup.R")
```

**Expected Output**:
```
=== Checking Installation ===
R version: R version 4.3.2 (2023-10-31)

Checking packages:
  ✓ DESeq2
  ✓ piano
  ✓ corrplot
  ✓ Hmisc
  ✓ reshape2
  ✓ ggplot2

✓✓✓ ALL PACKAGES INSTALLED SUCCESSFULLY! ✓✓✓
You are ready for the workshop!
```

**If you see any ✗ symbols**:
- Try installing that specific package manually
- See troubleshooting section below

---

## 📥 Step 5: Download Workshop Materials

### Option 1: Download ZIP (Easiest)

1. **Go to**: [GitHub repository URL will be provided]
2. Click green **"Code"** button
3. Click **"Download ZIP"**
4. **Extract** the ZIP file:
   - **Windows**: Right-click → "Extract All"
   - **Mac**: Double-click the ZIP file
5. **Remember the location!** You'll need this path in the workshop

**Recommended locations**:
- Windows: `C:\Users\YourName\Documents\workshop_materials`
- Mac: `/Users/YourName/Documents/workshop_materials`

---

### Option 2: Git Clone (If You Have Git)

```bash
# In Terminal (Mac) or Command Prompt (Windows):
cd Documents
git clone [repository-url]
cd workshop_materials
```

---

## 🔧 Troubleshooting

### General Issues

**Problem**: "Package X failed to install"

**Solutions** (try in order):
1. **Restart RStudio** and try again
2. **Install packages one by one**:
   ```r
   BiocManager::install("DESeq2")
   BiocManager::install("piano")
   install.packages("ggplot2")
   # etc.
   ```
3. **Check internet connection**
4. **Try different CRAN mirror**:
   ```r
   options(repos = c(CRAN = "https://cloud.r-project.org"))
   ```

---

**Problem**: "Cannot load package even though installed"

**Solution**:
```r
# Remove and reinstall
remove.packages("package-name")
install.packages("package-name")
```

---

**Problem**: "R session aborted" during installation

**Solutions**:
- **Increase memory** (if possible)
- **Install fewer packages at once**
- **Close other applications**

---

### Platform-Specific Help

#### Windows Additional Help

**Check Windows Version**:
- Compatible with Windows 10 or 11
- Windows 7/8: May need older R version

**Free Disk Space**:
- Need at least **2 GB free** for packages

**User Account**:
- Must have **admin rights** or ability to install software

---

#### Mac Additional Help

**Check macOS Version**:
- Minimum: macOS 10.13 (High Sierra)
- Recommended: macOS 11+ (Big Sur or newer)

**Security Settings**:
- If "app can't be opened": System Preferences → Security & Privacy → Allow

**Rosetta 2** (for Apple Silicon Macs):
- Some packages may require Rosetta 2
- Install if prompted: `softwareupdate --install-rosetta`

---

## 📧 Getting Help

### Before the Workshop

**If installation fails**:

1. **Take a screenshot** of the error message
2. **Email instructor** with:
   - Your operating system (Windows/Mac + version)
   - R version (`R.version.string`)
   - Error message screenshot
   - Which step failed

**Contact**: [your.email@institution.edu]

**Response time**: Within 24 hours

---

### During Setup

**Common Questions**:

**Q**: "How do I know if installation is complete?"
**A**: When you see `>` prompt in RStudio console and no more text is scrolling

**Q**: "Installation taking forever - is this normal?"
**A**: Yes! Can take 30-40 minutes. Don't close RStudio.

**Q**: "Red text in console - is this an error?"
**A**: Not always! Only if it says "ERROR" or "installation failed"

**Q**: "Can I use R without RStudio?"
**A**: Yes, but RStudio is HIGHLY recommended for the workshop

---

## ✅ Pre-Workshop Checklist (Final)

**24 Hours Before Workshop**, confirm:

- [ ] ✓ R installed (version ≥ 4.0)
- [ ] ✓ RStudio installed and launches
- [ ] ✓ All packages installed (ran verification script)
- [ ] ✓ Workshop materials downloaded
- [ ] ✓ Know the path to workshop materials folder
- [ ] ✓ Computer charged (or bring charger!)
- [ ] ✓ Stable internet connection available

**If all checked**: You're ready! 🎉

**If any unchecked**: Contact instructor ASAP

---

## 🎯 What to Bring to Workshop

- [ ] **Laptop** (fully charged or with charger)
- [ ] **R and RStudio** (pre-installed and verified)
- [ ] **Workshop materials** (downloaded and extracted)
- [ ] **Notebook** for notes (optional)
- [ ] **Questions** from setup process

---

## 📅 Workshop Day

**Arrive 10 minutes early** to:
- Test internet connection
- Open RStudio
- Load workshop materials
- Ask last-minute setup questions

**Workshop will start promptly** - no time for installation!

---

## 💾 Backup Plan

**If your installation fails completely**:

1. **Email instructor immediately** (at least 48h before)
2. **Options**:
   - Remote desktop to working setup
   - Pair with another student
   - Use cloud-based R (Posit Cloud)

**Do NOT wait until workshop day to report problems!**

---

## 📚 Additional Resources

### Video Tutorials

**Installing R and RStudio**:
- Windows: https://www.youtube.com/watch?v=YrEe2TLr3MI
- Mac: https://www.youtube.com/watch?v=5r2QZ_swM68

### Documentation

- R Installation: https://cran.r-project.org/
- RStudio Guide: https://posit.co/download/rstudio-desktop/
- Bioconductor: https://bioconductor.org/install/

---

## ⏰ Timeline Summary

**2-3 Days Before Workshop**:
- Complete Steps 1-5 (60 minutes)
- Verify all packages installed
- Email instructor if issues

**1 Day Before**:
- Re-run verification script
- Confirm workshop materials location
- Prepare questions

**Workshop Day**:
- Arrive early
- Open RStudio
- Ready to learn!

---

## ❓ FAQ

**Q: I use Linux. Can I attend?**
A: Yes! Installation is usually easier on Linux. Follow general R installation for your distribution.

**Q: Can I use an iPad/Chromebook?**
A: No - need Windows/Mac/Linux computer with R installed locally.

**Q: What if I have an older R version?**
A: Must upgrade to ≥ 4.0. Older versions won't work with required packages.

**Q: I'm not comfortable with R. Will I be lost?**
A: Basic R knowledge helpful, but we'll explain code step-by-step.

**Q: Can I install packages during the workshop?**
A: NO - takes too long. Must complete before workshop.

---

**Setup complete?** See you at the workshop! 🚀

**Questions?** Email: [instructor email]

**Next**: On workshop day, open `WORKSHOP_GUIDE.md` for the 3-hour schedule.

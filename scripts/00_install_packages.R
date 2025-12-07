################################################################################
# Workshop Setup: Package Installation Script
# Integrative Multi-Omics Analysis Workshop
################################################################################

cat("========================================\n")
cat("Installing Required R Packages\n")
cat("========================================\n\n")

# Check R version
r_version <- as.numeric(paste0(R.version$major, ".", R.version$minor))
if (r_version < 4.0) {
  stop("This workshop requires R version 4.0 or higher. Please upgrade R.")
} else {
  cat("✓ R version", R.version$major, ".", R.version$minor, "detected\n\n")
}

################################################################################
# 1. Install BiocManager (if not already installed)
################################################################################

cat("Step 1: Installing BiocManager...\n")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  cat("✓ BiocManager installed\n\n")
} else {
  cat("✓ BiocManager already installed\n\n")
}

################################################################################
# 2. Install Bioconductor Packages
################################################################################

cat("Step 2: Installing Bioconductor packages...\n")
cat("This may take several minutes...\n\n")

bioc_packages <- c(
  "DESeq2",           # Differential expression analysis
  "piano",            # Gene set enrichment analysis
  "MOFA2"             # Multi-omics factor analysis
)

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
    cat("✓", pkg, "installed\n\n")
  } else {
    cat("✓", pkg, "already installed\n")
  }
}

################################################################################
# 3. Install CRAN Packages
################################################################################

cat("\nStep 3: Installing CRAN packages...\n")

cran_packages <- c(
  # Data manipulation
  "tidyverse",        # Data wrangling and ggplot2
  "data.table",       # Fast data manipulation
  "reshape2",         # Data reshaping

  # Visualization
  "ggplot2",          # Grammar of graphics plotting
  "pheatmap",         # Pretty heatmaps
  "corrplot",         # Correlation plots
  "ggrepel",          # Better text labels in ggplot
  "RColorBrewer",     # Color palettes
  "viridis",          # Colorblind-friendly palettes

  # Statistical analysis
  "Hmisc",            # Correlation with p-values

  # Utilities
  "knitr",            # Report generation
  "rmarkdown"         # R Markdown documents
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg, dependencies = TRUE, repos = "https://cloud.r-project.org/")
    cat("✓", pkg, "installed\n\n")
  } else {
    cat("✓", pkg, "already installed\n")
  }
}

################################################################################
# 4. Optional: Install reticulate for Python integration
################################################################################

cat("\nStep 4: Optional Python integration...\n")
cat("(You can skip this if you don't plan to use Python)\n")

install_python <- readline(prompt = "Install reticulate for Python integration? (y/n): ")

if (tolower(install_python) == "y") {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    install.packages("reticulate")
    cat("✓ reticulate installed\n")
    cat("Note: You'll need Python 3.8+ installed separately\n")
    cat("MOFA2 Python package: pip install mofapy2\n\n")
  } else {
    cat("✓ reticulate already installed\n\n")
  }
}

################################################################################
# 5. Installation Summary
################################################################################

cat("\n========================================\n")
cat("Installation Complete!\n")
cat("========================================\n\n")

cat("Installed packages:\n")
cat("Bioconductor:", length(bioc_packages), "packages\n")
cat("CRAN:", length(cran_packages), "packages\n\n")

cat("Next step: Run 00_check_setup.R to verify installation\n")

################################################################################
# 6. Session Info (for troubleshooting)
################################################################################

cat("\nSession Information (for troubleshooting):\n")
cat("==========================================\n")
print(sessionInfo())

cat("\n\nIf you encounter any errors, please:\n")
cat("1. Check your internet connection\n")
cat("2. Restart R/RStudio\n")
cat("3. Try installing failed packages individually\n")
cat("4. Check package documentation for dependencies\n")

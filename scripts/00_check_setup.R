################################################################################
# Workshop Setup: Package Verification Script
# Integrative Multi-Omics Analysis Workshop
################################################################################

cat("========================================\n")
cat("Checking Workshop Setup\n")
cat("========================================\n\n")

################################################################################
# Function to check if package is installed and loadable
################################################################################

check_package <- function(pkg_name) {
  tryCatch({
    suppressPackageStartupMessages(library(pkg_name, character.only = TRUE))
    version <- as.character(packageVersion(pkg_name))
    cat("✓", pkg_name, "(version", version, ")\n")
    return(TRUE)
  }, error = function(e) {
    cat("✗", pkg_name, "- NOT INSTALLED or cannot be loaded\n")
    return(FALSE)
  })
}

################################################################################
# 1. Check R Version
################################################################################

cat("Step 1: Checking R version...\n")
r_version <- as.numeric(paste0(R.version$major, ".", R.version$minor))
if (r_version >= 4.0) {
  cat("✓ R version", R.version.string, "\n\n")
} else {
  cat("✗ R version too old. Please upgrade to R >= 4.0\n\n")
}

################################################################################
# 2. Check Required Packages
################################################################################

cat("Step 2: Checking required packages...\n")
cat("--------------------------------------\n")

required_packages <- c(
  # Bioconductor packages
  "DESeq2",
  "piano",
  "MOFA2",

  # CRAN packages
  "tidyverse",
  "ggplot2",
  "pheatmap",
  "corrplot",
  "Hmisc",
  "reshape2",
  "data.table"
)

results <- sapply(required_packages, check_package)

################################################################################
# 3. Check Optional Packages
################################################################################

cat("\nStep 3: Checking optional packages...\n")
cat("--------------------------------------\n")

optional_packages <- c(
  "ggrepel",
  "RColorBrewer",
  "viridis",
  "knitr",
  "rmarkdown",
  "reticulate"
)

opt_results <- sapply(optional_packages, check_package)

################################################################################
# 4. Test Key Functions
################################################################################

cat("\nStep 4: Testing key functions...\n")
cat("--------------------------------------\n")

# Test DESeq2
test_deseq <- tryCatch({
  suppressPackageStartupMessages(library(DESeq2))
  # Create minimal test data
  counts <- matrix(rpois(100, lambda = 10), ncol = 10)
  coldata <- data.frame(condition = factor(rep(c("A", "B"), each = 5)))
  dds <- DESeqDataSetFromMatrix(counts, coldata, ~ condition)
  cat("✓ DESeq2 functions working\n")
  TRUE
}, error = function(e) {
  cat("✗ DESeq2 test failed:", e$message, "\n")
  FALSE
})

# Test MOFA2
test_mofa <- tryCatch({
  suppressPackageStartupMessages(library(MOFA2))
  # Just check if main functions are available
  if (exists("create_mofa") && exists("run_mofa")) {
    cat("✓ MOFA2 functions available\n")
    TRUE
  } else {
    cat("✗ MOFA2 functions not found\n")
    FALSE
  }
}, error = function(e) {
  cat("✗ MOFA2 test failed:", e$message, "\n")
  FALSE
})

# Test plotting
test_plot <- tryCatch({
  suppressPackageStartupMessages(library(ggplot2))
  # Create simple plot
  p <- ggplot(data.frame(x = 1:10, y = 1:10), aes(x, y)) + geom_point()
  cat("✓ ggplot2 plotting working\n")
  TRUE
}, error = function(e) {
  cat("✗ Plotting test failed:", e$message, "\n")
  FALSE
})

################################################################################
# 5. Check Working Directory and Data Files
################################################################################

cat("\nStep 5: Checking data files...\n")
cat("--------------------------------------\n")

cat("Current working directory:\n")
cat(getwd(), "\n\n")

# Check for data files (adjust path as needed)
data_files <- c(
  "Counts_selected.txt",
  "FPKM_selected.txt",
  "patients.txt",
  "Ensembl2gene.tsv",
  "c5.bp.v6.2.symbols.gmt",
  "Genes_selected.txt"
)

cat("Looking for data files in current directory...\n")
files_found <- 0
for (file in data_files) {
  if (file.exists(file)) {
    cat("✓", file, "\n")
    files_found <- files_found + 1
  } else {
    cat("✗", file, "- NOT FOUND\n")
  }
}

if (files_found == 0) {
  cat("\nNote: No data files found in current directory.\n")
  cat("Make sure to set your working directory to the Materials folder.\n")
  cat("Use: setwd('path/to/Materials')\n")
}

################################################################################
# 6. System Information
################################################################################

cat("\nStep 6: System information...\n")
cat("--------------------------------------\n")
cat("Platform:", R.version$platform, "\n")
cat("OS:", Sys.info()["sysname"], Sys.info()["release"], "\n")
cat("R version:", R.version.string, "\n")

# Check memory
if (.Platform$OS.type == "unix") {
  mem_info <- system("free -h", intern = TRUE)
  cat("\nMemory:\n", mem_info[2], "\n")
} else {
  cat("\nMemory check not available on this OS\n")
}

################################################################################
# 7. Final Summary
################################################################################

cat("\n========================================\n")
cat("Setup Summary\n")
cat("========================================\n\n")

all_required <- all(results)
all_optional <- all(opt_results)
all_tests <- test_deseq && test_mofa && test_plot

if (all_required && all_tests) {
  cat("✓✓✓ ALL REQUIRED PACKAGES INSTALLED AND WORKING! ✓✓✓\n\n")
  cat("You are ready to start the workshop!\n")
  cat("Proceed to tutorials/Module1_DifferentialExpression.md\n\n")
} else {
  cat("⚠️  SOME ISSUES DETECTED  ⚠️\n\n")

  if (!all_required) {
    cat("Missing required packages:\n")
    missing <- names(results)[!results]
    for (pkg in missing) {
      cat("  -", pkg, "\n")
    }
    cat("\nPlease install missing packages:\n")
    cat("  BiocManager::install(c('", paste(missing, collapse = "', '"), "'))\n\n")
  }

  if (!all_tests) {
    cat("Some function tests failed. Try:\n")
    cat("  1. Restart R/RStudio\n")
    cat("  2. Reload packages\n")
    cat("  3. Check error messages above\n\n")
  }
}

if (!all_optional) {
  cat("Optional packages not installed (workshop will still work):\n")
  missing_opt <- names(opt_results)[!opt_results]
  for (pkg in missing_opt) {
    cat("  -", pkg, "\n")
  }
  cat("\n")
}

if (files_found < length(data_files)) {
  cat("⚠️  Data files not found in current directory\n")
  cat("Make sure to:\n")
  cat("  1. Set working directory: setwd('path/to/Materials')\n")
  cat("  2. Verify data files are in the correct location\n\n")
}

################################################################################
# 8. Save Session Info
################################################################################

cat("Saving session info to 'session_info.txt'...\n")
sink("session_info.txt")
cat("Workshop Setup Check\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
print(sessionInfo())
sink()
cat("✓ Session info saved\n\n")

cat("Setup check complete!\n")
cat("If you need help, share 'session_info.txt' with the instructor.\n")

################################################################################
# Reporter Metabolites Algorithm (Patil & Nielsen, 2005)
#
# The Reporter Metabolites algorithm identifies metabolites around which
# transcriptional regulation is significantly enriched.
#
# Reference: Patil KR, Nielsen J. Uncovering transcriptional regulation of
# metabolism by using metabolic network topology. Proc Natl Acad Sci USA
# 2005;102:2685-2689.
#
# This is a pure R implementation based on the RAVEN MATLAB code.
################################################################################

#' Reporter Metabolites Algorithm
#'
#' @param model List containing metabolic model with:
#'   - mets: metabolite IDs
#'   - metNames: metabolite names (optional)
#'   - S: stoichiometric matrix (metabolites x reactions)
#'   - rxnGeneMat: reaction-gene association matrix (reactions x genes)
#'   - genes: gene IDs
#' @param genes Character vector of gene names
#' @param genePValues Numeric vector of P-values for differential expression
#' @param geneFoldChanges Numeric vector of log fold changes (optional)
#' @param printResults Logical, print top 20 to console (default FALSE)
#'
#' @return List with Reporter Metabolites results
reporterMetabolites <- function(model, genes, genePValues, geneFoldChanges = NULL, printResults = FALSE) {

  # Input validation
  if (length(genes) != length(genePValues)) {
    stop("The number of genes and the number of P-values must be the same")
  }

  if (!all(c("mets", "S", "rxnGeneMat", "genes") %in% names(model))) {
    stop("Model must have: mets, S, rxnGeneMat, genes")
  }

  # Remove NaN P-values
  if (!is.null(geneFoldChanges)) {
    geneFoldChanges <- geneFoldChanges[!is.na(genePValues)]
  }
  genes <- genes[!is.na(genePValues)]
  genePValues <- genePValues[!is.na(genePValues)]

  # Remove genes not in the model
  in_model <- genes %in% model$genes
  genes <- genes[in_model]
  genePValues <- genePValues[in_model]
  if (!is.null(geneFoldChanges)) {
    geneFoldChanges <- geneFoldChanges[in_model]
  }

  cat("Genes in model:", length(genes), "\n")

  # Convert P-values to Z-scores
  geneZScores <- qnorm(genePValues, lower.tail = FALSE)

  # Prevent errors if P-values are exactly 0 or 1
  geneZScores[geneZScores == -Inf] <- -15
  geneZScores[geneZScores == Inf] <- 15

  cat("Gene Z-scores calculated\n")

  # For each metabolite, calculate aggregate Z-score
  n_mets <- length(model$mets)
  metZScores <- rep(NA, n_mets)
  metNGenes <- rep(NA, n_mets)
  meanZ <- rep(NA, n_mets)
  stdZ <- rep(NA, n_mets)

  cat("Calculating metabolite Z-scores for", n_mets, "metabolites...\n")

  for (i in 1:n_mets) {
    # Get reactions involving this metabolite
    rxn_indices <- which(model$S[i, ] != 0)

    if (length(rxn_indices) > 0) {
      # Get genes involved in these reactions
      gene_indices <- which(colSums(model$rxnGeneMat[rxn_indices, , drop = FALSE]) > 0)

      # Find these genes in our gene list
      gene_matches <- which(genes %in% model$genes[gene_indices])

      if (length(gene_matches) > 0) {
        # Calculate aggregated Z-score for the metabolite
        metZScores[i] <- sum(geneZScores[gene_matches]) / sqrt(length(gene_matches))
        meanZ[i] <- mean(geneZScores[gene_matches])
        stdZ[i] <- sd(geneZScores[gene_matches])
        metNGenes[i] <- length(gene_matches)
      }
    }

    # Progress indicator
    if (i %% 100 == 0) {
      cat("Processed", i, "/", n_mets, "metabolites\n")
    }
  }

  # Remove metabolites with no Z-scores
  valid <- !is.na(metZScores)
  mets <- model$mets[valid]
  if ("metNames" %in% names(model)) {
    metNames <- model$metNames[valid]
  } else {
    metNames <- mets
  }
  metNGenes <- metNGenes[valid]
  meanZ <- meanZ[valid]
  stdZ <- stdZ[valid]
  metZScores <- metZScores[valid]

  cat("Valid metabolites with scores:", length(mets), "\n")

  # Correct for background by calculating mean Z-score for random sets
  cat("Performing background correction...\n")
  sizes <- unique(metNGenes)
  nGenes <- length(genes)

  for (size in sizes) {
    # Sample 100,000 random sets of this size
    n_samples <- 100000
    random_indices <- matrix(sample(1:nGenes, n_samples * size, replace = TRUE),
                             ncol = size)
    random_zscores <- matrix(geneZScores[random_indices], ncol = size)
    random_aggregated <- rowSums(random_zscores) / sqrt(size)

    # Calculate background mean and std
    bg_mean <- mean(random_aggregated, na.rm = TRUE)
    bg_std <- sd(random_aggregated, na.rm = TRUE)

    # Correct all metabolites of this size
    size_indices <- which(metNGenes == size)
    metZScores[size_indices] <- (metZScores[size_indices] - bg_mean) / bg_std

    cat("  Size", size, ": background mean =", round(bg_mean, 3),
        ", std =", round(bg_std, 3), "\n")
  }

  # Calculate P-values from corrected Z-scores
  metPValues <- pnorm(metZScores, lower.tail = FALSE)

  # Sort by Z-score (descending)
  sort_order <- order(metZScores, decreasing = TRUE)
  mets <- mets[sort_order]
  metNames <- metNames[sort_order]
  metZScores <- metZScores[sort_order]
  metPValues <- metPValues[sort_order]
  metNGenes <- metNGenes[sort_order]
  meanZ <- meanZ[sort_order]
  stdZ <- stdZ[sort_order]

  # Prepare output
  repMets <- list(
    test = "all",
    mets = mets,
    metNames = metNames,
    metZScores = metZScores,
    metPValues = metPValues,
    metNGenes = metNGenes,
    meanZ = meanZ,
    stdZ = stdZ
  )

  # Call recursively for up/down regulated genes if fold changes provided
  results <- list(repMets)

  if (!is.null(geneFoldChanges) && any(!is.na(geneFoldChanges))) {
    cat("\nCalculating reporter metabolites for UP-regulated genes...\n")
    up_genes <- geneFoldChanges >= 0
    if (sum(up_genes) > 0) {
      repMets_up <- reporterMetabolites(
        model,
        genes[up_genes],
        genePValues[up_genes],
        printResults = FALSE
      )
      repMets_up[[1]]$test <- "only up"
      results <- c(results, repMets_up)
    }

    cat("\nCalculating reporter metabolites for DOWN-regulated genes...\n")
    down_genes <- geneFoldChanges < 0
    if (sum(down_genes) > 0) {
      repMets_down <- reporterMetabolites(
        model,
        genes[down_genes],
        genePValues[down_genes],
        printResults = FALSE
      )
      repMets_down[[1]]$test <- "only down"
      results <- c(results, repMets_down)
    }
  }

  # Print top 20 if requested
  if (printResults) {
    for (result in results) {
      cat("\n===========================================\n")
      cat("TOP 20 REPORTER METABOLITES\n")
      cat("TEST TYPE:", result$test, "\n")
      cat("===========================================\n")
      cat(sprintf("%-15s %-40s %s\n", "ID", "NAME", "P-VALUE"))
      cat("-------------------------------------------\n")

      n_print <- min(20, length(result$mets))
      for (j in 1:n_print) {
        cat(sprintf("%-15s %-40s %.2e\n",
                    result$mets[j],
                    substr(result$metNames[j], 1, 40),
                    result$metPValues[j]))
      }
      cat("\n")
    }
  }

  return(results)
}


#' Save Reporter Metabolites Results to File
#'
#' @param repMets List of reporter metabolites results
#' @param outputFile Path to output file
saveReporterMetabolites <- function(repMets, outputFile) {
  sink(outputFile)

  for (result in repMets) {
    cat("REPORTER METABOLITES USING TEST TYPE:", result$test, "\n")
    cat("ID\tNAME\tZ-SCORE\tP-VALUE\tNUMBER OF NEIGHBOURS\tAVERAGE Z-SCORE\tSTD Z-SCORE\n")

    for (i in 1:length(result$mets)) {
      cat(sprintf("%s\t%s\t%.4f\t%.2e\t%d\t%.4f\t%.4f\n",
                  result$mets[i],
                  result$metNames[i],
                  result$metZScores[i],
                  result$metPValues[i],
                  result$metNGenes[i],
                  result$meanZ[i],
                  result$stdZ[i]))
    }
    cat("\n")
  }

  sink()
  cat("Results saved to:", outputFile, "\n")
}

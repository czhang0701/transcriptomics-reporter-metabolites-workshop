################################################################################
# Advanced Co-expression Network Analysis with Module Detection
#
# This implements:
# 1. Spearman correlation with FDR correction
# 2. Network filtering (top X% correlations)
# 3. Community detection (Louvain/Leiden algorithm)
# 4. Module enrichment analysis
# 5. Cytoscape export
#
# Based on igraph and community detection algorithms
################################################################################

#' Calculate Co-expression Network
#'
#' @param expr_matrix Expression matrix (genes x samples)
#' @param pval_threshold FDR threshold for significance (default 0.05)
#' @param pos_only Keep only positive correlations (default FALSE)
#' @return Data frame with gene pairs and correlations
calculateCoexpression <- function(expr_matrix, pval_threshold = 0.05, pos_only = FALSE) {

  library(Hmisc)

  cat("Calculating Spearman correlation...\n")

  # Calculate correlation and p-values
  corr_result <- rcorr(t(expr_matrix), type = "spearman")

  corr <- corr_result$r
  pval <- corr_result$P

  cat("Creating edge list...\n")

  # Keep only upper triangle (avoid duplicates)
  corr[lower.tri(corr, diag = TRUE)] <- NA
  pval[lower.tri(pval, diag = TRUE)] <- NA

  # Convert to long format
  corr_long <- reshape2::melt(corr, na.rm = TRUE)
  pval_long <- reshape2::melt(pval, na.rm = TRUE)

  # Combine
  result <- data.frame(
    Gene1 = corr_long$Var1,
    Gene2 = corr_long$Var2,
    Correlation = corr_long$value,
    Pvalue = pval_long$value,
    stringsAsFactors = FALSE
  )

  cat("Adjusting P-values with FDR...\n")
  result$FDR <- p.adjust(result$Pvalue, method = "fdr")

  # Filter by significance
  result <- result[result$FDR < pval_threshold, ]

  # Filter by direction
  if (pos_only) {
    result <- result[result$Correlation > 0, ]
  }

  result <- result[, c("Gene1", "Gene2", "Correlation", "Pvalue", "FDR")]

  cat("Done! Found", nrow(result), "significant correlations\n")

  return(result)
}


#' Detect Network Modules
#'
#' @param edge_list Data frame with Gene1, Gene2, Correlation
#' @param top_percent Top percentage of correlations to keep (default 0.1 = 10%)
#' @param min_module_size Minimum module size (default 30)
#' @param positive Whether to analyze positive (TRUE) or negative (FALSE) correlations
#' @return List with network and clustering results
detectModules <- function(edge_list, top_percent = 0.1, min_module_size = 30, positive = TRUE) {

  library(igraph)

  cat("Filtering network to top", top_percent * 100, "% of correlations...\n")

  # Filter by direction
  if (positive) {
    edges <- edge_list[edge_list$Correlation > 0, ]
    threshold <- quantile(edges$Correlation, 1 - top_percent)
    edges <- edges[edges$Correlation > threshold, ]
  } else {
    edges <- edge_list[edge_list$Correlation < 0, ]
    threshold <- quantile(edges$Correlation, top_percent)
    edges <- edges[edges$Correlation < threshold, ]
  }

  cat("Building network with", nrow(edges), "edges...\n")

  # Create igraph object
  g <- graph_from_data_frame(
    d = edges[, c("Gene1", "Gene2", "Correlation")],
    directed = FALSE
  )

  cat("Detecting modules using Louvain algorithm...\n")

  # Community detection (Louvain algorithm - similar to Leiden)
  # Using absolute correlation as edge weight
  E(g)$weight <- abs(E(g)$Correlation)
  communities <- cluster_louvain(g, weights = E(g)$weight)

  # Get module assignments
  modules <- data.frame(
    Gene = V(g)$name,
    Module = communities$membership,
    stringsAsFactors = FALSE
  )

  # Filter small modules
  module_sizes <- table(modules$Module)
  large_modules <- names(module_sizes)[module_sizes >= min_module_size]
  modules <- modules[modules$Module %in% large_modules, ]

  # Renumber modules
  module_mapping <- setNames(1:length(large_modules), large_modules)
  modules$Module <- module_mapping[as.character(modules$Module)]

  cat("Found", length(unique(modules$Module)), "modules (size >=", min_module_size, ")\n")
  cat("Module sizes:", paste(table(modules$Module), collapse = ", "), "\n")

  # Calculate modularity
  modularity_score <- modularity(g, communities$membership)

  cat("Modularity:", round(modularity_score, 3), "\n")

  return(list(
    network = g,
    modules = modules,
    modularity = modularity_score,
    edges = edges
  ))
}


#' Enrich Network Modules
#'
#' @param modules Data frame with Gene and Module columns
#' @param gsc Gene set collection (from piano loadGSC)
#' @param background Total number of genes
#' @return List of enrichment results per module
enrichModules <- function(modules, gsc, background = NULL) {

  library(piano)

  if (is.null(background)) {
    background <- length(unique(modules$Gene))
  }

  cat("Enriching", length(unique(modules$Module)), "modules...\n")

  enrichment_results <- list()

  for (mod_id in sort(unique(modules$Module))) {
    cat("  Module", mod_id, "...")

    # Get genes in this module
    mod_genes <- modules$Gene[modules$Module == mod_id]

    # Create gene stats vector (1 for genes in module, 0 for others)
    # For over-representation analysis
    all_genes <- unique(modules$Gene)
    gene_stats <- rep(0, length(all_genes))
    names(gene_stats) <- all_genes
    gene_stats[mod_genes] <- 1

    # Run enrichment (Fisher's exact test for over-representation)
    gsaRes <- runGSA(
      geneLevelStats = gene_stats,
      gsc = gsc,
      geneSetStat = "fisher",
      signifMethod = "geneSampling",
      nPerm = 1000,
      gsSizeLim = c(10, 500),
      adjMethod = "fdr",
      ncpus = 1
    )

    # Get significant results
    results <- GSAsummaryTable(gsaRes)

    if (nrow(results) > 0) {
      # Add FDR column
      results$FDR <- results$`p adj (dist.dir.up)`

      # Keep only significant
      results <- results[results$FDR < 0.05, ]

      if (nrow(results) > 0) {
        cat(" Found", nrow(results), "enriched pathways\n")
        enrichment_results[[as.character(mod_id)]] <- results
      } else {
        cat(" No significant enrichment\n")
      }
    } else {
      cat(" No significant enrichment\n")
    }
  }

  return(enrichment_results)
}


#' Export Network for Cytoscape
#'
#' @param modules_result Output from detectModules()
#' @param output_prefix File prefix for output files
exportCytoscape <- function(modules_result, output_prefix) {

  cat("Exporting network for Cytoscape...\n")

  # Node file
  node_file <- paste0(output_prefix, "_nodes.txt")
  node_data <- data.frame(
    Gene = modules_result$modules$Gene,
    Module = modules_result$modules$Module,
    stringsAsFactors = FALSE
  )

  # Add module sizes
  module_sizes <- table(node_data$Module)
  node_data$ModuleSize <- module_sizes[as.character(node_data$Module)]

  write.table(
    node_data,
    file = node_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  cat("Nodes saved to:", node_file, "\n")

  # Edge file
  edge_file <- paste0(output_prefix, "_edges.txt")
  edge_data <- modules_result$edges[, c("Gene1", "Gene2", "Correlation")]
  names(edge_data) <- c("Source", "Target", "Weight")

  write.table(
    edge_data,
    file = edge_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  cat("Edges saved to:", edge_file, "\n")

  cat("\nTo visualize in Cytoscape:\n")
  cat("1. Import network: File > Import > Network from File\n")
  cat("2. Select", edge_file, "\n")
  cat("3. Import node attributes: File > Import > Table from File\n")
  cat("4. Select", node_file, "\n")
  cat("5. Style by Module column\n")
}


#' Save Module Enrichment Results
#'
#' @param enrichment_results List of enrichment results
#' @param output_file Excel file to save results
saveModuleEnrichment <- function(enrichment_results, output_file) {

  library(openxlsx)

  if (length(enrichment_results) == 0) {
    cat("No enrichment results to save\n")
    return()
  }

  cat("Saving enrichment results to:", output_file, "\n")

  wb <- createWorkbook()

  for (mod_id in names(enrichment_results)) {
    sheet_name <- paste0("Module_", mod_id)
    addWorksheet(wb, sheet_name)

    results <- enrichment_results[[mod_id]]

    # Select key columns
    results_clean <- data.frame(
      Pathway = rownames(results),
      Pathway_Name = results$Name,
      FDR = results$FDR,
      Genes_in_pathway = results$`Genes (tot)`,
      stringsAsFactors = FALSE
    )

    writeData(wb, sheet_name, results_clean)
  }

  saveWorkbook(wb, output_file, overwrite = TRUE)

  cat("Enrichment results saved!\n")
}

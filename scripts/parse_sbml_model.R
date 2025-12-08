################################################################################
# Parse SBML Metabolic Model
#
# Extracts metabolite-reaction-gene associations from Reference_model.xml
# Creates a model structure compatible with Reporter Metabolites algorithm
################################################################################

library(XML)

#' Parse SBML Model for Reporter Metabolites
#'
#' @param sbml_file Path to SBML file
#' @return List with model structure (mets, S, rxnGeneMat, genes)
parseSBMLModel <- function(sbml_file) {

  cat("Parsing SBML file:", sbml_file, "\n")
  cat("This may take 2-5 minutes...\n\n")

  # Parse XML
  doc <- xmlParse(sbml_file)
  root <- xmlRoot(doc)

  # Get namespace
  ns <- c(sbml = "http://www.sbml.org/sbml/level2/version3")

  # Get model node
  model_node <- xmlChildren(root)[["model"]]

  # Extract species (metabolites)
  cat("Step 1: Extracting metabolites...\n")
  species_list <- model_node[["listOfSpecies"]]
  n_species <- xmlSize(species_list)

  mets <- character(n_species)
  metNames <- character(n_species)

  for (i in 1:n_species) {
    species <- species_list[[i]]
    mets[i] <- xmlGetAttr(species, "id")
    metNames[i] <- xmlGetAttr(species, "name", default = mets[i])
  }

  cat("  Found", n_species, "metabolites\n")

  # Extract reactions
  cat("Step 2: Extracting reactions...\n")
  reactions_list <- model_node[["listOfReactions"]]
  n_reactions <- xmlSize(reactions_list)

  rxns <- character(n_reactions)
  rxnNames <- character(n_reactions)
  reaction_metabolites <- vector("list", n_reactions)
  reaction_genes <- vector("list", n_reactions)

  for (i in 1:n_reactions) {
    reaction <- reactions_list[[i]]
    rxns[i] <- xmlGetAttr(reaction, "id")
    rxnNames[i] <- xmlGetAttr(reaction, "name", default = rxns[i])

    # Get reactants
    reactants_list <- reaction[["listOfReactants"]]
    reactants <- character(0)
    if (!is.null(reactants_list)) {
      for (j in 1:xmlSize(reactants_list)) {
        species_ref <- xmlGetAttr(reactants_list[[j]], "species")
        reactants <- c(reactants, species_ref)
      }
    }

    # Get products
    products_list <- reaction[["listOfProducts"]]
    products <- character(0)
    if (!is.null(products_list)) {
      for (j in 1:xmlSize(products_list)) {
        species_ref <- xmlGetAttr(products_list[[j]], "species")
        products <- c(products, species_ref)
      }
    }

    reaction_metabolites[[i]] <- unique(c(reactants, products))

    # Extract gene associations from notes
    notes <- reaction[["notes"]]
    genes <- character(0)

    if (!is.null(notes)) {
      notes_text <- xmlValue(notes)

      # Look for GENE_ASSOCIATION or GENE ASSOCIATION
      if (grepl("GENE[_ ]ASSOCIATION", notes_text, ignore.case = TRUE)) {
        # Extract gene association text
        gene_pattern <- "GENE[_ ]ASSOCIATION[:\\s]+([^<\n]+)"
        gene_match <- regexpr(gene_pattern, notes_text, ignore.case = TRUE, perl = TRUE)

        if (gene_match[1] > 0) {
          gene_text <- regmatches(notes_text, gene_match)
          gene_text <- sub("GENE[_ ]ASSOCIATION[:\\s]+", "", gene_text, ignore.case = TRUE)

          # Parse gene associations: remove operators, parentheses
          gene_text <- gsub("\\(|\\)", " ", gene_text)
          gene_text <- gsub(" and | or ", " ", gene_text, ignore.case = TRUE)
          genes <- unique(trimws(unlist(strsplit(gene_text, "\\s+"))))
          genes <- genes[genes != "" & genes != "and" & genes != "or"]
        }
      }
    }

    reaction_genes[[i]] <- genes

    # Progress
    if (i %% 500 == 0) {
      cat("  Processed", i, "/", n_reactions, "reactions\n")
    }
  }

  cat("  Found", n_reactions, "reactions\n")

  # Extract unique genes
  cat("Step 3: Creating gene list...\n")
  all_genes <- unique(unlist(reaction_genes))
  all_genes <- all_genes[all_genes != ""]
  genes <- sort(all_genes)

  cat("  Found", length(genes), "unique genes\n")

  # Build stoichiometric matrix S (metabolites x reactions)
  cat("Step 4: Building stoichiometric matrix...\n")
  S <- matrix(0, nrow = n_species, ncol = n_reactions)
  rownames(S) <- mets
  colnames(S) <- rxns

  for (i in 1:n_reactions) {
    met_indices <- which(mets %in% reaction_metabolites[[i]])
    S[met_indices, i] <- 1  # Simplified: just mark presence (not stoichiometry)
  }

  cat("  S matrix:", nrow(S), "x", ncol(S), "\n")

  # Build reaction-gene matrix (reactions x genes)
  cat("Step 5: Building reaction-gene matrix...\n")
  rxnGeneMat <- matrix(0, nrow = n_reactions, ncol = length(genes))
  rownames(rxnGeneMat) <- rxns
  colnames(rxnGeneMat) <- genes

  for (i in 1:n_reactions) {
    if (length(reaction_genes[[i]]) > 0) {
      gene_indices <- which(genes %in% reaction_genes[[i]])
      rxnGeneMat[i, gene_indices] <- 1
    }
  }

  cat("  rxnGeneMat matrix:", nrow(rxnGeneMat), "x", ncol(rxnGeneMat), "\n")

  # Create model structure
  model <- list(
    mets = mets,
    metNames = metNames,
    rxns = rxns,
    rxnNames = rxnNames,
    genes = genes,
    S = S,
    rxnGeneMat = rxnGeneMat
  )

  # Summary statistics
  cat("\n=== MODEL SUMMARY ===\n")
  cat("Metabolites:", length(model$mets), "\n")
  cat("Reactions:", length(model$rxns), "\n")
  cat("Genes:", length(model$genes), "\n")
  cat("Reactions with gene associations:", sum(rowSums(rxnGeneMat) > 0), "\n")
  cat("Metabolites with reactions:", sum(rowSums(S) > 0), "\n")
  cat("\n")

  return(model)
}


#' Save parsed model for reuse
#'
#' @param model Model structure
#' @param output_file Path to save RDS file
saveModel <- function(model, output_file) {
  saveRDS(model, output_file)
  cat("Model saved to:", output_file, "\n")
  cat("Load with: model <- readRDS('", output_file, "')\n", sep = "")
}


#' Load saved model
#'
#' @param model_file Path to RDS file
#' @return Model structure
loadModel <- function(model_file) {
  if (!file.exists(model_file)) {
    stop("Model file not found: ", model_file)
  }
  model <- readRDS(model_file)
  cat("Model loaded from:", model_file, "\n")
  return(model)
}

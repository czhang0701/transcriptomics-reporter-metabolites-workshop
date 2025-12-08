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

  # Separate metabolites (M_xxx) from enzymes (E_xxx)
  is_enzyme <- grepl("^E_", mets)
  enzyme_ids <- mets[is_enzyme]
  enzyme_names <- metNames[is_enzyme]

  # Keep only metabolites
  mets <- mets[!is_enzyme]
  metNames <- metNames[!is_enzyme]

  cat("  Found", length(mets), "metabolites\n")
  cat("  Found", length(enzyme_ids), "enzymes\n")

  # Create enzyme-to-gene mapping
  # Extract Ensembl IDs or gene symbols from enzyme names
  genes_from_enzymes <- character(length(enzyme_ids))
  for (i in 1:length(enzyme_ids)) {
    # Enzyme name is either Ensembl ID or gene symbol
    gene_id <- enzyme_names[i]

    # Try to extract Ensembl ID from annotation or notes
    # For now, use the name directly
    if (grepl("^ENSG", gene_id)) {
      genes_from_enzymes[i] <- gene_id
    } else {
      # It might be a gene symbol in the notes
      genes_from_enzymes[i] <- gene_id
    }
  }

  # Create enzyme to gene mapping
  enzyme_to_gene <- data.frame(
    enzyme_id = enzyme_ids,
    gene_id = genes_from_enzymes,
    stringsAsFactors = FALSE
  )

  cat("  Enzyme-gene mapping created\n")

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

    # Filter to keep only metabolites (not enzymes)
    metabolites_only <- c(reactants, products)[grepl("^M_", c(reactants, products))]
    reaction_metabolites[[i]] <- unique(metabolites_only)

    # Extract gene associations from modifiers (enzymes)
    modifiers_list <- reaction[["listOfModifiers"]]
    enzyme_refs <- character(0)

    if (!is.null(modifiers_list)) {
      for (j in 1:xmlSize(modifiers_list)) {
        modifier_ref <- xmlGetAttr(modifiers_list[[j]], "species")
        # Only keep enzyme references (E_xxx)
        if (grepl("^E_", modifier_ref)) {
          enzyme_refs <- c(enzyme_refs, modifier_ref)
        }
      }
    }

    # Map enzymes to genes
    genes <- character(0)
    if (length(enzyme_refs) > 0) {
      gene_indices <- match(enzyme_refs, enzyme_to_gene$enzyme_id)
      genes <- enzyme_to_gene$gene_id[gene_indices[!is.na(gene_indices)]]
      genes <- unique(genes[!is.na(genes)])
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
  n_metabolites <- length(mets)
  S <- matrix(0, nrow = n_metabolites, ncol = n_reactions)
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

################################################################################
# Extract Metabolite-Gene Associations from SBML Model
#
# This script parses the Reference_model.xml SBML file to extract
# metabolite-gene associations, bypassing the rsbml unit validation issues.
#
# Output: metabolite_genes.gmt (GMT format for piano loadGSC)
################################################################################

# Load required library
library(XML)

cat("Extracting metabolite-gene associations from SBML model...\n")

# Parse SBML file as XML
sbml_file <- "data/Reference_model.xml"
doc <- xmlParse(sbml_file)
root <- xmlRoot(doc)

# Get the model node
model <- root[["model"]]

# Initialize lists to store associations
metabolite_genes <- list()

# Extract species (metabolites)
species_list <- model[["listOfSpecies"]]
species_df <- data.frame(
  id = character(),
  name = character(),
  stringsAsFactors = FALSE
)

if (!is.null(species_list)) {
  for (i in 1:xmlSize(species_list)) {
    species <- species_list[[i]]
    species_id <- xmlGetAttr(species, "id")
    species_name <- xmlGetAttr(species, "name")

    if (is.null(species_name)) {
      species_name <- species_id
    }

    species_df <- rbind(species_df, data.frame(
      id = species_id,
      name = species_name,
      stringsAsFactors = FALSE
    ))
  }
}

cat("Found", nrow(species_df), "metabolites\n")

# Extract reactions and gene associations
reactions_list <- model[["listOfReactions"]]
reaction_count <- 0
association_count <- 0

if (!is.null(reactions_list)) {
  for (i in 1:xmlSize(reactions_list)) {
    reaction <- reactions_list[[i]]
    reaction_id <- xmlGetAttr(reaction, "id")
    reaction_count <- reaction_count + 1

    # Get reactants and products (metabolites involved)
    reactants <- reaction[["listOfReactants"]]
    products <- reaction[["listOfProducts"]]

    metabolites <- c()

    # Extract reactant metabolites
    if (!is.null(reactants)) {
      for (j in 1:xmlSize(reactants)) {
        reactant <- reactants[[j]]
        species_ref <- xmlGetAttr(reactant, "species")
        if (!is.null(species_ref)) {
          metabolites <- c(metabolites, species_ref)
        }
      }
    }

    # Extract product metabolites
    if (!is.null(products)) {
      for (j in 1:xmlSize(products)) {
        product <- products[[j]]
        species_ref <- xmlGetAttr(product, "species")
        if (!is.null(species_ref)) {
          metabolites <- c(metabolites, species_ref)
        }
      }
    }

    # Get gene association (using notes or annotation)
    # SBML models typically store gene associations in different ways
    # Let's try to extract from notes
    notes <- xmlChildren(reaction)[["notes"]]
    genes <- c()

    if (!is.null(notes)) {
      # Try to find gene associations in notes
      notes_text <- xmlValue(notes)

      # Look for GENE_ASSOCIATION or GENE ASSOCIATION pattern
      if (grepl("GENE[_ ]ASSOCIATION", notes_text, ignore.case = TRUE)) {
        # Extract gene names (pattern matching)
        gene_match <- regexpr("GENE[_ ]ASSOCIATION[: ]+([^<\n]+)", notes_text, ignore.case = TRUE, perl = TRUE)
        if (gene_match[1] > 0) {
          gene_text <- regmatches(notes_text, gene_match)
          # Parse genes (remove operators like 'and', 'or', parentheses)
          genes_raw <- gsub("GENE[_ ]ASSOCIATION[: ]+", "", gene_text, ignore.case = TRUE)
          genes_raw <- gsub("[()]", " ", genes_raw)
          genes_raw <- gsub(" and | or ", " ", genes_raw, ignore.case = TRUE)
          genes <- unique(trimws(unlist(strsplit(genes_raw, "\\s+"))))
          genes <- genes[genes != ""]
        }
      }
    }

    # If we found genes and metabolites, create associations
    if (length(genes) > 0 && length(metabolites) > 0) {
      metabolites <- unique(metabolites)

      for (met in metabolites) {
        if (!met %in% names(metabolite_genes)) {
          metabolite_genes[[met]] <- c()
        }
        metabolite_genes[[met]] <- unique(c(metabolite_genes[[met]], genes))
        association_count <- association_count + length(genes)
      }
    }
  }
}

cat("Processed", reaction_count, "reactions\n")
cat("Found", association_count, "metabolite-gene associations\n")
cat("Unique metabolites with gene associations:", length(metabolite_genes), "\n")

# Filter metabolites with at least 3 genes and at most 100 genes
filtered_metabolites <- metabolite_genes[sapply(metabolite_genes, length) >= 3 &
                                          sapply(metabolite_genes, length) <= 100]

cat("Metabolites with 3-100 genes:", length(filtered_metabolites), "\n")

# Write to GMT format
# GMT format: NAME\tDESCRIPTION\tGENE1\tGENE2\t...
gmt_file <- "data/metabolite_genes.gmt"
gmt_lines <- c()

for (met_id in names(filtered_metabolites)) {
  # Get metabolite name
  met_name <- species_df$name[species_df$id == met_id]
  if (length(met_name) == 0 || is.na(met_name)) {
    met_name <- met_id
  }

  genes <- filtered_metabolites[[met_id]]

  # Create GMT line
  gmt_line <- paste(c(met_name, met_id, genes), collapse = "\t")
  gmt_lines <- c(gmt_lines, gmt_line)
}

writeLines(gmt_lines, gmt_file)

cat("\nGMT file saved to:", gmt_file, "\n")
cat("Ready to use with: loadGSC('", gmt_file, "')\n", sep = "")

# Also save as a data frame for alternative loading
metabolite_gene_df <- data.frame(
  Metabolite = character(),
  Gene = character(),
  stringsAsFactors = FALSE
)

for (met_id in names(filtered_metabolites)) {
  met_name <- species_df$name[species_df$id == met_id]
  if (length(met_name) == 0 || is.na(met_name)) {
    met_name <- met_id
  }

  genes <- filtered_metabolites[[met_id]]

  for (gene in genes) {
    metabolite_gene_df <- rbind(metabolite_gene_df, data.frame(
      Metabolite = met_name,
      Gene = gene,
      stringsAsFactors = FALSE
    ))
  }
}

write.table(
  metabolite_gene_df,
  file = "data/metabolite_genes.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Data frame saved to: data/metabolite_genes.txt\n")
cat("\nExtraction complete!\n")

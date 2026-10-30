# Protein annotation shared by 04_differential.R and 05_differential_ludovic.R
#
# DIA-NN msstats.csv uses UniProt entry names as Protein (e.g. ALBU_HUMAN, or
# "A_HUMAN;B_HUMAN" for protein groups), while older FragPipe files used
# accessions (A8K2U0) or sp|A8K2U0|A2ML1_HUMAN. Lookup order, per protein:
#   1. lead entry name   -> FASTA Protein_entry
#   2. lead accession    -> FASTA Protein_ID
#   3. DIA-NN protein_map.csv (Protein.Group accession) -> FASTA Protein_ID
#   4. still missing: accession / gene taken from protein_map.csv itself
# The Protein column itself is never modified (it is the join key for
# ProteinLevelData in prepare_comparison_results()).
# Uses `fasta_file` and (optionally) `protein_map_file` set by run_msstats.Rmd.

library(tidyverse)
library(Biostrings)

# FASTA is parsed once per report, even when sourced for every comparison
if (!exists(".fasta_annotation")) .fasta_annotation <- NULL

build_fasta_annotation <- function(fasta_file) {
  message("Extracting annotations from FASTA database...")
  fasta   <- readAAStringSet(fasta_file)
  headers <- names(fasta)

  # sp|P02768|ALBU_HUMAN ... (optionally prefixed, e.g. contam_sp|...)
  ids <- str_match(headers, "(sp|tr)\\|([^|\\s]+)\\|(\\S+)")

  tibble(header = headers) %>%
    mutate(
      db            = ids[, 2],
      Protein_ID    = ids[, 3],
      Protein_entry = ids[, 4],
      Protein_name  = str_extract(header, "(?<= )[^|]+?(?= OS=)"),
      Organism      = str_extract(header, "(?<=OS=).+?(?= OX=)"),
      TaxID         = str_extract(header, "(?<=OX=)\\d+"),
      Gene_name     = str_extract(header, "(?<=GN=)\\S+"),
      PE_level      = str_extract(header, "(?<=PE=)\\d"),
      Seq_version   = str_extract(header, "(?<=SV=)\\d+"),
      Length        = width(fasta)
    ) %>%
    filter(!str_detect(header, "^(rev_|DECOY_)"), !is.na(Protein_ID)) %>%
    # prefer the regular entry over a contaminant copy of the same protein
    arrange(str_detect(header, "^contam_")) %>%
    distinct(Protein_ID, .keep_all = TRUE)
}

read_protein_map <- function() {
  if (!exists("protein_map_file") || is.null(protein_map_file) ||
      !file.exists(protein_map_file) || file.size(protein_map_file) == 0) {
    return(NULL)
  }
  read.csv(protein_map_file, stringsAsFactors = FALSE, na.strings = c("", "NA"))
}

first_member <- function(x) sub(";.*$", "", x)

# Yes/No flag for the user's proteins of interest (highlight_genes_file).
# An entry matches the gene name, entry name, accession or Protein id
# (case-insensitive; any member of a protein group counts).
mark_highlighted <- function(df, highlight_genes) {
  targets <- toupper(trimws(highlight_genes))
  ids     <- paste(df$Protein, df$Protein_ID, df$Protein_entry, df$Gene_name, sep = ";")
  hit     <- vapply(strsplit(toupper(ids), "[;|]"),
                    function(x) any(trimws(x) %in% targets), logical(1))
  factor(ifelse(hit, "Yes", "No"), levels = c("Yes", "No"))
}

annotate_protein_results <- function(results_list,
                                     protein_col = "Protein") {

  if (is.null(.fasta_annotation)) {
    .fasta_annotation <<- build_fasta_annotation(fasta_file)
  }
  annotation  <- .fasta_annotation
  protein_map <- read_protein_map()

  annotation_cols <- c("Protein_ID", "Protein_entry", "Protein_name",
                       "Gene_name", "Organism", "Length")

  annotated_list <- lapply(results_list, function(df) {
    df       <- as.data.frame(df)
    proteins <- as.character(df[[protein_col]])
    lead     <- first_member(proteins)

    # sp|ACC|ENTRY -> accession + entry name; otherwise the id itself
    piped      <- str_match(lead, "\\|([^|]+)\\|([^|\\s]+)")
    lead_entry <- ifelse(is.na(piped[, 3]), lead, piped[, 3])
    lead_acc   <- ifelse(is.na(piped[, 2]), lead, piped[, 2])

    idx <- match(lead_entry, annotation$Protein_entry)
    idx <- ifelse(is.na(idx), match(lead_acc, annotation$Protein_ID), idx)

    map_acc <- map_gene <- rep(NA_character_, length(proteins))
    if (!is.null(protein_map)) {
      m        <- match(proteins, protein_map$ProteinName)
      map_acc  <- first_member(protein_map$Protein_group[m])
      if ("Genes_diann" %in% names(protein_map)) {
        map_gene <- first_member(protein_map$Genes_diann[m])
      }
      idx <- ifelse(is.na(idx), match(map_acc, annotation$Protein_ID), idx)
    }

    ann <- as.data.frame(annotation[idx, annotation_cols])
    ann$Protein_ID <- ifelse(is.na(ann$Protein_ID), map_acc, ann$Protein_ID)
    ann$Gene_name  <- ifelse(is.na(ann$Gene_name),  map_gene, ann$Gene_name)

    message(sprintf("Annotated %d / %d proteins from FASTA (%d more via DIA-NN protein map)",
                    sum(!is.na(idx)), length(idx),
                    sum(is.na(idx) & !is.na(ann$Protein_ID))))

    hl <- if (exists("highlight_genes")) highlight_genes else character(0)
    ann$Highlighted <- mark_highlighted(cbind(Protein = proteins, ann), hl)
    if (length(hl) > 0) {
      message(sprintf("Highlighted %d proteins from %d listed genes",
                      sum(ann$Highlighted == "Yes"), length(hl)))
    }

    # Annotations right after Protein, then the MSstats result columns
    df_annotated <- cbind(df[protein_col], ann,
                          df[setdiff(names(df), c(protein_col, annotation_cols, "Highlighted"))])

    # Sort by adjusted p-value if available
    if ("adj.pvalue" %in% names(df_annotated)) {
      df_annotated <- df_annotated[order(df_annotated$adj.pvalue), ]
    }

    df_annotated
  })

  message("Annotation complete from FASTA")

  annotated_list
}

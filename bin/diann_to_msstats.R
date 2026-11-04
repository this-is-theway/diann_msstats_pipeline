#!/usr/bin/env Rscript
# =============================================================================
#  Convert a DIA-NN 2.x report.parquet into MSstats input, the way recommended
#  by the DIA-NN developers:
#
#    1. DIA-NN is run with --export-quant, so report.parquet carries the
#       per-fragment intensities in Fr.0.Quantity ... Fr.N.Quantity.
#    2. These are joined with ';' into Fragment.Quant.Corrected, the column
#       MSstats' DIA-NN importer expects (DIA-NN 1.8 layout).
#    3. MSstats' DIANNtoMSstatsFormat() does the filtering and conversion
#       (decoys, q-values, unique peptides, few-measurement features).
#
#  Outputs
#    msstats.csv      fragment-level MSstats input (output of DIANNtoMSstatsFormat)
#    protein_map.csv  ProteinName -> DIA-NN protein group / accessions / genes,
#                     used as an annotation fallback in the report
#
#  Usage
#    diann_to_msstats.R --report=report.parquet --design=design.csv \
#        [--out=msstats.csv] [--protein-map=protein_map.csv] [--mbr=TRUE] \
#        [--qvalue=0.01] [--pg-qvalue=0.01] [--global-qvalue=0.01] \
#        [--log=MSstats_import.log]
# =============================================================================

suppressPackageStartupMessages({
  library(arrow)
  library(MSstatsAlt)
  library(data.table)
})

# ---- arguments (--key=value) ------------------------------------------------
defaults <- list(
  report          = NULL,
  design          = NULL,
  out             = "msstats.csv",
  "protein-map"   = "protein_map.csv",
  mbr             = "TRUE",
  qvalue          = "0.01",
  "pg-qvalue"     = "0.01",
  "global-qvalue" = "0.01",
  log             = "MSstats_import.log"
)
args <- defaults
for (a in commandArgs(trailingOnly = TRUE)) {
  kv <- regmatches(a, regexec("^--([^=]+)=(.*)$", a))[[1]]
  if (length(kv) != 3 || !kv[2] %in% names(defaults)) stop("Unknown argument: ", a)
  args[[kv[2]]] <- kv[3]
}
if (is.null(args$report) || is.null(args$design)) stop("--report and --design are required")
mbr <- as.logical(args$mbr)

ms_extensions <- c(".raw", ".mzml", ".mzxml", ".d", ".wiff", ".wiff2", ".dia")

# File name without directory and MS extension, as DIA-NN reports Run
run_name <- function(path) {
  name  <- basename(sub("[/\\\\]+$", "", trimws(path)))
  lower <- tolower(name)
  for (ext in ms_extensions) {
    if (endsWith(lower, ext)) return(substr(name, 1, nchar(name) - nchar(ext)))
  }
  name
}

# ---- annotation from the experimental design -------------------------------
read_design <- function(path) {
  design <- read.csv(path, colClasses = "character", na.strings = character(0),
                     check.names = FALSE)
  names(design) <- tolower(trimws(names(design)))
  for (col in c("file", "experiment")) {
    if (!col %in% names(design)) stop("experimental design is missing column '", col, "'")
  }
  if (!"bioreplicate" %in% names(design)) design$bioreplicate <- ""

  run  <- vapply(design$file, run_name, character(1), USE.NAMES = FALSE)
  cond <- trimws(design$experiment)
  rep  <- trimws(design$bioreplicate)
  # BioReplicate must be unique per sample (otherwise MSstats treats the
  # design as repeated measures and the QC pivot collapses samples).
  annot <- data.frame(
    Run          = run,
    Condition    = cond,
    BioReplicate = ifelse(rep != "", paste(cond, rep, sep = "_"), run),
    stringsAsFactors = FALSE
  )
  if (anyDuplicated(annot$Run)) {
    stop("duplicated runs in experimental design: ",
         paste(unique(annot$Run[duplicated(annot$Run)]), collapse = ", "))
  }
  annot
}

# ---- read report.parquet ----------------------------------------------------
input <- as.data.frame(read_parquet(args$report))
message(sprintf("Read %s precursor rows from %s", format(nrow(input), big.mark = ","), args$report))

# Fragment columns written by --export-quant: Fr.0.Quantity, Fr.1.Quantity, ...
c.frags <- grep("^Fr\\.[0-9]+\\.Quantity$", names(input), value = TRUE)
if (length(c.frags) == 0) {
  stop("no Fr.N.Quantity columns in ", args$report,
       ": run DIA-NN with --export-quant to export fragment intensities")
}
c.frags <- c.frags[order(as.integer(sub("^Fr\\.([0-9]+)\\.Quantity$", "\\1", c.frags)))]
message(sprintf("Fragment columns: %s ... %s (%d)", c.frags[1], tail(c.frags, 1), length(c.frags)))

# Same as the DIA-NN developers' recipe:
#   apply(input[, c.frags], 1, function(row) paste(as.character(row), collapse = ";"))
# vectorised (missing fragments become "NA" and are dropped by the importer)
input$Fragment.Quant.Corrected <- do.call(paste, c(lapply(input[c.frags], as.character), sep = ";"))
input[c.frags] <- NULL

# DIA-NN 2.x reports Run without extension; strip it anyway in case of paths
input$Run <- vapply(as.character(input$Run), run_name, character(1), USE.NAMES = FALSE)

annot <- read_design(args$design)
report_runs <- unique(input$Run)
extra   <- setdiff(report_runs, annot$Run)
missing <- setdiff(annot$Run, report_runs)
if (length(extra))   warning("runs not in experimental design (dropped): ", paste(sort(extra), collapse = ", "))
if (length(missing)) warning("design runs absent from DIA-NN report: ", paste(sort(missing), collapse = ", "))
input <- input[input$Run %in% annot$Run, , drop = FALSE]
annot <- annot[annot$Run %in% report_runs, , drop = FALSE]
if (nrow(input) == 0) {
  stop("no rows left after matching DIA-NN runs to the experimental design ",
       "(check the 'file' column against DIA-NN Run names)")
}

# ---- MSstats conversion ------------------------------------------------------
# Note MSstatsConvert's argument naming: global_qvalue_cutoff is applied to
# Q.Value (run-specific precursor FDR), qvalue_cutoff to Lib.Q.Value (MBR) or
# Global.Q.Value, pg_qvalue_cutoff to Lib.PG.Q.Value (MBR) or Global.PG.Q.Value.
msstats <- DIANNtoMSstatsFormat(
  input,
  annotation                 = annot,
  MBR                        = mbr,
  global_qvalue_cutoff       = as.numeric(args$qvalue),
  qvalue_cutoff              = as.numeric(args[["global-qvalue"]]),
  pg_qvalue_cutoff           = as.numeric(args[["pg-qvalue"]]),
  useUniquePeptide           = TRUE,
  removeFewMeasurements      = TRUE,
  removeProtein_with1Feature = FALSE,
  use_log_file               = TRUE,
  append                     = FALSE,
  verbose                    = TRUE,
  log_file_path              = args$log
)
msstats <- as.data.frame(msstats)

fwrite(msstats, args$out)
message(sprintf("Wrote %s: %s rows, %s proteins, %d runs", args$out,
                format(nrow(msstats), big.mark = ","),
                format(length(unique(msstats$ProteinName)), big.mark = ","),
                length(unique(msstats$Run))))

# ---- protein map (annotation fallback) --------------------------------------
map_cols <- c(ProteinName = "Protein.Names", Protein_group = "Protein.Group",
              Protein_ids = "Protein.Ids", Genes_diann = "Genes")
map_cols <- map_cols[map_cols %in% names(input)]
pmap <- input[, map_cols, drop = FALSE]
names(pmap) <- names(map_cols)
pmap <- pmap[!is.na(pmap$ProteinName) & pmap$ProteinName != "", , drop = FALSE]
pmap <- pmap[!duplicated(pmap$ProteinName), , drop = FALSE]
pmap <- pmap[pmap$ProteinName %in% msstats$ProteinName, , drop = FALSE]
write.csv(pmap, args[["protein-map"]], row.names = FALSE)
message(sprintf("Wrote %s: %s proteins", args[["protein-map"]], format(nrow(pmap), big.mark = ",")))

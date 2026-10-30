
# Load data
# Uses `input_file` and `log_dir` variable set by run_msstats.Rmd from params
# input_file is msstats.csv built from the DIA-NN report.parquet (DIANN_TO_MSSTATS).
# It has the standard 10-column MSstats layout, so FragPipetoMSstatsFormat() is used
# as a generic importer (unique peptides, few-measurement filter, duplicate handling).

msstats_format <- FragPipetoMSstatsFormat(
  input = input_file,
  useUniquePeptide = TRUE,
  removeFewMeasurements = TRUE,
  removeProtein_with1Feature = FALSE,
  summaryforMultipleRows = max,
  use_log_file = TRUE,
  append = FALSE,
  verbose = TRUE,
  log_file_path = file.path(log_dir, "MSstats_import.log")
)
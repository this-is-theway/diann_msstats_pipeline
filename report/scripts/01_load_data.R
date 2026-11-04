# Load data
# Uses `input_file` variable set by run_msstats.Rmd from params
# input_file is msstats.csv written by DIANN_TO_MSSTATS: the DIA-NN report.parquet
# (with --export-quant fragment intensities) already converted by MSstats'
# DIANNtoMSstatsFormat() (q-value filters, unique peptides, few-measurement filter).
# It is read as-is and passed straight to dataProcess().

msstats_format <- as.data.frame(data.table::fread(
  input_file,
  colClasses = list(character = c("ProteinName", "PeptideSequence", "FragmentIon",
                                  "IsotopeLabelType", "Condition", "BioReplicate", "Run"))
))

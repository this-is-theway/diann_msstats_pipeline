
# Data preprocessing
library(openxlsx)

summarized = dataProcess(
  msstats_format,
  logTrans = 2,
  normalization = "equalizeMedians",
  featureSubset = "all",
  n_top_feature = 10,
  min_feature_count = 2,
  remove_uninformative_feature_outlier = FALSE,
  summaryMethod = "TMP",
  equalFeatureVar = TRUE,
  censoredInt = "NA",
  imputation_method = "AFT",
  MBimpute = TRUE,
  numberOfCores = 4,
  use_log_file = TRUE,
  append = FALSE,
  log_file_path = file.path(log_dir, "MSstats_preprocessing.log")
)

# Export preprocessed results at both peptide and protein levels
protein_level_data <- summarized$ProteinLevelData
write.xlsx(protein_level_data, "protein_level_data.xlsx")

# Fragment-level data can exceed Excel's 1,048,576-row limit: write .csv.gz then
pepetide_level_data <- summarized$FeatureLevelData
if (nrow(pepetide_level_data) < 1e6) {
  write.xlsx(pepetide_level_data, "pepetide_level_data.xlsx")
} else {
  data.table::fwrite(pepetide_level_data, "pepetide_level_data.csv.gz")
}

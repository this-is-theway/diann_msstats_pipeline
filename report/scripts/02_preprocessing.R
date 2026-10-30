
# Data preprocessing
library(openxlsx)

summarized = dataProcess(
  msstats_format,
  logTrans = 2,
  normalization = "equalizeMedians",
  featureSubset = "topN",
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

pepetide_level_data <- summarized$FeatureLevelData
write.xlsx(pepetide_level_data, "pepetide_level_data.xlsx")

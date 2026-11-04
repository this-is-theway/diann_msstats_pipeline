# Data preprocessing no normalisation

summarized = dataProcess(
  msstats_format,
  logTrans = 2,
  normalization = FALSE,
  featureSubset = "all",
  n_top_feature = 10,
  min_feature_count = 2,
  remove_uninformative_feature_outlier = FALSE,
  summaryMethod = "TMP",
  equalFeatureVar = TRUE,
  censoredInt = "NA",
  MBimpute = FALSE,
  numberOfCores = 4
)
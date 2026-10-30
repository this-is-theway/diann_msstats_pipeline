# Data preprocessing no normalisation

summarized = dataProcess(
  msstats_format,
  logTrans = 2,
  normalization = FALSE,
  featureSubset = "topN",
  summaryMethod = "TMP",
  equalFeatureVar = TRUE,
  censoredInt = "NA",
  MBimpute = FALSE,
  numberOfCores = 4
)
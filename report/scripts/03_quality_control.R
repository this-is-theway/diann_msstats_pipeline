library(PCAtools)
library(tidyverse)

# From MSstats processed data
feature_data <- summarized$FeatureLevelData

# Distribution plot
dist_plot <- ggplot(feature_data, aes(x = SUBJECT, y = ABUNDANCE, fill = GROUP)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    x = "Sample",
    y = "Log2 Intensity",
    fill = "Condition",
    title = "Intensity Distribution per Sample (Post-Normalization)"
  )

# Protein matrix
protein_matrix <- summarized$ProteinLevelData %>%
  dplyr::select(Protein, SUBJECT, LogIntensities) %>%
  pivot_wider(names_from = SUBJECT, values_from = LogIntensities) %>%
  column_to_rownames("Protein")

message(sprintf("Total proteins: %d", nrow(protein_matrix)))

# Check missingness distribution
missing_pct <- rowMeans(is.na(protein_matrix))

# Filter proteins with >30% missing
protein_filtered <- protein_matrix[missing_pct < 0.3, ]
message(sprintf("After filtering: %d proteins", nrow(protein_filtered)))

# Impute remaining NA with row minimum
protein_final <- t(apply(protein_filtered, 1, function(x) {
  x[is.na(x)] <- min(x, na.rm = TRUE)
  x
}))

# Metadata
metadata <- summarized$ProteinLevelData %>%
  dplyr::select(SUBJECT, GROUP, SUBJECT) %>%
  dplyr::distinct() %>%
  dplyr::rename(Condition = GROUP, BioReplicate = SUBJECT) %>%
  as.data.frame()

rownames(metadata) <- metadata$BioReplicate
metadata <- metadata[colnames(protein_final), ]

# PCA
pca_result <- pca(
  mat = protein_final,
  metadata = metadata,
  removeVar = 0.1
)

pca_plot <- biplot(pca_result, colby = "Condition", legendPosition = "right", title = "PCA - Protein Level")
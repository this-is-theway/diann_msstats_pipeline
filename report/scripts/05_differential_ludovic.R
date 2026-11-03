library(openxlsx)
library(ggplot2)
library(dplyr)
library(plotly)
library(crosstalk)

# contrast matrix for the single comparison
contrast_mat <- matrix(c(-1, 1), nrow = 1)
colnames(contrast_mat) <- conditions
rownames(contrast_mat) <- comparison_name

# Fix direction: condition_A (first) gets +1, condition_B (second) gets -1
cond_A <- conditions[1]
cond_B <- conditions[2]
contrast_mat[1, cond_A] <- 1
contrast_mat[1, cond_B] <- -1

test_results <- groupComparison(contrast.matrix = contrast_mat,
                                data = summarized)





# Protein annotation (entry name / accession / DIA-NN protein map -> FASTA)
source(file.path(scripts_dir, "annotation_utils.R"))


# Function to split comparison results into list
split_comparison_results <- function(test_results) {
  
  # Extract comparison results
  results_df <- test_results$ComparisonResult
  
  # Get unique comparison labels
  comparison_labels <- unique(results_df$Label)
  
  # Split results into list of dataframes
  results_list <- list()
  
  for (label in comparison_labels) {
    results_list[[label]] <- results_df[results_df$Label == label, ]
  }
  
  message("Split results into ", length(comparison_labels), " comparisons")
  message("Comparison names: ", paste(names(results_list), collapse = ", "))
  
  return(results_list)
}


# Split results into list
results_list <- split_comparison_results(test_results)

# Annotate from FASTA (+ DIA-NN protein map fallback)
annotated_results <- annotate_protein_results(
  results_list,
  protein_col = "Protein"
)


# Export to excel file
library(openxlsx)
write.xlsx(annotated_results, paste("comparison_results_", comparison_name, ".xlsx"))





# Function to prepare comparison data and interactive plots linked via crosstalk
prepare_comparison_results <- function(comparison_data, comparison_name,
                                       pvalue_cutoff = 0.05, log2fc_cutoff = 1,
                                       highlight_genes = c()) {

  # Significance categories
  comparison_data$Significance <- "NS"
  comparison_data$Significance[comparison_data$adj.pvalue < pvalue_cutoff &
                                 comparison_data$log2FC > log2fc_cutoff] <- "Up"
  comparison_data$Significance[comparison_data$adj.pvalue < pvalue_cutoff &
                                 comparison_data$log2FC < -log2fc_cutoff] <- "Down"

  # Volcano y-axis: -log10(adj.pvalue), cap Inf produced by p==0
  comparison_data$neg_log_pvalue <- -log10(comparison_data$adj.pvalue)

  # MA x-axis: meanLog2Intensity
  protein_means_df <- summarized$ProteinLevelData %>%
              dplyr::select(Protein, LogIntensities) %>%
              dplyr::group_by(Protein) %>%
              dplyr::summarise(meanLog2Intensity = mean(LogIntensities, na.rm = TRUE))
  comparison_data <- comparison_data %>% dplyr::inner_join(protein_means_df, by = "Protein")


  # Hover tooltip (fall back to Protein ID when Gene_name is absent)
  gene_label <- ifelse(
    is.na(comparison_data$Gene_name) | comparison_data$Gene_name == "",
    comparison_data$Protein, comparison_data$Gene_name
  )
  comparison_data$hover_text <- paste0(
    "Gene: ", gene_label,
    "<br>Protein: ", comparison_data$Protein,
    "<br>log2FC: ", round(comparison_data$log2FC, 3),
    "<br>adj.p: ", signif(comparison_data$adj.pvalue, 3),
    "<br>Status: ", comparison_data$Significance,
    "<br>Highlighted: ", comparison_data$Highlighted
  )

  # Build the SharedData dataframe: visible display cols first, then plot-only cols
  display_cols <- c("Protein", "Protein_ID", "Gene_name", 
                    "log2FC", "pvalue", "adj.pvalue", "Highlighted", "Significance")
  plot_cols    <- c("neg_log_pvalue", "meanLog2Intensity", "hover_text")
  keep_cols    <- c(intersect(display_cols, names(comparison_data)), plot_cols)
  sd_data      <- comparison_data[, keep_cols, drop = FALSE]

  # 0-based column indices to hide in the DT (plot-only cols)
  hide_indices <- which(names(sd_data) %in% plot_cols) - 1L

  # Shared crosstalk key linking plots <-> table
  sd <- SharedData$new(sd_data, key = ~Protein)

  sig_colors <- c("Up" = "#d62728", "Down" = "#1f77b4", "NS" = "gray70")


  volcano_plot <- plot_ly(sd,
      x = ~log2FC, y = ~neg_log_pvalue,
      type = "scatter", mode = "markers",
      color = ~Significance, colors = sig_colors,
      text = ~hover_text, hoverinfo = "text",
      marker = list(size = 6, opacity = 0.7,
                    line = list(width = 0.5, color = "white"))
    ) %>%
    layout(
      title = list(
        text = paste0(
          "<b>Volcano Plot: ", comparison_name, "</b><br>",
          "<span style='font-size:11px;color:gray40'>",
          "Cutoffs: adj.p-value < ", pvalue_cutoff,
          ", |log2FC| > ", log2fc_cutoff, "</span>"
        ),
        x = 0.5, xanchor = "center",
        font = list(size = 14)
      ),
      xaxis = list(
        title = "log2(Fold Change)",
        zeroline = FALSE,
        showgrid = TRUE, gridcolor = "rgb(235,235,235)",
        showline = TRUE, linecolor = "black", mirror = TRUE
      ),
      yaxis = list(
        title = "-log10(Adjusted P-value)",
        zeroline = FALSE,
        showgrid = TRUE, gridcolor = "rgb(235,235,235)",
        showline = TRUE, linecolor = "black", mirror = TRUE
      ),
      plot_bgcolor  = "white",
      paper_bgcolor = "white",
      legend = list(
        title = list(text = "Regulation"),
        orientation = "v", x = 1.02, y = 0.5
      ),
      shapes = list(
        list(type = "line", x0 = log2fc_cutoff, x1 = log2fc_cutoff,
            y0 = 0, y1 = 1, yref = "paper",
            line = list(dash = "dash", color = "gray40", width = 1)),
        list(type = "line", x0 = -log2fc_cutoff, x1 = -log2fc_cutoff,
            y0 = 0, y1 = 1, yref = "paper",
            line = list(dash = "dash", color = "gray40", width = 1)),
        list(type = "line", x0 = 0, x1 = 1, xref = "paper",
            y0 = -log10(pvalue_cutoff), y1 = -log10(pvalue_cutoff),
            line = list(dash = "dash", color = "gray40", width = 1))
      ),
      margin = list(t = 80, r = 120)
    ) %>%
    highlight(on = "plotly_selected", off = "plotly_deselect")


    ma_plot <- plot_ly(sd,
      x = ~meanLog2Intensity, y = ~log2FC,
      type = "scatter", mode = "markers",
      color = ~Significance, colors = sig_colors,
      text = ~hover_text, hoverinfo = "text",
      marker = list(size = 6, opacity = 0.7,
                    line = list(width = 0.5, color = "white"))
    ) %>%
    layout(
      title = list(
        text = paste0(
          "<b>MA Plot: ", comparison_name, "</b><br>",
          "<span style='font-size:11px;color:gray40'>",
          "Mean expression vs. fold change (log2FC cutoff: ", log2fc_cutoff, ")",
          "</span>"
        ),
        x = 0.5, xanchor = "center",
        font = list(size = 14)
      ),
      xaxis = list(
        title = "log2(Mean Expression)",
        zeroline = FALSE,
        showgrid = TRUE, gridcolor = "rgb(235,235,235)",
        showline = TRUE, linecolor = "black", mirror = TRUE
      ),
      yaxis = list(
        title = "log2(Fold Change)",
        zeroline = FALSE,
        showgrid = TRUE, gridcolor = "rgb(235,235,235)",
        showline = TRUE, linecolor = "black", mirror = TRUE
      ),
      plot_bgcolor  = "white",
      paper_bgcolor = "white",
      legend = list(
        title = list(text = "Regulation"),
        orientation = "v", x = 1.02, y = 0.5
      ),
      shapes = list(
        list(type = "line", x0 = 0, x1 = 1, xref = "paper",
            y0 = 0, y1 = 0,
            line = list(dash = "dash", color = "black", width = 1)),
        list(type = "line", x0 = 0, x1 = 1, xref = "paper",
            y0 = log2fc_cutoff, y1 = log2fc_cutoff,
            line = list(dash = "dash", color = "gray40", width = 1)),
        list(type = "line", x0 = 0, x1 = 1, xref = "paper",
            y0 = -log2fc_cutoff, y1 = -log2fc_cutoff,
            line = list(dash = "dash", color = "gray40", width = 1))
      ),
      margin = list(t = 80, r = 120)
    ) %>%
    highlight(on = "plotly_selected", off = "plotly_deselect")


  n_up    <- sum(comparison_data$Significance == "Up",   na.rm = TRUE)
  n_down  <- sum(comparison_data$Significance == "Down", na.rm = TRUE)
  n_total <- nrow(comparison_data)
  n_sig   <- n_up + n_down

  list(
    data          = comparison_data,
    volcano_plot  = volcano_plot,
    ma_plot       = ma_plot,
    shared_data   = sd,
    hide_indices  = hide_indices,
    n_up          = n_up,
    n_down        = n_down,
    n_total       = n_total,
    n_sig         = n_sig,
    pvalue_cutoff = pvalue_cutoff,
    log2fc_cutoff = log2fc_cutoff
  )
}
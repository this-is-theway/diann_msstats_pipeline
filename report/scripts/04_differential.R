library(openxlsx)
library(ggplot2)
library(dplyr)
library(plotly)
library(crosstalk)



# Read the comparison file
# Uses `comparisons_file` variable set by run_msstats.Rmd from params
comparisons_df <- read.csv(comparisons_file, stringsAsFactors = FALSE, sep = ",")

# Function to create contrast matrix from CSV file
create_contrast_from_csv <- function(comparisons_df, conditions = NULL, 
                                     annotation = NULL) {
      
      # Check for required columns
      if (!all(c("condition_A", "condition_B") %in% colnames(comparisons_df))) {
        stop("CSV file must contain 'condition_A' and 'condition_B' columns")
      }
      
      # Clean up condition names (remove extra whitespace)
      comparisons_df$condition_A <- trimws(comparisons_df$condition_A)
      comparisons_df$condition_B <- trimws(comparisons_df$condition_B)
      
      # Get conditions from annotation if not provided
      if (is.null(conditions)) {
        if (is.null(annotation)) {
          # Extract all unique conditions from comparison file
          conditions <- unique(c(comparisons_df$condition_A, 
                                comparisons_df$condition_B))
        } else {
          # Get from annotation
          conditions <- unique(annotation$Condition)
        }
      }
      
      # Sort conditions for consistency
      conditions <- sort(conditions)
      
      # Initialize contrast matrix
      n_comparisons <- nrow(comparisons_df)
      n_conditions <- length(conditions)
      contrast_matrix <- matrix(0, nrow = n_comparisons, ncol = n_conditions)
      
      # Set column names
      colnames(contrast_matrix) <- conditions
      
      # Create row names and fill matrix
      row_names <- character(n_comparisons)
      
      for (i in 1:n_comparisons) {
        cond_A <- comparisons_df$condition_A[i]
        cond_B <- comparisons_df$condition_B[i]
        
        # Check if conditions exist
        if (!(cond_A %in% conditions)) {
          stop(paste("Condition", cond_A, "not found in available conditions"))
        }
        if (!(cond_B %in% conditions)) {
          stop(paste("Condition", cond_B, "not found in available conditions"))
        }
        
        # Create comparison name: ConditionA_vs_ConditionB
        row_names[i] <- paste0(cond_A, "_vs_", cond_B)
        
        # Set contrast values (A is numerator, B is denominator)
        contrast_matrix[i, cond_B] <- -1
        contrast_matrix[i, cond_A] <- 1
      }
      
      # Set row names
      rownames(contrast_matrix) <- row_names
      
      return(contrast_matrix)
}


contrast_mat <- create_contrast_from_csv(comparisons_df, 
                                         annotation = metadata)

test_results <- groupComparison(contrast.matrix = contrast_mat,
                                data = summarized,
                                use_log_file = TRUE,
                                append = FALSE,
                                log_file_path = file.path(log_dir, "MSstats_groupComparison.log"))



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


# Function to export results list to Excel
export_results_to_excel <- function(results_list, output_file = "comparison_results.xlsx") {
        
        # Load required package
        if (!require("openxlsx")) {
          stop("Package 'openxlsx' is required. Install with: install.packages('openxlsx')")
        }
        
        # Create workbook
        wb <- createWorkbook()
        
        # Add separate sheet for each comparison
        for (label in names(results_list)) {
          # Clean sheet name (Excel has restrictions)
          sheet_name <- gsub("[:\\\\/*?\\[\\]]", "_", label)
          sheet_name <- substr(sheet_name, 1, 31)  # Excel sheet name limit
          
          addWorksheet(wb, sheet_name)
          writeData(wb, sheet_name, results_list[[label]])
          
          # Format header
          headerStyle <- createStyle(fontSize = 11, fontColour = "#FFFFFF", 
                                     halign = "center", fgFill = "#4F81BD",
                                     border = "TopBottom", borderColour = "#4F81BD",
                                     textDecoration = "bold")
          addStyle(wb, sheet_name, headerStyle, rows = 1, 
                   cols = 1:ncol(results_list[[label]]), gridExpand = TRUE)
          
          # Auto-size columns
          setColWidths(wb, sheet_name, cols = 1:ncol(results_list[[label]]), widths = "auto")
          
          # Freeze top row
          freezePane(wb, sheet_name, firstRow = TRUE)
        }
        
        # Save workbook
        saveWorkbook(wb, output_file, overwrite = TRUE)
        
        message("Results exported to: ", output_file)
        message("Number of sheets: ", length(results_list))
        
        return(invisible(output_file))
}

# Split results into list
results_list <- split_comparison_results(test_results)

# Annotate from FASTA (+ DIA-NN protein map fallback)
annotated_results <- annotate_protein_results(
  results_list,
  protein_col = "Protein"
)

# Export to Excel
export_results_to_excel(
  annotated_results,
  output_file = "comparison_results.xlsx"
)





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
  display_cols <- c("Gene_name", "Highlighted", "Protein", "Protein_ID", "Protein_entry", "Protein_name",
                    "log2FC", "SE", "pvalue", "adj.pvalue", "Significance")
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
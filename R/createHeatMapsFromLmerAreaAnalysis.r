#' @title Generate Heatmaps and CSV Files from Model Results
#'
#' @description
#' This function automates the post-processing of linear mixed-effects model
#' results. It calculates q-values for multiple comparisons correction, sorts
#' brain regions anatomically if specified, and then generates both summary CSV
#' files and publication-ready heatmaps to visualize the statistical effects.
#'
#' @param results_df The dataframe of model results, typically in long format.
#' @param area_col The name of the column specifying the brain area (e.g., "ROI").
#' @param group_col The column to group by for q-value calculation and table generation.
#'   If `NULL`, the entire dataset is treated as a single group.
#' @param pvalue_col The column containing the p-values to be corrected.
#' @param area_file_type A string ("ROI" or "Channel") used to determine if
#'   anatomical sorting should be applied. Default: `"Channel"`.
#' @param output_prefix A string used as the base name for all output files.
#'   Default: `NULL`.
#' @param significance_threshold The significance threshold (based on q-values)
#'   for creating discrete effect tables and for asterisks on heatmaps. Default: 0.05.
#' @param scale_type (Optional) Specifies the color scale for the continuous heatmap.
#'   Use `"log"` (default) to emphasize values near zero or `"linear"` for a standard scale.
#' @param log_power (Optional) A numeric value that controls the intensity of the
#'   `"log"` scale transformation. Default: 0.2.
#' @param asterisk_sz The size (in mm) of the significance asterisks on the
#'   continuous heatmap. Default: 10.
#'
#' @return Nothing. The function is called for its side effect of writing files
#'   (CSVs and PDFs) to the current working directory.
#'
#' @seealso \code{\link{CreateSigBetaTable}}, \code{\link{getQvalues}}
#'
#' @importFrom utils write.table
#'
#' @export
createHeatMapsFromLmerAreaAnalysis <- function(results_df, area_col, group_col, pvalue_col, area_file_type = "Channel", output_prefix = NULL, significance_threshold = 0.05, scale_type = "log", log_power = 0.2, asterisk_sz = 10) {

  # This call correctly handles a NULL group_col.
  # Assumes getQvalues() is a function within your package.
  results_df <- getQvalues(results_df, group_col, pvalue_col)

  if (area_file_type == "ROI") {
    areas <- unique(results_df[[area_col]])
    # Assumes sort_rois_anatomically() is a function within your package.
    desired_anatomical_order <- sort_rois_anatomically(areas)
    results_df$Area_ordered <- factor(results_df[[area_col]], levels = desired_anatomical_order)
    results_df <- results_df[order(results_df$Area_ordered), ]
    results_df$Area_ordered <- NULL
  }
  
  # write.table is imported from 'utils' via the @importFrom tag
  write.table(results_df, paste0(output_prefix, "_FullOutput.csv"), quote = FALSE, sep = ",", row.names = FALSE)
  
  # Assumes CreateSigBetaTable() is a function within your package.
  out_list_discrete <- CreateSigBetaTable(results_df, area_col, "qVal", group_col = group_col, discrete = TRUE, sigThresh = significance_threshold)
  for (ty in names(out_list_discrete)) {
    # Assumes create_forced_color_heatmap() is a function within your package.
    create_forced_color_heatmap(
      out_list_discrete[[ty]], 3:ncol(out_list_discrete[[ty]]), 2,
      file_name = paste(output_prefix, ty, "Discrete_Effect_Table.pdf", sep = "_"),
      fig_width = 11, fig_height = 8.5
    )
  }

  out_list_continuous <- CreateSigBetaTable(results_df, area_col, "qVal", group_col = group_col, discrete = FALSE, sigThresh = significance_threshold)
  for (ty in names(out_list_continuous)) {
    # Assumes create_pvalue_heatmap() is a function within your package.
    create_pvalue_heatmap(
      out_list_continuous[[ty]], 3:ncol(out_list_continuous[[ty]]), 2,
      scale_type = scale_type,
      log_power = log_power,
      file_name = paste(output_prefix, ty, "Continuous_Effect_Table.pdf", sep = "_"),
      fig_width = 11, fig_height = 8.5,
      asterisk_size = asterisk_sz
    )
  }
}
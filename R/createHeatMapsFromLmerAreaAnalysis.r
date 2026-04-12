#' @title Generate Heatmaps and CSV Files from Model Results (Updated for Two-Stage Analysis)
#'
#' @description
#' Updated version that supports three heatmap types:
#'   - **"omnibus"**: Green gradient, unsigned p-values (for omnibus F-tests)
#'   - **"posthoc"**: Red/blue signed p-values with grey for untested ROIs
#'   - **"activation"**: Red/blue signed p-values (for baseline comparisons, no grey)
#'
#' @param results_df The dataframe of model results in long format.
#' @param area_col The name of the column specifying the brain area.
#' @param group_col The column to group by. If NULL, entire dataset is one group.
#' @param sigValue_col The column containing the corrected p-values.
#' @param dv_col The column containing the DV (for sign info). Default: "Beta".
#' @param area_file_type "ROI" or "Channel". Default: "Channel".
#' @param output_prefix Base name for output files. Default: NULL.
#' @param significance_threshold Threshold for asterisks. Default: 0.05.
#' @param heatmap_type One of "omnibus", "posthoc", or "activation". Default: "posthoc".
#' @param scale_type "log" or "linear". Default: "log".
#' @param log_power Controls log scale intensity. Default: 0.2.
#' @param asterisk_sz Asterisk size in mm. Default: 10.
#' @param removeIntercept Remove intercept rows? Default: TRUE.
#'
#' @return Nothing. Writes CSVs and PDFs to the current working directory.
#'
#' @importFrom utils write.table
#'
#' @export
createHeatMapsFromLmerAreaAnalysis <- function(results_df, area_col, group_col,
                                                sigValue_col, dv_col = "Beta",
                                                area_file_type = "Channel",
                                                output_prefix = NULL,
                                                significance_threshold = 0.05,
                                                heatmap_type = "posthoc",
                                                scale_type = "log",
                                                log_power = 0.2,
                                                asterisk_sz = 10,
                                                removeIntercept = TRUE) {

  # Sort ROIs anatomically if applicable
  if (area_file_type == "ROI") {
    areas <- unique(results_df[[area_col]])
    desired_anatomical_order <- sort_rois_anatomically(areas)
    results_df$Area_ordered <- factor(results_df[[area_col]], levels = desired_anatomical_order)
    results_df <- results_df[order(results_df$Area_ordered), ]
    results_df$Area_ordered <- NULL
  }

  # Save full output
  write.table(results_df, paste0(output_prefix, "_FullOutput.csv"),
              quote = FALSE, sep = ",", row.names = FALSE)

  # -----------------------------------------------------------------------
  # Determine output_mode based on heatmap_type
  # -----------------------------------------------------------------------
  if (heatmap_type == "omnibus") {
    output_mode_discrete <- "discrete"
    output_mode_continuous <- "continuous_unsigned"
  } else {
    # "posthoc" and "activation" both use signed values
    output_mode_discrete <- "discrete"
    output_mode_continuous <- "continuous_signed"
  }

  # -----------------------------------------------------------------------
  # DISCRETE heatmap (always uses forced color: -1, 0, 1)
  # For omnibus, sign comes from F-value (always positive), so only 0 and 1
  # -----------------------------------------------------------------------
  out_list_discrete <- CreateSigBetaTable(results_df, area_col, sigValue_col,
                                           dv_col = dv_col, group_col = group_col,
                                           output_mode = output_mode_discrete,
                                           sigThresh = significance_threshold,
                                           removeIntercept = removeIntercept)

  for (ty in names(out_list_discrete)) {
    if (heatmap_type == "omnibus") {
      # For omnibus discrete, use green/white/dark green mapping
      create_forced_color_heatmap(
        out_list_discrete[[ty]], 3:ncol(out_list_discrete[[ty]]), 2,
        color_mapping = c("0" = "white", "1" = "darkgreen"),
        file_name = paste(output_prefix, ty, "Discrete_Effect_Table.pdf", sep = "_"),
        fig_width = 11, fig_height = 8.5
      )
    } else {
      create_forced_color_heatmap(
        out_list_discrete[[ty]], 3:ncol(out_list_discrete[[ty]]), 2,
        file_name = paste(output_prefix, ty, "Discrete_Effect_Table.pdf", sep = "_"),
        fig_width = 11, fig_height = 8.5
      )
    }
  }

  # -----------------------------------------------------------------------
  # CONTINUOUS heatmap (dispatches to the correct function)
  # -----------------------------------------------------------------------
  out_list_continuous <- CreateSigBetaTable(results_df, area_col, sigValue_col,
                                             dv_col = dv_col, group_col = group_col,
                                             output_mode = output_mode_continuous,
                                             sigThresh = significance_threshold,
                                             removeIntercept = removeIntercept)

  for (ty in names(out_list_continuous)) {
    if (heatmap_type == "omnibus") {
      # Green gradient for unsigned omnibus p-values
      create_omnibus_heatmap(
        out_list_continuous[[ty]], 3:ncol(out_list_continuous[[ty]]), 2,
        scale_type = scale_type,
        log_power = log_power,
        file_name = paste(output_prefix, ty, "Continuous_Effect_Table.pdf", sep = "_"),
        fig_width = 11, fig_height = 8.5,
        asterisk_size = asterisk_sz
      )
    } else if (heatmap_type == "posthoc") {
      # Red/blue with grey for untested ROIs
      create_posthoc_heatmap(
        out_list_continuous[[ty]], 3:ncol(out_list_continuous[[ty]]), 2,
        scale_type = scale_type,
        log_power = log_power,
        file_name = paste(output_prefix, ty, "Continuous_Effect_Table.pdf", sep = "_"),
        fig_width = 11, fig_height = 8.5,
        asterisk_size = asterisk_sz
      )
    } else {
      # "activation" — standard red/blue, no grey (all cells are tested)
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
}

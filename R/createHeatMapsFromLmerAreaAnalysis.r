#' @title Generate Heatmaps and CSV Files from Model Results (Updated for Two-Stage Analysis)
#'
#' @description
#' Updated version that supports three heatmap types:
#'   - **"omnibus"**: Green gradient, unsigned p-values (for omnibus F-tests)
#'   - **"posthoc"**: Red/blue signed p-values with grey for untested ROIs
#'   - **"activation"**: Red/blue signed p-values (for baseline comparisons, no grey)
#'
#' When `split_by_type = TRUE` (default) the `_FullOutput.csv` and each
#' per-type heatmap PDF are written into a `<output_dir>/<tolower(type)>/`
#' subdirectory, one subdirectory per chromophore / measurement type.
#' When `FALSE`, all files are written directly into `output_dir`
#' (pre-split layout).
#'
#' @param results_df The dataframe of model results in long format.
#' @param area_col The name of the column specifying the brain area.
#' @param group_col The column to group by. If NULL, entire dataset is one group.
#' @param sigValue_col The column containing the corrected p-values.
#' @param dv_col The column containing the DV (for sign info). Default: "Beta".
#' @param area_order Optional character vector specifying the desired
#'   x-axis order of `area_col` levels. If `NULL` (default), areas
#'   appear in the order they occur in `results_df`. To sort
#'   anatomically, precompute the order with [sort_rois_anatomically()]
#'   (for ROI data) or [sort_fnirs_channels_by_name()] (for channel
#'   data) and pass the result here. Any sort function with signature
#'   `function(names, ...) -> ordered names` can be used.
#' @param output_prefix Base name for output files (no directory
#'   component). Default: NULL.
#' @param output_dir Parent directory for output files. When
#'   `split_by_type = TRUE`, per-type subdirectories are created under
#'   this directory. Default `"."` (current working directory).
#' @param split_by_type Logical. Split `_FullOutput.csv` and each
#'   per-type heatmap PDF into `<output_dir>/<tolower(type)>/`
#'   subdirectories (`TRUE`, default) or write everything flat into
#'   `output_dir` (`FALSE`). The split requires a `type` column in
#'   `results_df`.
#' @param type_col Name of the type column used for splitting. Default
#'   `"type"`.
#' @param significance_threshold Threshold for asterisks. Default: 0.05.
#' @param heatmap_type One of "omnibus", "posthoc", or "activation". Default: "posthoc".
#' @param scale_type "log" or "linear". Default: "log".
#' @param log_power Controls log scale intensity. Default: 0.2.
#' @param asterisk_sz Asterisk size in mm. Default: 10.
#' @param fig_width Figure width in inches, forwarded to every primitive
#'   heatmap call inside this wrapper. Default: 11.
#' @param fig_height Figure height in inches, forwarded to every
#'   primitive heatmap call inside this wrapper. Default: 8.5. Reduce
#'   when a heatmap has fewer rows than usual so each row keeps its
#'   typical visual height; a small helper such as
#'   `heatmap_height(df, full_n = 6)` in the calling script is one
#'   way to compute this.
#' @param removeIntercept Remove intercept rows? Default: TRUE.
#'
#' @return Nothing. Writes CSVs and PDFs under `output_dir`.
#'
#' @importFrom utils write.table
#'
#' @export
createHeatMapsFromLmerAreaAnalysis <- function(results_df, area_col, group_col,
                                                sigValue_col, dv_col = "Beta",
                                                area_order = NULL,
                                                output_prefix = NULL,
                                                output_dir    = ".",
                                                split_by_type = TRUE,
                                                type_col      = "type",
                                                significance_threshold = 0.05,
                                                heatmap_type = "posthoc",
                                                scale_type = "log",
                                                log_power = 0.2,
                                                asterisk_sz = 10,
                                                fig_width = 11,
                                                fig_height = 8.5,
                                                removeIntercept = TRUE) {

  # Apply caller-supplied area ordering, if any
  if (!is.null(area_order)) {
    missing_levels <- setdiff(unique(results_df[[area_col]]), area_order)
    if (length(missing_levels) > 0) {
      stop("area_order is missing ", length(missing_levels),
           " level(s) present in `", area_col, "`: ",
           paste(missing_levels, collapse = ", "), call. = FALSE)
    }
    results_df$Area_ordered <- factor(results_df[[area_col]], levels = area_order)
    results_df <- results_df[order(results_df$Area_ordered), ]
    results_df$Area_ordered <- NULL
  }

  # --- Path helpers ---------------------------------------------------------
  # Resolve the directory for a given type (type may be NULL if not splitting).
  resolve_dir <- function(ty) {
    if (split_by_type && !is.null(ty)) {
      d <- file.path(output_dir, tolower(as.character(ty)))
      if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
      d
    } else {
      if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      }
      output_dir
    }
  }

  # Save full output (split per-type if requested)
  if (split_by_type) {
    if (!(type_col %in% colnames(results_df))) {
      stop("`split_by_type = TRUE` requires a `", type_col,
           "` column in `results_df`.", call. = FALSE)
    }
    for (ty in unique(results_df[[type_col]])) {
      d <- resolve_dir(ty)
      df_sub <- results_df[results_df[[type_col]] == ty, , drop = FALSE]
      write.table(df_sub,
                  file.path(d, paste0(output_prefix, "_FullOutput.csv")),
                  quote = FALSE, sep = ",", row.names = FALSE)
    }
  } else {
    d <- resolve_dir(NULL)
    write.table(results_df,
                file.path(d, paste0(output_prefix, "_FullOutput.csv")),
                quote = FALSE, sep = ",", row.names = FALSE)
  }

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
    d <- resolve_dir(ty)
    fname <- file.path(d, paste(output_prefix, ty,
                                "Discrete_Effect_Table.pdf", sep = "_"))
    if (heatmap_type == "omnibus") {
      # For omnibus discrete, use green/white/dark green mapping
      create_forced_color_heatmap(
        out_list_discrete[[ty]], 3:ncol(out_list_discrete[[ty]]), 2,
        color_mapping = c("0" = "white", "1" = "darkgreen"),
        file_name = fname,
        fig_width = fig_width, fig_height = fig_height
      )
    } else {
      create_forced_color_heatmap(
        out_list_discrete[[ty]], 3:ncol(out_list_discrete[[ty]]), 2,
        file_name = fname,
        fig_width = fig_width, fig_height = fig_height
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
    d <- resolve_dir(ty)
    fname <- file.path(d, paste(output_prefix, ty,
                                "Continuous_Effect_Table.pdf", sep = "_"))
    if (heatmap_type == "omnibus") {
      # Green gradient for unsigned omnibus p-values
      create_omnibus_heatmap(
        out_list_continuous[[ty]], 3:ncol(out_list_continuous[[ty]]), 2,
        scale_type = scale_type,
        log_power = log_power,
        file_name = fname,
        fig_width = fig_width, fig_height = fig_height,
        asterisk_size = asterisk_sz
      )
    } else if (heatmap_type == "posthoc") {
      # Red/blue with grey for untested ROIs
      create_posthoc_heatmap(
        out_list_continuous[[ty]], 3:ncol(out_list_continuous[[ty]]), 2,
        scale_type = scale_type,
        log_power = log_power,
        file_name = fname,
        fig_width = fig_width, fig_height = fig_height,
        asterisk_size = asterisk_sz
      )
    } else {
      # "activation" - standard red/blue, no grey (all cells are tested)
      create_pvalue_heatmap(
        out_list_continuous[[ty]], 3:ncol(out_list_continuous[[ty]]), 2,
        scale_type = scale_type,
        log_power = log_power,
        file_name = fname,
        fig_width = fig_width, fig_height = fig_height,
        asterisk_size = asterisk_sz
      )
    }
  }
}

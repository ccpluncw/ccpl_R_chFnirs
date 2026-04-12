#' @title Iterate a per-cell analysis across type x area combinations
#'
#' @description
#' Shared iteration engine used internally by all Tier 2 area-wise
#' analysis functions (`run_omnibus_by_area`, `run_activation_tests`,
#' `run_gated_posthoc`, `run_slopes_vs_zero`, `run_gated_slopes`). It
#' loops over unique values of `type_col` and `area_col`, subsets the
#' data for each cell, fits a model via `fit_fn`, extracts a per-cell
#' result data frame via `extract_fn`, stamps `type` and `area` onto
#' every returned row, and `rbind`s the results.
#'
#' Model fitting and extraction errors for individual cells are caught
#' with `tryCatch` and logged via `message()`; the loop continues to
#' the next cell so that a single bad channel or region does not abort
#' an entire analysis. Cells with fewer than `min_rows` observations
#' are silently skipped.
#'
#' This helper exists to eliminate the ~15-line nested-loop boilerplate
#' that would otherwise be duplicated across every public Tier 2
#' function. Each public function builds short `fit_fn` and `extract_fn`
#' closures that capture its specific model type, emmeans call, and
#' column layout, and passes them here.
#'
#' This function is **internal**. It is not exported and is not part of
#' the chFnirs public API. Users should call the public wrappers
#' (`run_omnibus_by_area()` etc.), not this helper directly.
#'
#' @param data A long-format data.frame containing the raw (or
#'   preprocessed) observations.
#' @param area_col Character scalar naming the area column in `data`
#'   (e.g., `"ROI"`, `"Region"`, `"Channel"`).
#' @param type_col Character scalar naming the measurement-type column
#'   in `data` (defaults to `"type"` in all public callers; holds
#'   chromophore labels like `"hbo"`, `"hbr"`, `"cbsi"`, `"hbdiff"`).
#' @param min_rows Integer scalar. Cells with fewer observations than
#'   this are skipped without warning.
#' @param fit_fn A closure of the form `function(df_cell)` that takes
#'   a per-cell subset of `data` and returns a fitted model object
#'   (e.g., an `lm` or `lmerModLmerTest`).
#' @param extract_fn A closure of the form `function(fit, df_cell)`
#'   that takes a fitted model and the cell data, and returns a
#'   data.frame of per-cell results. May return multiple rows (one per
#'   effect, factor level, or contrast) or `NULL` / empty data.frame
#'   for cells that should be skipped silently.
#' @param fn_label Character scalar used in the tryCatch error message
#'   to identify the calling function for the user (e.g.,
#'   `"omnibus test"`, `"activation test"`).
#'
#' @return A data.frame containing every successful cell's results,
#'   with `type` and `area` as the first two columns. Returns `NULL`
#'   if no cell produced output.
#'
#' @keywords internal
iterate_by_area <- function(data, area_col, type_col, min_rows,
                            fit_fn, extract_fn, fn_label) {

  if (!area_col %in% colnames(data)) {
    stop("area_col '", area_col, "' not found in data.", call. = FALSE)
  }
  if (!type_col %in% colnames(data)) {
    stop("type_col '", type_col, "' not found in data.", call. = FALSE)
  }

  types <- unique(data[[type_col]])
  areas <- unique(data[[area_col]])
  df_out <- NULL

  for (ty in types) {
    for (ar in areas) {
      df_cell <- data[data[[area_col]] == ar & data[[type_col]] == ty, ,
                      drop = FALSE]
      if (nrow(df_cell) < min_rows) next

      tryCatch({
        fit <- fit_fn(df_cell)
        rows <- extract_fn(fit, df_cell)

        if (!is.null(rows) && is.data.frame(rows) && nrow(rows) > 0) {
          rows$type <- ty
          rows$area <- ar
          df_out <- rbind(df_out, rows)
        }
      }, error = function(e) {
        message("  [WARNING] ", fn_label, " failed for type = ", ty,
                ", area = ", ar, ": ", conditionMessage(e))
      })
    }
  }

  # Move type and area to the front for consistent output schema
  if (!is.null(df_out)) {
    id_cols    <- c("type", "area")
    other_cols <- setdiff(names(df_out), id_cols)
    df_out     <- df_out[, c(id_cols, other_cols), drop = FALSE]
    rownames(df_out) <- NULL
  }

  df_out
}

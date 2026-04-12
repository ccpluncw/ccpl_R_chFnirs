#' @title Create a Significance Table from Model Results
#'
#' @description
#' Reshapes a long-format data frame of model results (typically the
#' output of `run_omnibus_by_area()`, `run_gated_posthoc()`, or one of
#' the other Tier 2 analysis functions) into a wide-format
#' `predictor x area` table suitable for heatmap rendering. Supports
#' three output modes:
#' \itemize{
#'   \item \strong{discrete} — cells become -1, 0, or 1 based on the
#'     sign of `dv_col` and whether `sigColName` is below `sigThresh`.
#'   \item \strong{continuous_signed} — cells become the signed
#'     p-value (`sigColName * sign(dv_col)`), preserving direction
#'     information.
#'   \item \strong{continuous_unsigned} — cells become the absolute
#'     value of `sigColName`, for omnibus F-test results where there
#'     is no directional coefficient.
#' }
#' NA values in the input (e.g., from non-significant cells in a
#' gated post-hoc output) are preserved as NA in the output, so
#' downstream heatmap functions can render them as grey tiles.
#'
#' The function accepts both the Tier 2 snake_case schema
#' (`predictor`) and the pre-Tier-2 schema (`Predictor`). When
#' `predColName` is left at its default of `NULL`, the function
#' inspects `df.data` for a `"predictor"` column first and falls
#' back to `"Predictor"`; if neither is present, it errors and asks
#' the caller to pass `predColName` explicitly. Callers who want to
#' use a different column name can pass it directly. The output
#' column created inside each returned table uses the same name
#' that was detected or supplied, preserving the caller's schema.
#'
#' @param df.data A long-format data frame of model results. Must
#'   contain a predictor column (`"predictor"` or `"Predictor"`, or
#'   whatever name is passed via `predColName`) plus the columns
#'   named by `areaColName`, `sigColName`, and `dv_col`. Typically
#'   the output of `run_gated_posthoc()`, `run_activation_tests()`,
#'   or one of the other Tier 2 functions.
#' @param areaColName Character scalar naming the area column
#'   (e.g., `"area"`, `"ROI"`, `"Region"`). Default `"area"` matches
#'   the Tier 2 output schema.
#' @param sigColName Character scalar naming the significance column
#'   to use (e.g., `"p_bh"` for BH-adjusted omnibus p-values,
#'   `"p_bon"` for Bonferroni-adjusted post-hoc p-values,
#'   `"p_value"` for raw p-values). Required, no default — the caller
#'   must choose which p-value column is appropriate.
#' @param dv_col Character scalar naming the column whose sign gives
#'   the direction of the effect. Default `"estimate"` matches the
#'   Tier 2 output schema. For omnibus output use
#'   `output_mode = "continuous_unsigned"` and `dv_col` is unused.
#' @param predColName Character scalar naming the predictor column in
#'   `df.data`. Default `NULL` triggers auto-detection: `"predictor"`
#'   is tried first, `"Predictor"` second, and the function errors if
#'   neither is present. Pass explicitly to override.
#' @param group_col Optional character scalar naming a column to
#'   split on before reshaping, producing one wide table per group.
#'   Default `NULL` produces a single table under the name
#'   `"all_data"`.
#' @param output_mode One of `"discrete"`, `"continuous_signed"`, or
#'   `"continuous_unsigned"`. Default `"discrete"`. For backwards
#'   compatibility also accepts a logical `TRUE`/`FALSE` — `TRUE`
#'   maps to `"discrete"` and `FALSE` to `"continuous_signed"`.
#' @param sigThresh Numeric significance threshold used in
#'   `"discrete"` mode. Default `0.05`.
#' @param removeIntercept Logical. If `TRUE` (default), rows whose
#'   predictor begins with `"(Int"` are dropped before reshaping.
#'
#' @return A named list of data frames — one per level of `group_col`
#'   when provided, otherwise a single-element list with name
#'   `"all_data"`. Each element has one row per predictor (minus the
#'   intercept if removed) and one column per area, with the
#'   predictor column (and `group_col` if supplied) as the leading
#'   columns.
#'
#' @importFrom dplyr relocate
#' @importFrom rlang sym
#' @importFrom magrittr %>%
#'
#' @export
CreateSigBetaTable <- function(df.data, areaColName = "area", sigColName,
                               dv_col = "estimate", predColName = NULL,
                               group_col = NULL,
                               output_mode = "discrete",
                               sigThresh = 0.05, removeIntercept = TRUE) {

  # Backwards compatibility: accept discrete = TRUE/FALSE
  if (is.logical(output_mode)) {
    output_mode <- ifelse(output_mode, "discrete", "continuous_signed")
  }

  # Auto-detect the predictor column when the caller did not supply one.
  # Tier 2 data uses lowercase `predictor`; pre-Tier-2 data (e.g., ABC)
  # uses capitalized `Predictor`. The check is deterministic (lowercase
  # first, then capital) and errors loudly when neither is present, so
  # there is no silent wrong-column selection.
  if (is.null(predColName)) {
    if ("predictor" %in% colnames(df.data)) {
      predColName <- "predictor"
    } else if ("Predictor" %in% colnames(df.data)) {
      predColName <- "Predictor"
    } else {
      stop("Neither 'predictor' nor 'Predictor' column found in `df.data`. ",
           "Pass `predColName` explicitly.", call. = FALSE)
    }
  }

  masterTable <- list()

  df_list <- if (!is.null(group_col)) {
    split(df.data, df.data[[group_col]])
  } else {
    list(all_data = df.data)
  }

  for (group_name in names(df_list)) {
    df.tmp <- df_list[[group_name]]

    area <- unique(df.tmp[[areaColName]])
    predictors <- unique(df.tmp[[predColName]])

    # Initialize with NA (not 0) so untested cells stay NA
    tmpMatrix <- matrix(NA_real_, nrow = length(predictors), ncol = length(area))
    df.out <- as.data.frame(tmpMatrix)
    colnames(df.out) <- area
    df.out[[predColName]] <- predictors

    df.out <- df.out %>% relocate(!!sym(predColName))

    if (!is.null(group_col)) {
      df.out[[group_col]] <- group_name
      df.out <- df.out %>% relocate(!!sym(group_col))
    }

    if (removeIntercept) {
      df.out <- df.out[!grepl("\\(Int", df.out[[predColName]]), ]
    }

    # Populate the table
    for (ch in area) {
      df.tmp.2 <- df.tmp[df.tmp[[areaColName]] == ch, ]
      for (pr in predictors) {
        if (!pr %in% df.out[[predColName]]) next
        if (nrow(df.tmp.2[df.tmp.2[[predColName]] == pr, ]) > 0) {

          sig_val <- df.tmp.2[df.tmp.2[[predColName]] == pr, sigColName]
          dv_val  <- df.tmp.2[df.tmp.2[[predColName]] == pr, dv_col]

          if (length(sig_val) == 0 || is.na(sig_val)) {
            # Preserve NA for untested cells
            df.out[df.out[[predColName]] == pr, ch] <- NA_real_
          } else if (output_mode == "discrete") {
            df.out[df.out[[predColName]] == pr, ch] <- (sig_val < sigThresh) * sign(dv_val)
          } else if (output_mode == "continuous_signed") {
            df.out[df.out[[predColName]] == pr, ch] <- sig_val * sign(dv_val)
          } else if (output_mode == "continuous_unsigned") {
            df.out[df.out[[predColName]] == pr, ch] <- abs(sig_val)
          }
        }
      }
    }

    masterTable[[group_name]] <- df.out
  }

  return(masterTable)
}

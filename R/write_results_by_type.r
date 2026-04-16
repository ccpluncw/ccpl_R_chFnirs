#' @title Write a long-format results data.frame to disk, optionally split by type
#'
#' @description
#' Convenience wrapper around `utils::write.csv()` that encapsulates the
#' common second-level pattern of splitting a long-format results table
#' by its `type` column (e.g., chromophore: `hbo`, `hbr`, `hbdiff`,
#' `cbsi`) and writing one CSV per type into its own subdirectory.
#'
#' When `split_by_type = TRUE` (the default), the function creates
#' `<output_dir>/<tolower(type)>/` if needed and writes the filtered
#' rows to `<output_dir>/<tolower(type)>/<filename>` for every unique
#' value of `type_col`. When `split_by_type = FALSE`, the whole
#' data.frame is written to `<output_dir>/<filename>` in one piece
#' (the pre-split behavior).
#'
#' @param df A long-format data.frame containing a column named by
#'   `type_col` when `split_by_type = TRUE`.
#' @param filename Character scalar, the CSV file name (no directory
#'   component). Used verbatim inside each type subdirectory.
#' @param output_dir Character scalar, the parent directory that holds
#'   the per-type subdirectories (or the target for the combined file
#'   when `split_by_type = FALSE`). Default `"."` (current working
#'   directory).
#' @param split_by_type Logical. Split and write per-type
#'   subdirectories (`TRUE`, default) or write a single combined file
#'   (`FALSE`).
#' @param type_col Character scalar naming the type column used for
#'   splitting. Default `"type"`.
#' @param row.names Passed through to `utils::write.csv()`. Default
#'   `FALSE`.
#' @param ... Additional arguments passed through to
#'   `utils::write.csv()`.
#'
#' @return Invisibly returns `NULL`. Side effect: writes one or more
#'   CSV files to disk.
#'
#' @examples
#' \dontrun{
#' # After run_activation_tests() -> df has a `type` column
#' write_results_by_type(df, filename = "Activation_byGroup.csv")
#' # writes ./hbo/Activation_byGroup.csv, ./hbr/Activation_byGroup.csv, ...
#'
#' # Combined-file fallback
#' write_results_by_type(df, filename = "Activation_byGroup.csv",
#'                       split_by_type = FALSE)
#' }
#'
#' @importFrom utils write.csv
#' @export
write_results_by_type <- function(df,
                                  filename,
                                  output_dir    = ".",
                                  split_by_type = TRUE,
                                  type_col      = "type",
                                  row.names     = FALSE,
                                  ...) {
  if (!is.data.frame(df)) {
    stop("`df` must be a data.frame.", call. = FALSE)
  }
  if (!is.character(filename) || length(filename) != 1) {
    stop("`filename` must be a single character string.", call. = FALSE)
  }

  if (split_by_type) {
    if (!(type_col %in% colnames(df))) {
      stop("`split_by_type = TRUE` requires a `", type_col,
           "` column in `df`.", call. = FALSE)
    }
    for (ty in unique(df[[type_col]])) {
      subdir <- file.path(output_dir, tolower(as.character(ty)))
      if (!dir.exists(subdir)) {
        dir.create(subdir, recursive = TRUE, showWarnings = FALSE)
      }
      df_sub <- df[df[[type_col]] == ty, , drop = FALSE]
      utils::write.csv(df_sub, file.path(subdir, filename),
                       row.names = row.names, ...)
    }
  } else {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }
    utils::write.csv(df, file.path(output_dir, filename),
                     row.names = row.names, ...)
  }

  invisible(NULL)
}

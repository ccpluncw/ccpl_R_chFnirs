#' @title Check whether a results dataframe has any non-NA tested cells
#'
#' @description
#' Gate helper for conditional heatmap rendering. Returns `TRUE` when a
#' results dataframe contains at least one non-NA p-value, and `FALSE`
#' otherwise. Intended to be called before drawing a post-hoc or
#' effect-size heatmap so that empty figures are not produced when an
#' omnibus test was not significant (and therefore all post-hoc cells
#' were gated to NA).
#'
#' The function checks for a `p_bon` column first (used by Bonferroni-
#' corrected post-hocs) and falls back to `p_value` (used by omnibus or
#' unadjusted activation tests). If neither column is present, it
#' returns `FALSE`.
#'
#' @param df A data.frame of statistical results, or `NULL`.
#'
#' @return A single logical value. `TRUE` if `df` is non-null, non-empty,
#'   contains either a `p_bon` or `p_value` column, and that column has at
#'   least one non-NA value. `FALSE` otherwise.
#'
#' @examples
#' has_tested_results(NULL)                                 # FALSE
#' has_tested_results(data.frame(p_value = c(NA, NA)))         # FALSE
#' has_tested_results(data.frame(p_value = c(NA, 0.03)))       # TRUE
#' has_tested_results(data.frame(p_bon = 0.01, p_value = NA))   # TRUE
#'
#' @export
has_tested_results <- function(df) {
  if (is.null(df)) return(FALSE)
  if (nrow(df) == 0) return(FALSE)
  p_col <- if ("p_bon" %in% colnames(df)) {
    "p_bon"
  } else if ("p_value" %in% colnames(df)) {
    "p_value"
  } else {
    return(FALSE)
  }
  any(!is.na(df[[p_col]]))
}

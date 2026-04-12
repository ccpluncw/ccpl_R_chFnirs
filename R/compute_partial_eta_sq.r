#' @title Compute partial eta-squared from an ANOVA table
#'
#' @description
#' Computes partial eta-squared for every row of an ANOVA table produced
#' by [stats::anova()] or `car::Anova()`. Partial eta-squared is
#' defined as \eqn{SS_{effect} / (SS_{effect} + SS_{residual})} and
#' represents the proportion of variance in the outcome attributable to
#' each effect after partialing out the other effects in the model.
#'
#' The function auto-detects whether the sum-of-squares column is named
#' `"Sum Sq"` (base R / car) or `"SS"` (some alternative tables). If
#' neither column is present, or if no `Residuals` row is found, a
#' vector of `NA_real_` values of the correct length is returned so
#' that downstream code can safely append the result to the ANOVA table
#' without length mismatches.
#'
#' @param anova_table An ANOVA table (data.frame or coercible) containing
#'   a sum-of-squares column and a row named `"Residuals"`.
#'
#' @return A numeric vector of partial eta-squared values, one per row of
#'   `anova_table`. The `Residuals` row itself will be 0.5 by
#'   construction and should typically be dropped by the caller.
#'
#' @examples
#' \dontrun{
#' fit <- lm(mpg ~ cyl + hp, data = mtcars)
#' aov_tab <- car::Anova(fit, type = 3)
#' compute_partial_eta_sq(aov_tab)
#' }
#'
#' @export
compute_partial_eta_sq <- function(anova_table) {
  if ("Sum Sq" %in% colnames(anova_table)) {
    ss_col <- "Sum Sq"
  } else if ("SS" %in% colnames(anova_table)) {
    ss_col <- "SS"
  } else {
    return(rep(NA_real_, nrow(anova_table)))
  }
  ss_resid <- anova_table[rownames(anova_table) == "Residuals", ss_col]
  if (length(ss_resid) == 0) return(rep(NA_real_, nrow(anova_table)))
  anova_table[[ss_col]] / (anova_table[[ss_col]] + ss_resid)
}

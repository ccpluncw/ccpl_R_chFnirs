#' @title Compute Cohen's d from a fitted model and an estimate vector
#'
#' @description
#' Computes Cohen's d by dividing a vector of estimates (e.g., `emmeans`
#' marginal means, pairwise contrast estimates, or slope coefficients)
#' by the residual standard deviation of the fitted model from which
#' they were obtained.
#'
#' This is the standardized-mean-difference interpretation of d: the
#' effect is expressed in units of residual SD of the model, which is
#' appropriate when the estimates come from the same model whose error
#' variance reflects the relevant within-cell variability.
#'
#' If the residual SD is zero, missing, or undefined, a vector of
#' `NA_real_` values of the same length as `estimate_col` is returned.
#'
#' @param model A fitted model (e.g., from [stats::lm()] or
#'   `lmerTest::lmer()`) that supports [stats::sigma()].
#' @param estimate_col A numeric vector of effect estimates (e.g.,
#'   marginal means, contrast estimates, or slope coefficients) on the
#'   response scale of `model`.
#'
#' @return A numeric vector of Cohen's d values the same length as
#'   `estimate_col`.
#'
#' @examples
#' fit <- lm(mpg ~ cyl, data = mtcars)
#' compute_cohens_d(fit, coef(fit))
#'
#' @importFrom stats sigma
#' @export
compute_cohens_d <- function(model, estimate_col) {
  resid_sd <- sigma(model)
  if (is.null(resid_sd) || is.na(resid_sd) || resid_sd == 0) {
    return(rep(NA_real_, length(estimate_col)))
  }
  estimate_col / resid_sd
}

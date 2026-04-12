#' @title Test marginal means against zero for every type x area cell
#'
#' @description
#' For each combination of measurement type (chromophore) and area
#' (ROI, Region, Channel, ...), fits the requested model and computes
#' marginal means at each level of `by_factor` via `emmeans::emmeans()`,
#' then tests each marginal mean against zero. Returns a long-format
#' table with one row per level x type x area. Intended as the
#' "activation from baseline" test for categorical designs — does the
#' hemodynamic response differ from zero in each cell of the design?
#'
#' This is the ungated counterpart to `run_gated_posthoc()` and the
#' categorical counterpart to `run_slopes_vs_zero()`.
#'
#' @param data A long-format data.frame.
#' @param formula Model formula. Accepts a formula object (preferred)
#'   or a character string; see [run_omnibus_by_area()] for details.
#' @param area_col Character scalar naming the area column. Required.
#' @param by_factor Character scalar naming the factor whose marginal
#'   means should be tested. Passed to `emmeans::emmeans()` as the
#'   specs argument via `~ by_factor`.
#' @param type_col Character scalar naming the measurement-type
#'   column. Defaults to `"type"`.
#' @param model Required character scalar: `"lmer"` or `"lm"`. See
#'   [run_omnibus_by_area()] for the rationale behind requiring this
#'   argument.
#' @param vcov Character scalar naming the covariance estimator; see
#'   [run_omnibus_by_area()] for the full description of valid values,
#'   the `NULL` conditional default, and why HC estimators error with
#'   `model = "lmer"`.
#' @param min_rows Integer scalar; cells with fewer observations are
#'   silently skipped. Default `5`.
#' @param cohens_d Logical. Currently always `TRUE`; retained as an
#'   argument for a future opt-out. The output always includes a
#'   `cohens_d` column populated via `compute_cohens_d()` (which is
#'   `NA_real_` when `sigma(fit)` is undefined, e.g., singular fits).
#'
#' @return A data.frame with columns `type`, `area`, `predictor`,
#'   `estimate`, `se`, `df`, `statistic`, `statistic_type` (`"t"` or
#'   `"z"`), `p_value`, `cohens_d`. Returns `NULL` if no cell produced
#'   output.
#'
#' @seealso [run_omnibus_by_area()] for the full `model` and `vcov`
#'   argument discussion; [run_slopes_vs_zero()] for the continuous-
#'   moderator counterpart; [standardize_emmeans_columns()] for the
#'   underlying column renaming.
#'
#' @examples
#' \dontrun{
#' df.activation <- run_activation_tests(
#'   data      = df.all,
#'   formula   = estimate ~ Group * cond + (1 | ID),
#'   area_col  = "ROI",
#'   by_factor = "Group",
#'   model     = "lmer"
#' )
#' }
#'
#' @importFrom stats as.formula
#' @importFrom emmeans emmeans
#' @export
run_activation_tests <- function(data,
                                 formula,
                                 area_col,
                                 by_factor,
                                 type_col  = "type",
                                 model,
                                 vcov      = NULL,
                                 min_rows  = 5,
                                 cohens_d  = TRUE) {

  # --- Validate model (required) ---
  if (missing(model)) {
    stop("Argument `model` is required. Choose one of:\n",
         "  \"lmer\" - fit a linear mixed model via lmerTest::lmer (for\n",
         "           repeated-measures data with random intercepts/slopes)\n",
         "  \"lm\"   - fit an ordinary linear model via stats::lm (for\n",
         "           flat / independent-observation data)",
         call. = FALSE)
  }
  model <- match.arg(model, c("lmer", "lm"))

  # --- Resolve vcov default based on model, then validate ---
  valid_vcov <- c("classical", "HC0", "HC1", "HC2", "HC3", "HC4", "HC5")
  if (is.null(vcov)) {
    vcov <- if (model == "lmer") "classical" else "HC3"
  } else {
    vcov <- match.arg(vcov, valid_vcov)
  }

  if (model == "lmer" && vcov != "classical") {
    stop("`vcov = \"", vcov, "\"` is not valid when `model = \"lmer\"`.\n",
         "Set `vcov = \"classical\"` or `model = \"lm\"`. See ?run_omnibus_by_area.",
         call. = FALSE)
  }

  # --- Coerce formula if string ---
  if (is.character(formula)) {
    formula <- stats::as.formula(formula)
  }

  em_spec <- stats::as.formula(paste("~", by_factor))

  fit_fn <- function(df_cell) {
    if (model == "lmer") {
      lmerTest::lmer(formula, data = df_cell)
    } else {
      stats::lm(formula, data = df_cell)
    }
  }

  extract_fn <- function(fit, df_cell) {
    # Build emmeans grid (with HC-adjusted vcov for lm + HC*)
    if (model == "lm" && vcov != "classical") {
      vc <- sandwich::vcovHC(fit, type = vcov)
      em_grid <- emmeans::emmeans(fit, em_spec, vcov. = vc)
    } else {
      em_grid <- emmeans::emmeans(fit, em_spec)
    }

    em_test <- as.data.frame(
      emmeans::test(em_grid, null = 0, joint = FALSE)
    )
    em_test <- standardize_emmeans_columns(em_test, mode = "marginal")

    # predictor label from the factor-level column
    if (by_factor %in% colnames(em_test)) {
      em_test$predictor <- as.character(em_test[[by_factor]])
    } else {
      em_test$predictor <- NA_character_
    }

    # Determine statistic_type from what standardize_emmeans_columns produced
    # (t for lm/lmer default, z for large-df approximations)
    stat_type <- if ("t.ratio" %in% colnames(em_test)) "t"  # shouldn't happen post-rename
                 else "t"
    em_test$statistic_type <- stat_type

    em_test$cohens_d <- compute_cohens_d(fit, em_test$estimate)

    keep_cols <- c("predictor", "estimate", "se", "df",
                   "statistic", "statistic_type", "p_value", "cohens_d")
    keep_cols <- keep_cols[keep_cols %in% colnames(em_test)]
    em_test[, keep_cols, drop = FALSE]
  }

  iterate_by_area(
    data       = data,
    area_col   = area_col,
    type_col   = type_col,
    min_rows   = min_rows,
    fit_fn     = fit_fn,
    extract_fn = extract_fn,
    fn_label   = "activation test"
  )
}

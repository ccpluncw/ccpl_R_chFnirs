#' @title Test continuous-moderator slopes against zero for every cell
#'
#' @description
#' For each type x area cell, fits the requested model and uses
#' `emmeans::emtrends()` to compute the slope of `trend_var` (a
#' continuous moderator such as a hormone level or a personality
#' score), then tests each slope against zero. This is the
#' continuous-moderator analog of [run_activation_tests()]: ungated,
#' one row per slope per cell.
#'
#' When `by_factor` is `NULL`, a single overall slope is computed per
#' cell. When `by_factor` names a categorical factor, slopes are
#' computed within each level of that factor (e.g., one slope per
#' Group), yielding multiple rows per cell.
#'
#' @param data A long-format data.frame.
#' @param formula Model formula. Must include `trend_var` as a
#'   predictor for the slope to be identifiable. Accepts either a
#'   formula object or a character string.
#' @param area_col Character scalar naming the area column. Required.
#' @param trend_var Character scalar naming the continuous moderator
#'   whose slope is being tested (e.g., `"cortisol"`, `"zAveCort"`).
#'   Must appear as a predictor in `formula`.
#' @param by_factor Optional character scalar. If `NULL` (default),
#'   tests a single overall slope per cell. If a factor name, tests a
#'   slope within each level of that factor.
#' @param type_col Character scalar naming the measurement-type column.
#'   Defaults to `"type"`.
#' @param model Required character scalar: `"lmer"` or `"lm"`. See
#'   [run_omnibus_by_area()].
#' @param vcov Character scalar naming the covariance estimator. See
#'   [run_omnibus_by_area()] for valid values and the conditional
#'   default.
#' @param min_rows Integer scalar; default `5`.
#' @param cohens_d Logical; always `TRUE` for now.
#'
#' @return A data.frame with columns `type`, `area`, `predictor`,
#'   `estimate`, `se`, `df`, `statistic`, `statistic_type` (`"t"`),
#'   `p_value`, `cohens_d`. The `predictor` column is labeled
#'   `"<trend_var>_Slope"` for overall slopes and
#'   `"<level>_<trend_var>_Slope"` for slopes within `by_factor`
#'   levels.
#'
#' @seealso [run_omnibus_by_area()] for the full `model`/`vcov`
#'   discussion; [run_activation_tests()] for the categorical
#'   counterpart; [run_gated_slopes()] for the gated pairwise
#'   slope-comparison counterpart.
#'
#' @examples
#' \dontrun{
#' df.slopes <- run_slopes_vs_zero(
#'   data      = df.all,
#'   formula   = estimate ~ Group * cond * cortisol + (1 | ID),
#'   area_col  = "ROI",
#'   trend_var = "cortisol",
#'   by_factor = "Group",
#'   model     = "lmer"
#' )
#' }
#'
#' @importFrom stats as.formula
#' @importFrom emmeans emtrends
#' @export
run_slopes_vs_zero <- function(data,
                               formula,
                               area_col,
                               trend_var,
                               by_factor = NULL,
                               type_col  = "type",
                               model,
                               vcov      = NULL,
                               min_rows  = 5,
                               cohens_d  = TRUE) {

  # --- Validate model ---
  if (missing(model)) {
    stop("Argument `model` is required. Choose one of:\n",
         "  \"lmer\" - fit a linear mixed model via lmerTest::lmer (for\n",
         "           repeated-measures data with random intercepts/slopes)\n",
         "  \"lm\"   - fit an ordinary linear model via stats::lm (for\n",
         "           flat / independent-observation data)",
         call. = FALSE)
  }
  model <- match.arg(model, c("lmer", "lm"))

  # --- Resolve vcov ---
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

  if (is.character(formula)) formula <- stats::as.formula(formula)

  em_spec <- if (is.null(by_factor)) {
    stats::as.formula("~ 1")
  } else {
    stats::as.formula(paste("~", by_factor))
  }

  fit_fn <- function(df_cell) {
    if (model == "lmer") {
      lmerTest::lmer(formula, data = df_cell)
    } else {
      stats::lm(formula, data = df_cell)
    }
  }

  extract_fn <- function(fit, df_cell) {
    if (model == "lm" && vcov != "classical") {
      vc  <- sandwich::vcovHC(fit, type = vcov)
      emt <- emmeans::emtrends(fit, em_spec, var = trend_var, vcov. = vc)
    } else {
      emt <- emmeans::emtrends(fit, em_spec, var = trend_var)
    }

    em_df <- as.data.frame(
      summary(emt, infer = c(TRUE, TRUE), null = 0)
    )
    em_df <- standardize_emmeans_columns(em_df, mode = "trend",
                                         trend_var = trend_var)

    # Build predictor label
    if (!is.null(by_factor) && by_factor %in% colnames(em_df)) {
      em_df$predictor <- paste0(em_df[[by_factor]], "_", trend_var, "_Slope")
    } else {
      em_df$predictor <- paste0(trend_var, "_Slope")
    }

    em_df$statistic_type <- "t"
    em_df$cohens_d <- compute_cohens_d(fit, em_df$estimate)

    keep_cols <- c("predictor", "estimate", "se", "df",
                   "statistic", "statistic_type", "p_value", "cohens_d")
    keep_cols <- keep_cols[keep_cols %in% colnames(em_df)]
    em_df[, keep_cols, drop = FALSE]
  }

  iterate_by_area(
    data       = data,
    area_col   = area_col,
    type_col   = type_col,
    min_rows   = min_rows,
    fit_fn     = fit_fn,
    extract_fn = extract_fn,
    fn_label   = "slopes-vs-zero test"
  )
}

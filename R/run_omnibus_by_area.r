#' @title Run an area-wise omnibus F-test for every type x area cell
#'
#' @description
#' Fits a linear or linear mixed model to every unique combination of
#' measurement type (chromophore) and area (ROI, Region, Channel, ...)
#' in the input data, runs a Type III ANOVA on each fit, and returns a
#' long-format table with one row per effect x type x area combination.
#' This replaces the three near-identical `run_omnibus_by_area()` /
#' `run_omnibus_by_region()` helpers in the second-level analysis
#' scripts with a single parameterized function that supports both
#' `stats::lm()` and `lmerTest::lmer()` fits, optional
#' heteroscedasticity-consistent covariance, and automatic
#' partial eta-squared computation.
#'
#' @param data A long-format data.frame with one row per subject x
#'   type x area x condition observation.
#' @param formula Model formula. Accepts either a native R formula
#'   (preferred, e.g., `estimate ~ Group * cond + (1 | ID)`) or a
#'   character string for backwards compatibility with scripts that
#'   build formulas programmatically via `paste0()` / `sprintf()`.
#'   Strings are silently coerced via `stats::as.formula()`.
#' @param area_col Character scalar naming the area column in `data`
#'   (e.g., `"ROI"`, `"Region"`, `"Channel"`). Required.
#' @param type_col Character scalar naming the measurement-type column
#'   in `data` (chromophore labels). Defaults to `"type"`.
#' @param model Required character scalar, one of `"lmer"` or `"lm"`.
#'   There is no default — model choice is a statistical decision that
#'   must be made explicitly at every call. Omitting this argument
#'   raises an informative error. `"lmer"` dispatches to
#'   `lmerTest::lmer()`; `"lm"` dispatches to `stats::lm()`.
#' @param vcov Character scalar naming the covariance estimator. One of
#'   `"classical"`, `"HC0"`, `"HC1"`, `"HC2"`, `"HC3"`, `"HC4"`, or
#'   `"HC5"`. Defaults to `NULL`, which resolves to `"classical"` when
#'   `model = "lmer"` and to `"HC3"` when `model = "lm"`. The
#'   combination `model = "lmer"` with any HC value raises an error;
#'   see \strong{Details}.
#' @param min_rows Integer scalar. Cells with fewer observations than
#'   this are silently skipped. Default `5`.
#' @param effect_size Logical. Currently always `TRUE`; retained as an
#'   argument for a future opt-out but the output schema always
#'   includes a `partial_eta_sq` column. When `model = "lmer"` the
#'   column is populated with `NA_real_` because mixed-model ANOVA
#'   tables do not expose sums of squares.
#' @param refit_singular Logical. Only consulted when `model = "lmer"`.
#'   When `TRUE` (the default), any per-cell lmer fit that comes back
#'   singular (random-effect variance collapsed to the zero boundary,
#'   per `lme4::isSingular()`) is refit with a weak variance prior via
#'   `blme::blmer()` and coerced back through `lmerTest` so the ANOVA
#'   still reports Satterthwaite denominator df. This lifts the
#'   variance off the boundary and restores subject-level df for the
#'   between-subjects effects, which a singular fit otherwise inflates
#'   to the OLS residual (observation-level) df. Set to `FALSE` to keep
#'   the raw singular fit unchanged (the pre-remedy behaviour). The
#'   `singular` output column records, per cell, whether the original
#'   lmer fit was singular regardless of this setting. Requires the
#'   `blme` package (Suggests) when a refit is actually triggered.
#'
#' @return A data.frame with one row per effect x type x area, columns
#'   in order:
#'   \itemize{
#'     \item `type` — value from `type_col`
#'     \item `area` — value from `area_col`
#'     \item `effect` — effect name from the ANOVA table
#'     \item `statistic` — the F statistic
#'     \item `statistic_type` — always `"F"` for this function
#'     \item `df_num` — numerator degrees of freedom
#'     \item `df_den` — denominator degrees of freedom, or `NA` when
#'       the covariance adjustment does not report one
#'     \item `p_value` — raw omnibus p-value. Call `getBHvalues()` or
#'       use `run_omnibus_by_area_bh()` to add a `p_bh` column.
#'     \item `partial_eta_sq` — partial eta-squared from classical
#'       sums of squares. `NA_real_` when `model = "lmer"`.
#'     \item `singular` — logical; `TRUE` when the original lmer fit for
#'       that cell was singular (and was therefore refit via blme when
#'       `refit_singular = TRUE`). `NA` when `model = "lm"`.
#'   }
#'   Returns `NULL` if no cell produced output.
#'
#' @details
#' \strong{On the choice of HC estimator.} HC0–HC5 are a family of
#' heteroscedasticity-consistent (sandwich) covariance estimators that
#' differ in how residuals are reweighted when building the robust
#' variance. HC0 is the original White (1980) estimator; HC1 adds a
#' small-sample degrees-of-freedom correction; HC2 divides each
#' squared residual by \eqn{1 - h_{ii}}, where \eqn{h_{ii}} is the
#' leverage; HC3 divides by \eqn{(1 - h_{ii})^2} and is approximately
#' a jackknife estimator. HC4 and HC5 apply progressively stronger
#' corrections for high-leverage observations. Long and Ervin (2000)
#' showed that HC3 has the best small-sample properties among HC0–HC3
#' and is the recommended default for \eqn{N < 250}, which is why it
#' is the resolved default when `model = "lm"` and `vcov = NULL`.
#'
#' \strong{Why HC estimators do not apply to lmer models.}
#' Heteroscedasticity-consistent estimators are defined against the
#' OLS hat matrix and residual structure. Linear mixed models fit by
#' REML/ML do not have a hat matrix in the same sense — their
#' covariance arises from the estimation of variance components, not
#' from reweighting residuals — so HC0–HC5 are not mathematically
#' defined for `lmer` fits. Requesting `vcov = "HC3"` with
#' `model = "lmer"` therefore errors rather than silently fall back.
#' If you need robust inference for a mixed model, consider
#' cluster-robust sandwich estimators via `clubSandwich` or a
#' resampling procedure (bootstrap); neither is currently supported
#' by this package.
#'
#' \strong{Partial eta-squared with HC covariance.} When `model = "lm"`
#' and `vcov` is an HC estimator, the function runs
#' `car::Anova(fit, type = 3)` twice: once without the HC adjustment
#' to recover classical sums of squares for the `partial_eta_sq`
#' computation, and once with `white.adjust` set to the requested HC
#' type to produce the F statistics and p-values. This matches the
#' approach used in the HCvsNC script and ensures that effect-size
#' reporting is available even when inference uses robust covariance.
#'
#' \strong{References}
#' \itemize{
#'   \item White, H. (1980). A heteroskedasticity-consistent covariance
#'     matrix estimator and a direct test for heteroskedasticity.
#'     \emph{Econometrica}, 48(4), 817–838.
#'   \item MacKinnon, J. G., & White, H. (1985). Some
#'     heteroskedasticity-consistent covariance matrix estimators with
#'     improved finite sample properties. \emph{Journal of
#'     Econometrics}, 29(3), 305–325.
#'   \item Long, J. S., & Ervin, L. H. (2000). Using heteroscedasticity
#'     consistent standard errors in the linear regression model.
#'     \emph{The American Statistician}, 54(3), 217–224.
#' }
#'
#' @examples
#' \dontrun{
#' # Mixed model (repeated measures within subject)
#' df.omnibus <- run_omnibus_by_area(
#'   data     = df.all,
#'   formula  = estimate ~ Group * cond + (1 | ID),
#'   area_col = "ROI",
#'   model    = "lmer"
#' )
#'
#' # Linear model with HC3 robust covariance (flat / independent data)
#' df.omnibus <- run_omnibus_by_area(
#'   data     = df.flat,
#'   formula  = estimate ~ Group * sex,
#'   area_col = "Region",
#'   model    = "lm"
#'   # vcov defaults to "HC3" when model = "lm"
#' )
#' }
#'
#' @importFrom stats as.formula anova lm sigma
#' @importFrom lme4 isSingular
#' @importFrom lmerTest as_lmerModLmerTest
#' @export
run_omnibus_by_area <- function(data,
                                formula,
                                area_col,
                                type_col    = "type",
                                model,
                                vcov        = NULL,
                                min_rows    = 5,
                                effect_size = TRUE,
                                refit_singular = TRUE) {

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
         "Heteroscedasticity-consistent covariance estimators (HC0-HC5) are\n",
         "defined only for OLS residuals and are not applicable to linear\n",
         "mixed models.\n\n",
         "Set `vcov = \"classical\"` to use lmer's native covariance, or set\n",
         "`model = \"lm\"` to use the ", vcov, " estimator.",
         call. = FALSE)
  }

  # --- Coerce formula if the user passed a string ---
  if (is.character(formula)) {
    formula <- stats::as.formula(formula)
  }

  # --- Build fit closure ---
  # Returns list(fit, singular) so the extractor can both run the ANOVA
  # and stamp whether the *original* lmer fit was singular. When an lmer
  # fit is singular and refit_singular = TRUE, the same model is refit
  # with blme's weak default variance prior and coerced back through
  # lmerTest, which lifts the random-effect variance off the zero
  # boundary and restores subject-level Satterthwaite df for the
  # between-subjects effects (Chung et al. 2013).
  fit_fn <- function(df_cell) {
    if (model == "lmer") {
      fit  <- lmerTest::lmer(formula, data = df_cell)
      sing <- lme4::isSingular(fit)
      if (sing && refit_singular) {
        if (!requireNamespace("blme", quietly = TRUE)) {
          stop("`refit_singular = TRUE` needs the 'blme' package for the\n",
               "singular-fit remedy. Install it with install.packages(\"blme\"),\n",
               "or set `refit_singular = FALSE` to keep the singular lmer fit.",
               call. = FALSE)
        }
        fit <- lmerTest::as_lmerModLmerTest(
          blme::blmer(formula, data = df_cell)
        )
      }
      list(fit = fit, singular = sing)
    } else {
      list(fit = stats::lm(formula, data = df_cell), singular = NA)
    }
  }

  # --- Build extract closure ---
  extract_fn <- function(fit, df_cell) {
    # fit_fn hands us list(fit, singular); unpack before use.
    singular_flag <- fit$singular
    fit           <- fit$fit
    if (model == "lmer") {
      aov_tab <- as.data.frame(stats::anova(fit))
      partial_eta <- rep(NA_real_, nrow(aov_tab))
    } else {
      # Classical anova for partial eta-squared (needs Sum Sq column)
      aov_classical <- as.data.frame(car::Anova(fit, type = 3))
      partial_eta_all <- compute_partial_eta_sq(aov_classical)
      classical_keep <- !rownames(aov_classical) %in% c("Residuals", "(Intercept)")
      partial_eta <- partial_eta_all[classical_keep]

      # F-statistic table: classical or HC-adjusted
      if (vcov == "classical") {
        aov_tab <- aov_classical
      } else {
        aov_tab <- as.data.frame(
          car::Anova(fit, type = 3, white.adjust = tolower(vcov))
        )
      }
    }

    # Strip (Intercept) and Residuals rows
    keep_mask <- !rownames(aov_tab) %in% c("Residuals", "(Intercept)")
    effect_names <- rownames(aov_tab)[keep_mask]
    aov_tab <- aov_tab[keep_mask, , drop = FALSE]

    if (model == "lmer") {
      partial_eta <- partial_eta[keep_mask]
    }

    if (nrow(aov_tab) == 0) return(NULL)

    # Resolve column names across lmer / classical-lm / HC-lm variants
    # lmer:          "NumDF", "DenDF", "F value", "Pr(>F)"
    # classical lm:  "Df",             "F value", "Pr(>F)"
    # HC-lm:         "Df",             "F",       "Pr(>F)"
    cn <- colnames(aov_tab)

    statistic_col <- if ("F value" %in% cn) "F value"
                     else if ("F" %in% cn) "F"
                     else NA_character_

    p_col <- if ("Pr(>F)" %in% cn) "Pr(>F)" else NA_character_

    if ("NumDF" %in% cn) {
      df_num <- aov_tab[["NumDF"]]
      df_den <- aov_tab[["DenDF"]]
    } else if ("Df" %in% cn) {
      df_num <- aov_tab[["Df"]]
      df_den <- rep(NA_real_, nrow(aov_tab))
    } else {
      df_num <- rep(NA_real_, nrow(aov_tab))
      df_den <- rep(NA_real_, nrow(aov_tab))
    }

    data.frame(
      effect         = effect_names,
      statistic      = if (is.na(statistic_col)) NA_real_ else aov_tab[[statistic_col]],
      statistic_type = "F",
      df_num         = as.numeric(df_num),
      df_den         = as.numeric(df_den),
      p_value        = if (is.na(p_col)) NA_real_ else aov_tab[[p_col]],
      partial_eta_sq = as.numeric(partial_eta),
      singular       = singular_flag,
      stringsAsFactors = FALSE
    )
  }

  iterate_by_area(
    data       = data,
    area_col   = area_col,
    type_col   = type_col,
    min_rows   = min_rows,
    fit_fn     = fit_fn,
    extract_fn = extract_fn,
    fn_label   = "omnibus test"
  )
}

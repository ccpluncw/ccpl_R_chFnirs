#' @title Run pairwise categorical post-hoc contrasts gated on an omnibus
#'
#' @description
#' For each type x area cell, checks whether a named omnibus effect
#' cleared the significance gate (BH-adjusted p-value below
#' `p_threshold`) and, if so, runs pairwise Bonferroni-adjusted
#' contrasts over `contrast_factor` via `emmeans::contrast()`. Cells
#' that fail the gate optionally receive NA-filled placeholder rows so
#' that downstream heatmap functions render them as grey tiles
#' (preserving the predictor x area grid).
#'
#' This is the unified replacement for the three near-duplicate
#' `run_gated_posthoc()` functions in the second-level analysis
#' scripts. It accepts any `model` / `vcov` combination supported by
#' the rest of the Tier 2 family and applies the same input / output
#' schema conventions.
#'
#' @param omnibus_df A data.frame returned by [run_omnibus_by_area_bh()]
#'   or the equivalent sequence of [run_omnibus_by_area()] plus
#'   [getBHvalues()]. Must contain a column named by `p_col`. If the
#'   column is missing, the function errors with a message naming the
#'   missing column.
#' @param effect_name Character scalar naming the single omnibus effect
#'   to gate on (e.g., `"Group"`, `"cond"`, `"Group:cond"`). Rows in
#'   `omnibus_df` whose `effect` column equals this value supply the
#'   per-cell gate p-values.
#' @param data The long-format data.frame used to fit the omnibus.
#' @param formula Model formula used to fit each cell. Must match the
#'   formula used for the omnibus so that the contrast grid is
#'   well-defined. Accepts either a formula object or a character
#'   string.
#' @param area_col Character scalar naming the area column. Required.
#' @param contrast_factor Character scalar naming the factor to
#'   contrast over. For a simple `Group` main-effect post-hoc this is
#'   `"Group"`; for a `Group:cond` interaction post-hoc decomposed
#'   within each condition, pass `contrast_factor = "Group"` and
#'   `by_var = "cond"`.
#' @param by_var Optional character scalar naming a second factor to
#'   nest the pairwise contrasts within (passed to `emmeans::contrast()`
#'   as its `by=` argument). Default `NULL` runs a single pairwise
#'   family without nesting.
#' @param type_col Character scalar naming the measurement-type column.
#'   Defaults to `"type"`.
#' @param model Required character scalar: `"lmer"` or `"lm"`. See
#'   [run_omnibus_by_area()].
#' @param vcov Character scalar naming the covariance estimator. See
#'   [run_omnibus_by_area()] for valid values and the conditional
#'   default.
#' @param p_col Character scalar naming the column in `omnibus_df` to
#'   gate on. Default `"p_bh"`. If the column is absent, the function
#'   errors.
#' @param p_threshold Numeric scalar. Cells with `p_col < p_threshold`
#'   (and non-NA) are considered significant. Default `0.05`.
#' @param fill_empty Logical. If `TRUE` (default), cells that fail the
#'   omnibus gate still produce rows in the output with `NA_real_`
#'   values in the numeric columns, so heatmap functions can render
#'   them as grey. If `FALSE`, non-significant cells are dropped
#'   entirely from the output.
#' @param min_rows Integer scalar; default `5`.
#' @param cohens_d Logical; always `TRUE` for now. Output schema
#'   always includes a `cohens_d` column.
#' @param adjust Character scalar passed to `emmeans::contrast()` as
#'   its `adjust=` argument. Default `"bonferroni"`. Other valid
#'   values include `"tukey"`, `"sidak"`, `"holm"`, `"none"`, etc.;
#'   whichever is used, the resulting p-value column is named `p_bon`
#'   in the output regardless (column name reflects the "post-hoc
#'   pairwise p-value" slot, not the specific method).
#'
#' @return A data.frame with columns `type`, `area`, `contrast`,
#'   `predictor`, `estimate`, `se`, `df`, `statistic`, `statistic_type`
#'   (`"t"`), `p_bon`, `cohens_d`. Non-significant cells have NA
#'   numeric columns when `fill_empty = TRUE`; they are absent
#'   otherwise.
#'
#' @seealso [run_omnibus_by_area()] for the full `model`/`vcov`
#'   discussion; [run_omnibus_by_area_bh()] for the common upstream
#'   pattern; [run_gated_slopes()] for the continuous-moderator
#'   counterpart.
#'
#' @examples
#' \dontrun{
#' df.omnibus <- run_omnibus_by_area_bh(
#'   data = df.all, formula = estimate ~ Group * cond + (1 | ID),
#'   area_col = "ROI", model = "lmer"
#' )
#' df.ph.group <- run_gated_posthoc(
#'   omnibus_df      = df.omnibus,
#'   effect_name     = "Group",
#'   data            = df.all,
#'   formula         = estimate ~ Group * cond + (1 | ID),
#'   area_col        = "ROI",
#'   contrast_factor = "Group",
#'   model           = "lmer"
#' )
#' }
#'
#' @importFrom stats as.formula
#' @importFrom emmeans emmeans contrast
#' @export
run_gated_posthoc <- function(omnibus_df,
                              effect_name,
                              data,
                              formula,
                              area_col,
                              contrast_factor,
                              by_var       = NULL,
                              type_col     = "type",
                              model,
                              vcov         = NULL,
                              p_col        = "p_bh",
                              p_threshold  = 0.05,
                              fill_empty   = TRUE,
                              min_rows     = 5,
                              cohens_d     = TRUE,
                              adjust       = "bonferroni") {

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

  # --- Validate omnibus_df has the gate column ---
  if (!is.data.frame(omnibus_df)) {
    stop("`omnibus_df` must be a data.frame.", call. = FALSE)
  }
  if (!(p_col %in% colnames(omnibus_df))) {
    stop("`omnibus_df` does not contain the gate column `", p_col, "`.\n",
         "Did you forget to BH-correct the omnibus table before calling this\n",
         "function? Use `run_omnibus_by_area_bh()` or call `getBHvalues()` on\n",
         "the result of `run_omnibus_by_area()` first.",
         call. = FALSE)
  }
  if (!("effect" %in% colnames(omnibus_df))) {
    stop("`omnibus_df` does not contain an `effect` column. It should be the ",
         "output of `run_omnibus_by_area()` or `run_omnibus_by_area_bh()`.",
         call. = FALSE)
  }

  # --- Coerce formula ---
  if (is.character(formula)) formula <- stats::as.formula(formula)

  # --- Pre-filter omnibus to the effect of interest for faster lookup ---
  effect_rows <- omnibus_df[omnibus_df$effect == effect_name, , drop = FALSE]

  # --- emmeans spec formula ---
  em_spec <- if (is.null(by_var)) {
    stats::as.formula(paste("~", contrast_factor))
  } else {
    stats::as.formula(paste("~", contrast_factor, "*", by_var))
  }

  fit_fn <- function(df_cell) {
    if (model == "lmer") {
      lmerTest::lmer(formula, data = df_cell)
    } else {
      stats::lm(formula, data = df_cell)
    }
  }

  extract_fn <- function(fit, df_cell) {
    # Find this cell's row in effect_rows — iterate_by_area already stamps
    # type/area onto the output, so we need to look up by (type, area) here.
    # type_col / area_col values come from df_cell's first row.
    ty <- df_cell[[type_col]][1]
    ar <- df_cell[[area_col]][1]

    row_match <- effect_rows[effect_rows$type == ty & effect_rows$area == ar, ,
                             drop = FALSE]
    if (nrow(row_match) == 0) return(NULL)

    gate_p <- row_match[[p_col]][1]
    is_significant <- isTRUE(!is.na(gate_p) && gate_p < p_threshold)

    if (!is_significant && !fill_empty) return(NULL)

    # --- Build emmeans grid ---
    if (model == "lm" && vcov != "classical") {
      vc <- sandwich::vcovHC(fit, type = vcov)
      em_grid <- emmeans::emmeans(fit, em_spec, vcov. = vc)
    } else {
      em_grid <- emmeans::emmeans(fit, em_spec)
    }

    em_contrasts <- if (is.null(by_var)) {
      emmeans::contrast(em_grid, "pairwise", adjust = adjust)
    } else {
      emmeans::contrast(em_grid, "pairwise", by = by_var, adjust = adjust)
    }

    em_df <- as.data.frame(em_contrasts)
    em_df <- standardize_emmeans_columns(em_df, mode = "contrast")

    # Build predictor label: nested "byvar_contrast" when nested, else contrast alone
    if (!is.null(by_var) && by_var %in% colnames(em_df)) {
      em_df$predictor <- paste(em_df[[by_var]], em_df$contrast, sep = "_")
    } else {
      em_df$predictor <- em_df$contrast
    }

    em_df$statistic_type <- "t"
    em_df$cohens_d <- compute_cohens_d(fit, em_df$estimate)

    # If not significant, NA-fill numeric columns (keep structural labels)
    if (!is_significant) {
      em_df$estimate <- NA_real_
      em_df$se       <- NA_real_
      em_df$df       <- NA_real_
      em_df$statistic <- NA_real_
      em_df$p_bon    <- NA_real_
      em_df$cohens_d <- NA_real_
    }

    keep_cols <- c("contrast", "predictor", "estimate", "se", "df",
                   "statistic", "statistic_type", "p_bon", "cohens_d")
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
    fn_label   = "gated post-hoc"
  )
}

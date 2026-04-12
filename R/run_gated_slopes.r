#' @title Run gated pairwise slope comparisons for a continuous moderator
#'
#' @description
#' For each type x area cell, checks the gated significance of every
#' omnibus effect that involves a continuous moderator `trend_var` and,
#' for those that pass the gate, runs pairwise Bonferroni-adjusted
#' slope comparisons via `emmeans::emtrends()` + `emmeans::contrast()`.
#' This is the continuous-moderator counterpart of
#' [run_gated_posthoc()] — it tests "are slopes different between
#' groups / conditions?" rather than "are marginal means different?"
#'
#' The function automatically decomposes interaction effects into the
#' appropriate pairwise-slice specifications based on which of the
#' factors named in `categorical_factors` appear in each effect name.
#' A main-effect-only omnibus row (e.g., just `trend_var`) becomes a
#' single "overall slope vs zero" test. A two-way interaction with one
#' categorical factor (e.g., `Group:trend_var`) becomes a pairwise
#' slope comparison across that factor. A three-way interaction with
#' two categorical factors (e.g., `Group:cond:trend_var`) is
#' decomposed into two conditional slicings (`~ Group | cond` and
#' `~ cond | Group`), producing two groups of results tagged with
#' `posthoc_type`.
#'
#' The decomposition is driven entirely by the user-supplied
#' `categorical_factors` vector, so no factor names are hardcoded in
#' the package. See \strong{Details} for the decomposition rules.
#'
#' @param omnibus_df A data.frame from [run_omnibus_by_area_bh()] or
#'   the equivalent sequence. The function internally filters this to
#'   rows whose `effect` column matches `effect_pattern`, so passing
#'   either a pre-filtered slope-effect table (e.g., after
#'   `omnibus_df[grepl(trend_var, omnibus_df$effect), ]`) or the full
#'   omnibus table both work. Must contain the column named by
#'   `p_col` or the function errors.
#' @param data The long-format data.frame used to fit the omnibus.
#' @param formula Model formula used to fit each cell. Must include
#'   `trend_var` as a predictor. Accepts either a formula object or a
#'   character string.
#' @param area_col Character scalar naming the area column. Required.
#' @param trend_var Character scalar naming the continuous moderator
#'   whose slope is being compared.
#' @param categorical_factors Character vector naming the categorical
#'   factors in the design (e.g., `c("Group", "cond")`). The function
#'   uses these to detect which factors appear in each omnibus effect
#'   name and to build the appropriate pairwise decomposition. Order
#'   is preserved in the decomposition labels (the first-listed factor
#'   appears first in the `<A>By<B>` label).
#' @param effect_pattern Character scalar used in `grepl()` to select
#'   which rows of `omnibus_df` correspond to slope-containing
#'   effects. Defaults to `trend_var`, which is usually correct —
#'   effect names that contain the moderator's name are the ones you
#'   want to gate on.
#' @param type_col Character scalar naming the measurement-type
#'   column. Defaults to `"type"`.
#' @param model Required character scalar: `"lmer"` or `"lm"`. See
#'   [run_omnibus_by_area()].
#' @param vcov Character scalar; see [run_omnibus_by_area()].
#' @param p_col Character scalar naming the column in `omnibus_df` to
#'   gate on. Default `"p_bh"`.
#' @param p_threshold Numeric scalar; default `0.05`.
#' @param fill_empty Logical; if `TRUE` (default) non-significant
#'   cells receive NA-filled placeholder rows.
#' @param min_rows Integer scalar; default `5`.
#' @param cohens_d Logical; always `TRUE` for now.
#' @param adjust Character scalar passed to `emmeans::contrast()`'s
#'   `adjust=` argument. Default `"bonferroni"`.
#' @param suffix Character scalar appended to pairwise-contrast
#'   predictor labels to distinguish slope comparisons from mean
#'   comparisons in combined output. Default `"SlopeDiff"`.
#'
#' @return A data.frame with columns `type`, `area`, `omnibus_effect`,
#'   `posthoc_type`, `contrast`, `predictor`, `estimate`, `se`, `df`,
#'   `statistic`, `statistic_type` (`"t"`), `p_bon`, `cohens_d`.
#'   When an omnibus effect has no matching rows for a cell, that
#'   cell contributes nothing for that effect.
#'
#' @details
#' \strong{Decomposition rules (Design Y').} For each omnibus effect
#' matching `effect_pattern`, the function inspects the effect name
#' with `grepl()` against each entry of `categorical_factors` to
#' determine which categorical factors are involved. Based on the
#' count of factors present:
#' \itemize{
#'   \item \strong{Zero factors} (e.g., a main effect of just
#'     `trend_var`): a single spec is produced that tests the overall
#'     slope against zero via `emmeans::emtrends(fit, ~ 1, var =
#'     trend_var)` and `summary(..., null = 0)`. The `posthoc_type`
#'     column for this row is `"Overall"`.
#'   \item \strong{One factor} (e.g., `Group:trend_var`): a single
#'     spec computes slopes at each level of that factor via
#'     `emmeans::emtrends(fit, ~ Group, var = trend_var)` and runs
#'     `contrast(..., "pairwise")` over the levels. The
#'     `posthoc_type` column is the factor name.
#'   \item \strong{Two factors} (e.g., `Group:cond:trend_var`): two
#'     specs are produced, one for each conditional slice — slopes
#'     across the first factor within each level of the second, and
#'     slopes across the second factor within each level of the
#'     first. Both are returned, distinguished by their
#'     `posthoc_type` labels (`"<A>By<B>"` and `"<B>By<A>"`).
#'   \item \strong{Three or more factors:} currently handled by
#'     producing a single spec that conditions on all factors except
#'     the first, with a warning. Support for full 3+-way
#'     decomposition is deferred until a concrete study needs it.
#' }
#' The decomposition happens inside an internal helper
#' `get_posthoc_specs()` that is closed over the user-supplied
#' `categorical_factors` vector — there are no hardcoded factor
#' names anywhere in this function or in the package.
#'
#' \strong{Script-side filter pattern.} In the two-model workflow
#' used throughout the second-level scripts, the user typically
#' calls `run_omnibus_by_area_bh()` on the full slope model,
#' filters the omnibus table via
#' `df.omnibus.h <- df.omnibus.h.full[grepl(trend_var,
#' df.omnibus.h.full$effect), ]` to drop the categorical-only
#' rows that were already tested in the categorical phase, and
#' then passes the filtered table to this function. That filter
#' step is intentionally kept at the script level (rather than
#' inside this function) because it makes the two-model structure
#' visible where the analyst reads the code. You can pass either
#' the pre-filtered or the full omnibus table here — the internal
#' `effect_pattern` grep handles both.
#'
#' @seealso [run_omnibus_by_area()] for the `model`/`vcov`
#'   discussion; [run_slopes_vs_zero()] for the ungated counterpart;
#'   [run_gated_posthoc()] for the categorical counterpart.
#'
#' @examples
#' \dontrun{
#' # Phase 2 slopes analysis (two-model workflow)
#' df.omnibus.h.full <- run_omnibus_by_area_bh(
#'   data     = df.all,
#'   formula  = estimate ~ Group * cond * cortisol + (1 | ID),
#'   area_col = "ROI",
#'   model    = "lmer"
#' )
#' # Script-side filter: drop categorical-only rows
#' df.omnibus.h <- df.omnibus.h.full[
#'   grepl("cortisol", df.omnibus.h.full$effect),
#' ]
#' df.h.posthoc <- run_gated_slopes(
#'   omnibus_df          = df.omnibus.h,
#'   data                = df.all,
#'   formula             = estimate ~ Group * cond * cortisol + (1 | ID),
#'   area_col            = "ROI",
#'   trend_var           = "cortisol",
#'   categorical_factors = c("Group", "cond"),
#'   model               = "lmer"
#' )
#' }
#'
#' @importFrom stats as.formula
#' @importFrom emmeans emmeans emtrends contrast
#' @export
run_gated_slopes <- function(omnibus_df,
                             data,
                             formula,
                             area_col,
                             trend_var,
                             categorical_factors,
                             effect_pattern = trend_var,
                             type_col       = "type",
                             model,
                             vcov           = NULL,
                             p_col          = "p_bh",
                             p_threshold    = 0.05,
                             fill_empty     = TRUE,
                             min_rows       = 5,
                             cohens_d       = TRUE,
                             adjust         = "bonferroni",
                             suffix         = "SlopeDiff") {

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

  # --- Validate omnibus_df ---
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
    stop("`omnibus_df` does not contain an `effect` column.", call. = FALSE)
  }
  if (!is.character(categorical_factors) || length(categorical_factors) == 0) {
    stop("`categorical_factors` must be a non-empty character vector naming ",
         "the categorical factors in the design (e.g., c(\"Group\", \"cond\")).",
         call. = FALSE)
  }

  if (is.character(formula)) formula <- stats::as.formula(formula)

  # --- Identify matching effects in the omnibus table ---
  matching_effects <- unique(
    omnibus_df$effect[grepl(effect_pattern, omnibus_df$effect)]
  )
  if (length(matching_effects) == 0) return(NULL)

  # --- Decomposition helper (closes over trend_var + categorical_factors) ---
  # Label-prefix strings in `posthoc_type` use Title Case — the first
  # letter of each factor name is capitalized, even when the factor
  # itself is lowercase in the data. This produces readable composite
  # labels like "GroupByCond" (rather than "GroupBycond") when factor
  # names mix capital and lowercase, and it matches the labeling
  # convention used in the second-level analysis scripts' inline
  # `run_gated_emtrends()` helpers, so migration produces identical
  # `posthoc_type` columns.
  cap_first <- function(s) {
    if (is.null(s) || length(s) == 0 || nchar(s) == 0) return(s)
    paste0(toupper(substr(s, 1, 1)), substr(s, 2, nchar(s)))
  }

  get_posthoc_specs <- function(effect_name) {
    present <- categorical_factors[
      vapply(categorical_factors,
             function(f) grepl(f, effect_name),
             logical(1))
    ]
    n_present <- length(present)

    if (n_present == 0) {
      list(list(
        em_formula   = "~ 1",
        do_contrast  = FALSE,
        contrast_by  = NULL,
        label_prefix = "Overall"
      ))
    } else if (n_present == 1) {
      f <- present[1]
      list(list(
        em_formula   = paste("~", f),
        do_contrast  = TRUE,
        contrast_by  = NULL,
        label_prefix = cap_first(f)
      ))
    } else if (n_present == 2) {
      f1 <- present[1]
      f2 <- present[2]
      list(
        list(em_formula   = paste("~", f1, "|", f2),
             do_contrast  = TRUE,
             contrast_by  = f2,
             label_prefix = paste0(cap_first(f1), "By", cap_first(f2))),
        list(em_formula   = paste("~", f2, "|", f1),
             do_contrast  = TRUE,
             contrast_by  = f1,
             label_prefix = paste0(cap_first(f2), "By", cap_first(f1)))
      )
    } else {
      warning("More than 2 categorical factors detected in effect `",
              effect_name, "`. Full 3+-way decomposition is not yet ",
              "implemented; producing a single conditional spec.",
              call. = FALSE)
      f1 <- present[1]
      rest <- paste(present[-1], collapse = "*")
      list(list(
        em_formula   = paste("~", f1, "|", rest),
        do_contrast  = TRUE,
        contrast_by  = rest,
        label_prefix = paste(vapply(present, cap_first, character(1)),
                             collapse = "By")
      ))
    }
  }

  # --- Per-cell builder: fits model once, runs every (effect, spec) combo ---
  fit_fn <- function(df_cell) {
    if (model == "lmer") {
      lmerTest::lmer(formula, data = df_cell)
    } else {
      stats::lm(formula, data = df_cell)
    }
  }

  extract_fn <- function(fit, df_cell) {
    ty <- df_cell[[type_col]][1]
    ar <- df_cell[[area_col]][1]

    cell_rows <- NULL

    # Robust covariance matrix (computed once per cell when applicable)
    vc <- if (model == "lm" && vcov != "classical") {
      sandwich::vcovHC(fit, type = vcov)
    } else {
      NULL
    }

    for (eff in matching_effects) {
      row_match <- omnibus_df[omnibus_df$effect == eff &
                              omnibus_df$type == ty &
                              omnibus_df$area == ar, , drop = FALSE]
      if (nrow(row_match) == 0) next

      gate_p <- row_match[[p_col]][1]
      is_significant <- isTRUE(!is.na(gate_p) && gate_p < p_threshold)
      if (!is_significant && !fill_empty) next

      specs <- get_posthoc_specs(eff)

      for (spec in specs) {
        spec_df <- tryCatch({
          em_spec <- stats::as.formula(spec$em_formula)

          if (spec$do_contrast) {
            emt <- if (!is.null(vc)) {
              emmeans::emtrends(fit, em_spec, var = trend_var, vcov. = vc)
            } else {
              emmeans::emtrends(fit, em_spec, var = trend_var)
            }
            em_contrasts <- if (!is.null(spec$contrast_by)) {
              emmeans::contrast(emt, "pairwise",
                                by = spec$contrast_by, adjust = adjust)
            } else {
              emmeans::contrast(emt, "pairwise", adjust = adjust)
            }
            em_df <- as.data.frame(em_contrasts)
            em_df <- standardize_emmeans_columns(em_df, mode = "contrast")

            if (!is.null(spec$contrast_by) &&
                spec$contrast_by %in% colnames(em_df)) {
              em_df$predictor <- paste(em_df[[spec$contrast_by]],
                                       em_df$contrast, suffix, sep = "_")
            } else {
              em_df$predictor <- paste(em_df$contrast, suffix, sep = "_")
            }
          } else {
            # No contrast — just test slope against zero
            emt <- if (!is.null(vc)) {
              emmeans::emtrends(fit, em_spec, var = trend_var, vcov. = vc)
            } else {
              emmeans::emtrends(fit, em_spec, var = trend_var)
            }
            em_df <- as.data.frame(
              summary(emt, infer = c(TRUE, TRUE), null = 0)
            )
            em_df <- standardize_emmeans_columns(em_df, mode = "trend",
                                                 trend_var = trend_var)
            em_df$predictor <- paste(trend_var, "Slope_vsZero", sep = "_")
            # This branch uses p_value not p_bon — promote to p_bon
            # for output-schema uniformity, since this row still comes
            # out of a gated post-hoc call.
            if ("p_value" %in% colnames(em_df)) {
              colnames(em_df)[colnames(em_df) == "p_value"] <- "p_bon"
            }
            if (!"contrast" %in% colnames(em_df)) {
              em_df$contrast <- NA_character_
            }
          }

          em_df$statistic_type <- "t"
          em_df$cohens_d <- compute_cohens_d(fit, em_df$estimate)

          if (!is_significant) {
            em_df$estimate  <- NA_real_
            em_df$se        <- NA_real_
            em_df$df        <- NA_real_
            em_df$statistic <- NA_real_
            em_df$p_bon     <- NA_real_
            em_df$cohens_d  <- NA_real_
          }

          em_df$omnibus_effect <- eff
          em_df$posthoc_type   <- spec$label_prefix

          keep_cols <- c("omnibus_effect", "posthoc_type", "contrast",
                         "predictor", "estimate", "se", "df",
                         "statistic", "statistic_type", "p_bon", "cohens_d")
          keep_cols <- keep_cols[keep_cols %in% colnames(em_df)]
          em_df[, keep_cols, drop = FALSE]
        }, error = function(e) {
          message("  [WARNING] gated slopes failed for type = ", ty,
                  ", area = ", ar, ", effect = ", eff,
                  ", spec = ", spec$label_prefix, ": ",
                  conditionMessage(e))
          NULL
        })

        if (!is.null(spec_df) && nrow(spec_df) > 0) {
          cell_rows <- rbind(cell_rows, spec_df)
        }
      }
    }

    cell_rows
  }

  iterate_by_area(
    data       = data,
    area_col   = area_col,
    type_col   = type_col,
    min_rows   = min_rows,
    fit_fn     = fit_fn,
    extract_fn = extract_fn,
    fn_label   = "gated slopes"
  )
}

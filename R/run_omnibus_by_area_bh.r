#' @title Run area-wise omnibus and BH-correct the p-values in one call
#'
#' @description
#' Convenience wrapper that chains [run_omnibus_by_area()] with
#' [getBHvalues()] so that the common "fit omnibus, then BH-correct
#' for subsequent gated post-hocs" sequence is a single call at the
#' script level. The result is identical to calling the two functions
#' in sequence — this wrapper exists purely for ergonomics.
#'
#' The returned table has the same schema as `run_omnibus_by_area()`
#' plus an added `p_bh` column carrying Benjamini-Hochberg-adjusted
#' p-values. BH correction is applied in groups defined by
#' `bh_group_col`, which defaults to `c("type", "effect")` —
#' matching the convention used in the existing second-level
#' scripts, where each chromophore x omnibus-effect combination is
#' treated as an independent family of tests. Pass `"type"` for
#' single-column grouping or `NULL` to correct across the entire
#' table as one family.
#'
#' Downstream gated functions (`run_gated_posthoc()`,
#' `run_gated_slopes()`) default to gating on the `p_bh` column, so
#' passing the result of `run_omnibus_by_area_bh()` directly into a
#' gated function "just works" without any intermediate call.
#'
#' All arguments except the two `bh_*` arguments are forwarded to
#' `run_omnibus_by_area()`; see that function's documentation for
#' parameter semantics, the HC estimator discussion, and the rationale
#' for `model` being a required argument.
#'
#' @inheritParams run_omnibus_by_area
#' @param bh_group_col Character vector naming one or more columns
#'   in the omnibus table to group by when applying BH correction.
#'   Defaults to `c("type", "effect")` — each chromophore x omnibus
#'   effect combination is an independent family. This matches the
#'   family structure used in the pre-Tier-2 second-level analysis
#'   scripts, which applied BH correction grouped by chromophore
#'   (`type`) and effect name (`effect`), just under different
#'   column names in the old schema. Pass a single column name
#'   (e.g., `"type"`) for a coarser family, or `NULL` to correct
#'   across the entire table as one family.
#' @param bh_p_val_col Character scalar naming the raw p-value column
#'   in the omnibus table. Defaults to `"p_value"`, which is the name
#'   `run_omnibus_by_area()` produces.
#'
#' @return The omnibus table from `run_omnibus_by_area()` with an
#'   added `p_bh` column. Returns `NULL` if the underlying omnibus
#'   call produced no output.
#'
#' @seealso [run_omnibus_by_area()] for the underlying fit and the
#'   full discussion of model/vcov arguments;
#'   [getBHvalues()] for the BH correction semantics.
#'
#' @examples
#' \dontrun{
#' # Two-line "omnibus + BH + gated post-hoc" pattern:
#' df.omnibus <- run_omnibus_by_area_bh(
#'   data     = df.all,
#'   formula  = estimate ~ Group * cond + (1 | ID),
#'   area_col = "ROI",
#'   model    = "lmer"
#' )
#' df.ph.group <- run_gated_posthoc(
#'   omnibus_df       = df.omnibus,
#'   effect_name      = "Group",
#'   data             = df.all,
#'   formula          = estimate ~ Group * cond + (1 | ID),
#'   area_col         = "ROI",
#'   contrast_factor  = "Group",
#'   model            = "lmer"
#' )
#' }
#'
#' @export
run_omnibus_by_area_bh <- function(data,
                                   formula,
                                   area_col,
                                   type_col     = "type",
                                   model,
                                   vcov         = NULL,
                                   min_rows     = 5,
                                   effect_size  = TRUE,
                                   refit_singular = TRUE,
                                   bh_group_col = c("type", "effect"),
                                   bh_p_val_col = "p_value") {

  omnibus_df <- run_omnibus_by_area(
    data        = data,
    formula     = formula,
    area_col    = area_col,
    type_col    = type_col,
    model       = model,
    vcov        = vcov,
    min_rows    = min_rows,
    effect_size = effect_size,
    refit_singular = refit_singular
  )

  if (is.null(omnibus_df) || nrow(omnibus_df) == 0) return(omnibus_df)

  # getBHvalues() retains its legacy `pBH` output column name for
  # backwards compatibility with pre-Tier-2 chFnirs callers (e.g., the
  # ABC analysis scripts). Rename to the canonical snake_case `p_bh`
  # here so the Tier 2 output schema is consistent for new callers of
  # run_omnibus_by_area_bh() and downstream Tier 2 functions that
  # default to gating on `p_bh`.
  omnibus_bh <- getBHvalues(
    df_data   = omnibus_df,
    group_col = bh_group_col,
    p_val_col = bh_p_val_col
  )
  colnames(omnibus_bh)[colnames(omnibus_bh) == "pBH"] <- "p_bh"
  omnibus_bh
}

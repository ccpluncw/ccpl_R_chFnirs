#' @title Compute a row-aware heatmap figure height
#'
#' @description
#' Returns a `fig_height` value (in inches) suitable for passing to
#' [createHeatMapsFromLmerAreaAnalysis()], [create_effect_size_heatmap()],
#' or any of the primitive `create_*_heatmap()` functions, scaled so
#' each row of the heatmap keeps roughly the same visual height
#' regardless of how many predictors the data contains, and floored at
#' a minimum height that ensures the side legend (when present) is
#' fully visible.
#'
#' @details
#' Heatmaps reserve a fixed band of vertical space for the title,
#' rotated x-axis labels, and plot padding -- the chart strip itself
#' is whatever is left of `fig_height`. If the figure height is
#' simply scaled linearly with the row count, the chart strip becomes
#' margin-dominated at small `n` and rows look squished. This helper
#' instead holds the chart strip per row constant for `n` between 1
#' and `full_n`, accounting for `margin` overhead, and caps the
#' figure at `full_height` once `n >= full_n`. Above the cap, rows
#' shrink to fit the capped figure.
#'
#' At very small `n` the row-scaled figure can still be shorter than
#' the side-legend stack itself, which produces a clipped legend
#' (the colorbar headers or end labels disappear off the top or
#' bottom of the figure). `min_height` enforces a floor below which
#' the figure does not shrink, so the legend always renders fully.
#' At small `n` the rows are slightly taller than at `full_n` -- an
#' acceptable trade for a complete legend.
#'
#' When `show_legend = FALSE` the figure has no side-legend constraint
#' (the heatmap is being rendered without its legend, e.g. for a
#' publication panel that uses a single shared legend exported via
#' [create_heatmap_legend()]). In that case the floor drops from the
#' default 5 inches to 3 inches -- still enough for the title +
#' rotated x-axis labels + a couple of rows.
#'
#' Predictors whose name contains `"Intercept"` are excluded from the
#' row count by default, matching the chFnirs heatmap pipelines'
#' default `removeIntercept = TRUE` behaviour.
#'
#' @param df A long-format data.frame whose `pred_col` column gives
#'   the predictors that will become heatmap rows.
#' @param full_height Figure height in inches at and above the cap.
#'   Default `8.5`.
#' @param full_n Row count at which the figure reaches `full_height`.
#'   Default `10`. Increase to make rows smaller; decrease to make
#'   them larger.
#' @param margin Inches reserved for non-chart elements (title,
#'   rotated x-axis labels, plot padding). Default `2`. Tune if rows
#'   look too cramped (decrease) or too generous (increase) at small
#'   `n`.
#' @param show_legend Logical. Whether the heatmap this height feeds
#'   will be rendered with its side legend. Used to pick the default
#'   `min_height`. Default `TRUE` (matches the heatmap functions'
#'   default).
#' @param min_height Numeric floor on the returned height, in inches.
#'   When `NULL` (default), resolves to `5` if `show_legend` is `TRUE`
#'   (enough to fit the typical positive + negative colorbar legend
#'   stack with the not-tested swatch) or `3` if `show_legend` is
#'   `FALSE`. Pass an explicit value to override.
#' @param pred_col Name of the predictor column in `df`. Default
#'   `"predictor"`, matching the canonical chFnirs output schema.
#' @param removeIntercept Logical. If `TRUE` (default), predictors
#'   matching `"Intercept"` are dropped from the row count.
#'
#' @return A numeric scalar in inches, suitable for `fig_height`.
#'
#' @examples
#' \dontrun{
#' # Typical use: compute fig_height per heatmap call
#' createHeatMapsFromLmerAreaAnalysis(
#'   results_df = df.posthoc, area_col = "area", group_col = "type",
#'   sigValue_col = "p_bon", dv_col = "estimate",
#'   area_order = sort_rois_anatomically(unique(df.posthoc$area)),
#'   fig_height = heatmap_height(df.posthoc),
#'   output_prefix = "PostHoc"
#' )
#'
#' # Publication panel: heatmap without its legend, smaller floor
#' fh <- heatmap_height(df.posthoc, show_legend = FALSE)
#' create_posthoc_heatmap(..., show_legend = FALSE, fig_height = fh)
#'
#' # Explicit floor override
#' heatmap_height(df, min_height = 6)
#' }
#'
#' @export
heatmap_height <- function(df, full_height = 8.5, full_n = 10,
                           margin = 2,
                           show_legend = TRUE,
                           min_height = NULL,
                           pred_col = "predictor",
                           removeIntercept = TRUE) {
  preds <- unique(df[[pred_col]])
  if (removeIntercept) preds <- preds[!grepl("Intercept", preds)]
  n <- length(preds)
  per_row <- (full_height - margin) / full_n
  raw <- if (n >= full_n) full_height else margin + n * per_row

  if (is.null(min_height)) {
    min_height <- if (show_legend) 5 else 3
  }

  max(raw, min_height)
}

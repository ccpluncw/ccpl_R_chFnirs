#' @title Standardize emmeans output column names to the chFnirs schema
#'
#' @description
#' Renames the columns of a data frame produced by `emmeans::emmeans()`,
#' `emmeans::contrast()`, or `emmeans::emtrends()` to the canonical
#' chFnirs output schema (lowercase snake_case, `estimate` as the
#' generic value column, `statistic` as the generic test-statistic
#' column). Used internally by the area-wise wrappers but also exported so
#' that callers building their own emmeans-based analyses can produce
#' chFnirs-compatible output tables.
#'
#' The three `mode` values cover the three shapes of emmeans output:
#' \describe{
#'   \item{`"marginal"`}{Output of `summary(em_grid, ...)` or
#'     `test(em_grid, null = 0, joint = FALSE)`, where the value column
#'     is named `emmean` and the p-value is a raw unadjusted t-test vs
#'     zero. Renames `emmean` → `estimate`, `p.value` → `p_value`.}
#'   \item{`"contrast"`}{Output of `contrast(em_grid, "pairwise", ...)`,
#'     where the value column is already named `estimate` and the
#'     p-value column carries whatever adjustment was applied during
#'     the `contrast()` call (typically Bonferroni in chFnirs). Renames
#'     `p.value` → `p_bon`.}
#'   \item{`"trend"`}{Output of `summary(emtrends_object, ...)`, where
#'     the value column is named `<trend_var>.trend` and the p-value is
#'     raw. Renames `<trend_var>.trend` → `estimate`,
#'     `p.value` → `p_value`.}
#' }
#'
#' Common renames applied in all three modes: `SE` → `se`,
#' `t.ratio` → `statistic`, `z.ratio` → `statistic`. The `df` column
#' keeps its name. Any columns present in `df` that are not in the
#' rename map are passed through untouched, so emmeans factor-level
#' columns, confidence-interval columns, and `contrast` column labels
#' survive the call.
#'
#' @param df A data.frame produced by `as.data.frame()` on an emmeans,
#'   contrast, or emtrends result.
#' @param mode One of `"marginal"`, `"contrast"`, or `"trend"`. Picks
#'   the rename map appropriate to the shape of the input.
#' @param trend_var Character scalar naming the continuous moderator
#'   whose slope appears in the table when `mode = "trend"`. The column
#'   named `<trend_var>.trend` is renamed to `estimate`. Required when
#'   `mode = "trend"`; ignored otherwise.
#'
#' @return The input data frame with standardized column names.
#'
#' @examples
#' \dontrun{
#' # Marginal means
#' em <- emmeans::emmeans(fit, ~ Group)
#' df <- as.data.frame(emmeans::test(em, null = 0, joint = FALSE))
#' standardize_emmeans_columns(df, mode = "marginal")
#'
#' # Pairwise contrasts
#' em <- emmeans::emmeans(fit, ~ Group)
#' df <- as.data.frame(
#'   emmeans::contrast(em, "pairwise", adjust = "bonferroni")
#' )
#' standardize_emmeans_columns(df, mode = "contrast")
#'
#' # Slope from emtrends
#' em <- emmeans::emtrends(fit, ~ Group, var = "hormone")
#' df <- as.data.frame(summary(em, infer = c(TRUE, TRUE), null = 0))
#' standardize_emmeans_columns(df, mode = "trend", trend_var = "hormone")
#' }
#'
#' @export
standardize_emmeans_columns <- function(df,
                                        mode = c("marginal", "contrast", "trend"),
                                        trend_var = NULL) {
  mode <- match.arg(mode)

  # Dynamic rename first (trend-var column name varies)
  if (mode == "trend") {
    if (is.null(trend_var)) {
      stop("`trend_var` is required when mode = \"trend\".", call. = FALSE)
    }
    trend_col <- paste0(trend_var, ".trend")
    if (trend_col %in% colnames(df)) {
      colnames(df)[colnames(df) == trend_col] <- "estimate"
    }
  }

  # Static renames common across modes
  col_map <- c(
    "SE"      = "se",
    "t.ratio" = "statistic",
    "z.ratio" = "statistic"
  )

  if (mode == "marginal") {
    col_map[["emmean"]]  <- "estimate"
    col_map[["p.value"]] <- "p_value"
  } else if (mode == "contrast") {
    # emmeans contrast output already names the value column "estimate"
    col_map[["p.value"]] <- "p_bon"
  } else {
    # trend
    col_map[["p.value"]] <- "p_value"
  }

  for (old_name in names(col_map)) {
    if (old_name %in% colnames(df)) {
      colnames(df)[colnames(df) == old_name] <- unname(col_map[old_name])
    }
  }

  df
}

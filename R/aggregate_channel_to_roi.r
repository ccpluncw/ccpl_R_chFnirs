#' @title Aggregate per-channel response values to per-ROI weighted means
#'
#' @description
#' Given a long-format table of per-channel response values (e.g., GLM
#' betas, AUC amplitudes, or any other per-channel response measure)
#' together with a channel-to-ROI link table (a sparse weight matrix
#' produced by atlas projection, e.g., fOLD-style), returns a long-format
#' table of per-`(group_cols, ROI)` weighted means and a `coverage_frac`
#' diagnostic that captures what fraction of each ROI's total anatomical
#' weight survived per cell.
#'
#' The weighted mean is `sum(value * weight) / sum(weight)` across the
#' channels that contribute to each ROI, with `NA` values in `value_col`
#' tolerated and excluded from both numerator and denominator
#' (per-cell missingness from first-level QC or per-cell outlier
#' trimming). `coverage_frac` is `surviving_weight / total_weight` where
#' `total_weight` is the sum of weights for the ROI in the link table
#' and `surviving_weight` is the sum of weights for channels that
#' contributed a non-NA value to the cell.
#'
#' This function is pure compute - it does not apply a coverage gate or
#' emit any diagnostic figures. The caller filters on `coverage_frac`
#' (and writes any coverage histograms) after the function returns, so
#' that the aggregation and the gating decision stay independent.
#'
#' @details
#' **Channel-name convention.** `auc_df[[channel_col]]` and
#' `link_df[[channel_col]]` must use identical channel identifiers (e.g.,
#' `"Ch_11_1"` for source 11 / detector 1). The function joins on this
#' column verbatim - it does not construct channel names from
#' source/detector pairs. Caller is responsible for any pre-processing.
#'
#' **Many-to-many semantics.** The join between `auc_df` and `link_df`
#' is genuinely many-to-many (each channel contributes to multiple ROIs
#' with non-zero weight; each subject-cell has multiple channels). The
#' join uses `relationship = "many-to-many"` to suppress the dplyr
#' warning.
#'
#' **NA handling asymmetry.** NA in `value_col` is expected (per-channel
#' first-level QC dropouts or per-cell outlier trimming) and is
#' tolerated: the channel is excluded from both numerator and
#' denominator of the weighted mean and from `surviving_weight`. NA in
#' `weight_col` is NOT tolerated - atlas projections do not legitimately
#' produce NA weights (channels with no projection get weight `0`, not
#' `NA`); NA in weights indicates a corrupted link table and the
#' function errors out.
#'
#' **Empty cells.** A `(group_cols, ROI)` cell where every contributing
#' channel was NA-masked or had weight 0 yields `<value_col> = NaN` and
#' `coverage_frac = 0`. The row stays in the output; the caller's
#' coverage gate filters it.
#'
#' **Orphan channels.** Channels present in `auc_df` but absent from
#' `link_df` are dropped from the aggregation; the function emits a
#' warning listing them. (Common legitimate cause: short-separation
#' channels intentionally excluded from atlas projection.)
#'
#' @param auc_df Long-format data.frame of per-channel response values.
#'   Must contain `group_cols`, `channel_col`, and `value_col`.
#' @param link_df Channel-to-ROI link table. Must contain `channel_col`,
#'   `name_col`, and `weight_col`. Each `(channel, name)` pair must
#'   appear at most once. Weights must be non-NA numerics (typically
#'   in `[0, 1]`).
#' @param group_cols Character vector naming the columns in `auc_df`
#'   that define one subject-cell (e.g., subject + condition + type).
#'   The function groups by `c(group_cols, name_col)` for the aggregation.
#'   Default `c("ID", "cond", "type")` (SPT convention; ABC overrides
#'   with `c("ID", "Condition", "Type")`).
#' @param channel_col Character scalar naming the channel-identifier
#'   column in both `auc_df` and `link_df`. Default `"Channel"`.
#' @param value_col Character scalar naming the per-channel response
#'   column in `auc_df` (the value being aggregated). Could be GLM
#'   betas, AUC amplitudes, or any other per-channel scalar measure.
#'   Default `"beta"`.
#' @param weight_col Character scalar naming the atlas-weight column
#'   in `link_df`. Default `"weight"`.
#' @param name_col Character scalar naming the ROI-identifier column
#'   in `link_df`. Default `"Name"`.
#' @param out_roi_col Character scalar, the name the output ROI column
#'   should take. Default `"ROI"`.
#'
#' @return A data.frame with one row per `(group_cols, ROI)` cell.
#'   Columns, in order: `group_cols`, `out_roi_col`, `value_col`
#'   (carrying the input column name), `coverage_frac`. Sorted by
#'   `(group_cols..., out_roi_col)`. Empty `auc_df` returns an empty
#'   data.frame with this schema.
#'
#' @seealso [apply_roi_merge()] and [describe_merge_candidates()] for
#'   building / validating the link table used here.
#'
#' @examples
#' \dontrun{
#' library(readxl)
#' link <- read_excel("ROI_aal.xlsx")
#' link$Channel <- paste("Ch", link$source, link$detector, sep = "_")
#'
#' auc <- read.csv("Final_AUC_Results.csv")
#'
#' # SPT (lowercase group cols, value_col = "AUC"):
#' roi_df <- aggregate_channel_to_roi(
#'   auc, link,
#'   group_cols = c("ID", "cond", "type"),
#'   value_col  = "AUC"
#' )
#'
#' # ABC (uppercase group cols, same value_col override):
#' roi_df_abc <- aggregate_channel_to_roi(
#'   auc, link,
#'   group_cols = c("ID", "Condition", "Type"),
#'   value_col  = "AUC"
#' )
#'
#' # Caller applies the coverage gate downstream:
#' roi_df <- roi_df[roi_df$coverage_frac >= 0.5, ]
#' }
#'
#' @importFrom dplyr left_join group_by summarise mutate select all_of
#' @importFrom rlang sym !! := .data
#' @export
aggregate_channel_to_roi <- function(
  auc_df,
  link_df,
  group_cols  = c("ID", "cond", "type"),
  channel_col = "Channel",
  value_col   = "beta",
  weight_col  = "weight",
  name_col    = "Name",
  out_roi_col = "ROI"
) {
  # ---- input validation -----------------------------------------------------
  if (!is.data.frame(auc_df))  stop("`auc_df` must be a data.frame.",  call. = FALSE)
  if (!is.data.frame(link_df)) stop("`link_df` must be a data.frame.", call. = FALSE)

  for (arg in c("channel_col", "value_col", "weight_col",
                "name_col", "out_roi_col")) {
    v <- get(arg)
    if (!is.character(v) || length(v) != 1L || is.na(v) || nchar(v) == 0L) {
      stop("`", arg, "` must be a single non-empty character string.",
           call. = FALSE)
    }
  }
  if (!is.character(group_cols) || length(group_cols) < 1L) {
    stop("`group_cols` must be a non-empty character vector.", call. = FALSE)
  }

  need_auc  <- c(group_cols, channel_col, value_col)
  need_link <- c(channel_col, name_col, weight_col)
  miss_auc  <- setdiff(need_auc,  colnames(auc_df))
  miss_link <- setdiff(need_link, colnames(link_df))
  if (length(miss_auc) > 0L) {
    stop("`auc_df` is missing required columns: ",
         paste(miss_auc, collapse = ", "), ".", call. = FALSE)
  }
  if (length(miss_link) > 0L) {
    stop("`link_df` is missing required columns: ",
         paste(miss_link, collapse = ", "), ".", call. = FALSE)
  }

  if (nrow(link_df) == 0L) {
    stop("`link_df` is empty - cannot aggregate to zero ROIs.",
         call. = FALSE)
  }
  if (anyNA(link_df[[weight_col]])) {
    stop("NA values found in `link_df[[\"", weight_col, "\"]]`. Atlas ",
         "projections do not legitimately produce NA weights; channels ",
         "with no projection should have weight 0. Fix the link table ",
         "and retry.", call. = FALSE)
  }
  dup_key <- paste(link_df[[channel_col]], link_df[[name_col]], sep = "\r")
  if (anyDuplicated(dup_key)) {
    bad <- unique(dup_key[duplicated(dup_key)])
    stop("`link_df` has duplicate (", channel_col, ", ", name_col,
         ") rows (",  length(bad), " duplicate pair(s)): ",
         paste(utils::head(bad, 5), collapse = "; "),
         if (length(bad) > 5L) ", ..." else "", ".",
         call. = FALSE)
  }

  # Output schema for the empty-input shortcut.
  out_schema <- c(group_cols, out_roi_col, value_col, "coverage_frac")

  if (nrow(auc_df) == 0L) {
    empty <- as.data.frame(matrix(nrow = 0L, ncol = length(out_schema)))
    colnames(empty) <- out_schema
    return(empty)
  }

  # ---- orphan-channel warning ----------------------------------------------
  auc_channels  <- unique(auc_df[[channel_col]])
  link_channels <- unique(link_df[[channel_col]])
  orphans <- setdiff(auc_channels, link_channels)
  if (length(orphans) > 0L) {
    warning("`auc_df` contains ", length(orphans), " channel(s) with no ",
            "entry in `link_df` (dropped from aggregation): ",
            paste(orphans, collapse = ", "), ".", call. = FALSE)
  }

  # ---- ROI total weight (independent of auc_df) ----------------------------
  total_w <- dplyr::group_by(link_df, !!rlang::sym(name_col))
  total_w <- dplyr::summarise(
    total_w,
    .roi_total_weight = sum(.data[[weight_col]]),
    .groups           = "drop"
  )

  # ---- join and aggregate ---------------------------------------------------
  joined <- dplyr::left_join(
    auc_df, link_df,
    by           = channel_col,
    relationship = "many-to-many"
  )

  # Filter to rows that actually mapped to an ROI (orphan auc rows have NA Name
  # after the left join and shouldn't contribute to weighted sums).
  joined <- joined[!is.na(joined[[name_col]]), , drop = FALSE]

  grouped <- dplyr::group_by(
    joined,
    dplyr::across(dplyr::all_of(c(group_cols, name_col)))
  )

  # Treat NA in value_col as "missing" - exclude from numerator AND denominator.
  agg <- dplyr::summarise(
    grouped,
    .num             = sum(.data[[value_col]] * .data[[weight_col]],
                           na.rm = TRUE),
    .denom           = sum(.data[[weight_col]] *
                             as.numeric(!is.na(.data[[value_col]]))),
    .groups          = "drop"
  )

  out <- dplyr::left_join(agg, total_w, by = name_col)
  out <- dplyr::mutate(out,
    !!value_col      := ifelse(.data$.denom > 0,
                               .data$.num / .data$.denom,
                               NaN),
    coverage_frac    = .data$.denom / .data$.roi_total_weight
  )
  out <- dplyr::select(out, -dplyr::all_of(c(".num", ".denom",
                                             ".roi_total_weight")))

  # Rename name_col -> out_roi_col only if they differ.
  if (name_col != out_roi_col) {
    names(out)[names(out) == name_col] <- out_roi_col
  }

  # Reorder columns to the locked output schema, sort deterministically.
  out <- as.data.frame(out)
  out <- out[, out_schema, drop = FALSE]
  out <- out[do.call(order, out[, c(group_cols, out_roi_col), drop = FALSE]),
             , drop = FALSE]
  rownames(out) <- NULL
  out
}

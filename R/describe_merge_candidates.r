#' @title Describe candidate ROI merges from a channel-weight table
#'
#' @description
#' Diagnostic-only companion to [apply_roi_merge()]. Given a long-format
#' channel-weight table (one row per `(source, detector, ROI)` triple with a
#' `weight` column), reports, for every unordered ROI pair, the geometric
#' diagnostics that describe whether the fNIRS probe can resolve the two
#' regions: atlas-cosine similarity of the per-channel weight vectors,
#' dominant channel per ROI, whether the dominant channels coincide, and a
#' connected-component label derived from the per-pair geometric decision.
#'
#' Optionally also reports the Pearson correlation between the two ROIs'
#' AUC values across observations when a subject-level long-format AUC
#' table is supplied via `empirical_data`. Empirical r is reported for
#' methods-section disclosure; it is *not* used to gate the merge decision
#' (see Details).
#'
#' This function does not modify or merge anything — it produces the
#' descriptive table that informs the caller's hand-written `name_map`
#' passed to [apply_roi_merge()].
#'
#' @details
#' **Geometric rule composition.** Two diagnostics are combined according
#' to `geometric_rule`:
#' \itemize{
#'   \item `"any"` (default): pair passes if atlas-cosine similarity
#'     `>= atlas_cos_threshold` OR (when `use_dominant_channel = TRUE`)
#'     the two ROIs share their dominant channel. Reproduces the v3
#'     cap-resolution rule.
#'   \item `"all"`: both criteria must hold.
#'   \item `"cos_only"`: ignore dominant-channel sharing entirely;
#'     pair passes iff atlas cosine clears the threshold.
#' }
#'
#' **Why no adjacency gate.** For any reasonable single-array probe,
#' geometric resolvability already entails anatomical adjacency, because
#' fNIRS channels have finite spatial extent (~3 cm sensitivity banana).
#' If two ROIs share channels with non-trivial weight, they are physically
#' close in the head by construction. v3's "AND common Brodmann zone"
#' clause was therefore redundant and is not exposed in this API.
#'
#' **Why no empirical-r gate.** Empirical between-ROI correlation
#' confounds shared neural source with shared systemic physiology, shared
#' task-locked HRF, shared subject-level random effects, and shared scalp
#' coupling. None of those justify a *measurement-resolvability* merge
#' decision. Empirical r is reported here for disclosure purposes only
#' (e.g., to populate the v3 precentral_l / left_lpfc r=0.77 sentence in
#' the methods section).
#'
#' **Chain handling.** When a transitive chain arises in the per-pair
#' decision matrix (e.g., A-B passes, B-C passes, A-C fails), this
#' function reports each pair's literal `geometric_decision` and assigns
#' all three ROIs the same `component_id` from union-find. The caller is
#' the final arbiter — chains are visible in the table, and the
#' caller's `name_map` (passed to [apply_roi_merge()]) decides what
#' actually merges.
#'
#' @param weights_df Long-format data.frame with one row per
#'   `(source, detector, ROI)` triple. Must contain numeric columns named
#'   exactly `source`, `detector`, `weight` and a character column
#'   identifying the ROI (default `"Name"`; override via `name_col`).
#' @param atlas_cos_threshold Numeric scalar in `[0, 1]`. A pair passes
#'   the atlas-cosine criterion when `cos >= atlas_cos_threshold`. Set
#'   to `Inf` (or `> 1`) to effectively disable. Default `0.5`.
#' @param geometric_rule Character scalar in `c("any", "all", "cos_only")`.
#'   How atlas cosine and dominant-channel sharing combine. Default
#'   `"any"`.
#' @param use_dominant_channel Logical. Whether to incorporate
#'   dominant-channel sharing into the decision. Ignored when
#'   `geometric_rule = "cos_only"`. Default `TRUE`.
#' @param empirical_data Optional long-format data.frame with at least
#'   `<empirical_obs_cols>`, an ROI column (matching `name_col` in
#'   `weights_df`), and an AUC column (default `"AUC"`; override via
#'   `auc_col`). Passing this enables the `empirical_r` column in the
#'   output. Default `NULL`.
#' @param empirical_obs_cols Character vector naming the columns of
#'   `empirical_data` that define one observation. Default `"ID"`
#'   (subject-level r). For per-cell r, pass e.g. `c("ID", "cond", "type")`.
#' @param name_col Character scalar, name of the ROI-identifier column in
#'   `weights_df` and `empirical_data`. Default `"Name"`.
#' @param auc_col Character scalar, name of the AUC column in
#'   `empirical_data`. Default `"AUC"`.
#' @param include_zero_overlap Logical. Whether to retain rows for ROI
#'   pairs whose channel-weight vectors share no support
#'   (`n_shared_channels == 0`, equivalent to `atlas_cos == 0`). These
#'   pairs are trivially not merge candidates; filtering them out
#'   collapses the output to a scannable table of pairs with non-zero
#'   probe-level overlap. Default `FALSE`. Set `TRUE` to recover the
#'   full N-choose-2 enumeration (useful for sanity-checking probe
#'   coverage or for test fixtures).
#'
#' @return A data.frame (tibble) with one row per unordered ROI pair,
#'   columns: `roi_i`, `roi_j`, `atlas_cos`, `n_shared_channels`,
#'   `dom_channel_i`, `dom_channel_j`, `shared_dom_channel`,
#'   `empirical_r`, `geometric_decision`, `component_id`.
#'   Rows are sorted by `(roi_i, roi_j)` alphabetically.
#'
#' @seealso [apply_roi_merge()] for the executor that consumes the
#'   caller's hand-written `name_map`.
#'
#' @examples
#' \dontrun{
#' library(readxl)
#' weights <- read_excel("ROI_aal.xlsx")
#' audit <- describe_merge_candidates(
#'   weights,
#'   atlas_cos_threshold  = 0.5,
#'   geometric_rule       = "any",
#'   use_dominant_channel = TRUE
#' )
#' # Subset to the merge-eligible pairs:
#' audit[audit$geometric_decision, ]
#'
#' # With empirical r from spt_byROI.csv (long format):
#' roi_data <- read.csv("spt_byROI.csv")
#' audit2 <- describe_merge_candidates(
#'   weights, empirical_data = roi_data,
#'   empirical_obs_cols = c("ID", "cond", "type"),
#'   name_col = "Name", auc_col = "AUC"
#' )
#' }
#'
#' @importFrom stats cor
#' @importFrom utils combn
#' @export
describe_merge_candidates <- function(
  weights_df,
  atlas_cos_threshold   = 0.5,
  geometric_rule        = c("any", "all", "cos_only"),
  use_dominant_channel  = TRUE,
  empirical_data        = NULL,
  empirical_obs_cols    = "ID",
  name_col              = "Name",
  auc_col               = "AUC",
  include_zero_overlap  = FALSE
) {
  geometric_rule <- match.arg(geometric_rule)

  if (!is.data.frame(weights_df)) {
    stop("`weights_df` must be a data.frame.", call. = FALSE)
  }
  required_cols <- c("source", "detector", "weight", name_col)
  missing_cols  <- setdiff(required_cols, colnames(weights_df))
  if (length(missing_cols) > 0L) {
    stop("`weights_df` is missing required columns: ",
         paste(missing_cols, collapse = ", "), ".", call. = FALSE)
  }
  if (!is.numeric(atlas_cos_threshold) || length(atlas_cos_threshold) != 1L) {
    stop("`atlas_cos_threshold` must be a single numeric value.", call. = FALSE)
  }

  rois <- sort(unique(as.character(weights_df[[name_col]])))
  n    <- length(rois)
  if (n < 2L) {
    stop("`weights_df` must contain at least two distinct ROIs in column `",
         name_col, "`.", call. = FALSE)
  }

  # Wide channel x ROI weight matrix. Key (source, detector) -> single channel.
  channel_key <- paste0("S", weights_df$source, "-D", weights_df$detector)
  channels    <- sort(unique(channel_key))
  W           <- matrix(0, nrow = length(channels), ncol = n,
                        dimnames = list(channels, rois))
  for (i in seq_len(nrow(weights_df))) {
    ch  <- channel_key[i]
    roi <- as.character(weights_df[[name_col]][i])
    W[ch, roi] <- W[ch, roi] + weights_df$weight[i]
  }

  # Atlas-cosine similarity matrix.
  col_norms <- sqrt(colSums(W^2))
  zero_cols <- col_norms == 0
  if (any(zero_cols)) {
    warning("ROIs with all-zero weights (no probe coverage) detected: ",
            paste(rois[zero_cols], collapse = ", "),
            ". Their cosine similarity to any other ROI is undefined (NaN).",
            call. = FALSE)
  }
  # Denominator outer product, avoid 0/0 NaN by setting NaN cells later.
  denom    <- outer(col_norms, col_norms)
  cos_mat  <- crossprod(W) / denom
  cos_mat[!is.finite(cos_mat)] <- NA_real_

  # Dominant channel per ROI (NA if no positive weight).
  dom_idx <- apply(W, 2, function(v) {
    if (all(v == 0)) NA_integer_ else which.max(v)
  })
  dom_channel <- ifelse(is.na(dom_idx), NA_character_, channels[dom_idx])
  names(dom_channel) <- rois

  # Shared-channels count: # of channels both ROIs touch with non-zero weight.
  has_weight <- W > 0
  shared_n   <- crossprod(has_weight)

  # Empirical Pearson r between ROI columns, if data supplied.
  emp_r_mat <- matrix(NA_real_, nrow = n, ncol = n, dimnames = list(rois, rois))
  if (!is.null(empirical_data)) {
    if (!is.data.frame(empirical_data)) {
      stop("`empirical_data` must be a data.frame or NULL.", call. = FALSE)
    }
    needed <- c(empirical_obs_cols, name_col, auc_col)
    missing_emp <- setdiff(needed, colnames(empirical_data))
    if (length(missing_emp) > 0L) {
      stop("`empirical_data` is missing required columns: ",
           paste(missing_emp, collapse = ", "), ".", call. = FALSE)
    }
    obs_key <- do.call(paste, c(empirical_data[empirical_obs_cols], sep = "\r"))
    wide    <- tapply(empirical_data[[auc_col]],
                      list(obs_key, as.character(empirical_data[[name_col]])),
                      FUN = mean, na.rm = TRUE)
    wide[!is.finite(wide)] <- NA_real_
    present_rois <- intersect(colnames(wide), rois)
    if (length(present_rois) >= 2L) {
      r_sub <- stats::cor(wide[, present_rois, drop = FALSE],
                          use = "pairwise.complete.obs")
      emp_r_mat[present_rois, present_rois] <- r_sub
    }
  }

  # Per-pair table over the upper triangle.
  pair_idx <- utils::combn(n, 2)
  i_idx    <- pair_idx[1, ]
  j_idx    <- pair_idx[2, ]

  out <- data.frame(
    roi_i              = rois[i_idx],
    roi_j              = rois[j_idx],
    atlas_cos          = cos_mat[cbind(i_idx, j_idx)],
    n_shared_channels  = as.integer(shared_n[cbind(i_idx, j_idx)]),
    dom_channel_i      = dom_channel[i_idx],
    dom_channel_j      = dom_channel[j_idx],
    shared_dom_channel = !is.na(dom_channel[i_idx]) &
                         !is.na(dom_channel[j_idx]) &
                         (dom_channel[i_idx] == dom_channel[j_idx]),
    empirical_r        = emp_r_mat[cbind(i_idx, j_idx)],
    stringsAsFactors   = FALSE
  )

  # Geometric decision per pair under the chosen rule.
  cos_pass <- !is.na(out$atlas_cos) & out$atlas_cos >= atlas_cos_threshold
  dom_pass <- out$shared_dom_channel
  out$geometric_decision <- switch(
    geometric_rule,
    any       = if (use_dominant_channel) cos_pass | dom_pass else cos_pass,
    all       = if (use_dominant_channel) cos_pass & dom_pass else cos_pass,
    cos_only  = cos_pass
  )

  # Union-find connected components on the geometric decisions.
  parent <- seq_len(n)
  find_root <- function(x) {
    while (parent[x] != x) {
      parent[x] <<- parent[parent[x]]
      x <- parent[x]
    }
    x
  }
  union <- function(a, b) {
    ra <- find_root(a); rb <- find_root(b)
    if (ra != rb) parent[ra] <<- rb
  }
  passing <- which(out$geometric_decision)
  for (k in passing) union(i_idx[k], j_idx[k])
  roots <- vapply(seq_len(n), find_root, integer(1L))
  # Renumber to dense 1..K
  comp_label <- as.integer(factor(roots, levels = unique(roots)))
  names(comp_label) <- rois

  # component_id is the joint component label when both ROIs sit in the same
  # union-find group (whether by direct pass or transitive chain); NA when
  # they're in different components. The integer label itself is arbitrary
  # (renumbered 1..K above) — what matters is which rows share a value.
  comp_i <- comp_label[out$roi_i]
  comp_j <- comp_label[out$roi_j]
  out$component_id <- ifelse(comp_i == comp_j, comp_i, NA_integer_)

  if (!include_zero_overlap) {
    out <- out[out$n_shared_channels > 0L, , drop = FALSE]
  }

  out <- out[order(out$roi_i, out$roi_j), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' @title Apply a hand-written ROI merge map to a channel-weight table
#'
#' @description
#' Executor companion to [describe_merge_candidates()]. Given a
#' long-format channel-weight table and a caller-supplied named list
#' mapping merged-ROI names to their constituent ROI names, returns the
#' table with the constituent rows collapsed by `(source, detector,
#' merged_name)` via algebraic sum of weights (valid under AAL parcel
#' disjointness). ROIs not referenced anywhere in `name_map` are carried
#' through unchanged.
#'
#' This function makes no merge decisions of its own — it executes the
#' caller's `name_map` verbatim. The descriptive table from
#' [describe_merge_candidates()] informs the human who writes the
#' `name_map`, and any anatomical/theoretical carve-outs (e.g., retaining
#' precentral_l as standalone despite geometric similarity with left LPFC)
#' are expressed by *omitting* those ROIs from `name_map`.
#'
#' @details
#' **Why algebraic sum is correct.** Channel weights from atlas
#' projection (e.g., fOLD-style) are per-AAL-parcel sensitivities of one
#' source-detector channel. AAL parcels are spatially disjoint, so the
#' sensitivity of a single channel to the union of two parcels equals the
#' sum of its sensitivities to each parcel individually. Summing weights
#' is therefore exact, not an approximation.
#'
#' **Validation.** The function errors if a `name_map` key references an
#' ROI not present in `weights_df`, or if any ROI appears in more than one
#' merged group (which would create an ambiguous mapping).
#'
#' **Writing to disk.** When `save_dir` is non-`NULL`, the function writes
#' `file.path(save_dir, file_name)` for each element of `save_dir`. A
#' vector of paths handles the common pattern of writing the same output
#' to multiple canonical directories.
#'
#' @param weights_df Long-format data.frame with one row per
#'   `(source, detector, ROI)` triple. Must contain numeric columns
#'   `source`, `detector`, `weight` and a character ROI column (default
#'   `"Name"`; override via `name_col`).
#' @param name_map A named list where each name is the desired merged-ROI
#'   name and each value is a character vector of constituent ROI names.
#'   Constituent names must exist in `weights_df[[name_col]]`; no
#'   constituent may appear in more than one entry. Pass `list()` to
#'   leave the table unchanged.
#' @param save_dir Optional character vector of directory paths. When
#'   non-`NULL`, the merged table is written as an `.xlsx` file to
#'   `file.path(save_dir[i], file_name)` for every `i`. Directories are
#'   created if they do not exist. Default `NULL` (no write; return only).
#' @param file_name Character scalar, the output file name used in each
#'   `save_dir`. Default `"ROI_aal_merged.xlsx"`.
#' @param name_col Character scalar naming the ROI column in
#'   `weights_df`. Default `"Name"`.
#'
#' @return A data.frame (long-format, same schema as `weights_df`) with
#'   constituent ROIs collapsed to their merged names. Sorted by
#'   `(Name, source, detector)`.
#'
#' @seealso [describe_merge_candidates()] for producing the diagnostic
#'   table that informs the caller's `name_map`.
#'
#' @examples
#' \dontrun{
#' library(readxl)
#' weights <- read_excel("ROI_aal.xlsx")
#'
#' # The v3 12-ROI cap-resolution scheme. precentral_l is intentionally
#' # omitted (theoretical-centrality carve-out documented in the methods
#' # supplement).
#' name_map <- list(
#'   left_lpfc  = c("frontal_inf_tri_l", "frontal_mid_l"),
#'   right_lpfc = c("frontal_inf_tri_r", "frontal_mid_r"),
#'   right_sfg  = c("frontal_sup_r",     "frontal_sup_medial_r"),
#'   right_tpj  = c("angular_r",         "parietal_inf_r")
#' )
#'
#' merged <- apply_roi_merge(
#'   weights, name_map,
#'   save_dir  = c("path/to/study1", "path/to/study2"),
#'   file_name = "ROI_aal_merged.xlsx"
#' )
#' }
#'
#' @importFrom writexl write_xlsx
#' @export
apply_roi_merge <- function(
  weights_df,
  name_map,
  save_dir  = NULL,
  file_name = "ROI_aal_merged.xlsx",
  name_col  = "Name"
) {
  if (!is.data.frame(weights_df)) {
    stop("`weights_df` must be a data.frame.", call. = FALSE)
  }
  required_cols <- c("source", "detector", "weight", name_col)
  missing_cols  <- setdiff(required_cols, colnames(weights_df))
  if (length(missing_cols) > 0L) {
    stop("`weights_df` is missing required columns: ",
         paste(missing_cols, collapse = ", "), ".", call. = FALSE)
  }
  if (!is.list(name_map)) {
    stop("`name_map` must be a named list.", call. = FALSE)
  }
  if (length(name_map) > 0L &&
      (is.null(names(name_map)) || any(names(name_map) == ""))) {
    stop("Every element of `name_map` must be named.", call. = FALSE)
  }

  # Validate constituents.
  all_rois         <- unique(as.character(weights_df[[name_col]]))
  all_constituents <- unlist(name_map, use.names = FALSE)
  unknown <- setdiff(all_constituents, all_rois)
  if (length(unknown) > 0L) {
    stop("`name_map` references ROIs not present in `weights_df[[\"",
         name_col, "\"]]`: ",
         paste(unknown, collapse = ", "), ".", call. = FALSE)
  }
  dup <- all_constituents[duplicated(all_constituents)]
  if (length(dup) > 0L) {
    stop("ROI(s) appear in more than one `name_map` entry: ",
         paste(unique(dup), collapse = ", "), ".", call. = FALSE)
  }
  collision <- intersect(names(name_map), setdiff(all_rois, all_constituents))
  if (length(collision) > 0L) {
    stop("Merged name(s) collide with un-mapped ROI(s) already in ",
         "`weights_df`: ", paste(collision, collapse = ", "), ".",
         call. = FALSE)
  }

  # Build the constituent -> merged-name lookup.
  lookup <- setNames(
    rep(names(name_map), lengths(name_map)),
    all_constituents
  )

  # Apply the mapping; un-mapped ROIs stay as-is.
  out <- weights_df
  orig_names <- as.character(out[[name_col]])
  remap_idx <- orig_names %in% names(lookup)
  out[[name_col]][remap_idx] <- lookup[orig_names[remap_idx]]

  # Group-by-sum on (source, detector, name_col). Keep base R to avoid
  # introducing dplyr as a hard dependency for this single op.
  key <- paste(out$source, out$detector, out[[name_col]], sep = "\r")
  agg <- aggregate(
    out$weight,
    by  = list(key = key),
    FUN = sum
  )
  # Recover the keyed columns by splitting; preserve the column order from
  # the original frame.
  split_key <- do.call(rbind, strsplit(agg$key, "\r", fixed = TRUE))
  merged_df <- data.frame(
    source           = as.integer(split_key[, 1]),
    detector         = as.integer(split_key[, 2]),
    weight           = agg$x,
    stringsAsFactors = FALSE
  )
  merged_df[[name_col]] <- split_key[, 3]

  # Reorder columns to match the input frame's column ordering.
  in_cols <- intersect(colnames(weights_df), colnames(merged_df))
  merged_df <- merged_df[, in_cols, drop = FALSE]

  merged_df <- merged_df[
    order(merged_df[[name_col]], merged_df$source, merged_df$detector),
    , drop = FALSE
  ]
  rownames(merged_df) <- NULL

  if (!is.null(save_dir)) {
    if (!is.character(save_dir)) {
      stop("`save_dir` must be NULL or a character vector of directories.",
           call. = FALSE)
    }
    if (!is.character(file_name) || length(file_name) != 1L) {
      stop("`file_name` must be a single character string.", call. = FALSE)
    }
    for (d in save_dir) {
      if (!dir.exists(d)) {
        dir.create(d, recursive = TRUE, showWarnings = FALSE)
      }
      out_path <- file.path(d, file_name)
      writexl::write_xlsx(merged_df, out_path)
      message("Wrote: ", out_path)
    }
  }

  merged_df
}

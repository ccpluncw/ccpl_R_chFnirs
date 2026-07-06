#' @title Sort fNIRS channel names by anatomical position
#'
#' @description
#' Companion to [sort_rois_anatomically()] for fNIRS channel data. Takes
#' a character vector of channel names (e.g., `"Ch_13_10"`), parses out
#' the source and detector indices, looks up their 3D positions in a
#' `digpts.txt` file, and returns the names in anatomically-banded
#' plotting order (anterior to posterior, left to right within bands).
#'
#' Intended for use as the `area_order` argument of chFnirs heatmap
#' functions when areas are channels. Both this function and
#' [sort_rois_anatomically()] satisfy the chFnirs pluggable-sorter
#' contract `function(names, ...) -> ordered names`, so they are
#' interchangeable behind any heatmap that accepts `area_order`.
#'
#' @details
#' The default parser expects names of the form
#' `"<prefix>_<source>_<detector>"` where `<prefix>` is `Ch` or `ch`.
#' Names that do not match are an error -- supply a custom
#' `channel_parser` if the data uses a different convention. Original
#' names are preserved exactly in the output (no case normalisation).
#'
#' Banding, fiducial validation, and coordinate flipping are handled by
#' the underlying engine [sort_fnirs_channels()]; see that function for
#' details on the digpts format and the AP/LR orientation check.
#'
#' @param channel_names Character vector of channel names. Must parse
#'   to integer (source, detector) pairs via `channel_parser`.
#' @param digpts_path Path to a digpts.txt file (see
#'   [sort_fnirs_channels()] for format).
#' @param n_bands Number of AP bands. Default 2.
#' @param y_breaks Optional explicit band boundaries. Overrides
#'   `n_bands`.
#' @param channel_parser Function `function(names) -> matrix` that
#'   returns a 2-column integer matrix (source, detector) with one row
#'   per input name. Default parses names matching
#'   `^[Cc]h_(\\d+)_(\\d+)$`.
#' @param verbose Logical. Default `FALSE`.
#'
#' @return Character vector of input channel names in plotting order
#'   (anterior to posterior, left to right within each band). The
#'   per-channel diagnostics data.frame is attached as
#'   `attr(result, "diagnostics")` (see [sort_fnirs_channels()]).
#'
#' @examples
#' \dontrun{
#' channels <- unique(df$area)  # e.g., c("Ch_1_1", "Ch_2_3", ...)
#' ordered  <- sort_fnirs_channels_by_name(channels, "digpts.txt")
#' df$area  <- factor(df$area, levels = ordered)
#'
#' # Pass directly to a heatmap as `area_order`:
#' createHeatMapsFromLmerAreaAnalysis(
#'   results_df, area_col = "area", ...,
#'   area_order = sort_fnirs_channels_by_name(unique(df$area), "digpts.txt")
#' )
#'
#' # Custom parser for a different naming convention. The parser takes the
#' # character vector of names and returns a 2-column integer matrix of
#' # (source, detector). Examples for common formats:
#'
#' # "S1-D1", "S12-D3" (NIRStar / Aurora style):
#' p_dash <- function(names) {
#'   m <- regmatches(names, regexec("^S(\\d+)-D(\\d+)$", names))
#'   do.call(rbind, lapply(m, function(x) as.integer(x[2:3])))
#' }
#' sort_fnirs_channels_by_name(c("S1-D1", "S2-D3"), "digpts.txt",
#'                             channel_parser = p_dash)
#'
#' # "Source1_Detector1" with a trailing chromophore token, e.g.
#' # "Source1_Detector1 hbo":
#' p_long <- function(names) {
#'   m <- regmatches(names, regexec("Source(\\d+)_Detector(\\d+)", names))
#'   do.call(rbind, lapply(m, function(x) as.integer(x[2:3])))
#' }
#' }
#'
#' @seealso [sort_fnirs_channels()] (the matrix engine),
#'   [sort_rois_anatomically()] (the ROI counterpart).
#'
#' @export
sort_fnirs_channels_by_name <- function(channel_names,
                                        digpts_path,
                                        n_bands = 2,
                                        y_breaks = NULL,
                                        channel_parser = NULL,
                                        verbose = FALSE) {

  if (!is.character(channel_names) || length(channel_names) == 0) {
    stop("channel_names must be a non-empty character vector.", call. = FALSE)
  }
  if (anyDuplicated(channel_names) > 0) {
    stop("channel_names must be unique; got duplicates: ",
         paste(unique(channel_names[duplicated(channel_names)]), collapse = ", "),
         call. = FALSE)
  }

  if (is.null(channel_parser)) {
    channel_parser <- function(names) {
      m <- regmatches(names, regexec("^[Cc]h_(\\d+)_(\\d+)$", names))
      lens <- vapply(m, length, integer(1))
      bad <- names[lens != 3]
      if (length(bad) > 0) {
        stop("Could not parse channel names (expected 'Ch_<source>_<detector>'): ",
             paste(bad, collapse = ", "),
             ". Supply `channel_parser` for non-default formats.",
             call. = FALSE)
      }
      sd <- t(vapply(m, function(x) as.integer(x[2:3]), integer(2)))
      colnames(sd) <- c("source", "detector")
      sd
    }
  }

  sd <- channel_parser(channel_names)
  if (!is.matrix(sd) || ncol(sd) != 2) {
    stop("channel_parser must return a 2-column matrix (source, detector).",
         call. = FALSE)
  }
  if (nrow(sd) != length(channel_names)) {
    stop("channel_parser returned ", nrow(sd), " rows for ",
         length(channel_names), " channel names.", call. = FALSE)
  }

  # Closure: engine uses this to emit names; look up the exact input
  # name for each (source, detector) pair so capitalisation and any
  # parser-specific quirks are preserved unchanged in the output.
  name_lookup <- setNames(channel_names,
                          paste(sd[, 1], sd[, 2], sep = "_"))
  preserve_input_name <- function(s, d) {
    name_lookup[[paste(s, d, sep = "_")]]
  }

  sort_fnirs_channels(
    sd_matrix     = sd,
    digpts_path   = digpts_path,
    n_bands       = n_bands,
    y_breaks      = y_breaks,
    channel_namer = preserve_input_name,
    verbose       = verbose
  )
}

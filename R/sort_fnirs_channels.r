#' @title Sort fNIRS channels by anatomical position (matrix engine)
#'
#' @description
#' Internal engine that takes raw source/detector indices and a digpts
#' file and returns the channel names in plotting order. Most callers
#' should use [sort_fnirs_channels_by_name()], which parses channel
#' names like `"Ch_<source>_<detector>"` and delegates here.
#'
#' @details
#' Channels are placed at the midpoint of their source-detector pair
#' (read from `digpts_path`), grouped into anterior-posterior bands
#' (either via user-supplied `y_breaks` or by detecting the largest
#' gaps in the channel midpoint y-distribution), then sorted left-to-
#' right within each band.
#'
#' The z coordinate of each midpoint is parsed and reported in the
#' diagnostics attribute but is not used in the sort -- a curved-scalp
#' anatomical sort that respects superior/inferior would require
#' additional information beyond a flat AP/LR ordering.
#'
#' Fiducials (`nz`, `iz`, `ar`, `al`) are used to verify the AP/LR axes
#' are oriented as expected (+y anterior, +x right). When either axis
#' is reversed, the relevant coordinate is flipped on the fly so the
#' sort is correct regardless of input convention. With
#' `verbose = TRUE` a `message()` reports which axes were flipped.
#' If fiducials are not all present a `message()` notes that the
#' default convention is being assumed.
#'
#' @param sd_matrix A 2-column matrix or data.frame with one row per
#'   channel. Column 1 = source index (positive integer); column 2 =
#'   detector index (positive integer). Coerced to integer.
#' @param digpts_path Path to a digpts.txt file. Lines have the form
#'   `"<key>: <x> <y> <z>"`, with sources `s1, s2, ...`, detectors
#'   `d1, d2, ...`, and (ideally) fiducials `nz`, `iz`, `ar`, `al`.
#' @param n_bands Number of AP bands. Default 2. `1` produces a pure
#'   left-to-right sort. Ignored when `y_breaks` is supplied.
#' @param y_breaks Optional numeric vector of y thresholds (in digpts
#'   units) defining band boundaries explicitly. Overrides `n_bands`.
#' @param channel_namer Function `function(source, detector) -> string`.
#'   Default `function(s, d) sprintf("Ch_%d_%d", s, d)`, matching the
#'   convention used in SPT and ABC fNIRS data. Override for sites that
#'   use a different naming scheme.
#' @param verbose Logical. If `TRUE`, print orientation-check and band
#'   summary to the console. Default `FALSE`.
#'
#' @return A character vector of channel names in plotting order
#'   (anterior to posterior, left to right within each band). A
#'   data.frame of per-channel diagnostics is attached as
#'   `attr(result, "diagnostics")` with columns: `channel`, `source`,
#'   `detector`, `x`, `y`, `z`, `band`, `band_label`, `order_idx`.
#'
#' @keywords internal
sort_fnirs_channels <- function(sd_matrix,
                                digpts_path,
                                n_bands = 2,
                                y_breaks = NULL,
                                channel_namer = function(s, d) sprintf("Ch_%d_%d", s, d),
                                verbose = FALSE) {

  # ---- Validate sd_matrix ----
  sd_raw <- as.matrix(sd_matrix)
  if (ncol(sd_raw) != 2) {
    stop("sd_matrix must have 2 columns (source, detector); got ",
         ncol(sd_raw), ".", call. = FALSE)
  }
  if (!is.numeric(sd_raw)) {
    stop("sd_matrix must be numeric.", call. = FALSE)
  }
  if (any(is.na(sd_raw))) {
    stop("sd_matrix contains NA; every row must have a source and a detector.",
         call. = FALSE)
  }
  if (any(sd_raw != as.integer(sd_raw))) {
    stop("sd_matrix entries must be integers (got non-integer values).",
         call. = FALSE)
  }
  if (any(sd_raw < 1)) {
    stop("sd_matrix entries must be positive integers (>= 1).", call. = FALSE)
  }
  sd <- matrix(as.integer(sd_raw), ncol = 2)
  n  <- nrow(sd)

  # ---- Parse digpts ----
  if (!file.exists(digpts_path)) {
    stop("digpts file not found: ", digpts_path, call. = FALSE)
  }
  lines <- readLines(digpts_path)
  lines <- trimws(lines); lines <- lines[nzchar(lines)]
  pts <- list()
  for (ln in lines) {
    m <- regmatches(ln, regexec(
      "^(\\w+):\\s*(-?[0-9.]+)\\s+(-?[0-9.]+)\\s+(-?[0-9.]+)", ln))[[1]]
    if (length(m) == 5) pts[[m[2]]] <- as.numeric(m[3:5])
  }
  if (length(pts) == 0) {
    stop("No coordinates parsed from digpts file: ", digpts_path,
         ". Expected lines like 's1: 12.3 45.6 78.9'.", call. = FALSE)
  }

  # ---- Verify and (if needed) flip orientation from fiducials ----
  if (all(c("nz", "iz", "ar", "al") %in% names(pts))) {
    ap_flipped <- pts$nz[2] <= pts$iz[2]
    lr_flipped <- pts$ar[1] <= pts$al[1]
    if (ap_flipped || lr_flipped) {
      for (key in names(pts)) {
        if (ap_flipped) pts[[key]][2] <- -pts[[key]][2]
        if (lr_flipped) pts[[key]][1] <- -pts[[key]][1]
      }
      if (verbose) {
        flipped <- c(if (ap_flipped) "AP (y)", if (lr_flipped) "LR (x)")
        message("Flipped coordinate axes: ", paste(flipped, collapse = " and "),
                " (so +y = anterior, +x = right).")
      }
    } else if (verbose) {
      message(sprintf(
        "Coordinate orientation OK: nz_y=%.1f, iz_y=%.1f, ar_x=%.1f, al_x=%.1f",
        pts$nz[2], pts$iz[2], pts$ar[1], pts$al[1]))
    }
  } else if (verbose) {
    message("Note: fiducials (nz/iz/ar/al) not all present; ",
            "assuming +y = anterior, +x = right.")
  }

  # ---- Verify all needed source/detector keys exist before allocating ----
  needed <- unique(c(paste0("s", sd[, 1]), paste0("d", sd[, 2])))
  missing_keys <- setdiff(needed, names(pts))
  if (length(missing_keys) > 0) {
    stop("Missing positions in digpts for: ",
         paste(missing_keys, collapse = ", "), call. = FALSE)
  }

  # ---- Compute channel midpoints ----
  channels <- data.frame(
    channel  = character(n),
    source   = integer(n),
    detector = integer(n),
    x        = numeric(n),
    y        = numeric(n),
    z        = numeric(n),
    stringsAsFactors = FALSE
  )
  for (i in seq_len(n)) {
    s_i <- sd[i, 1]; d_i <- sd[i, 2]
    mid <- (pts[[paste0("s", s_i)]] + pts[[paste0("d", d_i)]]) / 2
    channels$channel[i]  <- channel_namer(s_i, d_i)
    channels$source[i]   <- s_i
    channels$detector[i] <- d_i
    channels$x[i] <- mid[1]; channels$y[i] <- mid[2]; channels$z[i] <- mid[3]
  }

  # ---- Decide bands ----
  if (!is.null(y_breaks)) {
    breaks <- sort(unique(c(-Inf, y_breaks, Inf)))
    n_bands_actual <- length(breaks) - 1
  } else if (n_bands <= 1) {
    breaks <- c(-Inf, Inf)
    n_bands_actual <- 1
  } else {
    sy   <- sort(channels$y, decreasing = TRUE)  # anterior first
    gaps <- -diff(sy)                            # all positive
    n_breaks_needed <- min(n_bands - 1, length(gaps))
    if (n_breaks_needed < 1) {
      breaks <- c(-Inf, Inf)
      n_bands_actual <- 1
    } else {
      gap_idx <- order(gaps, decreasing = TRUE)[seq_len(n_breaks_needed)]
      bdry    <- vapply(gap_idx, function(k) (sy[k] + sy[k + 1]) / 2, numeric(1))
      breaks  <- sort(c(-Inf, bdry, Inf))
      n_bands_actual <- length(breaks) - 1
    }
  }

  # cut() returns 1 for lowest y, n for highest. We want band 1 = highest y.
  channels$band <- (n_bands_actual + 1L) -
    as.integer(cut(channels$y, breaks = breaks,
                   labels = FALSE, include.lowest = TRUE))

  band_y_ranges <- vapply(seq_len(n_bands_actual), function(b) {
    yb <- channels$y[channels$band == b]
    if (length(yb)) sprintf("y=[%.1f, %.1f]", min(yb), max(yb)) else "empty"
  }, character(1))
  channels$band_label <- sprintf("Band %d (%s)", channels$band,
                                 band_y_ranges[channels$band])

  # ---- Sort: band primary, x ascending ----
  channels <- channels[order(channels$band, channels$x), ]
  channels$order_idx <- seq_len(nrow(channels))
  rownames(channels) <- NULL

  if (verbose) {
    message("Sorted ", n, " channels into ", n_bands_actual, " band(s):")
    for (b in seq_len(n_bands_actual)) {
      yb <- channels$y[channels$band == b]
      message(sprintf("  Band %d: %2d channels, y range [%6.2f, %6.2f]",
                      b, length(yb), min(yb), max(yb)))
    }
  }

  ordered_names <- channels$channel
  attr(ordered_names, "diagnostics") <- channels
  ordered_names
}

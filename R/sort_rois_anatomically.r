#' @title Sort ROI Names Based on a Neuroanatomical Hierarchy
#'
#' @description
#' This function takes a character vector of ROI (Region of Interest) names
#' and reorders them according to a predefined anatomical logic. The sorting
#' prioritizes lobes (anterior to posterior), then sub-regions (superior to inferior,
#' medial to lateral), and finally hemisphere (left before right).
#'
#' Both AAL-style names (`<region>_<sub>_<l|r>`) and merged-ROI names with a
#' `left_` / `right_` prefix (e.g., `left_lpfc`, `right_tpj`) are recognized.
#'
#' @details
#' The anatomical sorting logic is hardcoded within the function, defining a
#' specific order for:
#' \itemize{
#'   \item{\strong{Primary Regions/Lobes:}} Frontal -> Lateral PFC (merged token \code{lpfc}) -> Precentral -> Parietal -> TPJ (merged token \code{tpj}) -> Angular -> Temporal (anterior to posterior).
#'   \item{\strong{Sub-Regions:}} Superior Medial -> Superior -> Middle -> Inferior Triangularis -> Inferior (superior to inferior/medial to lateral).
#'   \item{\strong{Hemisphere:}} Left ('_l' or 'left_') before Right ('_r' or 'right_').
#' }
#' Warnings will be issued if any ROI name cannot be parsed or matched to the
#' defined anatomical categories, and such ROIs will be placed at the end of the
#' sorted list. This function is designed for a specific naming scheme and
#' may require modification for different conventions.
#'
#' @param roi_vector A character vector containing the names of the ROIs to be sorted.
#'   Expected formats:
#'   \itemize{
#'     \item AAL-style suffix: \code{"<lobe>_<subregion>_<hemisphere>"}, e.g.,
#'       \code{"frontal_sup_medial_l"}, \code{"angular_r"}, \code{"temporal_mid_l"}.
#'     \item Merged-ROI prefix: \code{"<hemisphere>_<token>"}, e.g.,
#'       \code{"left_lpfc"}, \code{"right_tpj"}.
#'   }
#'
#' @return A character vector containing the ROI names sorted according to the
#'   specified neuroanatomical hierarchy.
#'
#' @seealso \code{\link[base]{order}}
#'
#' @examples
#' my_rois <- c("angular_l", "angular_r", "frontal_inf_tri_l", "frontal_inf_tri_r",
#'              "frontal_mid_l", "frontal_mid_r", "frontal_sup_l",
#'              "frontal_sup_medial_l", "frontal_sup_medial_r", "frontal_sup_r",
#'              "parietal_inf_l", "parietal_inf_r", "precentral_l", "precentral_r",
#'              "temporal_mid_l", "temporal_mid_r",
#'              "left_lpfc", "right_lpfc", "left_tpj", "right_tpj")
#' sorted_rois <- sort_rois_anatomically(my_rois)
#' print(sorted_rois)
#'
#' @export
sort_rois_anatomically <- function(roi_vector) {

  # 1. Define Custom Anatomical Orders (manual ranking)
  #
  # `lpfc` (lateral prefrontal cortex) is a merged token absorbing
  # frontal_inf_tri + frontal_mid, so it sorts between frontal_sup and
  # precentral (rank 1.5).
  #
  # `sfg` (superior frontal gyrus) is a merged token absorbing
  # frontal_sup + frontal_sup_medial, so it sorts within the frontal
  # region (rank 1.2) — anatomically more superior than lpfc but
  # below the lateral/medial frontal_sup pair when both are still
  # standalone in the same vector.
  #
  # `tpj` (temporoparietal junction) is a merged token absorbing
  # angular + parietal_inf, so it sorts between parietal and temporal
  # (rank 3.5; same rank as angular which it largely overlaps).
  lobe_order <- c(
    "frontal" = 1, "sfg" = 1.2, "lpfc" = 1.5, "precentral" = 2,
    "parietal" = 3, "tpj" = 3.5, "angular" = 4, "temporal" = 5
  )
  sub_region_order <- c(
    "sup_medial" = 1, "sup" = 2, "mid" = 3, "inf_tri" = 4, "inf" = 5,
    "precentral" = 1, "angular" = 1, "none" = 99,
    "lpfc" = 1, "sfg" = 1, "tpj" = 1
  )
  hemisphere_order <- c("l" = 0, "r" = 1)

  # 2. Parse and Rank Each ROI
  roi_data <- data.frame(
    original_name = roi_vector, lobe_rank = 999,
    sub_region_rank = 999, hemisphere_rank = 999,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(roi_vector)) {
    name <- roi_vector[i]

    # Hemisphere detection — supports both naming conventions.
    #   AAL-style:    `<region>_<l|r>`         (e.g., "frontal_sup_l")
    #   Merged-style: `<left|right>_<token>`   (e.g., "left_lpfc")
    hemisphere <- NA
    base_name  <- name
    if (startsWith(name, "left_")) {
      hemisphere <- "l"
      base_name  <- sub("^left_", "", name)
    } else if (startsWith(name, "right_")) {
      hemisphere <- "r"
      base_name  <- sub("^right_", "", name)
    } else if (endsWith(name, "_l")) {
      hemisphere <- "l"
      base_name  <- sub("_l$", "", name)
    } else if (endsWith(name, "_r")) {
      hemisphere <- "r"
      base_name  <- sub("_r$", "", name)
    }

    if (is.na(hemisphere)) {
      warning(paste("Could not determine hemisphere for:", name))
      next
    }
    roi_data$hemisphere_rank[i] <- hemisphere_order[hemisphere]

    found_match <- FALSE
    for (lobe in names(lobe_order)) {
      if (startsWith(base_name, lobe)) {
        roi_data$lobe_rank[i] <- lobe_order[lobe]
        sub_name <- sub(paste0("^", lobe, "_?"), "", base_name)

        if (sub_name %in% names(sub_region_order)) {
          roi_data$sub_region_rank[i] <- sub_region_order[sub_name]
        } else if (sub_name == "") {
          # No sub-region: solo lobes (precentral, angular) and merged
          # tokens (lpfc, tpj) get a sub-region rank from their own
          # entry in sub_region_order; everything else gets "none".
          if (lobe %in% c("precentral", "angular", "lpfc", "sfg", "tpj")) {
            roi_data$sub_region_rank[i] <- sub_region_order[lobe]
          } else {
            roi_data$sub_region_rank[i] <- sub_region_order["none"]
          }
        } else {
          warning(paste("Could not determine sub-region for:", name,
                        " using sub-name:", sub_name))
        }
        found_match <- TRUE
        break
      }
    }
    if (!found_match) {
      warning(paste("Could not determine primary region/lobe for:", name))
    }
  }

  # 3. Sort the data frame by the ranks
  sorted_roi_data <- roi_data[order(
    roi_data$lobe_rank, roi_data$sub_region_rank, roi_data$hemisphere_rank
  ), ]

  # 4. Return the sorted original names
  return(sorted_roi_data$original_name)
}

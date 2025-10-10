#' @title Sort ROI Names Based on a Neuroanatomical Hierarchy
#'
#' @description
#' This function takes a character vector of ROI (Region of Interest) names
#' and reorders them according to a predefined anatomical logic. The sorting
#' prioritizes lobes (anterior to posterior), then sub-regions (superior to inferior,
#' medial to lateral), and finally hemisphere (left before right).
#'
#' @details
#' The anatomical sorting logic is hardcoded within the function, defining a
#' specific order for:
#' \itemize{
#'   \item{\strong{Primary Regions/Lobes:}} Frontal -> Precentral -> Parietal -> Angular -> Temporal (anterior to posterior).
#'   \item{\strong{Sub-Regions:}} Superior Medial -> Superior -> Middle -> Inferior Triangularis -> Inferior (superior to inferior/medial to lateral).
#'   \item{\strong{Hemisphere:}} Left ('_l') before Right ('_r').
#' }
#' Warnings will be issued if any ROI name cannot be parsed or matched to the
#' defined anatomical categories, and such ROIs will be placed at the end of the
#' sorted list. This function is designed for a specific naming scheme and
#' may require modification for different conventions.
#'
#' @param roi_vector A character vector containing the names of the ROIs to be sorted.
#'   Expected format: `"[lobe/region]_[subregion/gyrus]_[hemisphere]"`.
#'   Examples: "frontal_sup_medial_l", "angular_r", "temporal_mid_l".
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
#'              "temporal_mid_l", "temporal_mid_r", "unparseable_roi_name")
#' sorted_rois <- sort_rois_anatomically(my_rois)
#' print(sorted_rois)
#'
#' @export
sort_rois_anatomically <- function(roi_vector) {

  # 1. Define Custom Anatomical Orders (manual ranking)
  lobe_order <- c(
    "frontal" = 1, "precentral" = 2, "parietal" = 3, "angular" = 4, "temporal" = 5
  )
  sub_region_order <- c(
    "sup_medial" = 1, "sup" = 2, "mid" = 3, "inf_tri" = 4, "inf" = 5,
    "precentral" = 1, "angular" = 1, "none" = 99
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
    hemisphere <- ifelse(endsWith(name, "_l"), "l", ifelse(endsWith(name, "_r"), "r", NA))

    if (is.na(hemisphere)) {
      warning(paste("Could not determine hemisphere for:", name))
      next
    }
    roi_data$hemisphere_rank[i] <- hemisphere_order[hemisphere]
    base_name <- sub("_[lr]$", "", name)

    found_match <- FALSE
    for (lobe in names(lobe_order)) {
      if (startsWith(base_name, lobe)) {
        roi_data$lobe_rank[i] <- lobe_order[lobe]
        sub_name <- sub(paste0(lobe, "_?"), "", base_name)
        
        if (sub_name %in% names(sub_region_order)) {
          roi_data$sub_region_rank[i] <- sub_region_order[sub_name]
        } else if (sub_name == "") {
          roi_data$sub_region_rank[i] <- sub_region_order["none"]
        } else if (lobe %in% c("precentral", "angular")) {
           roi_data$sub_region_rank[i] <- sub_region_order[lobe]
        } else {
          warning(paste("Could not determine sub-region for:", name, " using sub-name:", sub_name))
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
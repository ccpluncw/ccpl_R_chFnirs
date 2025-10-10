#' @title Create a Significance Table from Model Results
#'
#' @description
#' Reshapes a long-format data frame of model results into a wide-format table.
#' The table shows the significance and direction of predictors across different
#' areas, and can optionally be grouped by another variable.
#'
#' @param df.data A data frame in long format containing model results. It should
#'   include columns for predictors, areas, beta values, and significance values.
#' @param areaColName The name of the column that specifies the area (e.g., "ROI").
#' @param sigColName The name of the column containing the significance value (e.g., "pVal").
#' @param group_col The name of an optional column to group the results by.
#'   If provided, a separate table is created for each level of this factor. Default: `NULL`.
#' @param discrete Logical. If `TRUE`, the output will be -1, 0, or 1 based on the
#'   significance threshold. If `FALSE`, it will be the signed p-value. Default: `TRUE`.
#' @param sigThresh The numeric significance threshold (e.g., 0.05). Default: 0.05.
#' @param removeIntercept Logical. If `TRUE`, rows corresponding to the model
#'   intercept will be removed from the output table. Default: `TRUE`.
#'
#' @return A list of data frames. If `group_col` is `NULL`, the list will contain a
#'   single data frame. Otherwise, it will contain a data frame for each group.
#'
#' @importFrom dplyr relocate
#' @importFrom rlang sym
#' @importFrom magrittr %>%
#'
#' @examples
#' # Create sample model output data
#' model_data <- data.frame(
#'   Predictor = rep(c("(Intercept)", "Age", "Sex"), times = 4),
#'   Beta = rnorm(12),
#'   pVal = runif(12, 0, 0.1),
#'   ROI = rep(c("Area1", "Area2"), each = 6),
#'   Group = rep(c("A", "B"), each = 3, times = 2)
#' )
#'
#' # --- Example 1: Ungrouped analysis ---
#' result_list_ungrouped <- CreateSigBetaTable(
#'   df.data = model_data,
#'   areaColName = "ROI",
#'   sigColName = "pVal"
#' )
#' print(result_list_ungrouped$all_data)
#'
#' # --- Example 2: Grouped analysis ---
#' result_list_grouped <- CreateSigBetaTable(
#'   df.data = model_data,
#'   areaColName = "ROI",
#'   sigColName = "pVal",
#'   group_col = "Group"
#' )
#' print(result_list_grouped)
#'
#' @export
CreateSigBetaTable2 <- function(df.data, areaColName, sigColName, group_col = NULL, discrete = TRUE, sigThresh = 0.05, removeIntercept = TRUE) {

  masterTable <- list()

  # 1. Create a list of data frames to process.
  # This handles both the grouped and ungrouped cases elegantly.
  df_list <- if (!is.null(group_col)) {
    # If group_col is provided, split the data frame by that column.
    split(df.data, df.data[[group_col]])
  } else {
    # If group_col is NULL, create a list containing the single, whole data frame.
    list(all_data = df.data)
  }

  # 2. Loop through the list of data frames (will be one or many).
  for (group_name in names(df_list)) {
    df.tmp <- df_list[[group_name]]

    area <- unique(df.tmp[[areaColName]])
    predictors <- unique(df.tmp$Predictor)

    # Create output table
    tmpMatrix <- matrix(0, nrow = length(predictors), ncol = length(area))
    df.out <- as.data.frame(tmpMatrix)
    colnames(df.out) <- area
    df.out$Predictor <- predictors
    
    # The pipe (%>%), relocate(), and sym() are imported via @importFrom
    df.out <- df.out %>% relocate(Predictor)

    # Conditionally add the group column back if it exists
    if (!is.null(group_col)) {
      df.out[[group_col]] <- group_name
      df.out <- df.out %>% relocate(!!sym(group_col))
    }

    if (removeIntercept) {
      df.out <- df.out[!grepl("\\(Int", df.out$Predictor), ]
    }

    # Populate the table with significance values
    for (ch in area) {
      df.tmp.2 <- df.tmp[df.tmp[[areaColName]] == ch, ]
      for (pr in predictors) {
        # Skip if the predictor doesn't exist in the output after removing intercept
        if (!pr %in% df.out$Predictor) next

        # Check if there is data to pull from
        if (nrow(df.tmp.2[df.tmp.2$Predictor == pr, ]) > 0) {
            if (discrete == TRUE) {
              df.out[df.out$Predictor == pr, ch] <- (df.tmp.2[df.tmp.2$Predictor == pr, sigColName] < sigThresh) * sign(df.tmp.2[df.tmp.2$Predictor == pr, "Beta"])
            } else {
              df.out[df.out$Predictor == pr, ch] <- df.tmp.2[df.tmp.2$Predictor == pr, sigColName] * sign(df.tmp.2[df.tmp.2$Predictor == pr, "Beta"])
            }
        }
      }
    }

    masterTable[[group_name]] <- df.out
  }

  return(masterTable)
}
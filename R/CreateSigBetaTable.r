#' @title Create Grouped Significance Tables from Model Results
#'
#' @description
#' Reshapes a long-format data frame of model results into a wide-format table.
#' This function iterates through a required 'Type' column in the input data,
#' creating a separate significance table for each level of 'Type'.
#'
#' @param df.data A data frame in long format containing model results. It must
#'   include a column named 'Type' that will be used for grouping.
#' @param areaColName The name of the column that specifies the area (e.g., "ROI").
#' @param sigColName The name of the column containing the significance value (e.g., "pVal").
#' @param discrete Logical. If `TRUE`, the output will be -1, 0, or 1 based on the
#'   significance threshold. If `FALSE`, it will be the signed p-value. Default: `TRUE`.
#' @param sigThresh The numeric significance threshold (e.g., 0.05). Default: 0.05.
#' @param removeIntercept Logical. If `TRUE`, rows corresponding to the model
#'   intercept will be removed from the output table. Default: `TRUE`.
#'
#' @return A list of data frames, with each data frame corresponding to a unique
#'   value from the 'Type' column.
#'
#' @importFrom dplyr relocate
#' @importFrom magrittr %>%
#'
#' @examples
#' # Create sample model output data with a 'Type' column
#' model_data_typed <- data.frame(
#'   Predictor = rep(c("(Intercept)", "Age"), times = 4),
#'   Beta = rnorm(8),
#'   pVal = runif(8, 0, 0.1),
#'   ROI = rep(c("Area1", "Area2"), each = 4),
#'   Type = rep(c("Fast", "Slow"), each = 2, times = 2)
#' )
#'
#' # Generate the list of significance tables
#' result_list <- CreateSigBetaTable(
#'   df.data = model_data_typed,
#'   areaColName = "ROI",
#'   sigColName = "pVal"
#' )
#' print(result_list)
#'
#' @export
CreateSigBetaTable <- function(df.data, areaColName, sigColName, discrete = TRUE, sigThresh = 0.05, removeIntercept = TRUE) {

  masterTable <- list()
  types <- unique(df.data$Type)

  for (tp in types) {
    type <- tp
    df.tmp <- df.data[df.data$Type == tp, ]

    area <- unique(df.tmp[[areaColName]])
    predictors <- unique(df.tmp$Predictor)

    # Create output table
    tmpMatrix <- matrix(0, nrow = length(predictors), ncol = length(area))
    df.out <- as.data.frame(tmpMatrix)
    colnames(df.out) <- area
    df.out$Predictor <- predictors
    
    # move predictor and type to the front
    df.out <- df.out %>% relocate(Predictor)
    df.out$Type <- type
    df.out <- df.out %>% relocate(Type)
    
    if (removeIntercept) {
      df.out <- df.out[!grepl("\\(Int", df.out$Predictor), ]
    }

    # Populate the table with significance values
    for (ch in area) {
      df.tmp.2 <- df.tmp[df.tmp[[areaColName]] == ch, ]
      for (pr in predictors) {
        # Skip if this predictor was removed (e.g., the intercept)
        if (!pr %in% df.out$Predictor) next
        
        # Check that a result exists for this predictor/area combination
        df.row <- df.tmp.2[df.tmp.2$Predictor == pr, ]
        if (nrow(df.row) > 0) {
          if (discrete == TRUE) {
            df.out[df.out$Predictor == pr, ch] <- (df.row[[sigColName]] < sigThresh) * sign(df.row[["Beta"]])
          } else {
            df.out[df.out$Predictor == pr, ch] <- df.row[[sigColName]] * sign(df.row[["Beta"]])
          }
        }
      }
    }
    masterTable[[type]] <- df.out
  }
  return(masterTable)
}
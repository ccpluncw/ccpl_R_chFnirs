#' @title Calculate False Discovery Rate q-values
#'
#' @description
#' A wrapper for the `qvalue::qvalue` function to calculate q-values from
#' p-values. This function's main utility is its ability to perform the
#' calculation on subsets of data defined by a grouping column.
#'
#' @param df_data The input dataframe.
#' @param group_col The name of the column to group by. If `NULL`, the calculation
#'   is performed on the entire dataframe at once. Default: `NULL`.
#' @param p_val_col The name of the column that contains the p-values.
#'
#' @return The original dataframe with an added column named "qVal".
#'
#' @importFrom qvalue qvalue
#'
#' @examples
#' # Create sample data
#' set.seed(123)
#' my_data <- data.frame(
#'   p_value = runif(100),
#'   group = rep(c("A", "B"), each = 50)
#' )
#'
#' # --- Example 1: Calculate q-values for the whole dataset ---
#' data_with_qvals <- getQvalues(
#'   df_data = my_data,
#'   group_col = NULL,
#'   p_val_col = "p_value"
#' )
#' head(data_with_qvals)
#'
#' # --- Example 2: Calculate q-values separately for each group ---
#' data_with_grouped_qvals <- getQvalues(
#'   df_data = my_data,
#'   group_col = "group",
#'   p_val_col = "p_value"
#' )
#' head(data_with_grouped_qvals)
#'
#' @export
getQvalues <- function(df_data, group_col = NULL, p_val_col) {

  # 1. Create a list of data frames to process.
  df_list <- if (!is.null(group_col)) {
    # If group_col is provided, split the data frame by group into a list.
    split(df_data, df_data[[group_col]])
  } else {
    # If group_col is NULL, create a list containing the single, whole data frame.
    list(all_data = df_data)
  }

  # 2. Apply the q-value calculation to each data frame in the list.
  processed_list <- lapply(df_list, function(df_chunk) {
    if (nrow(df_chunk) > 1) {
      # The tryCatch block handles a common error in qvalue() when pi0 cannot be estimated.
      # In such cases, it retries by setting pi0 = 1.
      df_chunk$qVal <- tryCatch({
        qvalue(df_chunk[[p_val_col]])$qvalues
      }, error = function(e) {
        qvalue(df_chunk[[p_val_col]], pi0 = 1)$qvalues
      })
    } else {
      # q-value is the same as p-value when there is only one test.
      df_chunk$qVal <- df_chunk[[p_val_col]]
    }
    return(df_chunk)
  })

  # 3. Combine the processed data frames back into one.
  df_data_out <- do.call(rbind, processed_list)

  # Restore original row order by removing the names set by split()
  rownames(df_data_out) <- NULL

  return(df_data_out)
}
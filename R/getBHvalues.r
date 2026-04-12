#' @title Calculate Benjamini-Hochberg FDR-Adjusted p-values
#'
#' @description
#' Calculates Benjamini-Hochberg (1995) FDR-adjusted p-values using base R's
#' `p.adjust` function. This function's main utility is its ability to perform
#' the correction on subsets of data defined by one or more grouping columns,
#' so that separate "families" of tests can be FDR-controlled independently
#' (e.g., each chromophore x omnibus effect combination gets its own family).
#'
#' Unlike Storey's q-value method, BH requires no estimation of pi0 and is
#' reliable regardless of the number of tests — making it appropriate for
#' fNIRS channel-level analysis where test counts are typically 20-50.
#'
#' @param df_data The input dataframe.
#' @param group_col The name(s) of the column(s) to group by. Accepts either
#'   a single character string (e.g., `"type"`) or a character vector of
#'   column names (e.g., `c("type", "effect")`). When multiple columns are
#'   supplied, the data frame is split on the unique combinations of their
#'   values so that each combination forms an independent family of tests.
#'   If `NULL`, the correction is performed on the entire dataframe at once.
#'   Default: `NULL`.
#' @param p_val_col The name of the column that contains the p-values.
#'
#' @return The original dataframe with an added column named `"pBH"`. Row
#'   order may be reordered into grouping order when `group_col` is
#'   supplied; downstream code should join or merge by the grouping columns
#'   rather than rely on positional order.
#'
#' @references
#' Benjamini, Y. & Hochberg, Y. (1995). Controlling the false discovery rate:
#'   a practical and powerful approach to multiple testing. \emph{Journal of the
#'   Royal Statistical Society Series B}, 57(1), 289-300.
#'
#' Singh, A.K. & Dan, I. (2006). Exploring the false discovery rate in
#'   multichannel NIRS. \emph{NeuroImage}, 33(2), 542-549.
#'
#' Yucel, M.A. et al. (2021). Best practices for fNIRS publications.
#'   \emph{Neurophotonics}, 8(1), 012101.
#'
#' @examples
#' # Create sample data
#' set.seed(123)
#' my_data <- data.frame(
#'   p_value = runif(100),
#'   type    = rep(c("hbo", "hbr"), each = 50),
#'   effect  = rep(c("Group", "cond"), times = 50)
#' )
#'
#' # --- Example 1: Correct across the whole dataset ---
#' data_with_bh <- getBHvalues(
#'   df_data   = my_data,
#'   group_col = NULL,
#'   p_val_col = "p_value"
#' )
#' head(data_with_bh)
#'
#' # --- Example 2: Correct separately for each type ---
#' data_with_grouped_bh <- getBHvalues(
#'   df_data   = my_data,
#'   group_col = "type",
#'   p_val_col = "p_value"
#' )
#' head(data_with_grouped_bh)
#'
#' # --- Example 3: Correct separately for each type x effect family ---
#' data_with_nested_bh <- getBHvalues(
#'   df_data   = my_data,
#'   group_col = c("type", "effect"),
#'   p_val_col = "p_value"
#' )
#' head(data_with_nested_bh)
#'
#' @importFrom stats p.adjust
#' @export
getBHvalues <- function(df_data, group_col = NULL, p_val_col) {

  # 1. Validate input columns.
  if (!is.null(group_col)) {
    missing_cols <- setdiff(group_col, colnames(df_data))
    if (length(missing_cols) > 0) {
      stop("`group_col` contains columns not found in `df_data`: ",
           paste(missing_cols, collapse = ", "), call. = FALSE)
    }
  }
  if (!(p_val_col %in% colnames(df_data))) {
    stop("`p_val_col` '", p_val_col, "' not found in `df_data`.",
         call. = FALSE)
  }

  # 2. Build the grouping key and split.
  df_list <- if (!is.null(group_col)) {
    # interaction() produces a single factor whose levels are the unique
    # combinations of the supplied columns. Works uniformly for one or
    # more group_col entries; drop = TRUE removes empty level combinations.
    grp_key <- interaction(df_data[, group_col, drop = FALSE],
                           drop = TRUE, sep = ".")
    split(df_data, grp_key)
  } else {
    list(all_data = df_data)
  }

  # 3. Apply the BH correction to each family in the list.
  processed_list <- lapply(df_list, function(df_chunk) {
    if (nrow(df_chunk) > 1) {
      df_chunk$pBH <- p.adjust(df_chunk[[p_val_col]], method = "BH")
    } else {
      # With a single test there is nothing to correct.
      df_chunk$pBH <- df_chunk[[p_val_col]]
    }
    df_chunk
  })

  # 4. Combine the processed chunks back into one data frame.
  df_data_out <- do.call(rbind, processed_list)
  rownames(df_data_out) <- NULL

  df_data_out
}

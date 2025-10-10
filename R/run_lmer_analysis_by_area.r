#' @title Run Linear Mixed-Effects Models Iteratively by Area
#'
#' @description
#' Iterates through unique values in an 'area' column (and optionally a
#' 'group' column) of a dataframe, running an LME model on each subset of data.
#'
#' @details
#' This function provides two main execution paths:
#' 1. If `group_col` is `NULL`, it runs one model for each unique area specified in `area_col`.
#' 2. If a `group_col` is provided, it runs a model for each unique combination of area and group.
#'
#' It's a wrapper around `lme4::lmer` designed to automate repetitive modeling tasks.
#'
#' @param input_df The dataframe containing the data to be modeled.
#' @param area_col The name of the column specifying the brain area (e.g., "ROI").
#' @param group_col The column to group by. If `NULL`, the analysis runs per area on the whole dataset. Default: `NULL`.
#' @param formula_str A string containing the formula for the lmer model (e.g., `"Value ~ Condition + (1|Subject)"`).
#' @param run_contrasts A logical value indicating whether to compute emmeans contrasts. Default: `FALSE`.
#' @param contrast_variable The name of the variable to use for pairwise contrasts (e.g., "Group" or "cond").
#' @param sinkFileName The name of the sink file for the output. If `NULL`, output prints to the console. Default: `NULL`.
#' @param sinkAppend A boolean specifying whether the sink file should append to an existing file. Default: `TRUE`.
#'
#' @return A list containing two dataframes: `model_results` and `contrast_results`.
#'
#' @seealso \code{\link[lme4]{lmer}}, \code{\link[emmeans]{emmeans}}, \code{\link[lmerTest]{anova}}
#'
#' @importFrom lme4 lmer
#' @importFrom emmeans emmeans
#'
#' @export
run_lmer_analysis_by_area <- function(input_df, area_col, group_col = NULL, formula_str, run_contrasts = FALSE, contrast_variable = NULL, sinkFileName = NULL, sinkAppend = TRUE) {

  df_lmr_out <- NULL
  df_contrast_out <- NULL
  areas <- unique(input_df[[area_col]])

  if(!is.null(sinkFileName)) {
    # sink() is a base R function, no namespace needed.
    sink(sinkFileName, append = sinkAppend)
  }

  if (!is.null(group_col)) {
    # --- PATH 1: A group_col IS provided ---
    groups <- unique(input_df[[group_col]])

    for (grp in groups) {
      for (ar in areas) {
        df_tmp <- input_df[input_df[[area_col]] == ar & input_df[[group_col]] == grp, ]
        if(nrow(df_tmp) > 0) {
            # lmer() is imported via @importFrom lme4 lmer
            lmr_tmp <- lmer(formula_str, data = df_tmp)
            df_lmr_summary <- as.data.frame(summary(lmr_tmp)$coefficient)

            df_lmr_tmp <- data.frame(Predictor = rownames(df_lmr_summary), df_lmr_summary)
            colnames(df_lmr_tmp) <- c("Predictor", "Beta", "SE", "df", "tVal", "pVal")
            df_lmr_tmp[[group_col]] <- grp
            df_lmr_tmp[[area_col]] <- ar
            df_lmr_out <- rbind(df_lmr_out, df_lmr_tmp)

            if (run_contrasts) {
              if (is.null(contrast_variable)) stop("Error: Must provide a 'contrast_variable' when 'run_contrasts' is TRUE.")
              contrast_formula <- as.formula(paste("pairwise ~", contrast_variable))
              # emmeans() is imported via @importFrom emmeans emmeans
              emm_pairwise <- emmeans(lmr_tmp, contrast_formula, adjust = "none")
              em_tmp <- as.data.frame(emm_pairwise$contrasts)
              colnames(em_tmp) <- c("Predictor", "Beta", "SE", "df", "tVal", "pVal")

              em_tmp[[group_col]] <- grp
              em_tmp[[area_col]] <- ar
              df_contrast_out <- rbind(df_contrast_out, em_tmp)
            }

            cat("\n\n ************* Group = ", grp, " Area = ", ar, "****************\n\n")
            # Using lmerTest::anova to get p-values for the lmer model.
            print(lmerTest::anova(lmr_tmp))
            print(summary(lmr_tmp))
            if(run_contrasts) {
              cat("\n\n ************* Contrast \n")
              print(emm_pairwise)
            }
        }
      }
    }
  } else {
    # --- PATH 2: group_col IS NULL ---
    for (ar in areas) {
      df_tmp <- input_df[input_df[[area_col]] == ar, ]
      if(nrow(df_tmp) > 0) {
          # lmer() is imported via @importFrom lme4 lmer
          lmr_tmp <- lmer(formula_str, data = df_tmp)
          df_lmr_summary <- as.data.frame(summary(lmr_tmp)$coefficient)

          df_lmr_tmp <- data.frame(Predictor = rownames(df_lmr_summary), df_lmr_summary)
          colnames(df_lmr_tmp) <- c("Predictor", "Beta", "SE", "df", "tVal", "pVal")
          df_lmr_tmp[[area_col]] <- ar
          df_lmr_out <- rbind(df_lmr_out, df_lmr_tmp)

          if (run_contrasts) {
            if (is.null(contrast_variable)) stop("Error: Must provide a 'contrast_variable' when 'run_contrasts' is TRUE.")
            contrast_formula <- as.formula(paste("pairwise ~", contrast_variable))
            # emmeans() is imported via @importFrom emmeans emmeans
            emm_pairwise <- emmeans(lmr_tmp, contrast_formula, adjust = "none")
            em_tmp <- as.data.frame(emm_pairwise$contrasts)
            colnames(em_tmp) <- c("Predictor", "Beta", "SE", "df", "tVal", "pVal")

            em_tmp[[area_col]] <- ar
            df_contrast_out <- rbind(df_contrast_out, em_tmp)
          }

          cat("\n\n ************* Area = ", ar, "****************\n\n")
          # Using lmerTest::anova to get p-values for the lmer model.
          print(lmerTest::anova(lmr_tmp))
          print(summary(lmr_tmp))
          if(run_contrasts) {
            cat("\n\n ************* Contrast \n")
            print(emm_pairwise)
          }
      }
    }
  }

  if(!is.null(sinkFileName)) {
    sink(NULL)
  }

  return(list(model_results = df_lmr_out, contrast_results = df_contrast_out))
}
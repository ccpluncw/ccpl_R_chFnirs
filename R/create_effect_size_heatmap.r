#' @title Create an effect-size heatmap, optionally split by group
#'
#' @description
#' Draws one heatmap of effect-size values (e.g., Cohen's d, partial
#' eta-squared, slope magnitudes) per level of an optional grouping
#' column, with predictors on the y-axis and ROIs/areas on the x-axis.
#' When the effect-size column contains any negative values, a diverging
#' palette centered at zero is used; otherwise a single-tailed
#' sequential palette is used. Cell values are labeled in-tile and each
#' group is saved to its own PDF.
#'
#' This fills the gap between [create_pvalue_heatmap()] and the beta /
#' omnibus heatmap functions for the case where the quantity being
#' visualized is an effect size rather than a p-value or raw coefficient.
#'
#' @param results_df A long-format data.frame containing at least the
#'   columns named by `area_col`, `predictor`, `es_col`, and (optionally)
#'   `group_col`.
#' @param area_col Name of the column containing ROI/area labels. These
#'   become the x-axis of the heatmap.
#' @param es_col Name of the column containing the effect-size values to
#'   plot (e.g., `"cohens_d"`, `"partial_eta_sq"`).
#' @param group_col (Optional) Name of a column to split on before
#'   plotting. One heatmap is produced per unique value. If `NULL`, a
#'   single heatmap of the entire dataset is produced. Default `NULL`.
#' @param title (Optional) Plot title. If `NULL`, defaults to
#'   `paste("Effect Size:", es_col)`. Default `NULL`.
#' @param low_color,mid_color,high_color Colors for the diverging
#'   palette used when `es_col` contains negative values. Defaults
#'   `"blue"`, `"white"`, `"red"`.
#' @param seq_low_color,seq_high_color Colors for the sequential palette
#'   used when `es_col` contains only non-negative values. Defaults
#'   `"white"`, `"darkorange"`.
#' @param neg_label_digits,pos_label_digits Number of decimal digits
#'   shown in-cell for the diverging and sequential palettes
#'   respectively. Defaults 2 and 3.
#' @param save_dir Directory to write PDFs to. Default `"."`.
#' @param file_prefix File-name prefix. The final file name is
#'   `<file_prefix>_<group>_<es_col>_Heatmap.pdf` when a `group_col` is
#'   supplied and `<file_prefix>_<es_col>_Heatmap.pdf` otherwise.
#'   Default `"EffectSize"`.
#' @param fig_width,fig_height Figure dimensions in inches. Defaults
#'   11 and 8.5.
#' @param landscape If `TRUE`, `fig_width` and `fig_height` are swapped
#'   before saving. Default `FALSE`.
#'
#' @return Invisibly returns `NULL`. Side effect: writes one PDF per
#'   group level (or a single PDF if `group_col` is `NULL`) into
#'   `save_dir`, and prints a message for each file written.
#'
#' @examples
#' \dontrun{
#' # results_df from run_activation_tests() or run_gated_posthoc()
#' create_effect_size_heatmap(
#'   results_df  = my_results,
#'   area_col    = "area",
#'   es_col      = "cohens_d",
#'   group_col   = "type",
#'   save_dir    = "figures/effect_sizes",
#'   file_prefix = "Activation"
#' )
#' }
#'
#' @importFrom dplyr relocate
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes geom_tile geom_text scale_fill_gradient
#' @importFrom ggplot2 scale_fill_gradient2 labs theme_minimal theme
#' @importFrom ggplot2 element_text element_blank ggsave
#' @export
create_effect_size_heatmap <- function(
    results_df,
    area_col,
    es_col,
    group_col      = NULL,
    title          = NULL,
    low_color      = "blue",
    mid_color      = "white",
    high_color     = "red",
    seq_low_color  = "white",
    seq_high_color = "darkorange",
    neg_label_digits = 2,
    pos_label_digits = 3,
    save_dir       = ".",
    file_prefix    = "EffectSize",
    fig_width      = 11,
    fig_height     = 8.5,
    landscape      = FALSE
) {
  if (!area_col %in% colnames(results_df)) {
    stop("area_col '", area_col, "' not found in results_df.")
  }
  if (!es_col %in% colnames(results_df)) {
    stop("es_col '", es_col, "' not found in results_df.")
  }
  if (!"predictor" %in% colnames(results_df)) {
    stop("results_df must contain a 'predictor' column.")
  }
  if (!is.null(group_col) && !group_col %in% colnames(results_df)) {
    stop("group_col '", group_col, "' not found in results_df.")
  }

  if (landscape) {
    tmp <- fig_width
    fig_width <- fig_height
    fig_height <- tmp
  }

  df_list <- if (!is.null(group_col)) {
    split(results_df, results_df[[group_col]])
  } else {
    list(all = results_df)
  }

  plot_title <- if (is.null(title)) paste("Effect Size:", es_col) else title

  for (grp_name in names(df_list)) {
    df_tmp <- df_list[[grp_name]]
    areas      <- unique(df_tmp[[area_col]])
    predictors <- unique(df_tmp$predictor)

    # Build wide predictor x area table of effect sizes
    mat <- matrix(NA_real_, nrow = length(predictors), ncol = length(areas))
    df_wide <- as.data.frame(mat)
    colnames(df_wide) <- areas
    df_wide$predictor <- predictors

    for (ar in areas) {
      for (pr in predictors) {
        row <- df_tmp[df_tmp[[area_col]] == ar & df_tmp$predictor == pr, ]
        if (nrow(row) > 0 && !is.na(row[[es_col]][1])) {
          df_wide[df_wide$predictor == pr, ar] <- row[[es_col]][1]
        }
      }
    }

    df_wide <- dplyr::relocate(df_wide, "predictor")
    df_long <- reshape2::melt(df_wide, id.vars = "predictor")
    colnames(df_long) <- c("y_label", "x_label", "value")
    df_long$x_label <- factor(df_long$x_label, levels = areas)

    has_negative <- any(df_long$value < 0, na.rm = TRUE)

    if (has_negative) {
      fmt <- paste0("%.", neg_label_digits, "f")
      p <- ggplot(df_long, aes(x = x_label, y = y_label, fill = value)) +
        geom_tile(color = "grey50") +
        scale_fill_gradient2(
          low = low_color, mid = mid_color, high = high_color,
          midpoint = 0, na.value = "grey80", name = es_col
        ) +
        geom_text(
          aes(label = ifelse(!is.na(value), sprintf(fmt, value), "")),
          size = 3.5, color = "black"
        )
    } else {
      fmt <- paste0("%.", pos_label_digits, "f")
      p <- ggplot(df_long, aes(x = x_label, y = y_label, fill = value)) +
        geom_tile(color = "grey50") +
        scale_fill_gradient(
          low = seq_low_color, high = seq_high_color,
          na.value = "grey80", name = es_col
        ) +
        geom_text(
          aes(label = ifelse(!is.na(value), sprintf(fmt, value), "")),
          size = 3.5, color = "black"
        )
    }

    p <- p +
      labs(title = plot_title, x = "", y = "") +
      theme_minimal() +
      theme(
        axis.text.x      = element_text(angle = 45, vjust = 1, hjust = 1, size = 12),
        axis.text.y      = element_text(size = 12),
        plot.title       = element_text(hjust = 0.5, size = 16, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )

    fname <- if (!is.null(group_col)) {
      paste0(file_prefix, "_", grp_name, "_", es_col, "_Heatmap.pdf")
    } else {
      paste0(file_prefix, "_", es_col, "_Heatmap.pdf")
    }
    out_path <- file.path(save_dir, fname)
    ggsave(filename = out_path, plot = p,
           width = fig_width, height = fig_height, units = "in")
    message("Effect size heatmap saved to: ", out_path)
  }

  invisible(NULL)
}

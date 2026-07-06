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
#' @param p_col (Optional) Name of a column containing p-values. When
#'   non-NULL, cells whose p-value is at or below `significance_threshold`
#'   are marked with an asterisk. When `NULL` (default), no significance
#'   markers are drawn. Matches the `significance_threshold` /
#'   `asterisk_size` convention used by [create_omnibus_heatmap()] and
#'   [create_posthoc_heatmap()].
#' @param significance_threshold (Optional) P-value threshold for
#'   asterisks. Only used when `p_col` is non-NULL. Default: 0.05.
#' @param asterisk_size (Optional) Size of asterisks in mm. Default: 6.
#' @param show_values Logical. If `TRUE` (default), the numeric
#'   effect-size value is printed in each cell. Combined with `p_col`,
#'   the four label modes are: values only (`show_values = TRUE`,
#'   `p_col = NULL`), both (`TRUE`, non-NULL), asterisks only (`FALSE`,
#'   non-NULL), no text (`FALSE`, `NULL`).
#' @param area_order Optional character vector specifying the desired
#'   x-axis order of `area_col` levels. If `NULL` (default), areas
#'   appear in the order they occur in `results_df`. To sort
#'   anatomically, precompute the order with [sort_rois_anatomically()]
#'   (for ROI data) or [sort_fnirs_channels_by_name()] (for channel
#'   data) and pass the result here. Any sort function with signature
#'   `function(names, ...) -> ordered names` can be used.
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
#' @param split_by_type Logical. When `TRUE` (default) and `group_col`
#'   is non-NULL, each group's PDF is written into
#'   `<save_dir>/<tolower(group_level)>/` instead of flat in
#'   `save_dir`. Ignored when `group_col` is `NULL`.
#' @param show_legend Logical. If `TRUE` (default), the fill-scale
#'   legend is rendered. If `FALSE`, the legend is suppressed --
#'   useful for publication panels that share a single legend
#'   exported via [create_heatmap_legend()].
#' @param value_limits Numeric length-2 vector giving the colorbar's
#'   numeric domain. For diverging palettes (effect-size data
#'   contains negatives), pass e.g. `c(-1, 1)` for Cohen's d -- the
#'   positive colorbar then spans `c(0, 1)` and the negative spans
#'   `c(-1, 0)`. For sequential palettes (all non-negative), pass
#'   e.g. `c(0, 0.2)` for partial η². When `NULL` (default), limits
#'   are auto-computed from the data range. There is no universal
#'   default -- pick limits that match the effect-size scale you
#'   are plotting.
#' @param legend_breaks Numeric vector of values to label on the
#'   colorbar (e.g., `c(-0.8, -0.5, -0.2, 0.2, 0.5, 0.8)` for
#'   Cohen's d, or `c(0.01, 0.06, 0.14)` for partial η²). For
#'   diverging palettes the function automatically routes
#'   non-negative values to the positive colorbar and negative
#'   values to the negative colorbar. When `NULL` (default), ggplot
#'   chooses break positions automatically.
#' @param pos_legend_title Title for the positive colorbar
#'   (diverging case) or for the sole colorbar (sequential case).
#'   When `NULL` (default), uses `es_col`.
#' @param neg_legend_title Title for the negative colorbar.
#'   Diverging case only. When `NULL` (default), uses `es_col`.
#'
#' @return Invisibly returns `NULL`. Side effect: writes one PDF per
#'   group level (or a single PDF if `group_col` is `NULL`) into
#'   `save_dir`, and prints a message for each file written.
#'
#' @examples
#' \dontrun{
#' # Values only (default) — no significance markers
#' create_effect_size_heatmap(
#'   results_df  = my_results,
#'   area_col    = "area",
#'   es_col      = "cohens_d",
#'   group_col   = "type",
#'   save_dir    = "figures/effect_sizes",
#'   file_prefix = "Activation"
#' )
#'
#' # Both values and asterisks for significant cells
#' create_effect_size_heatmap(
#'   results_df  = my_results,
#'   area_col    = "area",
#'   es_col      = "cohens_d",
#'   p_col       = "p_bh",
#'   group_col   = "type",
#'   save_dir    = "figures/effect_sizes",
#'   file_prefix = "Activation"
#' )
#'
#' # Asterisks only — no numeric labels
#' create_effect_size_heatmap(
#'   results_df  = my_results,
#'   area_col    = "area",
#'   es_col      = "cohens_d",
#'   p_col       = "p_bh",
#'   show_values = FALSE,
#'   group_col   = "type",
#'   save_dir    = "figures/effect_sizes",
#'   file_prefix = "Activation"
#' )
#' }
#'
#' @importFrom dplyr relocate
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes geom_tile geom_text scale_fill_gradient
#' @importFrom ggplot2 scale_x_discrete guide_colorbar labs theme_minimal theme waiver
#' @importFrom ggplot2 element_text element_blank ggsave
#' @importFrom ggnewscale new_scale_fill
#' @export
create_effect_size_heatmap <- function(
    results_df,
    area_col,
    es_col,
    p_col          = NULL,
    significance_threshold = 0.05,
    asterisk_size  = 6,
    show_values    = TRUE,
    area_order     = NULL,
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
    landscape      = FALSE,
    split_by_type  = TRUE,
    show_legend    = TRUE,
    value_limits   = NULL,
    legend_breaks  = NULL,
    pos_legend_title = NULL,
    neg_legend_title = NULL
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
  if (!is.null(p_col) && !p_col %in% colnames(results_df)) {
    stop("p_col '", p_col, "' not found in results_df.")
  }

  # Apply caller-supplied area ordering, if any
  if (!is.null(area_order)) {
    missing_levels <- setdiff(unique(results_df[[area_col]]), area_order)
    if (length(missing_levels) > 0) {
      stop("area_order is missing ", length(missing_levels),
           " level(s) present in `", area_col, "`: ",
           paste(missing_levels, collapse = ", "), call. = FALSE)
    }
    results_df$Area_ordered <- factor(results_df[[area_col]],
                                      levels = area_order)
    results_df <- results_df[order(results_df$Area_ordered), ]
    results_df$Area_ordered <- NULL
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
  use_sig <- !is.null(p_col) && significance_threshold > 0

  for (grp_name in names(df_list)) {
    df_tmp <- df_list[[grp_name]]
    areas      <- unique(df_tmp[[area_col]])
    predictors <- unique(df_tmp$predictor)

    # Build wide predictor x area table of effect sizes
    mat <- matrix(NA_real_, nrow = length(predictors), ncol = length(areas))
    df_wide <- as.data.frame(mat)
    colnames(df_wide) <- areas
    df_wide$predictor <- predictors

    if (use_sig) {
      p_mat <- matrix(NA_real_, nrow = length(predictors), ncol = length(areas))
      p_wide <- as.data.frame(p_mat)
      colnames(p_wide) <- areas
      p_wide$predictor <- predictors
    }

    for (ar in areas) {
      for (pr in predictors) {
        row <- df_tmp[df_tmp[[area_col]] == ar & df_tmp$predictor == pr, ]
        if (nrow(row) > 0 && !is.na(row[[es_col]][1])) {
          df_wide[df_wide$predictor == pr, ar] <- row[[es_col]][1]
        }
        if (use_sig && nrow(row) > 0 && !is.na(row[[p_col]][1])) {
          p_wide[p_wide$predictor == pr, ar] <- row[[p_col]][1]
        }
      }
    }

    df_wide <- dplyr::relocate(df_wide, "predictor")
    df_long <- reshape2::melt(df_wide, id.vars = "predictor")
    colnames(df_long) <- c("y_label", "x_label", "value")
    df_long$x_label <- factor(df_long$x_label, levels = areas)

    # -- build cell labels --
    has_negative <- any(df_long$value < 0, na.rm = TRUE)
    fmt <- if (has_negative) {
      paste0("%.", neg_label_digits, "f")
    } else {
      paste0("%.", pos_label_digits, "f")
    }

    if (use_sig) {
      p_wide <- dplyr::relocate(p_wide, "predictor")
      p_long <- reshape2::melt(p_wide, id.vars = "predictor")
      df_long$p_val <- p_long$value
      df_long$is_sig <- !is.na(df_long$p_val) &
                         df_long$p_val <= significance_threshold
    } else {
      df_long$is_sig <- FALSE
    }

    val_str <- if (show_values) sprintf(fmt, df_long$value) else
                 rep("", nrow(df_long))
    sig_str <- ifelse(df_long$is_sig, "*", "")
    df_long$cell_label <- ifelse(
      is.na(df_long$value), "",
      paste0(val_str, sig_str)
    )

    # -- legend titles default to es_col --
    pos_title <- if (!is.null(pos_legend_title)) pos_legend_title else es_col
    neg_title <- if (!is.null(neg_legend_title)) neg_legend_title else es_col

    # -- partition legend_breaks by sign for diverging case --
    if (!is.null(legend_breaks)) {
      breaks_pos <- legend_breaks[legend_breaks >= 0]
      breaks_neg <- legend_breaks[legend_breaks <  0]
    } else {
      breaks_pos <- waiver()
      breaks_neg <- waiver()
    }

    # -- build plot --
    if (has_negative) {
      # Diverging case: split into positive and negative subsets, render
      # each through its own scale_fill_gradient joined via ggnewscale,
      # mirroring create_pvalue_heatmap / create_posthoc_heatmap.
      df_na  <- df_long[is.na(df_long$value), ]
      df_pos <- df_long[!is.na(df_long$value) & df_long$value >= 0, ]
      df_neg <- df_long[!is.na(df_long$value) & df_long$value <  0, ]

      # Resolve colorbar limits. When value_limits is NULL, use the
      # observed data range (matches the prior scale_fill_gradient2
      # auto-limits behaviour).
      if (!is.null(value_limits)) {
        lim_pos <- c(0, value_limits[2])
        lim_neg <- c(value_limits[1], 0)
      } else {
        max_abs <- max(abs(df_long$value), na.rm = TRUE)
        lim_pos <- c(0, max_abs)
        lim_neg <- c(-max_abs, 0)
      }

      p <- ggplot(mapping = aes(x = x_label, y = y_label))

      if (nrow(df_na) > 0) {
        p <- p + geom_tile(data = df_na, fill = "grey80", color = "grey50")
      }

      if (nrow(df_pos) > 0) {
        p <- p +
          geom_tile(data = df_pos, aes(fill = value), color = "grey50") +
          scale_fill_gradient(
            name   = pos_title,
            low    = mid_color, high = high_color,
            limits = lim_pos,
            breaks = breaks_pos,
            labels = if (inherits(breaks_pos, "waiver")) waiver() else breaks_pos,
            na.value = "grey80",
            guide  = guide_colorbar(order = 1)
          ) +
          new_scale_fill()
      }

      if (nrow(df_neg) > 0) {
        p <- p +
          geom_tile(data = df_neg, aes(fill = value), color = "grey50") +
          scale_fill_gradient(
            name   = neg_title,
            low    = low_color, high = mid_color,
            limits = lim_neg,
            breaks = breaks_neg,
            labels = if (inherits(breaks_neg, "waiver")) waiver() else breaks_neg,
            na.value = "grey80",
            guide  = guide_colorbar(order = 2)
          )
      }
    } else {
      # Sequential case: single scale, white -> seq_high_color.
      seq_limits <- if (!is.null(value_limits)) value_limits else
                      c(0, max(df_long$value, na.rm = TRUE))

      p <- ggplot(df_long, aes(x = x_label, y = y_label, fill = value)) +
        geom_tile(color = "grey50") +
        scale_fill_gradient(
          name   = pos_title,
          low    = seq_low_color, high = seq_high_color,
          limits = seq_limits,
          breaks = breaks_pos,
          labels = if (inherits(breaks_pos, "waiver")) waiver() else breaks_pos,
          na.value = "grey80"
        )
    }

    # Enforce x-axis order. Required because the diverging path uses a
    # base ggplot with no data attached (line ~314), feeding three
    # geom_tile layers with subsets (df_na, df_pos, df_neg) that may
    # have differing levels populated. Without an explicit
    # scale_x_discrete, ggplot can rebuild the discrete axis level
    # order from the first non-empty geom's data and append later
    # levels out-of-sequence — visible as misplaced ROIs (e.g.,
    # temporal_mid_l between precentral_l and precentral_r) when the
    # input contains NA cells. Anchoring to `areas` (which is the
    # caller-supplied area_order, factor-validated at line ~252) makes
    # the axis order deterministic in both diverging and sequential
    # paths. Mirrors the scale_x_discrete(limits = ...) safety net
    # already used in create_pvalue_heatmap and create_posthoc_heatmap.
    p <- p + ggplot2::scale_x_discrete(limits = areas)

    # Cell labels. Pass df_long explicitly because the diverging path's
    # base ggplot has no data attached (data is supplied per geom_tile)
    # and geom_text would otherwise inherit nothing.
    if (show_values) {
      p <- p + geom_text(
        data = df_long, aes(label = cell_label),
        size = 3.5, color = "black"
      )
    } else if (use_sig) {
      p <- p + geom_text(
        data = df_long, aes(label = cell_label),
        size = asterisk_size, color = "black"
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
        panel.grid.minor = element_blank(),
        legend.position  = if (show_legend) "right" else "none"
      )

    fname <- if (!is.null(group_col)) {
      paste0(file_prefix, "_", grp_name, "_", es_col, "_Heatmap.pdf")
    } else {
      paste0(file_prefix, "_", es_col, "_Heatmap.pdf")
    }
    out_subdir <- if (split_by_type && !is.null(group_col)) {
      file.path(save_dir, tolower(as.character(grp_name)))
    } else {
      save_dir
    }
    if (!dir.exists(out_subdir)) {
      dir.create(out_subdir, recursive = TRUE, showWarnings = FALSE)
    }
    out_path <- file.path(out_subdir, fname)
    ggsave(filename = out_path, plot = p,
           width = fig_width, height = fig_height, units = "in")
    message("Effect size heatmap saved to: ", out_path)
  }

  invisible(NULL)
}

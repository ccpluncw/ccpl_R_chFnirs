#' @title Save a standalone heatmap legend as a PDF
#'
#' @description
#' Renders just the legend that the chFnirs heatmap functions
#' produce, with no panel or axes, and saves it to its own PDF.
#' Designed for the publication-figure workflow where multiple
#' heatmaps are arranged into a panel and share a single shared
#' legend, rather than each carrying its own embedded copy.
#'
#' @details
#' The legend is composed of up to three independent pieces, each
#' toggled on or off:
#'
#' - **Positive colorbar** (`show_pos`): a gradient from `mid_color`
#'   to `pos_color`, labeled at `pos_legend_breaks`. Used by the
#'   activation and post-hoc heatmaps for positive effects, and by
#'   the omnibus heatmap for unsigned p-values (with `pos_color`
#'   set to e.g. `"darkgreen"`).
#'
#' - **Negative colorbar** (`show_neg`): a gradient from `neg_color`
#'   to `mid_color`, labeled at `neg_legend_breaks`. Used by the
#'   activation and post-hoc heatmaps for negative effects.
#'
#' - **Not-tested swatch** (`show_not_tested`): a single colored
#'   square labeled by `not_tested_label`. Used by the post-hoc
#'   heatmap to indicate ROIs whose omnibus test was not significant.
#'
#' At least one of the three flags must be `TRUE`. The function
#' builds a minimal `ggplot2` plot with the requested fill scales,
#' uses [cowplot::get_legend()] to extract just the legend grob,
#' and saves the result. This guarantees the standalone legend
#' renders identically to the embedded version produced by the
#' heatmap functions when called with matching parameters.
#'
#' Common patterns:
#'
#' ```r
#' # Activation legend (positive + negative, no not-tested)
#' create_heatmap_legend(show_pos = TRUE, show_neg = TRUE,
#'                       show_not_tested = FALSE)
#'
#' # Omnibus legend (single green colorbar)
#' create_heatmap_legend(show_pos = TRUE, show_neg = FALSE,
#'                       show_not_tested = FALSE,
#'                       pos_color = "darkgreen",
#'                       pos_legend_title = "p-value",
#'                       pos_legend_breaks = c(0.001, 0.01, 0.05, 0.10, 0.50))
#'
#' # Post-hoc legend (full stack: positive + negative + not-tested)
#' create_heatmap_legend(show_pos = TRUE, show_neg = TRUE,
#'                       show_not_tested = TRUE)
#'
#' # Effect-size legend, Cohen's d (diverging, no transformation)
#' create_heatmap_legend(show_pos = TRUE, show_neg = TRUE,
#'                       scale_type = "none",
#'                       pos_value_limits = c(0, 1),
#'                       neg_value_limits = c(-1, 0),
#'                       pos_legend_breaks = c(0.2, 0.5, 0.8),
#'                       neg_legend_breaks = c(-0.8, -0.5, -0.2),
#'                       pos_legend_title = "Cohen's d (+)",
#'                       neg_legend_title = "Cohen's d (-)")
#'
#' # Effect-size legend, partial eta-squared (sequential)
#' create_heatmap_legend(show_pos = TRUE, show_neg = FALSE,
#'                       scale_type = "none",
#'                       pos_value_limits = c(0, 0.2),
#'                       pos_legend_breaks = c(0.01, 0.06, 0.14),
#'                       pos_legend_title = "partial eta-sq",
#'                       pos_color = "darkorange")
#' ```
#'
#' @param show_pos Logical. Include the positive colorbar.
#'   Default `TRUE`.
#' @param show_neg Logical. Include the negative colorbar.
#'   Default `TRUE`.
#' @param show_not_tested Logical. Include the not-tested swatch.
#'   Default `FALSE`.
#' @param scale_type One of `"linear"` (default), `"log"`, or
#'   `"none"`. Must match the `scale_type` used when rendering the
#'   corresponding heatmaps; otherwise the colorbar break positions
#'   will not line up with the heatmap cells.
#'   - `"linear"` and `"log"` apply the p-value transformation
#'     `1 - |p|^log_power` (or `1 - |p|`) to map breaks into the
#'     `[0, 1]` color space used by `create_pvalue_heatmap`,
#'     `create_posthoc_heatmap`, and `create_omnibus_heatmap`.
#'   - `"none"` applies no transformation: `pos_legend_breaks` and
#'     `neg_legend_breaks` are used as the colorbar tick positions
#'     directly. Use this for `create_effect_size_heatmap` where the
#'     plotted value IS the color position (e.g., Cohen's d, partial
#'     η²). Pair with `pos_value_limits` / `neg_value_limits`.
#' @param log_power Numeric. Power applied when `scale_type = "log"`.
#'   Default `0.2`. Must match the heatmap's `log_power`.
#' @param pos_legend_breaks Numeric vector of values shown on the
#'   positive colorbar. For p-value modes, these are the p-values
#'   to label (default `c(0.001, 0.01, 0.05)`). For
#'   `scale_type = "none"`, these are the effect-size values to
#'   label (e.g., `c(0.2, 0.5, 0.8)` for Cohen's d).
#' @param neg_legend_breaks Numeric vector of (negative) values shown
#'   on the negative colorbar. Default `c(-0.05, -0.01, -0.001)`.
#' @param pos_value_limits Numeric length-2 vector giving the positive
#'   colorbar's domain. Default `c(0, 1)` (correct for `"linear"` /
#'   `"log"` p-value modes, where the transformation always lands in
#'   `[0, 1]`). For `scale_type = "none"`, override to match the
#'   heatmap's `value_limits` -- e.g. `c(0, 1)` for Cohen's d or
#'   `c(0, 0.2)` for partial η².
#' @param neg_value_limits Numeric length-2 vector giving the negative
#'   colorbar's domain. Default `c(-1, 0)`. Override for
#'   `scale_type = "none"` to match the heatmap (e.g. `c(-1, 0)` for
#'   Cohen's d).
#' @param pos_legend_title Title for the positive colorbar. Default
#'   `"Positive Effect"`. For omnibus legends, pass `"p-value"`.
#' @param neg_legend_title Title for the negative colorbar. Default
#'   `"Negative Effect"`.
#' @param not_tested_label Label shown next to the not-tested swatch.
#'   Default `"Not tested\n(omnibus n.s.)"` (newline matches the
#'   embedded post-hoc heatmap legend exactly).
#' @param pos_color Color for the high (most-significant) end of the
#'   positive colorbar. Default `"red"`. Pass `"darkgreen"` for the
#'   omnibus convention.
#' @param neg_color Color for the low (most-significant negative) end
#'   of the negative colorbar. Default `"blue"`.
#' @param mid_color Color for the white midpoint of both colorbars.
#'   Default `"white"`.
#' @param not_tested_color Color for the not-tested swatch. Default
#'   `"grey80"`.
#' @param fig_width Width of the output PDF in inches. Default `3`.
#' @param fig_height Height of the output PDF in inches. Default `5`
#'   (enough to fit the typical full three-piece legend stack).
#' @param save_dir Directory to write the PDF to. Created if it does
#'   not exist. Default `"."`.
#' @param file_name Output PDF file name. Default
#'   `"Heatmap_Legend.pdf"`.
#'
#' @return Invisibly returns `NULL`. Side effect: writes one PDF and
#'   prints a confirmation message.
#'
#' @examples
#' \dontrun{
#' # In a publication-figure script: render N heatmaps without
#' # legends, then export one shared legend.
#'
#' fig_dir <- "paper/figures"
#'
#' for (chrom in c("hbo", "cbsi")) {
#'   df_sub <- subset(df.posthoc, type == chrom)
#'   create_posthoc_heatmap(
#'     df_sub, ..., show_legend = FALSE,
#'     fig_height = heatmap_height(df_sub, show_legend = FALSE),
#'     save_dir = fig_dir,
#'     file_name = sprintf("Panel_%s.pdf", chrom)
#'   )
#' }
#'
#' create_heatmap_legend(
#'   show_pos = TRUE, show_neg = TRUE, show_not_tested = TRUE,
#'   save_dir = fig_dir,
#'   file_name = "Shared_Legend_PostHoc.pdf"
#' )
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradient
#' @importFrom ggplot2 scale_fill_manual guide_colorbar guide_legend
#' @importFrom ggplot2 theme_void theme element_text ggsave
#' @importFrom ggnewscale new_scale_fill
#' @importFrom cowplot get_plot_component ggdraw draw_grob
#' @importFrom stats setNames
#' @export
create_heatmap_legend <- function(
    show_pos          = TRUE,
    show_neg          = TRUE,
    show_not_tested   = FALSE,
    scale_type        = "linear",
    log_power         = 0.2,
    pos_legend_breaks = c(0.001, 0.01, 0.05),
    neg_legend_breaks = c(-0.05, -0.01, -0.001),
    pos_value_limits  = c(0, 1),
    neg_value_limits  = c(-1, 0),
    pos_legend_title  = "Positive Effect",
    neg_legend_title  = "Negative Effect",
    not_tested_label  = "Not tested\n(omnibus n.s.)",
    pos_color         = "red",
    neg_color         = "blue",
    mid_color         = "white",
    not_tested_color  = "grey80",
    fig_width         = 3,
    fig_height        = 5,
    save_dir          = ".",
    file_name         = "Heatmap_Legend.pdf"
) {
  if (!any(show_pos, show_neg, show_not_tested)) {
    stop("At least one of `show_pos`, `show_neg`, or `show_not_tested` ",
         "must be TRUE.", call. = FALSE)
  }
  if (!scale_type %in% c("linear", "log", "none")) {
    stop("`scale_type` must be \"linear\", \"log\", or \"none\".",
         call. = FALSE)
  }

  # For "linear" / "log" (p-value modes), apply the same break
  # transformation as create_pvalue_heatmap / create_posthoc_heatmap so
  # the standalone legend's ticks line up with the embedded heatmap
  # legend's ticks. For "none" (effect-size mode), the plotted value
  # IS the color position; breaks are used as-is.
  if (scale_type == "log") {
    pos_breaks_t <- 1 - abs(pos_legend_breaks)^log_power
    neg_breaks_t <- -(1 - abs(neg_legend_breaks)^log_power)
    pos_labels   <- pos_legend_breaks
    neg_labels   <- abs(neg_legend_breaks)
  } else if (scale_type == "linear") {
    pos_breaks_t <- 1 - abs(pos_legend_breaks)
    neg_breaks_t <- -(1 - abs(neg_legend_breaks))
    pos_labels   <- pos_legend_breaks
    neg_labels   <- abs(neg_legend_breaks)
  } else {
    pos_breaks_t <- pos_legend_breaks
    neg_breaks_t <- neg_legend_breaks
    pos_labels   <- pos_legend_breaks
    neg_labels   <- neg_legend_breaks
  }

  # Each enabled scale needs at least one row of data inside its
  # limits to register; the tiles themselves are discarded with the
  # rest of the panel during legend extraction. Use the midpoint of
  # the user-supplied limits so the dummy is always in range.
  df_pos <- data.frame(x = 1, y = 1, value = mean(pos_value_limits))
  df_neg <- data.frame(x = 2, y = 1, value = mean(neg_value_limits))
  df_nt  <- data.frame(x = 3, y = 1, status = not_tested_label)

  p <- ggplot()

  if (show_pos) {
    p <- p +
      geom_tile(data = df_pos, aes(x = x, y = y, fill = value)) +
      scale_fill_gradient(
        name   = pos_legend_title,
        low    = mid_color, high = pos_color,
        limits = pos_value_limits,
        breaks = pos_breaks_t,
        labels = pos_labels,
        guide  = guide_colorbar(order = 1)
      ) +
      new_scale_fill()
  }

  if (show_neg) {
    p <- p +
      geom_tile(data = df_neg, aes(x = x, y = y, fill = value)) +
      scale_fill_gradient(
        name   = neg_legend_title,
        low    = neg_color, high = mid_color,
        limits = neg_value_limits,
        breaks = neg_breaks_t,
        labels = neg_labels,
        guide  = guide_colorbar(order = 2)
      ) +
      new_scale_fill()
  }

  if (show_not_tested) {
    p <- p +
      geom_tile(data = df_nt, aes(x = x, y = y, fill = status)) +
      scale_fill_manual(
        name   = "",
        values = setNames(not_tested_color, not_tested_label),
        guide  = guide_legend(order = 3)
      )
  }

  p <- p +
    theme_void() +
    theme(
      legend.title    = element_text(size = 12, face = "bold"),
      legend.text     = element_text(size = 10),
      legend.position = "right",
      legend.box      = "vertical"
    )

  # ggplot2 >= 4.0 splits guide boxes into one slot per legend.position
  # ("guide-box-right", "guide-box-left", etc.); cowplot::get_legend()
  # uses partial matching on "guide-box" and warns when multiple slots
  # are present. Match the right-side slot exactly to silence the
  # warning and pick the legend we actually built.
  legend_grob <- get_plot_component(p, "guide-box-right")
  legend_plot <- ggdraw() + draw_grob(legend_grob)

  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  }
  out_path <- file.path(save_dir, file_name)
  ggsave(filename = out_path, plot = legend_plot,
         width = fig_width, height = fig_height, units = "in")
  message("Heatmap legend saved to: ", out_path)
  invisible(NULL)
}

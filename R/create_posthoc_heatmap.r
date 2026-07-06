#' @title Create a P-Value Heatmap with Grey Fill for Untested ROIs
#'
#' @description
#' An extension of `create_pvalue_heatmap` that renders NA cells (ROIs where
#' post-hoc tests were not conducted because the omnibus test was not significant)
#' as grey. All other behavior is identical to the original function: signed
#' p-values are displayed on a red/blue diverging scale, with asterisks for
#' significant cells.
#'
#' @param table_data A data.frame containing the data to visualize.
#'   Cells with NA indicate ROIs that were not tested (omnibus gate not passed).
#' @param data_col_indices A numeric vector of indices for columns with continuous data.
#' @param descriptor_col_indices A numeric vector of indices for columns with y-axis labels.
#' @param scale_type (Optional) `"linear"` (default) or `"log"` to emphasize values near zero.
#' @param log_power (Optional) Controls the intensity of the `"log"` scale (default: 0.2).
#' @param legend_breaks (Optional) A vector of original data values for the legend.
#'   Default: `c(-0.05, -0.01, 0, 0.01, 0.05)`.
#' @param significance_threshold (Optional) Threshold for asterisks. Default: 0.05.
#' @param asterisk_size (Optional) Size of asterisks in mm. Default: 6.
#' @param pos_legend_title (Optional) Title for the positive value legend. Default: "Positive Effect".
#' @param neg_legend_title (Optional) Title for the negative value legend. Default: "Negative Effect".
#' @param low_color (Optional) Color for negative values. Default: "blue".
#' @param mid_color (Optional) Color for the midpoint (0). Default: "white".
#' @param high_color (Optional) Color for positive values. Default: "red".
#' @param untested_color (Optional) Color for NA (untested) cells. Default: "grey80".
#' @param save_dir (Optional) Folder for saving the PDF. Default: current directory.
#' @param file_name (Optional) Name of the output PDF. Default: "PostHoc_Heatmap.pdf".
#' @param fig_width (Optional) Figure width in inches. Default: 8.5.
#' @param fig_height (Optional) Figure height in inches. Default: 11.
#' @param title (Optional) Title of the heatmap. Default: "Post-Hoc P-Value Heatmap".
#' @param landscape (Optional) If TRUE, swaps width/height. Default: FALSE.
#' @param show_legend (Optional) If `TRUE` (default), the side legend
#'   is rendered. If `FALSE`, the legend is suppressed -- useful for
#'   publication panels that share a single legend exported via
#'   [create_heatmap_legend()].
#'
#' @return Saves a plot to a file and prints a confirmation message.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select all_of mutate filter
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradient scale_fill_manual
#' @importFrom ggplot2 guide_colorbar geom_text scale_x_discrete labs theme_minimal
#' @importFrom ggplot2 theme element_text element_blank ggsave
#' @importFrom ggnewscale new_scale_fill
#'
#' @export
create_posthoc_heatmap <- function(
    table_data,
    data_col_indices,
    descriptor_col_indices,
    scale_type = "linear",
    log_power = 0.2,
    legend_breaks = c(-0.05, -0.01, -0.001, 0.001, 0.01, 0.05),
    significance_threshold = 0.05,
    asterisk_size = 6,
    pos_legend_title = "Positive Effect",
    neg_legend_title = "Negative Effect",
    low_color = "blue",
    mid_color = "white",
    high_color = "red",
    untested_color = "grey80",
    save_dir = ".",
    file_name = "PostHoc_Heatmap.pdf",
    fig_width = 8.5,
    fig_height = 11,
    title = "Post-Hoc P-Value Heatmap",
    landscape = FALSE,
    show_legend = TRUE
) {
    # --- 1. Data Preparation ---
    original_x_order <- colnames(table_data)[data_col_indices]

    heatmap_data <- table_data %>%
        select(all_of(data_col_indices))

    y_labels <- table_data[[descriptor_col_indices[1]]]
    heatmap_data$y_labels <- y_labels

    df_long <- melt(heatmap_data, id.vars = "y_labels")
    colnames(df_long) <- c("y_label", "x_label", "value")

    df_long$x_label <- factor(df_long$x_label, levels = original_x_order)

    # Separate NA (untested) from tested cells
    df_untested <- df_long %>% filter(is.na(value))
    df_tested <- df_long %>% filter(!is.na(value))

    # --- 2. Color Transformation for tested cells ---
    if (scale_type == "log") {
        df_tested$color_value <- sign(df_tested$value) * (1 - abs(df_tested$value)^log_power)
        transformed_breaks <- sign(legend_breaks) * (1 - abs(legend_breaks)^log_power)
    } else {
        df_tested$color_value <- sign(df_tested$value) * (1 - abs(df_tested$value))
        transformed_breaks <- sign(legend_breaks) * (1 - abs(legend_breaks))
    }

    df_tested <- df_tested %>%
      mutate(
        significance_label = ifelse(abs(value) <= significance_threshold & value != 0, "*", NA_character_)
      )

    # --- 3. Legend Preparation ---
    pos_breaks <- transformed_breaks[legend_breaks >= 0]
    pos_labels <- abs(legend_breaks[legend_breaks >= 0])

    neg_breaks <- transformed_breaks[legend_breaks < 0]
    neg_labels <- abs(legend_breaks[legend_breaks < 0])

    # --- 4. Plotting ---
    p <- ggplot(mapping = aes(x = x_label, y = y_label))

    # --- GREY TILES for untested ROIs ---
    if (nrow(df_untested) > 0) {
        # Add a dummy column for the grey fill legend
        df_untested$status <- "Not tested\n(omnibus n.s.)"
        p <- p +
            geom_tile(data = df_untested, aes(fill = status), color = "grey50") +
            scale_fill_manual(
                name = "",
                values = c("Not tested\n(omnibus n.s.)" = untested_color),
                guide = ggplot2::guide_legend(order = 3)
            ) +
            new_scale_fill()
    }

    # --- POSITIVE SCALE ---
    df_pos <- df_tested %>% filter(value >= 0)
    if (nrow(df_pos) > 0) {
        p <- p +
            geom_tile(data = df_pos, aes(fill = color_value), color = "grey50") +
            scale_fill_gradient(
                name = pos_legend_title, low = mid_color, high = high_color,
                limits = c(0, 1), breaks = pos_breaks, labels = pos_labels,
                guide = guide_colorbar(order = 1)
            ) +
            new_scale_fill()
    }

    # --- NEGATIVE SCALE ---
    df_neg <- df_tested %>% filter(value < 0)
    if (nrow(df_neg) > 0) {
        p <- p +
            geom_tile(data = df_neg, aes(fill = color_value), color = "grey50") +
            scale_fill_gradient(
                name = neg_legend_title, low = low_color, high = mid_color,
                limits = c(-1, 0), breaks = neg_breaks, labels = neg_labels,
                guide = guide_colorbar(order = 2)
            )
    }

    # --- SIGNIFICANCE ASTERISKS ---
    p <- p +
        geom_text(data = df_tested, aes(label = significance_label),
                  color = "black", size = asterisk_size, na.rm = TRUE) +
        scale_x_discrete(limits = original_x_order) +
        labs(title = title, x = "", y = "") +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12),
            axis.text.y = element_text(size = 12),
            plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            legend.title = element_text(size = 12, face = "bold"),
            legend.text = element_text(size = 10),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = if (show_legend) "right" else "none",
            legend.box = "vertical"
        )

    # --- 5. Saving ---
    if (landscape) {
        temp <- fig_width
        fig_width <- fig_height
        fig_height <- temp
    }

    ggsave(
        filename = file.path(save_dir, file_name), plot = p,
        width = fig_width, height = fig_height, units = "in"
    )

    message(paste("Post-hoc heatmap saved to:", file.path(save_dir, file_name)))
}

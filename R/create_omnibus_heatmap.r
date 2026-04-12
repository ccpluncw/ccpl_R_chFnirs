#' @title Create an Unsigned P-Value Heatmap (Green Gradient) for Omnibus Results
#'
#' @description
#' Generates a heatmap using a single-color green gradient for unsigned p-values
#' (e.g., from omnibus F-tests). Darker green = smaller p-value = stronger evidence.
#' White = p-value of 1 (no evidence). Cells meeting the significance threshold
#' are marked with an asterisk. This function is designed for omnibus test results
#' where there is no directional sign.
#'
#' @param table_data A data.frame containing the data to visualize.
#' @param data_col_indices A numeric vector of indices for columns with p-value data.
#' @param descriptor_col_indices A numeric vector of indices for columns with y-axis labels.
#' @param scale_type (Optional) `"log"` (default) to emphasize values near zero or `"linear"`.
#' @param log_power (Optional) Controls the intensity of the log scale (default: 0.2).
#' @param legend_breaks (Optional) A vector of original p-values for the legend.
#'   Default: `c(0.001, 0.01, 0.05, 0.10, 0.50)`.
#' @param significance_threshold (Optional) Threshold for asterisks. Default: 0.05.
#' @param asterisk_size (Optional) Size of asterisks in mm. Default: 6.
#' @param low_color (Optional) Color for smallest p-values (strongest evidence). Default: "darkgreen".
#' @param high_color (Optional) Color for largest p-values (no evidence). Default: "white".
#' @param save_dir (Optional) Folder for saving the PDF. Default: current directory.
#' @param file_name (Optional) Name of the output PDF. Default: "Omnibus_Heatmap.pdf".
#' @param fig_width (Optional) Figure width in inches. Default: 8.5.
#' @param fig_height (Optional) Figure height in inches. Default: 11.
#' @param title (Optional) Title of the heatmap. Default: "Omnibus P-Value Heatmap".
#' @param landscape (Optional) If TRUE, swaps width/height. Default: FALSE.
#'
#' @return Saves a plot to a file and prints a confirmation message.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select all_of mutate
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradient guide_colorbar geom_text
#' @importFrom ggplot2 scale_x_discrete labs theme_minimal theme element_text element_blank ggsave
#'
#' @export
create_omnibus_heatmap <- function(
    table_data,
    data_col_indices,
    descriptor_col_indices,
    scale_type = "log",
    log_power = 0.2,
    legend_breaks = c(0.001, 0.01, 0.05, 0.10, 0.50),
    significance_threshold = 0.05,
    asterisk_size = 6,
    low_color = "darkgreen",
    high_color = "white",
    save_dir = ".",
    file_name = "Omnibus_Heatmap.pdf",
    fig_width = 8.5,
    fig_height = 11,
    title = "Omnibus P-Value Heatmap",
    landscape = FALSE
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

    # The raw values here are unsigned p-values (0 to 1).
    # We transform so that small p -> high color intensity (near 1),
    # large p -> low color intensity (near 0).
    if (scale_type == "log") {
        df_long$color_value <- (1 - abs(df_long$value)^log_power)
        transformed_breaks <- (1 - abs(legend_breaks)^log_power)
    } else {
        df_long$color_value <- (1 - abs(df_long$value))
        transformed_breaks <- (1 - abs(legend_breaks))
    }

    df_long <- df_long %>%
      mutate(
        significance_label = ifelse(abs(value) <= significance_threshold & value != 0, "*", NA_character_)
      )

    # --- 2. Plotting ---
    p <- ggplot(df_long, aes(x = x_label, y = y_label)) +
        geom_tile(aes(fill = color_value), color = "grey50") +
        scale_fill_gradient(
            name = "p-value",
            low = high_color,
            high = low_color,
            limits = c(0, 1),
            breaks = transformed_breaks,
            labels = legend_breaks,
            guide = guide_colorbar(order = 1)
        ) +
        # --- SIGNIFICANCE ASTERISKS ---
        geom_text(aes(label = significance_label), color = "black", size = asterisk_size, na.rm = TRUE) +
        scale_x_discrete(limits = original_x_order) +
        # --- GENERAL THEME AND LABELS ---
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
            legend.position = "right",
            legend.box = "vertical"
        )

    # --- 3. Saving the Plot ---
    if (landscape) {
        temp <- fig_width
        fig_width <- fig_height
        fig_height <- temp
    }

    ggsave(
        filename = file.path(save_dir, file_name), plot = p,
        width = fig_width, height = fig_height, units = "in"
    )

    message(paste("Omnibus heatmap saved to:", file.path(save_dir, file_name)))
}

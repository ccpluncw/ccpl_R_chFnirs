#' @title Create a P-Value-Focused Continuous Heatmap with Split Legends
#'
#' @description This function generates a heatmap for continuous data from -1 to 1. It uses
#' two separate, side-by-side color legends: one for positive values (e.g., red)
#' and one for negative values (e.g., blue). This is ideal for p-values where
#' the sign indicates direction and the magnitude indicates significance. Color
#' intensity is highest near zero and can be scaled linearly or with a log-style
#' transformation. The original x-axis column order is preserved. Cells with
#' values below a specified alpha threshold are marked with an asterisk.
#'
#' @param table_data A data.frame containing the data to visualize.
#' @param data_col_indices A numeric vector of indices for columns with continuous data.
#' @param descriptor_col_indices A numeric vector of indices for columns with y-axis labels.
#' @param scale_type (Optional) `"linear"` (default) or `"log"` to emphasize values near zero.
#' @param log_power (Optional) Controls the intensity of the `"log"` scale (default: 0.1).
#' @param legend_breaks (Optional) A vector of original data values for the legend.
#'   Default: `c(-0.05, -0.01, 0, 0.01, 0.05)`.
#' @param significance_threshold (Optional) An asterisk is placed in cells where
#'   the absolute value is <= this threshold (and not zero). Default: 0.05.
#' @param asterisk_size (Optional) The numeric size of the significance asterisk in mm.
#'   If set to 0, no asterisks will be plotted. Default: 6.
#' @param pos_legend_title (Optional) Title for the positive value legend.
#'   Default: `"Positive Effect"`.
#' @param neg_legend_title (Optional) Title for the negative value legend.
#'   Default: `"Negative Effect"`.
#' @param low_color (Optional) Color for negative values. Default: "blue".
#' @param mid_color (Optional) Color for the midpoint (0). Default: "white".
#' @param high_color (Optional) Color for positive values. Default: "red".
#' @param save_dir (Optional) Folder for saving the PDF. Default: current directory.
#' @param file_name (Optional) Name of the output PDF. Default: "P_Value_Heatmap.pdf".
#' @param fig_width (Optional) Figure width in inches. Default: 8.5.
#' @param fig_height (Optional) Figure height in inches. Default: 11.
#' @param title (Optional) The title of the heatmap. Default: "P-Value Heatmap".
#' @param landscape (Optional) If TRUE, swaps figure width/height. Default: FALSE.
#'
#' @return Saves a plot to a file and prints a confirmation message to the console.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select all_of mutate filter
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradient guide_colorbar geom_text
#' @importFrom ggplot2 scale_x_discrete labs theme_minimal theme element_text element_blank ggsave
#' @importFrom ggnewscale new_scale_fill
#'
#' @examples
#' # 1. Create sample data with a deliberate column order
#' set.seed(42)
#' p_value_data <- data.frame(
#'   Gene_Symbol = paste("GENE", LETTERS[1:15]),
#'   Z_Contrast = c(runif(5, -0.05, 0.05), runif(5, -1, -0.5), runif(5, 0.5, 1)),
#'   A_Contrast = sample(c(runif(5, -0.01, 0.01), runif(10, -1, 1))),
#'   B_Contrast = runif(15, -1, 1),
#'   stringsAsFactors = FALSE
#' )
#' p_value_data[,2:4] <- sapply(p_value_data[,2:4], function(x) pmax(-1, pmin(1, x)))
#'
#' \dontrun{
#' # 2. Generate the heatmap with a custom asterisk size
#' create_pvalue_heatmap(
#'   table_data = p_value_data,
#'   data_col_indices = 2:4,
#'   descriptor_col_indices = 1,
#'   title = "Heatmap with Large Asterisks",
#'   file_name = "Large_Asterisk_Heatmap.pdf",
#'   asterisk_size = 8 # Making the asterisks larger
#' )
#'
#' # 3. Generate the heatmap with no asterisks
#' create_pvalue_heatmap(
#'   table_data = p_value_data,
#'   data_col_indices = 2:4,
#'   descriptor_col_indices = 1,
#'   title = "Heatmap with No Asterisks",
#'   file_name = "No_Asterisk_Heatmap.pdf",
#'   asterisk_size = 0 # Disables asterisks
#' )
#' }
#' @export
create_pvalue_heatmap <- function(
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
    save_dir = ".",
    file_name = "Adjusted_P_Value_Heatmap.pdf",
    fig_width = 8.5,
    fig_height = 11,
    title = "Adjusted P-Value Heatmap",
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

    # --- 2. Key Transformation for Color Intensity ---
    if (scale_type == "log") {
        df_long$color_value <- sign(df_long$value) * (1 - abs(df_long$value)^log_power)
        transformed_breaks <- sign(legend_breaks) * (1 - abs(legend_breaks)^log_power)
    } else {
        df_long$color_value <- sign(df_long$value) * (1 - abs(df_long$value))
        transformed_breaks <- sign(legend_breaks) * (1 - abs(legend_breaks))
    }

    df_long <- df_long %>%
      mutate(
        significance_label = ifelse(abs(value) <= significance_threshold & value != 0, "*", NA_character_)
      )

    # --- 3. Legend Preparation for Split Scales ---
    pos_breaks <- transformed_breaks[legend_breaks >= 0]
    pos_labels <- abs(legend_breaks[legend_breaks >= 0])

    neg_breaks <- transformed_breaks[legend_breaks < 0]
    neg_labels <- abs(legend_breaks[legend_breaks < 0])

    # --- 4. Plotting with ggplot2 and ggnewscale ---
    p <- ggplot(df_long, aes(x = x_label, y = y_label)) +
        # --- POSITIVE SCALE ---
        geom_tile(data = . %>% filter(value >= 0), aes(fill = color_value)) +
        scale_fill_gradient(
            name = pos_legend_title, low = mid_color, high = high_color,
            limits = c(0, 1), breaks = pos_breaks, labels = pos_labels,
            guide = guide_colorbar(order = 1)
        ) +
        # --- NEW SCALE ---
        new_scale_fill() +
        # --- NEGATIVE SCALE ---
        geom_tile(data = . %>% filter(value < 0), aes(fill = color_value)) +
        scale_fill_gradient(
            name = neg_legend_title, low = low_color, high = mid_color,
            limits = c(-1, 0), breaks = neg_breaks, labels = neg_labels,
            guide = guide_colorbar(order = 2)
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

    # --- 5. Saving the Plot ---
    if (landscape) {
        temp <- fig_width
        fig_width <- fig_height
        fig_height <- temp
    }

    ggsave(
        filename = file.path(save_dir, file_name), plot = p,
        width = fig_width, height = fig_height, units = "in"
    )

    message(paste("Heatmap saved to:", file.path(save_dir, file_name)))
}
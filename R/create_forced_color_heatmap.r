#' @title Create and Save a Forced-Color Heatmap
#'
#' @description
#' Generates a heatmap from a data frame with discrete, forced colors for
#' specific values (e.g., -1, 0, 1). This is useful for visualizing categorical
#' or significance-coded data. The plot is created with `ggplot2` and saved as a PDF.
#'
#' @param table_data A data.frame containing the data to visualize.
#' @param data_col_indices A numeric vector of indices for the columns containing
#'   the numerical data for the heatmap.
#' @param descriptor_col_indices A numeric vector of indices for the columns
#'   containing descriptive strings for the y-axis labels. The first index
#'   is used for the labels.
#' @param color_mapping (Optional) A named vector of colors for the values.
#'   The names must match the unique values in the data columns.
#'   Default: `c("-1" = "blue", "0" = "white", "1" = "red")`.
#' @param save_dir (Optional) A character string specifying the folder for saving
#'   the PDF. Default: current working directory (".").
#' @param file_name (Optional) The name of the output PDF file. Default: "Effect_Table.pdf".
#' @param fig_width (Optional) Width of the figure in inches. Default: 8.5.
#' @param fig_height (Optional) Height of the figure in inches. Default: 11.
#' @param title (Optional) The title of the heatmap. Default: "Effect Table".
#' @param landscape (Optional) A logical value. If TRUE, swaps figure width and
#'   height for landscape orientation. Default: FALSE.
#'
#' @return Nothing. The function is called for its side effect of saving a plot to a file.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select all_of
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_manual labs theme_minimal
#' @importFrom ggplot2 theme element_text element_blank ggsave
#'
#' @examples
#' # 1. Create sample data
#' my_data_table <- data.frame(
#'   Gene_Name = paste("Gene", 1:10),
#'   Description = paste("Description for gene", 1:10),
#'   Sample_A = sample(c(-1, 0, 1), 10, replace = TRUE),
#'   Sample_B = sample(c(-1, 0, 1), 10, replace = TRUE),
#'   Sample_C = sample(c(-1, 0, 1), 10, replace = TRUE),
#'   stringsAsFactors = FALSE
#' )
#'
#' \dontrun{
#' # 2. Basic usage (portrait, default settings)
#' create_forced_color_heatmap(
#'   table_data = my_data_table,
#'   data_col_indices = 3:5,
#'   descriptor_col_indices = 1
#'  )
#'
#' # 3. Customizing colors and saving location
#' if (!dir.exists("results")) dir.create("results")
#' create_forced_color_heatmap(
#'   table_data = my_data_table,
#'   data_col_indices = 3:5,
#'   descriptor_col_indices = 1,
#'   color_mapping = c("-1" = "darkblue", "0" = "gray90", "1" = "darkred"),
#'   save_dir = "results",
#'   file_name = "Custom_Heatmap.pdf"
#' )
#' }
#' @export
create_forced_color_heatmap <- function(
    table_data,
    data_col_indices,
    descriptor_col_indices,
    color_mapping = c("-1" = "blue", "0" = "white", "1" = "red"),
    save_dir = ".",
    file_name = "Effect_Table.pdf",
    fig_width = 8.5,
    fig_height = 11,
    title = "Effect Table",
    landscape = FALSE
) {
    heatmap_data <- table_data %>%
        select(all_of(data_col_indices))

    y_labels <- table_data[[descriptor_col_indices[1]]]
    heatmap_data$y_labels <- y_labels

    df_long <- melt(heatmap_data, id.vars = "y_labels")
    colnames(df_long) <- c("y_label", "x_label", "value")

    df_long$value <- factor(df_long$value, levels = names(color_mapping))

    p <- ggplot(df_long, aes(x = x_label, y = y_label, fill = value)) +
        geom_tile(color = "black") +
        scale_fill_manual(values = color_mapping, drop = FALSE) +
        labs(title = title, x = "", y = "") +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.title = element_blank()
        )

    if (landscape) {
        temp <- fig_width
        fig_width <- fig_height
        fig_height <- temp
    }

    ggsave(
        filename = file.path(save_dir, file_name),
        plot = p,
        width = fig_width,
        height = fig_height,
        units = "in"
    )

    message(paste("Heatmap saved to:", file.path(save_dir, file_name)))
}
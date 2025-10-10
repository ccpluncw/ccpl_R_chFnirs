#' @title Create and Save a Colored Heatmap
#'
#' @description
#' Generates a heatmap from a data frame using the `ggplot2` package, designed to
#' be similar to MATLAB's default heatmap, and saves it as a PDF file.
#'
#' @param table_data A data.frame containing the data to visualize.
#' @param data_col_indices A numeric vector of indices for the columns containing
#'   the numerical data for the heatmap.
#' @param descriptor_col_indices A numeric vector of indices for the columns
#'   containing descriptive strings for the y-axis labels. The first index
#'   is used for the labels.
#' @param color_map (Optional) A character vector of three colors for the low,
#'   mid, and high points of the gradient. Default: `c("blue", "white", "red")`.
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
#' @importFrom dplyr select all_of mutate
#' @importFrom tidyr pivot_longer
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradient2 labs theme_minimal
#' @importFrom ggplot2 theme element_text margin element_blank element_rect ggsave
#'
#' @examples
#' # 1. Create sample data
#' my_data_table <- data.frame(
#'   Gene_Name = paste("Gene", 1:10),
#'   Description = paste("Description for gene", 1:10),
#'   Sample_A = rnorm(10),
#'   Sample_B = rnorm(10, mean = 1),
#'   Sample_C = rnorm(10, mean = -1),
#'   stringsAsFactors = FALSE
#' )
#'
#' \dontrun{
#' # 2. Basic usage (portrait, default settings)
#' create_colored_heatmap(
#'   table_data = my_data_table,
#'   data_col_indices = 3:5,
#'   descriptor_col_indices = 1
#' )
#'
#' # 3. Custom colormap and file path
#' my_colors <- c("darkgreen", "yellow", "red")
#' # Create a directory for the results if it doesn't exist
#' if (!dir.exists("results")) dir.create("results")
#' create_colored_heatmap(
#'   table_data = my_data_table,
#'   data_col_indices = 3:5,
#'   descriptor_col_indices = 1,
#'   color_map = my_colors,
#'   save_dir = "results",
#'   file_name = "my_heatmap.pdf"
#' )
#' }
#' @export
create_colored_heatmap <- function(table_data,
                                     data_col_indices,
                                     descriptor_col_indices,
                                     color_map = c("#0000FF", "#FFFFFF", "#FF0000"),
                                     save_dir = ".",
                                     file_name = "Effect_Table.pdf",
                                     fig_width = 8.5,
                                     fig_height = 11,
                                     title = "Effect Table",
                                     landscape = FALSE) {

  # --- 1. Input Validation and Setup ---
  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
    message("Created directory: ", save_dir)
  }

  if (length(color_map) != 3) {
    stop("`color_map` must be a vector of three colors (low, mid, high).")
  }

  # --- 2. Extract and Reshape Data ---
  if (!is.data.frame(table_data)) {
    table_data <- as.data.frame(table_data)
  }

  descriptor_col_name <- names(table_data)[descriptor_col_indices[1]]

  data_long <- table_data %>%
    select(all_of(descriptor_col_name), all_of(data_col_indices)) %>%
    pivot_longer(
      cols = -all_of(descriptor_col_name),
      names_to = "x_axis_labels",
      values_to = "heatmap_value"
    ) %>%
    mutate(x_axis_labels = factor(x_axis_labels, levels = names(table_data)[data_col_indices]))

  # --- 3. Create the Heatmap Plot using ggplot2 ---
  heatmap_plot <- ggplot(
      data_long,
      aes(x = x_axis_labels, y = .data[[descriptor_col_name]], fill = heatmap_value)
    ) +
    geom_tile(color = "grey85", size = 0.5) +
    scale_fill_gradient2(
      low = color_map[1], mid = color_map[2], high = color_map[3],
      midpoint = 0, guide = "colorbar"
    ) +
    labs(
      title = title, x = NULL, y = NULL, fill = "Value"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b = 20)),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "right" # Changed from "none" to show the legend by default
    )

  # --- 4. Save the Figure ---
  save_width <- if (landscape) fig_height else fig_width
  save_height <- if (landscape) fig_width else fig_height
  full_path <- file.path(save_dir, file_name)

  ggsave(
    filename = full_path, plot = heatmap_plot,
    width = save_width, height = save_height,
    units = "in", device = "pdf"
  )

  message("Heatmap saved to: ", full_path)
}
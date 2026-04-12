#' @title Plot simple slopes for continuous x continuous interactions
#'
#' @description
#' For each significant (or matching) type x area cell in a BH-corrected
#' omnibus table, fits the requested model and produces a pair of
#' simple-slopes plots: one with `var1` on the x-axis and `var2` as the
#' moderator (split at levels of the moderator expressed in SD units),
#' and one with the roles reversed. Intended as the canonical
#' visualization for a continuous-by-continuous interaction
#' (`var1:var2`) in an area-wise fNIRS analysis — typically used after
#' `run_omnibus_by_area_bh()` and `run_gated_slopes()` have identified
#' which cells show significant interactions.
#'
#' The function always returns a nested list of ggplot objects
#' organized by chromophore type, so callers can inspect or further
#' customize the plots interactively. If `pdf_prefix` is supplied,
#' the function also writes one PDF per type as a side effect, with
#' both pages printed per cell. Empty types (no matching cells) are
#' included in the returned list with value `list()` so the structure
#' is predictable.
#'
#' @details
#' \strong{Output format.} When `pdf_prefix` is supplied, the function
#' writes **PDF files only**. PDF is the only built-in output format
#' because it natively supports multi-page output: each chromophore
#' type produces a single PDF containing two pages per cell (page A:
#' `var1` on the x-axis, page B: `var2` on the x-axis), in the order
#' the cells appear in the filtered omnibus table. Raster formats like
#' JPEG, PNG, SVG, and TIFF do not have multi-page support, so adding
#' them would fan each analysis out to dozens of separate files and is
#' intentionally not implemented here.
#'
#' If you need raster output for presentations or publications, use
#' the **returned plot list**: each element is a
#' `list(page_A = <ggplot>, page_B = <ggplot>)` that you can pass
#' directly to `ggplot2::ggsave()` with the raster format and DPI of
#' your choice. This is the recommended workflow because it also lets
#' you customize individual plots (axis labels, annotations,
#' dimensions) before saving, which is usually what "I want this for a
#' slide" turns out to mean. See the examples for a concrete
#' `ggsave()` call.
#'
#' \strong{How the lines are constructed.} For each moderator level
#' \eqn{m}, the function uses `emmeans::emtrends()` to obtain the
#' slope of the focal variable at that moderator value, and
#' `emmeans::emmeans()` to obtain the fitted intercept at
#' `(focal = mean(focal), moderator = m)`. The predicted line at
#' moderator level \eqn{m} is then reconstructed analytically as
#' \eqn{\hat{y}(x, m) = b_0(m) + s(m) \cdot (x - \bar{x})}, where
#' \eqn{b_0(m)} is the fitted intercept, \eqn{s(m)} is the slope, and
#' \eqn{\bar{x}} is the focal variable's sample mean in the cell.
#' This approach is numerically equivalent to using
#' `ggeffects::ggpredict()` or `interactions::interact_plot()` on the
#' same model but avoids adding either as a dependency.
#'
#' \strong{Effect pattern matching.} The default `effect_pattern` is a
#' regex that matches any effect name containing both `var1` and `var2`
#' as whole words, in either order. This means cells gated on
#' `var1:var2`, `var2:var1`, `factor:var1:var2`, etc. all get picked up
#' by a single call. Pass a narrower pattern if you want to restrict
#' to a single specific effect name.
#'
#' \strong{Output layout.} Each type's PDF contains two consecutive
#' pages per cell: page A plots `var1` on the x-axis with `var2` as the
#' moderator, page B reverses the roles. Cells are ordered as they
#' appear in the filtered omnibus table. Non-significant cells are
#' skipped silently when `gate = TRUE` (the default); pass
#' `gate = FALSE` to plot every cell in the filtered omnibus
#' regardless of p-value.
#'
#' @param omnibus_df A BH-corrected omnibus table in Tier 2 schema
#'   (lowercase `type`, `area`, `effect`, `p_bh`, etc.), typically
#'   produced by `run_omnibus_by_area_bh()`. Must contain an `effect`
#'   column, and (when `gate = TRUE`) a column named by `p_col`.
#' @param data The long-format data.frame used to fit the omnibus.
#' @param formula Model formula used to fit each cell. Accepts either
#'   a formula object (preferred) or a character string. Must include
#'   `var1` and `var2` as predictors (with their interaction) for the
#'   slopes to be identifiable.
#' @param area_col Character scalar naming the area column in `data`.
#' @param var1 Character scalar naming the first continuous variable.
#' @param var2 Character scalar naming the second continuous variable.
#' @param type_col Character scalar naming the measurement-type column
#'   in `data`. Default `"type"`.
#' @param model Required character scalar: `"lmer"` or `"lm"`. See
#'   `run_omnibus_by_area()` for the rationale behind requiring this
#'   argument.
#' @param vcov Character scalar naming the covariance estimator. See
#'   `run_omnibus_by_area()` for valid values and the `NULL`
#'   conditional default.
#' @param effect_pattern Regex passed to `grepl()` against
#'   `omnibus_df$effect` to identify the interaction effect(s) of
#'   interest. Default matches any effect name containing `var1` and
#'   `var2` as whole words in either order.
#' @param p_col Character scalar naming the gate column in
#'   `omnibus_df`. Default `"p_bh"`. Ignored when `gate = FALSE`.
#' @param p_threshold Numeric scalar. Cells with `p_col < p_threshold`
#'   (and non-NA) pass the gate. Default `0.05`. Ignored when
#'   `gate = FALSE`.
#' @param gate Logical. When `TRUE` (default), only cells whose
#'   gate-column value falls below `p_threshold` are plotted. When
#'   `FALSE`, every cell whose effect name matches `effect_pattern` is
#'   plotted regardless of significance — useful for exploratory
#'   visualization before a study is powered or when inspecting
#'   non-significant trends.
#' @param min_rows Integer scalar; cells with fewer rows are skipped.
#'   Default `5`.
#' @param y_label Character scalar used as the y-axis label on every
#'   plot. Default `NULL` auto-detects the label from the left-hand
#'   side of `formula` (so a formula of `estimate ~ ...` produces
#'   `"estimate"`, `AUC ~ ...` produces `"AUC"`, etc.). Pass an
#'   explicit string to override.
#' @param mod_sd_levels Numeric vector of moderator levels to plot,
#'   expressed in SD units relative to the moderator's mean in the
#'   current cell. Default `c(-1, 0, 1)` produces three lines: Mean-1SD,
#'   Mean, Mean+1SD. Use `c(-1, 1)` to drop the mean line, `c(-2, 0, 2)`
#'   for wider fans, etc. Level labels are auto-generated from the
#'   values as `<moderator> +/-<n>SD` / `<moderator> Mean`.
#' @param line_colors Character vector of colors for the moderator
#'   levels. Recycled to match `length(mod_sd_levels)`. Default
#'   `c("grey65", "steelblue", "navyblue")`.
#' @param line_types Character vector of ggplot2 linetype names,
#'   recycled to match. Default `c("dashed", "solid", "longdash")`.
#' @param theme_base_size Numeric scalar passed to
#'   `ggplot2::theme_minimal()` as `base_size`. Default `14`.
#' @param pdf_prefix Character scalar naming the file-name prefix for
#'   the PDF outputs, or `NULL` (default) to skip file writing
#'   entirely. The function **always** returns the plot list
#'   regardless of this argument — `pdf_prefix` only controls whether
#'   PDFs are also written as a side effect. When supplied, one PDF
#'   per chromophore type is written to the current working directory
#'   as `<pdf_prefix>_<type>_SimpleSlopes.pdf`, each containing all
#'   page A + page B pairs for that type in order. PDF is the only
#'   supported output format; see \strong{Details} for why, and for
#'   how to produce raster output (JPEG, PNG, etc.) from the returned
#'   plot list via `ggplot2::ggsave()`.
#' @param pdf_width,pdf_height Numeric. PDF dimensions in inches.
#'   Defaults 10 and 7.
#'
#' @return Invisibly returns a named nested list with one element per
#'   chromophore type present in the filtered omnibus. Each element is
#'   itself a named list keyed by `"<area>|<effect>"`, whose values are
#'   two-element lists with names `page_A` and `page_B`, each holding a
#'   ggplot object. Types with no surviving cells contribute an empty
#'   `list()`. Callers can `print(result$hbo[[1]]$page_A)` interactively
#'   or iterate over the nested structure programmatically.
#'
#' @seealso [run_omnibus_by_area_bh()] for building the omnibus input;
#'   [run_gated_slopes()] for the statistical post-hoc companion.
#'
#' @examples
#' \dontrun{
#' # Typical two-phase workflow: omnibus on the full interaction model,
#' # filter to the slope effects, then plot the simple slopes for the
#' # significant cells.
#' df.omnibus <- run_omnibus_by_area_bh(
#'   data     = df.all,
#'   formula  = estimate ~ Group * cond * cortisol * testosterone + (1 | ID),
#'   area_col = "ROI",
#'   model    = "lmer"
#' )
#' # Optional script-side filter (keep only slope-containing effects)
#' df.omnibus.h <- df.omnibus[grepl("cortisol|testosterone",
#'                                  df.omnibus$effect), ]
#'
#' # Plot and write PDFs (one per chromophore type)
#' slopes <- plot_simple_slopes(
#'   omnibus_df    = df.omnibus.h,
#'   data          = df.all,
#'   formula       = estimate ~ Group * cond * cortisol * testosterone + (1 | ID),
#'   area_col      = "ROI",
#'   var1          = "cortisol",
#'   var2          = "testosterone",
#'   model         = "lmer",
#'   pdf_prefix    = "ROI_CortTest"
#' )
#'
#' # Inspect interactively without writing files
#' slopes_only <- plot_simple_slopes(
#'   omnibus_df = df.omnibus.h,
#'   data       = df.all,
#'   formula    = estimate ~ Group * cond * cortisol * testosterone + (1 | ID),
#'   area_col   = "ROI",
#'   var1       = "cortisol",
#'   var2       = "testosterone",
#'   model      = "lmer"
#' )
#' print(slopes_only$hbo[[1]]$page_A)
#'
#' # Save a specific page as a high-DPI JPEG for a slide deck
#' # (raster escape hatch — uses the returned plot list directly)
#' ggplot2::ggsave(
#'   filename = "hbo_dlpfc_cortisol_page.jpg",
#'   plot     = slopes_only$hbo[[1]]$page_A,
#'   width    = 8, height = 6, dpi = 300
#' )
#' }
#'
#' @importFrom stats as.formula sd
#' @importFrom emmeans emmeans emtrends
#' @importFrom ggplot2 ggplot aes geom_line scale_color_manual
#' @importFrom ggplot2 scale_linetype_manual labs theme_minimal theme
#' @importFrom ggplot2 element_text
#' @importFrom grDevices pdf dev.off
#' @export
plot_simple_slopes <- function(
    omnibus_df,
    data,
    formula,
    area_col,
    var1,
    var2,
    type_col         = "type",
    model,
    vcov             = NULL,
    effect_pattern   = paste0("\\b", var1, "\\b.*\\b", var2, "\\b"),
    p_col            = "p_bh",
    p_threshold      = 0.05,
    gate             = TRUE,
    min_rows         = 5,
    y_label          = NULL,
    mod_sd_levels    = c(-1, 0, 1),
    line_colors      = c("grey65", "steelblue", "navyblue"),
    line_types       = c("dashed", "solid", "longdash"),
    theme_base_size  = 14,
    pdf_prefix    = NULL,
    pdf_width        = 10,
    pdf_height       = 7
) {

  # --- Validate model (required) ---
  if (missing(model)) {
    stop("Argument `model` is required. Choose one of:\n",
         "  \"lmer\" - fit a linear mixed model via lmerTest::lmer (for\n",
         "           repeated-measures data with random intercepts/slopes)\n",
         "  \"lm\"   - fit an ordinary linear model via stats::lm (for\n",
         "           flat / independent-observation data)",
         call. = FALSE)
  }
  model <- match.arg(model, c("lmer", "lm"))

  # --- Resolve vcov ---
  valid_vcov <- c("classical", "HC0", "HC1", "HC2", "HC3", "HC4", "HC5")
  if (is.null(vcov)) {
    vcov <- if (model == "lmer") "classical" else "HC3"
  } else {
    vcov <- match.arg(vcov, valid_vcov)
  }
  if (model == "lmer" && vcov != "classical") {
    stop("`vcov = \"", vcov, "\"` is not valid when `model = \"lmer\"`.\n",
         "Set `vcov = \"classical\"` or `model = \"lm\"`. See ?run_omnibus_by_area.",
         call. = FALSE)
  }

  # --- Coerce formula ---
  if (is.character(formula)) formula <- stats::as.formula(formula)

  # --- Resolve y_label from formula LHS if NULL ---
  if (is.null(y_label)) {
    form_chars <- as.character(formula)
    y_label <- if (length(form_chars) >= 2) form_chars[2] else "Response"
  }

  # --- Validate inputs ---
  if (!is.data.frame(omnibus_df)) {
    stop("`omnibus_df` must be a data.frame.", call. = FALSE)
  }
  if (!("effect" %in% colnames(omnibus_df))) {
    stop("`omnibus_df` does not contain an `effect` column. Expected the ",
         "output of `run_omnibus_by_area_bh()` or `run_omnibus_by_area()`.",
         call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.", call. = FALSE)
  }
  if (!(area_col %in% colnames(data))) {
    stop("`area_col` '", area_col, "' not found in `data`.", call. = FALSE)
  }
  if (!(type_col %in% colnames(data))) {
    stop("`type_col` '", type_col, "' not found in `data`.", call. = FALSE)
  }
  if (!(var1 %in% colnames(data))) {
    stop("`var1` '", var1, "' not found in `data`.", call. = FALSE)
  }
  if (!(var2 %in% colnames(data))) {
    stop("`var2` '", var2, "' not found in `data`.", call. = FALSE)
  }
  if (!is.numeric(mod_sd_levels) || length(mod_sd_levels) == 0) {
    stop("`mod_sd_levels` must be a non-empty numeric vector.", call. = FALSE)
  }

  # --- Find matching effects ---
  matching_effects <- unique(
    omnibus_df$effect[grepl(effect_pattern, omnibus_df$effect)]
  )
  if (length(matching_effects) == 0) {
    message("No omnibus effects match pattern '", effect_pattern,
            "'. Returning empty result.")
    return(invisible(list()))
  }

  sig_rows <- omnibus_df[omnibus_df$effect %in% matching_effects, , drop = FALSE]

  # --- Apply gate if requested ---
  if (gate) {
    if (!(p_col %in% colnames(omnibus_df))) {
      stop("`omnibus_df` does not contain the gate column `", p_col, "`.\n",
           "Either pass a BH-corrected omnibus (see `run_omnibus_by_area_bh()`)\n",
           "or set `gate = FALSE` to skip the significance filter.",
           call. = FALSE)
    }
    sig_rows <- sig_rows[!is.na(sig_rows[[p_col]]) &
                           sig_rows[[p_col]] < p_threshold, , drop = FALSE]
  }

  if (nrow(sig_rows) == 0) {
    mode_label <- if (gate) "significant" else "matching"
    message("No ", mode_label, " cells to plot for pattern '",
            effect_pattern, "'.")
    return(invisible(list()))
  }

  # --- Per-cell fit closure ---
  fit_fn <- function(df_cell) {
    if (model == "lmer") {
      lmerTest::lmer(formula, data = df_cell)
    } else {
      stats::lm(formula, data = df_cell)
    }
  }

  # --- Inner helper: build one simple-slopes page ---
  # focal_var on the x-axis, moderator_var split at mod_sd_levels SD
  build_slope_page <- function(fit, df_sub, focal_var, moderator_var,
                               ty, ar, eff, pval, vc) {

    mod_mean  <- mean(df_sub[[moderator_var]], na.rm = TRUE)
    mod_sd_v  <- stats::sd(df_sub[[moderator_var]], na.rm = TRUE)
    mod_values <- mod_mean + mod_sd_levels * mod_sd_v

    # Level labels: "<mod> +1SD" / "<mod> Mean" / "<mod> -1SD" etc.
    make_label <- function(sd_mult) {
      if (sd_mult == 0) {
        paste0(moderator_var, " Mean")
      } else if (sd_mult > 0) {
        paste0(moderator_var, " +", format(sd_mult, trim = TRUE), "SD")
      } else {
        paste0(moderator_var, " ", format(sd_mult, trim = TRUE), "SD")
      }
    }
    lev_labels <- vapply(mod_sd_levels, make_label, character(1))

    focal_range <- seq(min(df_sub[[focal_var]], na.rm = TRUE),
                       max(df_sub[[focal_var]], na.rm = TRUE),
                       length.out = 50)
    focal_m <- mean(df_sub[[focal_var]], na.rm = TRUE)
    trend_col <- paste0(focal_var, ".trend")

    preds_list <- lapply(seq_along(mod_values), function(j) {
      mv    <- mod_values[j]
      label <- lev_labels[j]

      at_slope <- setNames(list(mv), moderator_var)
      at_int   <- setNames(list(focal_m, mv),
                           c(focal_var, moderator_var))

      if (!is.null(vc)) {
        emt <- emmeans::emtrends(fit, ~ 1, var = focal_var,
                                 at = at_slope, vcov. = vc)
        em  <- emmeans::emmeans(fit, ~ 1, at = at_int, vcov. = vc)
      } else {
        emt <- emmeans::emtrends(fit, ~ 1, var = focal_var, at = at_slope)
        em  <- emmeans::emmeans(fit, ~ 1, at = at_int)
      }

      slope     <- as.data.frame(summary(emt))[[trend_col]]
      intercept <- as.data.frame(summary(em))$emmean

      data.frame(
        focal    = focal_range,
        y_pred   = intercept + slope * (focal_range - focal_m),
        ModLevel = factor(label, levels = lev_labels),
        stringsAsFactors = FALSE
      )
    })
    preds <- do.call(rbind, preds_list)

    # Recycle colors and linetypes to match the number of levels
    n_levels   <- length(lev_labels)
    colors_used <- rep_len(line_colors, n_levels)
    ltypes_used <- rep_len(line_types, n_levels)
    names(colors_used) <- lev_labels
    names(ltypes_used) <- lev_labels

    subtitle <- paste0("Effect: ", eff)
    if (!is.na(pval)) {
      subtitle <- paste0(subtitle, "  (", p_col, " = ",
                         format(round(pval, 4), nsmall = 4), ")")
    }

    ggplot(preds, aes(x = focal, y = y_pred,
                      color = ModLevel, linetype = ModLevel)) +
      geom_line(linewidth = 1.2) +
      scale_color_manual(values = colors_used) +
      scale_linetype_manual(values = ltypes_used) +
      labs(
        title    = paste0(ty, " | ", ar),
        subtitle = subtitle,
        x        = focal_var,
        y        = y_label,
        color    = "Moderator",
        linetype = "Moderator"
      ) +
      theme_minimal(base_size = theme_base_size) +
      theme(
        plot.title      = element_text(face = "bold"),
        legend.position = "bottom"
      )
  }

  # --- Main loop: build pages per type ---
  types_in_sig <- unique(sig_rows$type)
  result <- setNames(vector("list", length(types_in_sig)), types_in_sig)
  for (ty in names(result)) result[[ty]] <- list()

  for (ty in types_in_sig) {
    sig_rows_ty <- sig_rows[sig_rows$type == ty, , drop = FALSE]
    type_pages <- list()

    for (i in seq_len(nrow(sig_rows_ty))) {
      ar   <- sig_rows_ty$area[i]
      eff  <- sig_rows_ty$effect[i]
      pval <- if (p_col %in% colnames(sig_rows_ty)) {
        sig_rows_ty[[p_col]][i]
      } else {
        NA_real_
      }

      df_cell <- data[data[[area_col]] == ar & data[[type_col]] == ty, ,
                      drop = FALSE]
      if (nrow(df_cell) < min_rows) next

      tryCatch({
        fit <- fit_fn(df_cell)
        vc <- if (model == "lm" && vcov != "classical") {
          sandwich::vcovHC(fit, type = vcov)
        } else {
          NULL
        }

        page_A <- build_slope_page(fit, df_cell,
                                   focal_var = var1, moderator_var = var2,
                                   ty = ty, ar = ar, eff = eff,
                                   pval = pval, vc = vc)
        page_B <- build_slope_page(fit, df_cell,
                                   focal_var = var2, moderator_var = var1,
                                   ty = ty, ar = ar, eff = eff,
                                   pval = pval, vc = vc)

        cell_key <- paste(ar, eff, sep = "|")
        type_pages[[cell_key]] <- list(page_A = page_A, page_B = page_B)
      }, error = function(e) {
        message("  [WARNING] plot_simple_slopes failed for type = ", ty,
                ", area = ", ar, ", effect = ", eff, ": ",
                conditionMessage(e))
      })
    }

    result[[ty]] <- type_pages
  }

  # --- Optionally write per-type PDFs ---
  if (!is.null(pdf_prefix)) {
    for (ty in names(result)) {
      type_pages <- result[[ty]]
      if (length(type_pages) == 0) next

      pdf_name <- paste0(pdf_prefix, "_", ty, "_SimpleSlopes.pdf")
      grDevices::pdf(pdf_name, width = pdf_width, height = pdf_height)
      for (cell_name in names(type_pages)) {
        cell <- type_pages[[cell_name]]
        print(cell$page_A)
        print(cell$page_B)
      }
      grDevices::dev.off()
      message("Wrote simple slopes plot: ", pdf_name)
    }
  }

  invisible(result)
}

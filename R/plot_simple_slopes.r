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
#' \strong{Confidence bands and colors.} When \code{show_ci = TRUE} (the
#' default), a ribbon behind each line shows the pointwise 95% confidence
#' interval from \code{emmeans::emmeans()} at the plotted x values, with
#' the moderator at each requested level and any other factors treated
#' the same way as for the line. Degrees of freedom are whatever
#' \code{emmeans} uses for that model. The y-axis range includes the
#' bands. With \code{line_colors = NULL} the lines and bands run from
#' blue through gray to red across the sorted levels of
#' \code{mod_sd_levels}.
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
#' @param omnibus_df A BH-corrected omnibus table in the canonical schema
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
#' @param facet_by Character vector of factor names to facet the
#'   resulting plot by, or `NULL` (default) for auto-detection. When
#'   the matched effect is a pure 2-way `var1:var2`, no faceting is
#'   needed and the function behaves identically to the unfaceted
#'   form. When the matched effect is higher-order (e.g.,
#'   `cond:var1:var2`, `Group:var1:var2`, or the 4-way
#'   `Group:cond:var1:var2`), the var1xvar2 surface differs across
#'   levels of the extra factor(s), and a single unfaceted slice is
#'   misleading. Auto-detection tokenizes the matched effect string
#'   on `":"`, removes `var1`/`var2`, and keeps tokens that are
#'   columns in `data` — so a `cond:var1:var2` effect auto-facets by
#'   `cond`. One facet variable -> `facet_wrap`; two -> `facet_grid`
#'   on the first two; more -> `facet_wrap` on the interaction (with
#'   an info message). Pass an explicit character vector to override
#'   the auto-detection.
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
#' @param line_colors Character vector of colors, one per value of
#'   \code{mod_sd_levels} and in the same order. A vector of any other
#'   length errors. For lines in a grey, blue and navy scheme use
#'   \code{line_colors = c("grey65", "steelblue", "navyblue")} with three
#'   levels. Default \code{NULL} interpolates blue, gray and red across
#'   the sorted levels, so the levels -1, 0 and 1 give \code{"steelblue"},
#'   \code{"grey45"} and \code{"tomato"}.
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
#' @param output_dir Character scalar, the parent directory for PDF
#'   output. When `split_by_type = TRUE`, per-type subdirectories are
#'   created under this directory. Default `"."` (current working
#'   directory).
#' @param split_by_type Logical. Write each type's PDF into
#'   `<output_dir>/<tolower(type)>/` (`TRUE`, default) or write all
#'   PDFs flat into `output_dir` (`FALSE`).
#' @param show_points Logical. Overlay raw subject-level points on the
#'   simple-slope lines. Default `FALSE`: in the continuous x
#'   continuous case the drawn lines are synthetic moderator slices
#'   (Mean +/- 1SD), not actual subjects, so overlaying raw points
#'   suggests a 1:1 correspondence that does not exist. Set to `TRUE`
#'   for exploratory inspection.
#' @param points_aggregate One of `"auto"` (default) or `"never"`.
#'   Controls whether overlaid raw-data points are aggregated to
#'   subject-level. With long-format input (one row per
#'   Subject x within-subject factor) and between-subject moderators
#'   (e.g., a hormone z-score), plotting rows as-is over-plots each
#'   subject once per within-subject level — every subject contributes
#'   multiple dots stacked at their single x-value. `"auto"` groups
#'   by `id_col` and plots one mean dot per subject so the dot count
#'   matches N. Only takes effect when `show_points = TRUE`.
#' @param id_col Character scalar naming the subject-ID column in
#'   `data`. Required when `points_aggregate = "auto"` and
#'   `show_points = TRUE`. Default `"ID"`.
#' @param point_alpha Numeric in (0, 1]. Alpha transparency for raw
#'   points when `show_points = TRUE`. Default `0.35`.
#' @param show_ci Logical. Draw a pointwise 95% confidence band behind
#'   each line, filled with the line color at alpha 0.15. \code{FALSE}
#'   draws the lines alone. For lines alone in a grey, blue and navy
#'   scheme use \code{show_ci = FALSE} with
#'   \code{line_colors = c("grey65", "steelblue", "navyblue")}. Default
#'   \code{TRUE}.
#' @param x_label Named character vector of x-axis labels keyed by the
#'   values of \code{var1} and \code{var2}, for example
#'   \code{c(varA = "Variable A", varB = "Variable B")}. Names are needed
#'   because the two variables swap axes between \code{page_A} and
#'   \code{page_B}. A variable that is not named keeps its variable name
#'   on the x-axis. Names matching neither variable error. Default
#'   \code{NULL}.
#' @param moderator_label Named character vector keyed like
#'   \code{x_label}. On a page whose moderator is named, the value titles
#'   the legend and the level labels drop the variable name
#'   (\code{"-1SD"}, \code{"Mean"}, \code{"+1SD"}). A moderator that is
#'   not named keeps the legend title \code{"Moderator"} and level labels
#'   built from its variable name. Names matching neither variable error.
#'   Default \code{NULL}.
#' @param panel_labels Named character vector of facet-strip labels keyed
#'   by factor level, for example
#'   \code{c(levelA = "Level A", levelB = "Level B")}. Unnamed levels keep
#'   their level names. Names that match no level of the facet factors for
#'   the matched effects error, as does supplying labels when no effect
#'   is faceted. Panel order does not change. Default \code{NULL}.
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
#' @importFrom stats as.formula sd aggregate
#' @importFrom emmeans emmeans emtrends
#' @importFrom ggplot2 ggplot aes geom_line scale_color_manual
#' @importFrom ggplot2 scale_linetype_manual labs theme_minimal theme
#' @importFrom ggplot2 element_text facet_wrap facet_grid
#' @importFrom grDevices pdf dev.off colorRampPalette
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
    facet_by         = NULL,
    p_col            = "p_bh",
    p_threshold      = 0.05,
    gate             = TRUE,
    min_rows         = 5,
    y_label          = NULL,
    mod_sd_levels    = c(-1, 0, 1),
    line_colors      = NULL,
    line_types       = c("dashed", "solid", "longdash"),
    theme_base_size  = 14,
    pdf_prefix    = NULL,
    pdf_width        = 10,
    pdf_height       = 7,
    output_dir       = ".",
    split_by_type    = TRUE,
    show_points      = FALSE,
    points_aggregate = c("auto", "never"),
    id_col           = "ID",
    point_alpha      = 0.35,
    show_ci          = TRUE,
    x_label          = NULL,
    moderator_label  = NULL,
    panel_labels     = NULL
) {

  points_aggregate <- match.arg(points_aggregate)

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
  if (!is.logical(show_ci) || length(show_ci) != 1 || is.na(show_ci)) {
    stop("`show_ci` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(line_colors) &&
      length(line_colors) != length(mod_sd_levels)) {
    stop("`line_colors` must have one color per value of `mod_sd_levels` (",
         length(mod_sd_levels), "); got ", length(line_colors), ".",
         call. = FALSE)
  }
  if (is.null(line_colors)) {
    #blue-gray-red across the sorted levels so colors do not repeat
    blue_gray_red <- grDevices::colorRampPalette(
      c("steelblue", "grey45", "tomato")
    )
    line_colors <- blue_gray_red(length(mod_sd_levels))[
      rank(mod_sd_levels, ties.method = "first")
    ]
  }
  for (lab_arg in c("x_label", "moderator_label")) {
    lab <- get(lab_arg)
    if (is.null(lab)) next
    if (!is.character(lab) || is.null(names(lab)) || !all(nzchar(names(lab)))) {
      stop("`", lab_arg, "` must be a named character vector keyed by `var1` ",
           "and `var2`, because the two variables swap roles between pages.",
           call. = FALSE)
    }
    unmatched <- setdiff(names(lab), c(var1, var2))
    if (length(unmatched) > 0) {
      stop("`", lab_arg, "` names match neither `var1` nor `var2`: ",
           paste(unmatched, collapse = ", "), call. = FALSE)
    }
  }
  if (show_points && points_aggregate == "auto" &&
      !(id_col %in% colnames(data))) {
    stop("`points_aggregate = \"auto\"` requires `id_col = \"", id_col,
         "\"` to exist in `data`. Pass `id_col` explicitly or set ",
         "`points_aggregate = \"never\"`.", call. = FALSE)
  }
  if (!is.null(facet_by)) {
    if (!is.character(facet_by) || length(facet_by) == 0) {
      stop("`facet_by` must be a non-empty character vector or NULL.",
           call. = FALSE)
    }
    missing_facets <- setdiff(facet_by, colnames(data))
    if (length(missing_facets) > 0) {
      stop("`facet_by` names not found in `data`: ",
           paste(missing_facets, collapse = ", "), call. = FALSE)
    }
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

  #facet factors for an effect, from facet_by or parsed from the effect name
  derive_facet_vars <- function(eff) {
    if (!is.null(facet_by)) return(facet_by)
    tokens <- strsplit(as.character(eff), ":", fixed = TRUE)[[1]]
    candidates <- setdiff(tokens, c(var1, var2))
    candidates[candidates %in% colnames(data)]
  }

  #check against every matched effect so the result does not depend on the gate
  if (!is.null(panel_labels)) {
    if (!is.character(panel_labels) || is.null(names(panel_labels)) ||
        !all(nzchar(names(panel_labels)))) {
      stop("`panel_labels` must be a named character vector keyed by ",
           "factor level.", call. = FALSE)
    }
    panel_vars <- unique(unlist(lapply(matching_effects, derive_facet_vars)))
    known_levels <- unique(unlist(lapply(panel_vars, function(fv) {
      as.character(data[[fv]])
    })))
    unmatched <- setdiff(names(panel_labels), known_levels)
    if (length(unmatched) > 0) {
      stop("`panel_labels` names match no facet level: ",
           paste(unmatched, collapse = ", "), ". Available levels: ",
           if (length(known_levels) > 0) {
             paste(known_levels, collapse = ", ")
           } else {
             "none (no matched effect is faceted)"
           }, ".", call. = FALSE)
    }
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

  #level names with any requested replacements
  strip_labeller <- if (is.null(panel_labels)) {
    ggplot2::label_value
  } else {
    ggplot2::as_labeller(function(lv) {
      lab <- unname(panel_labels[lv])
      ifelse(is.na(lab), lv, lab)
    })
  }

  # --- Inner helper: build one simple-slopes page ---
  # focal_var on the x-axis, moderator_var split at mod_sd_levels SD.
  # facet_vars (character vector) names categorical factors to facet
  # on when the matched effect is higher-order (e.g., cond:v1:v2 or
  # the 4-way Group:cond:v1:v2). Empty -> unfaceted (original 2-way
  # behavior preserved byte-for-byte).
  build_slope_page <- function(fit, df_sub, focal_var, moderator_var,
                               ty, ar, eff, pval, vc, facet_vars) {

    mod_mean  <- mean(df_sub[[moderator_var]], na.rm = TRUE)
    mod_sd_v  <- stats::sd(df_sub[[moderator_var]], na.rm = TRUE)
    mod_values <- mod_mean + mod_sd_levels * mod_sd_v

    #a named moderator titles the legend so its levels drop the variable name
    mod_named <- moderator_var %in% names(moderator_label)
    lev_prefix <- if (mod_named) "" else paste0(moderator_var, " ")

    # Level labels: "<mod> +1SD" / "<mod> Mean" / "<mod> -1SD" etc.
    make_label <- function(sd_mult) {
      if (sd_mult == 0) {
        paste0(lev_prefix, "Mean")
      } else if (sd_mult > 0) {
        paste0(lev_prefix, "+", format(sd_mult, trim = TRUE), "SD")
      } else {
        paste0(lev_prefix, format(sd_mult, trim = TRUE), "SD")
      }
    }
    lev_labels <- vapply(mod_sd_levels, make_label, character(1))

    focal_range <- seq(min(df_sub[[focal_var]], na.rm = TRUE),
                       max(df_sub[[focal_var]], na.rm = TRUE),
                       length.out = 50)
    focal_m <- mean(df_sub[[focal_var]], na.rm = TRUE)
    trend_col <- paste0(focal_var, ".trend")

    # Reference grid: ~ 1 collapses every other factor in the model
    # (correct ONLY for a pure 2-way var1:var2 effect). For
    # higher-order effects we condition on the residual factor(s) so
    # the rendered slopes match what the F-test said is varying.
    if (length(facet_vars) > 0) {
      em_form <- stats::as.formula(
        paste0("~ ", paste(facet_vars, collapse = " + "))
      )
    } else {
      em_form <- ~ 1
    }

    #band grid adds the focal variable so emmeans returns one row per x
    band_form <- stats::as.formula(
      paste0("~ ", paste(c(facet_vars, focal_var), collapse = " + "))
    )

    preds_list <- lapply(seq_along(mod_values), function(j) {
      mv    <- mod_values[j]
      label <- lev_labels[j]

      at_slope <- setNames(list(mv), moderator_var)
      at_int   <- setNames(list(focal_m, mv),
                           c(focal_var, moderator_var))

      if (!is.null(vc)) {
        emt <- emmeans::emtrends(fit, em_form, var = focal_var,
                                 at = at_slope, vcov. = vc)
        em  <- emmeans::emmeans(fit, em_form, at = at_int, vcov. = vc)
      } else {
        emt <- emmeans::emtrends(fit, em_form, var = focal_var, at = at_slope)
        em  <- emmeans::emmeans(fit, em_form, at = at_int)
      }

      emt_df <- as.data.frame(summary(emt))
      em_df  <- as.data.frame(summary(em))

      band_df <- NULL
      if (show_ci) {
        at_band <- setNames(list(focal_range, mv),
                            c(focal_var, moderator_var))
        em_band <- if (!is.null(vc)) {
          emmeans::emmeans(fit, band_form, at = at_band, vcov. = vc)
        } else {
          emmeans::emmeans(fit, band_form, at = at_band)
        }
        band_df <- as.data.frame(summary(em_band, infer = c(TRUE, FALSE)))
      }

      add_band <- function(out, rows) {
        if (is.null(rows)) return(out)
        rows <- rows[order(rows[[focal_var]]), , drop = FALSE]
        lcl_col <- if ("lower.CL" %in% colnames(rows)) "lower.CL" else "asymp.LCL"
        ucl_col <- if ("upper.CL" %in% colnames(rows)) "upper.CL" else "asymp.UCL"
        out$lcl <- rows[[lcl_col]]
        out$ucl <- rows[[ucl_col]]
        out
      }

      if (length(facet_vars) == 0) {
        # Unfaceted 2-way path — preserve original behavior.
        slope     <- emt_df[[trend_col]]
        intercept <- em_df$emmean
        out <- data.frame(
          focal    = focal_range,
          y_pred   = intercept + slope * (focal_range - focal_m),
          ModLevel = factor(label, levels = lev_labels),
          stringsAsFactors = FALSE
        )
        add_band(out, band_df)
      } else {
        # One slope/intercept per facet-level combination.
        do.call(rbind, lapply(seq_len(nrow(emt_df)), function(k) {
          slope     <- emt_df[k, trend_col]
          intercept <- em_df[k, "emmean"]
          out <- data.frame(
            focal    = focal_range,
            y_pred   = intercept + slope * (focal_range - focal_m),
            ModLevel = factor(label, levels = lev_labels),
            stringsAsFactors = FALSE
          )
          for (fv in facet_vars) out[[fv]] <- emt_df[k, fv]
          band_rows <- NULL
          if (show_ci) {
            in_panel <- Reduce(`&`, lapply(facet_vars, function(fv) {
              as.character(band_df[[fv]]) == as.character(emt_df[k, fv])
            }))
            band_rows <- band_df[in_panel, , drop = FALSE]
          }
          add_band(out, band_rows)
        }))
      }
    })
    preds <- do.call(rbind, preds_list)

    #line_colors has one color per level, linetypes recycle
    n_levels   <- length(lev_labels)
    colors_used <- line_colors
    ltypes_used <- rep_len(line_types, n_levels)
    names(colors_used) <- lev_labels
    names(ltypes_used) <- lev_labels

    subtitle <- paste0("Effect: ", eff)
    if (!is.na(pval)) {
      subtitle <- paste0(subtitle, "  (", p_col, " = ",
                         format(round(pval, 4), nsmall = 4), ")")
    }

    p <- ggplot(preds, aes(x = focal, y = y_pred,
                           color = ModLevel, linetype = ModLevel))
    if (show_ci) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = lcl, ymax = ucl, fill = ModLevel),
        alpha = 0.15, color = NA
      )
    }
    p <- p + geom_line(linewidth = 1.2)

    if (show_points) {
      y_col <- if (length(as.character(formula)) >= 2) {
        as.character(formula)[2]
      } else {
        NULL
      }
      if (!is.null(y_col) && y_col %in% colnames(df_sub)) {
        # df_sub is long-format (one row per Subject x within-subject
        # factor). With between-subject moderators (e.g., a hormone
        # z-score), plotting rows as-is stacks multiple dots at each
        # subject's single x. When points_aggregate = "auto":
        #   - unfaceted: aggregate by ID -> one mean dot per subject.
        #   - faceted:   aggregate by ID x facet_vars -> one dot per
        #     subject per facet panel (preserves the within-subject
        #     factor that defines the panel).
        if (points_aggregate == "auto") {
          by_list <- c(list(.id = df_sub[[id_col]]),
                       setNames(lapply(facet_vars,
                                       function(fv) df_sub[[fv]]),
                                facet_vars))
          agg_df <- stats::aggregate(
            df_sub[, c(focal_var, y_col), drop = FALSE],
            by  = by_list,
            FUN = mean, na.rm = TRUE
          )
          raw <- data.frame(
            focal  = agg_df[[focal_var]],
            y_pred = agg_df[[y_col]]
          )
          for (fv in facet_vars) raw[[fv]] <- agg_df[[fv]]
        } else {
          raw <- data.frame(
            focal  = df_sub[[focal_var]],
            y_pred = df_sub[[y_col]]
          )
          for (fv in facet_vars) raw[[fv]] <- df_sub[[fv]]
        }
        p <- p + ggplot2::geom_point(
          data = raw, inherit.aes = FALSE,
          ggplot2::aes(x = focal, y = y_pred),
          alpha = point_alpha, color = "grey30", size = 1.2
        )
      }
    }

    x_lab <- if (focal_var %in% names(x_label)) {
      x_label[[focal_var]]
    } else {
      focal_var
    }
    legend_title <- if (mod_named) moderator_label[[moderator_var]] else "Moderator"

    final_plot <- p +
      scale_color_manual(values = colors_used) +
      scale_linetype_manual(values = ltypes_used) +
      labs(
        title    = paste0(ty, " | ", ar),
        subtitle = subtitle,
        x        = x_lab,
        y        = y_label,
        color    = legend_title,
        linetype = legend_title
      )
    if (show_ci) {
      final_plot <- final_plot +
        ggplot2::scale_fill_manual(values = colors_used) +
        labs(fill = legend_title)
    }

    # Faceting layer. One facet -> facet_wrap. Two -> facet_grid.
    # More -> facet_wrap on the combined interaction (with a notice).
    if (length(facet_vars) == 1) {
      final_plot <- final_plot + ggplot2::facet_wrap(
        stats::as.formula(paste0("~ ", facet_vars[1])),
        labeller = strip_labeller
      )
    } else if (length(facet_vars) == 2) {
      final_plot <- final_plot + ggplot2::facet_grid(
        stats::as.formula(paste0(facet_vars[1], " ~ ", facet_vars[2])),
        labeller = strip_labeller
      )
    } else if (length(facet_vars) > 2) {
      message("  [INFO] plot_simple_slopes: ", length(facet_vars),
              " facet variables; using facet_wrap on the combined ",
              "interaction. Consider passing `facet_by` explicitly.")
      final_plot <- final_plot + ggplot2::facet_wrap(
        stats::as.formula(paste0("~ ", paste(facet_vars, collapse = " + "))),
        labeller = strip_labeller
      )
    }

    final_plot +
      theme_minimal(base_size = theme_base_size) +
      theme(
        plot.title      = element_text(face = "bold"),
        legend.position = "bottom"
      )
  }

  # --- Main loop: build pages per type ---
  types_in_sig <- unique(as.character(sig_rows[[type_col]]))
  result <- setNames(vector("list", length(types_in_sig)), types_in_sig)
  for (ty in names(result)) result[[ty]] <- list()

  for (ty in types_in_sig) {
    sig_rows_ty <- sig_rows[sig_rows[[type_col]] == ty, , drop = FALSE]
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

        facet_vars <- derive_facet_vars(eff)

        page_A <- build_slope_page(fit, df_cell,
                                   focal_var = var1, moderator_var = var2,
                                   ty = ty, ar = ar, eff = eff,
                                   pval = pval, vc = vc,
                                   facet_vars = facet_vars)
        page_B <- build_slope_page(fit, df_cell,
                                   focal_var = var2, moderator_var = var1,
                                   ty = ty, ar = ar, eff = eff,
                                   pval = pval, vc = vc,
                                   facet_vars = facet_vars)

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

      out_subdir <- if (split_by_type) {
        file.path(output_dir, tolower(as.character(ty)))
      } else {
        output_dir
      }
      if (!dir.exists(out_subdir)) {
        dir.create(out_subdir, recursive = TRUE, showWarnings = FALSE)
      }
      pdf_name <- file.path(out_subdir,
                            paste0(pdf_prefix, "_", ty, "_SimpleSlopes.pdf"))
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

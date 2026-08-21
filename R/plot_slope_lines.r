#' @title Plot significant slopes from activation or gated-slope tests
#'
#' @description
#' Companion plot function for [run_slopes_vs_zero()] and
#' [run_gated_slopes()]. For every **significant row** in the input
#' results table, draws one figure showing the fitted slope line(s)
#' involved in that specific test, with an optional CI ribbon and
#' optional raw-data overlay. Designed for the "categorical by
#' continuous" side of second-level analyses — the discrete counterpart
#' to [plot_simple_slopes()].
#'
#' Two input shapes are supported and auto-detected:
#'
#' * **Activation input** — output of [run_slopes_vs_zero()]. Each
#'   significant row is plotted as a single slope line on the
#'   continuous moderator (`trend_var`), optionally filtered to the
#'   level named in the row's `predictor` column.
#' * **Gated post-hoc input** — output of [run_gated_slopes()]. Each
#'   significant row is plotted as the two slope lines involved in the
#'   pairwise contrast recorded in that row (e.g., the two levels of
#'   `Group` at one level of `cond`, for a three-way
#'   Group:cond:trend_var effect).
#'
#' Mode is detected automatically from the presence of the
#' `posthoc_type` column. Per-row (not per-cell) plotting keeps each
#' figure focused on exactly the effect that made the row significant
#' — no extra lines that were not part of the test.
#'
#' @details
#' \strong{Three-way decomposition.} For rows with
#' `posthoc_type = "<A>By<B>"` (a three-way interaction decomposed as
#' "contrast A within each level of B"), the plot shows two lines —
#' the two levels of A named in the row's `contrast` column — both at
#' the specific level of B recorded implicitly in the row's `predictor`
#' column. The conditioning level is parsed out of the `predictor`
#' value by stripping the contrast text and the `suffix`. Example:
#' `predictor = "lvl1_grpA - grpB_SlopeDiff"`,
#' `contrast = "grpA - grpB"`, `suffix = "SlopeDiff"` →
#' conditioning level `"lvl1"`.
#'
#' \strong{Raw-data overlay.} When `show_points = TRUE` (default), the
#' subject-level data points belonging to the group(s) being plotted
#' are overlaid on the lines — filtered to match the level of each
#' line. For activation plots that means subjects matching the row's
#' `predictor` level (e.g., only Eval subjects on the "Eval slope"
#' plot). For two-way and three-way post-hoc plots that means subjects
#' matching the contrast level (and conditioning level, for three-way).
#' Outliers and sample size become visible; the slope CIs are
#' less easy to take on faith.
#'
#' \strong{CI ribbon.} When `show_ci = TRUE` (default), a ribbon
#' around each line shows the pointwise 95% CI for the predicted
#' response, obtained from `emmeans::emmeans(fit, ..., at = list(grid))`
#' with `infer = c(TRUE, TRUE)`. For mixed models this uses the
#' Kenward-Roger / Satterthwaite machinery emmeans is already
#' configured with; for `lm` + HC\* the CI uses the robust covariance
#' passed as `vcov.`.
#'
#' @param results_df A data.frame in the canonical schema. Either the output
#'   of [run_slopes_vs_zero()] (activation mode — no `posthoc_type`
#'   column) or the output of [run_gated_slopes()] (post-hoc mode —
#'   has `posthoc_type`).
#' @param data The long-format data.frame used to fit the analysis.
#' @param formula Model formula (formula object preferred, character
#'   string accepted for backwards compatibility). Must include
#'   `trend_var`, and (for post-hoc mode) the categorical factors.
#' @param area_col Character scalar naming the area column.
#' @param trend_var Character scalar naming the continuous moderator
#'   (x-axis of every plot).
#' @param type_col Character scalar naming the measurement-type
#'   column. Default `"type"`.
#' @param model Required character scalar: `"lmer"` or `"lm"`.
#' @param vcov Character scalar; see [run_omnibus_by_area()]. The
#'   `NULL` conditional default applies.
#' @param p_col Character scalar naming the gate column in `results_df`.
#'   Default `"p_bh"` for activation input, `"p_bon"` for post-hoc
#'   input. Pass explicitly to avoid ambiguity.
#' @param p_threshold Numeric scalar. Rows with `p_col < p_threshold`
#'   are plotted. Default `0.05`.
#' @param by_factor Character scalar naming the level-column for
#'   activation mode (the `by_factor` passed to `run_slopes_vs_zero()`).
#'   Ignored in post-hoc mode. Default `NULL` (treat the row as an
#'   overall slope with no level filtering).
#' @param categorical_factors Character vector of categorical factor
#'   names used in post-hoc mode to parse `posthoc_type` labels back
#'   to factor names. Ignored in activation mode. Required in
#'   post-hoc mode.
#' @param suffix Character scalar matching the `suffix` argument that
#'   was passed to [run_gated_slopes()] when the input was produced.
#'   Used to strip predictor-column suffixes during parsing. Default
#'   `"SlopeDiff"`.
#' @param show_ci Logical; overlay pointwise 95% CI ribbon per line.
#'   Default `TRUE`.
#' @param show_points Logical; overlay raw-data points matching each
#'   line's group. Default `TRUE`.
#' @param points_aggregate One of `"auto"` (default) or `"never"`.
#'   Controls whether raw-row points are aggregated to subject-level
#'   before plotting. With long-format input (one row per
#'   Subject x within-subject factor), plotting rows as-is over-plots
#'   each subject once per within-subject level whenever the spec
#'   does not condition on every within-subject factor (e.g., an
#'   "Overall" main-effect line, or a Group-only filter when the
#'   model also has a Condition factor). `"auto"` groups each spec's
#'   filtered rows by `id_col` and plots one mean dot per subject,
#'   so the dots are at the same level of analysis as the line.
#'   This is a no-op when the spec already conditions on every
#'   within-subject factor (one row per subject after filter). Set
#'   `"never"` to preserve the legacy long-format dot cloud.
#' @param id_col Character scalar naming the subject-ID column in
#'   `data`. Required when `points_aggregate = "auto"` and
#'   `show_points = TRUE`. Default `"ID"`.
#' @param point_alpha Numeric in (0, 1]. Point transparency. Default
#'   `0.4`.
#' @param n_grid Integer. Number of x-axis grid points for line +
#'   ribbon rendering. Default `40`.
#' @param y_label Character scalar used as the y-axis label. Default
#'   `NULL` auto-detects from the formula LHS.
#' @param min_rows Integer. Skip cells with fewer rows after
#'   filtering. Default `5`.
#' @param pdf_prefix Character scalar file-name prefix for the PDF
#'   outputs, or `NULL` to skip file writing. When supplied, one PDF
#'   per chromophore type is written containing all pages for that
#'   type: `<output_dir>/<type>/<pdf_prefix>_<type>_SlopeLines.pdf`.
#' @param output_dir Character scalar, parent directory for PDF
#'   output. Default `"."`.
#' @param split_by_type Logical. Route per-type PDFs to
#'   `<output_dir>/<tolower(type)>/` (`TRUE`, default) or flat into
#'   `output_dir`.
#' @param pdf_width,pdf_height Numeric inches. Defaults 10, 7.
#' @param theme_base_size Numeric base size for `theme_minimal()`.
#'   Default 14.
#'
#' @return Invisibly returns a nested list keyed by type and then by
#'   a per-row cell key; each element is a single ggplot. When
#'   `pdf_prefix` is supplied, PDFs are also written as a side effect.
#'
#' @seealso [plot_simple_slopes()] for continuous x continuous;
#'   [run_slopes_vs_zero()] and [run_gated_slopes()] for statistical
#'   companions.
#'
#' @examples
#' \dontrun{
#' # Activation: plot Eval-only and NoEval-only testosterone slopes
#' # at each significant (type, area) cell.
#' df.slope.byGroup <- run_slopes_vs_zero(
#'   data = df.all, formula = AUC ~ Group * cond * zAveTest + (1|ID),
#'   area_col = "ROI", trend_var = "zAveTest",
#'   by_factor = "Group", model = "lmer"
#' )
#' df.slope.byGroup$p_bh <- p.adjust(df.slope.byGroup$p_value, "BH")
#' plot_slope_lines(
#'   results_df = df.slope.byGroup,
#'   data = df.all, formula = AUC ~ Group * cond * zAveTest + (1|ID),
#'   area_col = "ROI", trend_var = "zAveTest",
#'   by_factor = "Group", model = "lmer",
#'   pdf_prefix = "ROI_zAve_Test_SlopeVsZero_byGroup"
#' )
#'
#' # Gated post-hoc: plot just the pairs that drove each significant
#' # three-way decomposition row.
#' plot_slope_lines(
#'   results_df = df.h.posthoc,
#'   data = df.all, formula = AUC ~ Group * cond * zAveTest + (1|ID),
#'   area_col = "ROI", trend_var = "zAveTest",
#'   categorical_factors = c("Group", "cond"),
#'   model = "lmer", p_col = "p_bon", suffix = "TestSlopeDiff",
#'   pdf_prefix = "ROI_zAve_Test_PostHoc_SlopeLines"
#' )
#' }
#'
#' @importFrom stats as.formula sd setNames aggregate
#' @importFrom emmeans emmeans
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon geom_point
#' @importFrom ggplot2 scale_color_manual scale_fill_manual
#' @importFrom ggplot2 labs theme_minimal theme element_text
#' @importFrom grDevices pdf dev.off
#' @export
plot_slope_lines <- function(
    results_df,
    data,
    formula,
    area_col,
    trend_var,
    type_col            = "type",
    model,
    vcov                = NULL,
    p_col               = NULL,
    p_threshold         = 0.05,
    by_factor           = NULL,
    categorical_factors = NULL,
    suffix              = "SlopeDiff",
    show_ci             = TRUE,
    show_points         = TRUE,
    points_aggregate    = c("auto", "never"),
    id_col              = "ID",
    point_alpha         = 0.4,
    n_grid              = 40,
    y_label             = NULL,
    min_rows            = 5,
    pdf_prefix          = NULL,
    output_dir          = ".",
    split_by_type       = TRUE,
    pdf_width           = 10,
    pdf_height          = 7,
    theme_base_size     = 14
) {

  # --- Validate model / vcov ---
  if (missing(model)) {
    stop("Argument `model` is required. Choose \"lmer\" or \"lm\".",
         call. = FALSE)
  }
  model <- match.arg(model, c("lmer", "lm"))
  points_aggregate <- match.arg(points_aggregate)
  if (show_points && points_aggregate == "auto" &&
      !(id_col %in% colnames(data))) {
    stop("`points_aggregate = \"auto\"` requires `id_col = \"", id_col,
         "\"` to exist in `data`. Pass `id_col` explicitly or set ",
         "`points_aggregate = \"never\"`.", call. = FALSE)
  }
  valid_vcov <- c("classical", "HC0", "HC1", "HC2", "HC3", "HC4", "HC5")
  if (is.null(vcov)) {
    vcov <- if (model == "lmer") "classical" else "HC3"
  } else {
    vcov <- match.arg(vcov, valid_vcov)
  }
  if (model == "lmer" && vcov != "classical") {
    stop("`vcov = \"", vcov, "\"` is not valid when `model = \"lmer\"`.",
         call. = FALSE)
  }

  if (is.character(formula)) formula <- stats::as.formula(formula)
  if (is.null(y_label)) {
    fc <- as.character(formula)
    y_label <- if (length(fc) >= 2) fc[2] else "Response"
  }

  # --- Detect mode from columns ---
  is_posthoc <- "posthoc_type" %in% colnames(results_df)
  if (is_posthoc && (is.null(categorical_factors) ||
                      length(categorical_factors) == 0)) {
    stop("`categorical_factors` is required when `results_df` has a ",
         "`posthoc_type` column (gated post-hoc input).", call. = FALSE)
  }
  if (is.null(p_col)) {
    p_col <- if (is_posthoc) "p_bon" else "p_bh"
  }
  if (!(p_col %in% colnames(results_df))) {
    stop("`results_df` does not contain gate column `", p_col, "`.",
         call. = FALSE)
  }

  # --- Filter to significant rows ---
  sig <- results_df[!is.na(results_df[[p_col]]) &
                      results_df[[p_col]] < p_threshold, , drop = FALSE]
  if (nrow(sig) == 0) {
    message("No significant rows to plot (", p_col, " < ", p_threshold, ").")
    return(invisible(list()))
  }

  cap_first <- function(s) {
    if (is.null(s) || !nchar(s)) return(s)
    paste0(toupper(substr(s, 1, 1)), substr(s, 2, nchar(s)))
  }

  # --- Parse a post-hoc row into plotting specs ---
  parse_posthoc_row <- function(row) {
    pt <- row$posthoc_type
    if (identical(pt, "Overall")) {
      return(list(kind = "overall"))
    }
    # Two-way: posthoc_type matches cap_first(factor) for exactly one factor
    single <- categorical_factors[
      vapply(categorical_factors,
             function(f) identical(pt, cap_first(f)),
             logical(1))
    ]
    if (length(single) == 1) {
      pair <- strsplit(as.character(row$contrast), " - ", fixed = TRUE)[[1]]
      if (length(pair) != 2) {
        stop("Could not split contrast '", row$contrast,
             "' into two levels.", call. = FALSE)
      }
      return(list(kind = "twoway", contrast_factor = single, levels = pair))
    }
    # Three-way: "<A>By<B>"
    m <- regmatches(pt, regexec("^(.+)By(.+)$", pt))[[1]]
    if (length(m) == 3) {
      a_fac <- categorical_factors[
        vapply(categorical_factors,
               function(f) identical(cap_first(f), m[2]),
               logical(1))
      ]
      b_fac <- categorical_factors[
        vapply(categorical_factors,
               function(f) identical(cap_first(f), m[3]),
               logical(1))
      ]
      if (length(a_fac) != 1 || length(b_fac) != 1) {
        stop("Could not map posthoc_type '", pt,
             "' to categorical_factors.", call. = FALSE)
      }
      pair <- strsplit(as.character(row$contrast), " - ", fixed = TRUE)[[1]]
      if (length(pair) != 2) {
        stop("Could not split contrast '", row$contrast,
             "' into two levels.", call. = FALSE)
      }
      # Parse conditioning level from predictor
      pred <- as.character(row$predictor)
      strip <- paste0("_", row$contrast, "_", suffix, "$")
      by_level <- sub(strip, "", pred)
      return(list(kind = "threeway",
                  contrast_factor = a_fac,
                  by_factor_name  = b_fac,
                  levels          = pair,
                  by_level        = by_level))
    }
    stop("Unrecognized posthoc_type: '", pt, "'.", call. = FALSE)
  }

  # --- Fit model helper ---
  fit_cell <- function(df_cell) {
    if (model == "lmer") {
      lmerTest::lmer(formula, data = df_cell)
    } else {
      stats::lm(formula, data = df_cell)
    }
  }

  # --- Build one page (ggplot) for one sig row ---
  build_page <- function(fit, df_cell, row, vc) {
    ty <- as.character(row[[type_col]])
    ar <- as.character(row$area)

    # Decide the specs: a list of (label, at_list, data_filter_idx)
    specs <- list()

    if (is_posthoc) {
      ph <- parse_posthoc_row(row)
      if (ph$kind == "overall") {
        specs[[1]] <- list(label = "Overall",
                           at = list(),
                           filter = rep(TRUE, nrow(df_cell)))
      } else if (ph$kind == "twoway") {
        for (lv in ph$levels) {
          at_l <- stats::setNames(list(lv), ph$contrast_factor)
          specs[[length(specs) + 1]] <- list(
            label  = paste0(ph$contrast_factor, " = ", lv),
            at     = at_l,
            filter = df_cell[[ph$contrast_factor]] == lv
          )
        }
      } else {  # threeway
        for (lv in ph$levels) {
          at_l <- stats::setNames(list(lv, ph$by_level),
                                  c(ph$contrast_factor, ph$by_factor_name))
          specs[[length(specs) + 1]] <- list(
            label  = paste0(ph$contrast_factor, " = ", lv,
                            ", ", ph$by_factor_name, " = ", ph$by_level),
            at     = at_l,
            filter = df_cell[[ph$contrast_factor]] == lv &
                     df_cell[[ph$by_factor_name]] == ph$by_level
          )
        }
      }
    } else {
      # Activation mode: one line, filtered (optionally) to predictor level.
      # run_slopes_vs_zero() labels predictor as "<level>_<trend_var>_Slope"
      # (when by_factor supplied) or "<trend_var>_Slope" (when not), so strip
      # the "_<trend_var>_Slope" suffix to recover the plain factor level.
      pred_raw <- as.character(row$predictor)
      pred_lvl <- sub(paste0("_", trend_var, "_Slope$"), "", pred_raw)
      if (identical(pred_lvl, pred_raw) && !is.null(by_factor)) {
        # Suffix didn't match; warn and fall back to raw value so user can see
        message("  [INFO] plot_slope_lines: predictor '", pred_raw,
                "' did not match '_", trend_var, "_Slope' suffix; using raw.")
      }
      if (!is.null(by_factor)) {
        at_l <- stats::setNames(list(pred_lvl), by_factor)
        flt <- df_cell[[by_factor]] == pred_lvl
        lbl <- paste0(by_factor, " = ", pred_lvl)
      } else {
        at_l <- list()
        flt <- rep(TRUE, nrow(df_cell))
        lbl <- "Overall"
      }
      specs[[1]] <- list(label = lbl, at = at_l, filter = flt)
    }

    # Build the x-grid from the relevant subjects' range
    all_filter <- Reduce(`|`, lapply(specs, function(s) s$filter))
    xvals <- df_cell[[trend_var]][all_filter]
    if (length(xvals) < min_rows) return(NULL)
    xgrid <- seq(min(xvals, na.rm = TRUE),
                 max(xvals, na.rm = TRUE), length.out = n_grid)

    # For each spec, call emmeans at (at_l + trend_var = xgrid)
    line_rows <- list()
    for (s in specs) {
      at_full <- c(s$at, stats::setNames(list(xgrid), trend_var))
      em_form <- stats::as.formula(paste0("~ ", trend_var))
      em <- tryCatch({
        if (!is.null(vc)) {
          emmeans::emmeans(fit, em_form, at = at_full, vcov. = vc)
        } else {
          emmeans::emmeans(fit, em_form, at = at_full)
        }
      }, error = function(e) {
        message("  [WARNING] emmeans failed for ", ty, " / ", ar,
                " / ", s$label, ": ", conditionMessage(e))
        NULL
      })
      if (is.null(em)) next
      sm <- as.data.frame(summary(em, infer = c(TRUE, TRUE)))
      xcol <- trend_var
      ycol <- if ("emmean" %in% colnames(sm)) "emmean" else "estimate"
      lcl  <- if ("lower.CL" %in% colnames(sm)) "lower.CL" else "asymp.LCL"
      ucl  <- if ("upper.CL" %in% colnames(sm)) "upper.CL" else "asymp.UCL"
      line_rows[[length(line_rows) + 1]] <- data.frame(
        x       = sm[[xcol]],
        y       = sm[[ycol]],
        lcl     = sm[[lcl]],
        ucl     = sm[[ucl]],
        label   = s$label,
        stringsAsFactors = FALSE
      )
    }
    if (length(line_rows) == 0) return(NULL)
    preds <- do.call(rbind, line_rows)
    preds$label <- factor(preds$label, levels = vapply(specs, `[[`,
                                                       character(1), "label"))

    # Points (filtered per spec). With long-format input (one row per
    # Subject x within-subject level), plotting rows as-is over-plots
    # each subject once per within-subject level whenever the spec
    # does not condition on every within-subject factor (Overall,
    # by_factor = "Group", twoway Group posthoc). When
    # points_aggregate = "auto", group by id_col within each spec and
    # plot one mean dot per subject — a no-op when the filter already
    # leaves one row per subject (twoway Cond, threeway).
    pts_rows <- list()
    if (show_points) {
      y_var <- as.character(formula)[2]
      if (y_var %in% colnames(df_cell)) {
        do_agg <- points_aggregate == "auto"
        for (s in specs) {
          sub <- df_cell[s$filter, , drop = FALSE]
          if (nrow(sub) == 0) next
          if (do_agg) {
            agg_df <- stats::aggregate(
              sub[, c(trend_var, y_var), drop = FALSE],
              by  = list(.id = sub[[id_col]]),
              FUN = mean, na.rm = TRUE
            )
            pts_rows[[length(pts_rows) + 1]] <- data.frame(
              x = agg_df[[trend_var]], y = agg_df[[y_var]],
              label = s$label, stringsAsFactors = FALSE
            )
          } else {
            pts_rows[[length(pts_rows) + 1]] <- data.frame(
              x = sub[[trend_var]], y = sub[[y_var]],
              label = s$label, stringsAsFactors = FALSE
            )
          }
        }
      }
    }
    pts <- if (length(pts_rows) > 0) do.call(rbind, pts_rows) else NULL
    if (!is.null(pts)) {
      pts$label <- factor(pts$label, levels = levels(preds$label))
    }

    # Build subtitle with the row's effect and p-value
    p_val <- row[[p_col]]
    sub_parts <- character(0)
    if (!is.null(row$omnibus_effect)) {
      sub_parts <- c(sub_parts,
                     paste0("Effect: ", as.character(row$omnibus_effect)))
    } else if (!is.null(row$predictor)) {
      sub_parts <- c(sub_parts,
                     paste0("Predictor: ", as.character(row$predictor)))
    }
    if (!is.na(p_val)) {
      sub_parts <- c(sub_parts, paste0(p_col, " = ",
                                       format(round(p_val, 4), nsmall = 4)))
    }
    subtitle <- paste(sub_parts, collapse = "  |  ")

    # Colors — reuse a small palette
    n_lines <- nlevels(preds$label)
    palette_fill <- c("steelblue", "tomato", "darkolivegreen",
                     "mediumorchid", "goldenrod", "slategrey")
    cols <- rep_len(palette_fill, n_lines)
    names(cols) <- levels(preds$label)

    p <- ggplot2::ggplot(preds, ggplot2::aes(x = x, y = y, color = label,
                                             fill = label))
    if (show_ci) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = lcl, ymax = ucl),
        alpha = 0.15, color = NA
      )
    }
    p <- p + ggplot2::geom_line(linewidth = 1.2)
    if (!is.null(pts)) {
      p <- p + ggplot2::geom_point(
        data = pts, ggplot2::aes(x = x, y = y, color = label),
        inherit.aes = FALSE, alpha = point_alpha, size = 1.3
      )
    }
    p +
      ggplot2::scale_color_manual(values = cols) +
      ggplot2::scale_fill_manual(values = cols) +
      ggplot2::labs(
        title    = paste0(ty, " | ", ar),
        subtitle = subtitle,
        x        = trend_var,
        y        = y_label,
        color    = NULL,
        fill     = NULL
      ) +
      ggplot2::theme_minimal(base_size = theme_base_size) +
      ggplot2::theme(
        plot.title      = ggplot2::element_text(face = "bold"),
        legend.position = "bottom"
      )
  }

  # --- Iterate rows, grouped by type ---
  # Use [[type_col]] for type-column lookups so callers with a
  # non-default type_col (e.g., a capitalized "Type" column) match.
  # Hardcode `$area` because iterate_by_area always stamps
  # lowercase `area` on outputs regardless of input area_col.
  types_in_sig <- unique(as.character(sig[[type_col]]))
  result <- stats::setNames(vector("list", length(types_in_sig)),
                            types_in_sig)
  for (ty in types_in_sig) result[[ty]] <- list()

  for (ty in types_in_sig) {
    sig_ty <- sig[sig[[type_col]] == ty, , drop = FALSE]
    for (i in seq_len(nrow(sig_ty))) {
      row <- sig_ty[i, , drop = FALSE]
      ar  <- as.character(row$area)
      df_cell <- data[data[[type_col]] == ty & data[[area_col]] == ar, ,
                      drop = FALSE]
      if (nrow(df_cell) < min_rows) next

      tryCatch({
        fit <- fit_cell(df_cell)
        vc <- if (model == "lm" && vcov != "classical") {
          sandwich::vcovHC(fit, type = vcov)
        } else {
          NULL
        }
        g <- build_page(fit, df_cell, row, vc)
        if (!is.null(g)) {
          key_parts <- c(ar,
                         as.character(row$omnibus_effect),
                         as.character(row$posthoc_type),
                         as.character(row$predictor),
                         as.character(row$contrast))
          key_parts <- key_parts[!is.na(key_parts) & nzchar(key_parts)]
          cell_key <- paste(key_parts, collapse = "|")
          result[[ty]][[cell_key]] <- g
        }
      }, error = function(e) {
        message("  [WARNING] plot_slope_lines failed for ", ty,
                " / ", ar, ": ", conditionMessage(e))
      })
    }
  }

  # --- Optionally write per-type PDFs ---
  if (!is.null(pdf_prefix)) {
    for (ty in names(result)) {
      pages <- result[[ty]]
      if (length(pages) == 0) next
      out_subdir <- if (split_by_type) {
        file.path(output_dir, tolower(ty))
      } else {
        output_dir
      }
      if (!dir.exists(out_subdir)) {
        dir.create(out_subdir, recursive = TRUE, showWarnings = FALSE)
      }
      pdf_name <- file.path(out_subdir,
                            paste0(pdf_prefix, "_", ty, "_SlopeLines.pdf"))
      grDevices::pdf(pdf_name, width = pdf_width, height = pdf_height)
      for (nm in names(pages)) print(pages[[nm]])
      grDevices::dev.off()
      message("Wrote slope-lines PDF: ", pdf_name)
    }
  }

  invisible(result)
}

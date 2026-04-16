# Package-level miscellany. By convention this file is named zzz.R so
# that it sorts last during package load and any `.onLoad`-style hooks
# run after every other file's top-level assignments are in place.

# ---------------------------------------------------------------------
# Silence "no visible binding for global variable" R CMD check NOTEs.
# ---------------------------------------------------------------------
# The heatmap functions in this package use ggplot2's `aes()` and
# dplyr's `mutate()` / `filter()` with bare column names that refer to
# columns in the data frame being plotted or transformed. From R's
# static analysis perspective these look like undefined globals, so
# `R CMD check` emits a NOTE for each one. The check is safe to
# ignore -- the code runs correctly -- but the NOTEs clutter the
# check output. Declaring the names here tells R CMD check that they
# are intentional and should not be reported.
#
# Any new heatmap function that uses bare column names should add
# them to this list rather than sprinkling `.data$` pronouns through
# the function body. The house style for chFnirs heatmaps is bare
# aes() column references (see create_pvalue_heatmap, etc.).
utils::globalVariables(c(
  # Long-format melt columns used by every heatmap in this package
  "x_label",
  "y_label",
  "value",
  "color_value",
  "significance_label",
  # Post-hoc untested-cell marker (create_posthoc_heatmap)
  "status",
  # magrittr pipe self-reference in `. %>% filter(...)` inside
  # create_pvalue_heatmap's geom_tile(data = . %>% filter(...)) calls.
  ".",
  # Local data-frame columns constructed inside plot_simple_slopes()
  # and passed to aes() via bare references.
  "focal",
  "y_pred",
  "ModLevel",
  # Cell label column built inside create_effect_size_heatmap()
  "cell_label"
))

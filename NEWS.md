# chFnirs 0.0.0.9000

## plot_simple_slopes()

* New argument `show_ci`, default `TRUE`, draws a pointwise 95% confidence
  band behind each line, computed with `emmeans::emmeans()` over the plotted
  x values. The lines themselves are unchanged. The y-axis now extends to
  include the bands.
* `line_colors` now defaults to `NULL`, which interpolates blue, gray and red
  across the sorted `mod_sd_levels` (the default levels give `steelblue`,
  `grey45` and `tomato`) instead of `c("grey65", "steelblue", "navyblue")`.
  A supplied `line_colors` must have one color per value of `mod_sd_levels`;
  it is no longer recycled.
* New optional arguments `x_label`, `moderator_label` and `panel_labels` set
  the x-axis label, legend title and facet-strip labels. `x_label` and
  `moderator_label` are named by `var1`/`var2`, since those variables swap
  roles between pages. `panel_labels` is named by factor level. Unmatched
  names error. Defaults (`NULL`) keep variable and level names.
* To get the earlier appearance, call with `show_ci = FALSE` and
  `line_colors = c("grey65", "steelblue", "navyblue")` (for other numbers of
  levels, supply one color per level, for example with `rep_len()`).

## Documentation

* `run_gated_posthoc()` and `run_gated_slopes()`: the `adjust` argument help
  now states the adjustment scope. With `by_var`, or for effects decomposed by
  a second categorical factor, the adjustment applies within each level of
  that factor and not across the area; the two directions of a two-factor
  slope decomposition are adjusted separately.
* Tutorial: the gated post-hoc section now says what gating protects against
  and that it does not control error across effects, models, or signals. The
  Dmitrienko and Tamhane (2007) reference is removed.

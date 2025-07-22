#' Generic Plot Method for 'gam_influence' Objects
#'
#' Dispatches to specific plotting functions based on the `type` argument. This
#' is the primary user-facing function for creating visualizations.
#'
#' @param x A `gam_influence` object that has been processed by `calculate_influence()`.
#' @param type The type of plot to generate. Can be one of:
#'   - `"stan"`: Standardisation plot.
#'   - `"step"`: Step plot.
#'   - `"influ"`: Influence plot.
#'   - `"cdi"`: Coefficient-Distribution-Influence plot (requires `term` argument).
#'   - `"all"`: A combined layout of the step and influence plots.
#' @param ... Additional arguments passed to specific plot functions (e.g., `term` for cdi).
#' @return A ggplot object or a patchwork object for combined plots.
#' @export
plot.gam_influence <- function(x, type = "all", ...) {
  if (is.null(x$calculated)) {
    stop("Calculations not performed. Please run `calculate_influence()` first.")
  }

  switch(type,
    stan = plot_standardisation(x, ...),
    step = plot_stepwise_index(x, ...),
    influ = plot_term_influence(x, ...),
    cdi = plot_cdi(x, ...),
    all = plot_step_and_influence(x, ...),
    stop("Invalid plot type. Choose from 'stan', 'step', 'influ', 'cdi', 'all'.")
  )
}

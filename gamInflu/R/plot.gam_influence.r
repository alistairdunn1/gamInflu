#' @title Generic Plot Method for gam_influence Objects
#' @description Dispatches to specific plotting functions based on the `type` argument. This
#' is the primary user-facing function for creating visualisations.
#' @param x A `gam_influence` object that has been processed by `calculate_influence()`.
#' @param type The type of plot to generate. Can be one of:
#'   - `"stan"`: Standardisation plot.
#'   - `"step"`: Step plot.
#'   - `"influ"`: Influence plot.
#'   - `"cdi"`: Coefficient-Distribution-Influence plot (requires `term` argument).
#'   - `"distribution"`: Data distribution plot for a specific term (requires `term` argument).
#'   - `"all"`: A combined layout of the step and influence plots.
#' @param ... Additional arguments passed to specific plot functions (e.g., `term` for cdi and distribution).
#' @return A ggplot object or a patchwork object for combined plots.
#' @export
plot.gam_influence <- function(x, type = "all", ...) {
  if (is.null(x$calculated)) {
    stop("Calculations not performed. Please run `calculate_influence()` first.", call. = FALSE)
  }

  switch(type,
    stan = plot_standardisation(x, ...),
    step = plot_stepwise_index(x, ...),
    influ = plot_term_influence(x, ...),
    cdi = plot_cdi(x, ...),
    distribution = plot_term_distribution(x, ...),
    all = plot_step_and_influence(x, ...),
    stop("Invalid plot type. Choose from 'stan', 'step', 'influ', 'cdi', 'distribution', 'all'.", call. = FALSE)
  )
}

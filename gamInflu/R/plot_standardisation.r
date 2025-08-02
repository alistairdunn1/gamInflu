#' @title Standardisation Plot
#' @description Creates a standardisation plot comparing the unstandardised (raw) index to the final
#' standardised index for the focus term. Works with all supported GLM families (Gaussian, binomial,
#' gamma, Poisson) and displays family-appropriate confidence intervals.
#' @param obj A `gam_influence` object containing calculated indices from `calculate_influence()`.
#' @param show_unstandardised Logical. Should the unstandardised index be displayed? Default is TRUE.
#'   When FALSE, only the standardised index is shown and the legend is hidden for a cleaner plot.
#' @return A ggplot object showing both unstandardised and standardised indices, with confidence ribbon
#'   around the standardised index. The plot includes:
#'   - Unstandardised index (grey line and points): Raw aggregated values by focus level (if show_unstandardised = TRUE)
#'   - Standardised index (blue line and points): Model-adjusted values accounting for all terms
#'   - Confidence ribbon (light blue): Uncertainty bounds around the standardised index
#'   - Reference line at y=1: Baseline for relative comparison
#' @details
#' The standardisation plot visualizes the effect of model standardisation on the focus term index.
#' The unstandardised index shows the raw relationship, while the standardised index accounts for
#' all other model terms and their interactions.
#'
#' **Family-specific behaviour:**
#' - **Gaussian**: Traditional CPUE-style indices with geometric/arithmetic mean aggregation
#' - **Binomial**: Probability-based indices for presence/absence data
#' - **Gamma**: Positive continuous indices using geometric mean methods
#' - **Poisson**: Count-based indices with appropriate statistical methods
#'
#' The confidence intervals are calculated using the specified confidence level from
#' `calculate_influence()` (default 95%) and reflect the uncertainty in the model predictions.
#' @examples
#' \dontrun{
#' # Basic usage - shows both standardised and unstandardised
#' gi <- gam_influence(your_model, focus = "year")
#' gi <- calculate_influence(gi)
#' plot_standardisation(gi)
#'
#' # Show only standardised index for cleaner presentation
#' plot_standardisation(gi, show_unstandardised = FALSE)
#'
#' # With different families
#' # Binomial model
#' gi_binom <- calculate_influence(gi_binomial_model)
#' plot_standardisation(gi_binom) # Shows probability indices
#'
#' # Gamma model - standardised only
#' gi_gamma <- calculate_influence(gi_gamma_model)
#' plot_standardisation(gi_gamma, show_unstandardised = FALSE) # Clean standardised plot
#' }
#' @importFrom ggplot2 ggplot aes geom_hline geom_line geom_point geom_ribbon labs scale_colour_manual scale_y_continuous theme
#' @importFrom rlang .data
#' @export
plot_standardisation <- function(obj, show_unstandardised = TRUE) {
  df <- obj$calculated$indices
  if (is.null(df)) {
    stop("No indices calculated. Please run `calculate_influence()` first.", call. = FALSE)
  }
  if (!("unstan" %in% names(df)) || !("standardised_index" %in% names(df))) {
    stop("Data frame must contain 'unstan' and 'standardised_index' columns.", call. = FALSE)
  }
  if (!("level" %in% names(df))) {
    stop("Data frame must contain a 'level' column for the focus term.", call. = FALSE)
  }
  if (!("stan_lower" %in% names(df)) || !("stan_upper" %in% names(df))) {
    stop("Data frame must contain 'stan_lower' and 'stan_upper' columns for the standardised index range.", call. = FALSE)
  }
  if (is.null(obj$focus)) {
    stop("Focus term is not set. Please set the focus term in the gam_influence object.", call. = FALSE)
  }

  # Convert level to numeric if possible
  if (is.factor(df$level) && all(!is.na(as.numeric(as.character(levels(df$level)))))) {
    df$level <- as.numeric(as.character(df$level))
  }

  # Start building the plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$level, group = 1)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = 1), linetype = "dashed", colour = "grey")
  
  # Conditionally add unstandardised index
  if (show_unstandardised) {
    p <- p +
      ggplot2::geom_line(ggplot2::aes(y = .data$unstan, colour = "Unstandardised")) +
      ggplot2::geom_point(ggplot2::aes(y = .data$unstan, colour = "Unstandardised"))
  }
  
  # Add standardised index and confidence ribbon
  p <- p +
    ggplot2::geom_line(ggplot2::aes(y = .data$standardised_index, colour = "Standardised")) +
    ggplot2::geom_point(ggplot2::aes(y = .data$standardised_index, colour = "Standardised")) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$stan_lower, ymax = .data$stan_upper), fill = "royalblue", alpha = 0.2) +
    ggplot2::labs(x = obj$focus, y = "Index") +
    ggplot2::scale_y_continuous(limits = c(0, NA))
  
  # Conditionally add legend and colour scale
  if (show_unstandardised) {
    p <- p +
      ggplot2::labs(colour = "") +
      ggplot2::scale_colour_manual(values = c("Unstandardised" = "grey40", "Standardised" = "royalblue"))
  } else {
    # No legend when only showing standardised
    p <- p +
      ggplot2::scale_colour_manual(values = c("Standardised" = "royalblue"), guide = "none")
  }

  return(p)
}

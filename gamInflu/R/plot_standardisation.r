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
#' The standardisation plot visualises the effect of model standardisation on the focus term index.
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

#' @title Unstandardised Index Plot with Boxplots or Violin Plots
#' @description Creates a plot showing the raw year effect as violin plots or boxplots by year, with an optional
#' overlay of the standardised index as a line with confidence interval ribbon. The standardised
#' index is rescaled to have the same mean as the unstandardised data for direct comparison.
#' If the response was fitted on the log scale (islog = TRUE), both the raw data and standardised
#' indices are exp-transformed back to the original scale before plotting.
#' @param obj A `gam_influence` object containing calculated indices from `calculate_influence()`.
#' @param show_standardised Logical. Should the standardised index be displayed as an overlay? Default is TRUE.
#'   When TRUE, adds a line and confidence ribbon showing the rescaled standardised index for comparison.
#' @param plot_type Character. Type of plot to display raw data distribution. Options are "violin" (default) or "boxplot".
#'   - "violin": Shows kernel density estimation of the distribution shape
#'   - "boxplot": Shows median, quartiles, whiskers, and outliers
#' @return A ggplot object showing:
#'   - Violin plots or boxplots: Distribution of raw observations by focus term level (typically year)
#'   - Rescaled standardised index line (blue): Model-adjusted values rescaled to match unstandardised mean (if show_standardised = TRUE)
#'   - Confidence ribbon (light blue): Uncertainty bounds around rescaled standardised index (if show_standardised = TRUE)
#'   - Automatic exp-transformation: Applied when islog = TRUE to show data on original scale
#' @details
#' This plot provides insight into the raw data distribution within each level of the focus term
#' (typically years) and how the model's standardised index relates to this underlying variability.
#'
#' **Plot types:**
#' - **Violin plots** show: Kernel density estimation revealing the full distribution shape, including multimodality and skewness (default)
#' - **Boxplots** show: Median (centre line), quartiles (box boundaries), whiskers (1.5 × IQR), and outliers (points beyond whiskers)
#'
#' **Log-scale handling:**
#' When the response was fitted on the log scale (obj$islog = TRUE), the function automatically:
#' - Applies exp() transformation to raw data before plotting for interpretability
#' - Uses standardised indices as-is (already transformed by calculate_influence())
#'
#' The standardised index overlay is rescaled to have the same mean as the unstandardised data,
#' allowing direct visual comparison between the raw data patterns and the model-derived trend.
#'
#' Rescaling formula: rescaled_index = standardised_index * (mean(unstandardised) / mean(standardised))
#' @examples
#' \dontrun{
#' # Basic usage - shows violin plots with rescaled standardised overlay
#' gi <- gam_influence(your_model, focus = "year")
#' gi <- calculate_influence(gi)
#' plot_unstandardised(gi)
#'
#' # Use boxplots for traditional quartile visualization
#' plot_unstandardised(gi, plot_type = "boxplot")
#'
#' # Show only violin plots without standardised overlay
#' plot_unstandardised(gi, show_standardised = FALSE)
#'
#' # For log-scale models, data is automatically exp-transformed
#' gi_log <- gam_influence(log_model, focus = "year") # islog will be TRUE
#' gi_log <- calculate_influence(gi_log)
#' plot_unstandardised(gi_log) # Shows exp-transformed data on original scale
#'
#' # Compare patterns between raw data and rescaled model trend
#' plot_unstandardised(gi, show_standardised = TRUE, plot_type = "violin")
#' }
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_violin geom_line geom_ribbon labs scale_y_continuous theme_minimal
#' @importFrom rlang .data
#' @export
plot_unstandardised <- function(obj, show_standardised = TRUE, plot_type = "violin") {
  # Validate inputs
  if (is.null(obj$data)) {
    stop("No data found in gam_influence object.", call. = FALSE)
  }
  if (is.null(obj$focus)) {
    stop("Focus term is not set in gam_influence object.", call. = FALSE)
  }
  if (is.null(obj$response)) {
    stop("Response variable is not set in gam_influence object.", call. = FALSE)
  }

  # Validate plot_type parameter
  if (!plot_type %in% c("boxplot", "violin")) {
    stop("plot_type must be either 'boxplot' or 'violin'.", call. = FALSE)
  }

  # Get the raw data for boxplots
  focus_var <- obj$focus
  response_var <- obj$response

  if (!(focus_var %in% names(obj$data))) {
    stop("Focus variable '", focus_var, "' not found in data.", call. = FALSE)
  }
  if (!(response_var %in% names(obj$data))) {
    stop("Response variable '", response_var, "' not found in data.", call. = FALSE)
  }

  # Get plot data
  plot_data <- obj$data

  # Check if response is on log scale and transform if necessary
  islog <- isTRUE(obj$islog)
  if (islog) {
    # Transform response to original scale using exp()
    plot_data[[response_var]] <- exp(plot_data[[response_var]])
    # Update y-axis label to reflect transformation
    y_label <- paste(response_var, "(exp-transformed)")
  } else {
    y_label <- "Index"
  }

  # Calculate the mean of unstandardised data for rescaling (after transformation if needed)
  unstd_mean <- mean(plot_data[[response_var]], na.rm = TRUE)

  # Convert focus variable to numeric if possible (same approach as plot_standardisation)
  if (is.factor(plot_data[[focus_var]]) && all(!is.na(as.numeric(as.character(levels(plot_data[[focus_var]])))))) {
    plot_data[[focus_var]] <- as.numeric(as.character(plot_data[[focus_var]]))
  }

  # Create the base plot with either boxplots or violin plots
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[[focus_var]], y = .data[[response_var]])) +
    ggplot2::labs(x = focus_var, y = y_label)

  # Add the appropriate plot type
  if (plot_type == "boxplot") {
    p <- p + ggplot2::geom_boxplot(aes(group = .data[[focus_var]]), fill = "grey90", colour = "grey60", alpha = 0.7)
  } else if (plot_type == "violin") {
    p <- p + ggplot2::geom_violin(aes(group = .data[[focus_var]]), fill = "grey90", colour = "grey60", alpha = 0.7)
  }

  # Add standardised index overlay if requested
  if (show_standardised) {
    # Get standardised index data
    indices_df <- obj$calculated$indices
    if (is.null(indices_df)) {
      warning("No calculated indices found. Run `calculate_influence()` first to show standardised overlay.")
    } else {
      # Check required columns exist
      required_cols <- c("level", "standardised_index", "stan_lower", "stan_upper")
      if (all(required_cols %in% names(indices_df))) {
        # Convert level to numeric if possible (same approach as plot_standardisation)
        if (is.factor(indices_df$level) && all(!is.na(as.numeric(as.character(levels(indices_df$level)))))) {
          indices_df$level <- as.numeric(as.character(indices_df$level))
        }

        # Note: standardised indices are already transformed by calculate_influence() when islog=TRUE
        # so no additional transformation needed for indices_df

        # Calculate rescaling factor to match unstandardised mean
        std_mean <- mean(indices_df$standardised_index, na.rm = TRUE)
        rescale_factor <- unstd_mean / std_mean

        # Rescale standardised index and confidence bounds
        indices_df$standardised_index_rescaled <- indices_df$standardised_index * rescale_factor
        indices_df$stan_lower_rescaled <- indices_df$stan_lower * rescale_factor
        indices_df$stan_upper_rescaled <- indices_df$stan_upper * rescale_factor

        # Add standardised index as line and ribbon
        p <- p +
          ggplot2::geom_ribbon(
            data = indices_df,
            ggplot2::aes(
              x = .data$level,
              ymin = .data$stan_lower_rescaled,
              ymax = .data$stan_upper_rescaled,
              group = 1
            ),
            fill = "royalblue", alpha = 0.3, inherit.aes = FALSE
          ) +
          ggplot2::geom_line(
            data = indices_df,
            ggplot2::aes(x = .data$level, y = .data$standardised_index_rescaled, group = 1),
            colour = "royalblue", inherit.aes = FALSE
          )
      } else {
        warning(
          "Required columns for standardised index not found. Missing: ",
          paste(setdiff(required_cols, names(indices_df)), collapse = ", ")
        )
      }
    }
  }

  return(p)
}

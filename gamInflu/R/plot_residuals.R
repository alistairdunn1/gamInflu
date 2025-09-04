#' @title Residuals Plot
#' @description Creates residual plots for GAM models with options for standard GAM residuals
#' and violin plots of residuals by focus term levels. Supports all GLM families and provides
#' comprehensive residual diagnostics.
#' @param obj A `gam_influence` object containing the fitted GAM model and data.
#' @param type Character. Type of residual plot to create:
#'   - `"standard"`: Standard GAM residual plots (4-panel diagnostic plot)
#'   - `"violin"`: Violin plot of residuals by focus term levels
#'   - `"both"`: Both standard and violin plots (combined layout)
#' @param residual_type Character. Type of residuals to calculate:
#'   - `"deviance"`: Deviance residuals (default, family-appropriate)
#'   - `"pearson"`: Pearson residuals
#'   - `"response"`: Response residuals (observed - fitted)
#'   - `"working"`: Working residuals

#' @param violin_trim Logical. Should violin plots be trimmed to data range? Default is TRUE.
#' @param add_boxplot Logical. Should boxplots be overlaid on violin plots? Default is TRUE.
#' @param add_smooth Logical. Should a loess smooth line be added to the plot? Default is TRUE.
#' @param by Character. Optional variable name for faceting plots. If provided, creates separate
#' panels for each level of this variable. Must be a factor or character variable (not numeric).
#' @param ... Additional arguments passed to ggplot2 functions.
#' @return A ggplot object or patchwork object (for combined plots) showing:
#'   **Standard plots:**
#'   - Residuals vs fitted values (or linear predictor for deviance residuals)
#'   - Q-Q plot of residuals
#'   - Scale-location plot (sqrt(|residuals|) vs fitted values or linear predictor)
#'   - Residuals vs leverage (Cook's distance contours)
#'
#'   **Violin plots:**
#'   - Distribution of residuals by focus term levels
#'   - Optional boxplots
#' @details
#' **Residual Types:**
#' - **Deviance residuals**: Default choice, appropriate for all GLM families. Plotted against linear predictor.
#' - **Pearson residuals**: Standardised residuals, useful for checking variance assumptions. Plotted against fitted values.
#' - **Response residuals**: Simple observed - fitted, interpretable but may not be appropriate for all families. Plotted against fitted values.
#' - **Working residuals**: Used in iterative fitting, less commonly needed for diagnostics. Plotted against fitted values.
#'
#' **Family-specific behaviour:**
#' - **Gaussian**: All residual types available, deviance = Pearson for identity link
#' - **Binomial**: Deviance residuals recommended for binary data
#' - **Gamma**: Deviance residuals account for variance structure
#' - **Poisson**: Deviance residuals appropriate for count data
#' - **Tweedie**: Deviance residuals handle compound Poisson-Gamma structure
#'
#' **Violin Plots:**
#' Show the full distribution of residuals within each level of the focus term,
#' revealing patterns that summary statistics might miss. If the focus term can be
#' converted to numeric (e.g., years stored as factors), the x-axis will display
#' as a numeric sequence for better readability.
#' @examples
#' \dontrun{
#' # Basic usage - standard residual plots
#' gi <- gam_influence(your_model, focus = "year")
#' plot_residuals(gi, type = "standard")
#'
#' # Violin plot of residuals by focus levels
#' plot_residuals(gi, type = "violin")
#'
#' # Both plots combined
#' plot_residuals(gi, type = "both")
#'
#' # Using different residual types
#' plot_residuals(gi, type = "violin", residual_type = "pearson")
#' plot_residuals(gi, type = "standard", residual_type = "response")
#'
#' # Customised violin plot
#' plot_residuals(gi,
#'   type = "violin",
#'   add_boxplot = FALSE
#' )
#'
#' # Numeric focus terms (e.g., years) will display as numeric sequence
#' gi_year <- gam_influence(your_model, focus = "year") # year as factor
#' plot_residuals(gi_year, type = "violin") # x-axis shows as 1990, 1991, 1992...
#'
#' # Faceted plots by another variable
#' plot_residuals(gi, type = "violin", by = "gear_type")
#' plot_residuals(gi, type = "both", by = "area")
#'
#' # For different model families
#' # Binomial model
#' plot_residuals(gi_binomial, type = "violin") # Shows residual patterns for binary data
#'
#' # Gamma model
#' plot_residuals(gi_gamma, type = "standard") # Family-appropriate diagnostics
#' }
#' @importFrom ggplot2 ggplot aes geom_point geom_violin geom_boxplot geom_smooth geom_abline geom_qq geom_qq_line labs theme element_text facet_wrap scale_fill_manual scale_x_continuous
#' @importFrom patchwork wrap_plots plot_annotation
#' @importFrom stats residuals fitted predict qqnorm qqline hatvalues cooks.distance quantile
#' @importFrom tools toTitleCase
#' @importFrom rlang .data
#' @export
plot_residuals <- function(obj,
                           type = c("standard", "violin", "both"),
                           residual_type = c("deviance", "pearson", "response", "working"),
                           violin_trim = TRUE,
                           add_boxplot = TRUE,
                           add_smooth = TRUE,
                           by = NULL,
                           ...) {
  # Validate inputs
  if (!inherits(obj, "gam_influence")) {
    stop("Object must be of class 'gam_influence'", call. = FALSE)
  }
  if (is.null(obj$model)) {
    stop("No model found in gam_influence object", call. = FALSE)
  }
  if (is.null(obj$data)) {
    stop("No data found in gam_influence object", call. = FALSE)
  }
  if (is.null(obj$focus)) {
    stop("Focus term not set in gam_influence object", call. = FALSE)
  }

  type <- match.arg(type)
  residual_type <- match.arg(residual_type)

  # Validate by parameter
  if (!is.null(by)) {
    if (!is.character(by) || length(by) != 1) {
      stop("'by' must be a single character string specifying a variable name", call. = FALSE)
    }
    if (!by %in% names(obj$data)) {
      stop("Variable '", by, "' not found in data", call. = FALSE)
    }
    if (is.numeric(obj$data[[by]])) {
      stop("Variable '", by, "' must not be numeric. Use factor or character variables for faceting.", call. = FALSE)
    }
  }

  # Calculate residuals and fitted values/linear predictor
  model <- obj$model
  fitted_vals <- fitted(model)

  # Get linear predictor (systematic component before link function)
  linear_predictor <- predict(model, type = "link")

  # Get residuals based on type
  residuals_vals <- switch(residual_type,
    "deviance" = residuals(model, type = "deviance"),
    "pearson" = residuals(model, type = "pearson"),
    "response" = residuals(model, type = "response"),
    "working" = residuals(model, type = "working")
  )

  # Create base data frame
  focus_values <- obj$data[[obj$focus]]

  # Try to convert focus term to numeric if possible (for proper x-axis sequencing)
  focus_numeric <- NULL
  if (is.factor(focus_values) || is.character(focus_values)) {
    # Attempt numeric conversion
    numeric_attempt <- suppressWarnings(as.numeric(as.character(focus_values)))
    if (!any(is.na(numeric_attempt)) && length(unique(numeric_attempt)) > 1) {
      focus_numeric <- numeric_attempt
    }
  } else if (is.numeric(focus_values)) {
    focus_numeric <- focus_values
  }

  resid_data <- data.frame(
    fitted = fitted_vals,
    linear_predictor = linear_predictor,
    residuals = residuals_vals,
    focus_level = focus_values,
    focus_numeric = focus_numeric,
    observation = seq_along(fitted_vals)
  )

  # Add by variable if specified
  if (!is.null(by)) {
    resid_data$by_var <- obj$data[[by]]
  }

  # Generate plots based on type
  if (type == "standard") {
    plot_result <- create_standard_residual_plots(resid_data, residual_type, model, by, ...)
  } else if (type == "violin") {
    plot_result <- create_violin_residual_plot(
      resid_data, obj$focus, residual_type,
      violin_trim, add_boxplot, by, ...
    )
  } else if (type == "both") {
    standard_plot <- create_standard_residual_plots(
      resid_data, residual_type, model, by, ...
    )
    violin_plot <- create_violin_residual_plot(
      resid_data, obj$focus, residual_type,
      violin_trim, add_boxplot, by, ...
    )

    plot_result <- patchwork::wrap_plots(standard_plot, violin_plot, ncol = 1, heights = c(2, 1))
  }

  return(plot_result)
}

#' @title Create Standard GAM Residual Plots
#' @description Internal function to create the standard 4-panel GAM diagnostic plots
#' @param resid_data Data frame containing residuals and fitted values
#' @param residual_type Character string indicating residual type
#' @param model The GAM model object
#' @param by Character. Optional variable name for faceting
#' @param ... Additional arguments
#' @return A patchwork object with 4 diagnostic panels
#' @noRd
create_standard_residual_plots <- function(resid_data, residual_type, model, by = NULL, add_smooth = TRUE, ...) {
  # Determine x-axis variable based on residual type
  # For deviance residuals, use linear predictor; for others, use fitted values
  if (residual_type == "deviance") {
    x_var <- resid_data$linear_predictor
    x_label <- "Linear Predictor"
  } else {
    x_var <- resid_data$fitted
    x_label <- "Fitted Values"
  }

  # Panel 1: Residuals vs Fitted/Linear Predictor
  p1 <- ggplot2::ggplot(resid_data, ggplot2::aes(x = x_var, y = .data$residuals)) +
    ggplot2::geom_point(alpha = 0.4, colour = "royalblue", size = 0.7) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::labs(x = x_label, y = paste(tools::toTitleCase(residual_type), "Residuals"))

  if (add_smooth) {
    p1 <- p1 + ggplot2::geom_smooth(method = "loess", se = FALSE, colour = "red")
  }

  if (!is.null(by)) {
    p1 <- p1 + ggplot2::facet_wrap(~ .data$by_var)
  }

  # Panel 2: Q-Q Plot
  p2 <- ggplot2::ggplot(resid_data, ggplot2::aes(sample = .data$residuals)) +
    ggplot2::geom_qq(alpha = 0.4, colour = "royalblue", size = 0.7) +
    ggplot2::geom_qq_line(colour = "black") +
    ggplot2::labs(x = "Theoretical Quantiles", y = "Sample Quantiles")

  if (!is.null(by)) {
    p2 <- p2 + ggplot2::facet_wrap(~ .data$by_var)
  }

  # Panel 3: Scale-Location (sqrt of absolute residuals vs fitted/linear predictor)
  resid_data$sqrt_abs_resid <- sqrt(abs(resid_data$residuals))
  p3 <- ggplot2::ggplot(resid_data, ggplot2::aes(x = x_var, y = .data$sqrt_abs_resid)) +
    ggplot2::geom_point(alpha = 0.4, colour = "royalblue", size = 0.7) +
    ggplot2::labs(x = x_label, y = expression(sqrt("|Residuals|")))

  if (add_smooth) {
    p3 <- p3 + ggplot2::geom_smooth(method = "loess", se = FALSE, colour = "red")
  }

  if (!is.null(by)) {
    p3 <- p3 + ggplot2::facet_wrap(~ .data$by_var)
  }

  # Panel 4: Residuals vs Leverage (with Cook's distance)
  # Calculate leverage and Cook's distance - handle GAM models safely
  p4 <- tryCatch(
    {
      leverage <- hatvalues(model)
      cooks_d <- cooks.distance(model)

      resid_data$leverage <- leverage
      resid_data$cooks_distance <- cooks_d

      # Identify high leverage and high Cook's distance points
      resid_data$high_influence <- (leverage > 2 * mean(leverage, na.rm = TRUE)) |
        (cooks_d > 4 / length(cooks_d))

      ggplot2::ggplot(resid_data, ggplot2::aes(x = .data$leverage, y = .data$residuals)) +
        ggplot2::geom_point(ggplot2::aes(colour = .data$high_influence, size = .data$cooks_distance), alpha = 0.6, size = 0.7) +
        {
          if (add_smooth) {
            ggplot2::geom_smooth(method = "loess", se = FALSE, colour = "red")
          } else {
            NULL
          }
        } +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
        ggplot2::scale_colour_manual(
          values = c("FALSE" = "black", "TRUE" = "red"),
          name = "High Influence"
        ) +
        ggplot2::scale_size_continuous(name = "Cook's Distance") +
        ggplot2::labs(x = "Leverage", y = paste(tools::toTitleCase(residual_type), "Residuals")) +
        {
          if (!is.null(by)) ggplot2::facet_wrap(~ .data$by_var) else NULL
        }
    },
    error = function(e) {
      # Fallback plot without leverage/Cook's distance if calculation fails
      ggplot2::ggplot(resid_data, ggplot2::aes(x = .data$observation, y = .data$residuals)) +
        ggplot2::geom_point(alpha = 0.6, colour = "royalblue", size = 0.7) +
        {
          if (add_smooth) ggplot2::geom_smooth(method = "loess", se = FALSE, colour = "red") else NULL
        } +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
        ggplot2::labs(x = "Observation Index", y = paste(tools::toTitleCase(residual_type), "Residuals")) +
        {
          if (!is.null(by)) ggplot2::facet_wrap(~ .data$by_var) else NULL
        }
    }
  )

  # Combine plots
  combined_plot <- patchwork::wrap_plots(p1, p2, p3, p4, ncol = 2)

  return(combined_plot)
}

#' @title Create Violin Plot of Residuals by Focus Levels
#' @description Internal function to create violin plot showing residual distribution by focus term levels
#' @param resid_data Data frame containing residuals and focus levels
#' @param focus_name Name of the focus term
#' @param residual_type Character string indicating residual type
#' @param violin_trim Logical indicating whether to trim violins to data range
#' @param add_boxplot Logical indicating whether to add boxplots
#' @param by Character. Optional variable name for faceting
#' @param ... Additional arguments
#' @return A ggplot object showing violin plot
#' @noRd
create_violin_residual_plot <- function(resid_data, focus_name, residual_type,
                                        violin_trim, add_boxplot, by = NULL, ...) {
  # Determine whether to use numeric or original focus values for x-axis
  use_numeric <- !is.null(resid_data$focus_numeric) && !any(is.na(resid_data$focus_numeric))

  # Initialise mapping variable
  unique_mappings <- NULL

  # For violin plots, we need discrete groups even if the underlying data is numeric
  # We'll use the original focus_level for violin creation but apply numeric scaling for the axis
  if (use_numeric) {
    # Create factor version for violin plotting but keep numeric for axis
    resid_data$focus_factor <- factor(resid_data$focus_level)
    # Create mapping from factor levels to numeric values for axis labeling
    level_to_numeric <- setNames(resid_data$focus_numeric, as.character(resid_data$focus_level))
    unique_mappings <- level_to_numeric[!duplicated(names(level_to_numeric))]

    # Create the base violin plot using factor levels
    p <- ggplot2::ggplot(resid_data, ggplot2::aes(x = .data$focus_factor, y = .data$residuals))
  } else {
    # Create the base violin plot with original factor/character levels
    p <- ggplot2::ggplot(resid_data, ggplot2::aes(x = .data$focus_level, y = .data$residuals))
  }

  # Add violin plot
  p <- p + ggplot2::geom_violin(
    fill = "royalblue",
    trim = violin_trim,
    alpha = 0.7,
    show.legend = FALSE
  )

  # Add boxplot if requested
  if (add_boxplot) {
    p <- p + ggplot2::geom_boxplot(width = 0.2, alpha = 0.8, show.legend = FALSE)
  }

  # Add reference line at y = 0
  p <- p + ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50")

  # Formatting
  p <- p +
    ggplot2::labs(
      x = focus_name,
      y = paste(tools::toTitleCase(residual_type), "Residuals")
    ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  # Add faceting if by variable is specified
  if (!is.null(by)) {
    p <- p + ggplot2::facet_wrap(~ .data$by_var)
  }

  # Apply appropriate x-axis scaling
  if (use_numeric) {
    # Use numeric values as axis labels while keeping factor structure for violins
    numeric_labels <- sort(unique(resid_data$focus_numeric))
    factor_labels <- names(sort(unique_mappings))
    p <- p + ggplot2::scale_x_discrete(labels = setNames(numeric_labels, factor_labels))
  }

  return(p)
}

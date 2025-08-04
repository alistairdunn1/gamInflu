#' @title Binomial Model Visualization
#' @description Creates a specialized plot for binomial models showing both raw proportions
#' and model predictions of positive outcomes by the focus term on the same plot.
#' @param obj A `gam_influence` object calculated for a binomial model.
#' @param confidence_level Numeric. Confidence level for confidence intervals (default 0.95).
#' @return A ggplot object showing both observed proportions and predicted probabilities.
#' @details
#' This function is specifically designed for binomial models and creates a single plot
#' with two lines:
#'
#' **Blue line (Observed)** - Shows the observed proportion of positive outcomes
#' (1s vs 0s) for each level of the focus term, with binomial confidence intervals.
#'
#' **Red line (Predicted)** - Shows the model-predicted probabilities for each
#' level of the focus term, with confidence intervals based on prediction standard errors.
#'
#' The function expects the response variable to be binary (0/1) or proportional (0-1).
#' For binary data, proportions are calculated automatically. The plot helps assess
#' how well the model captures the observed patterns in the data.
#' @examples
#' \dontrun{
#' # Binomial model example
#' data$presence <- rbinom(nrow(data), 1, 0.3) # Binary response
#' data$year <- factor(data$year)
#' mod <- mgcv::gam(presence ~ year + s(depth), data = data, family = binomial())
#' gi <- gam_influence(mod, focus = "year")
#' gi <- calculate_influence(gi, family = "binomial")
#' plot_binomial(gi)
#' }
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_ribbon scale_color_manual scale_fill_manual
#' @importFrom ggplot2 labs theme_minimal scale_y_continuous scale_x_continuous theme .data
#' @importFrom scales percent_format
#' @export
plot_binomial <- function(obj, confidence_level = 0.95) {
  # --- Validation ---
  if (!inherits(obj, "gam_influence")) {
    stop("Object must be of class 'gam_influence'.", call. = FALSE)
  }

  if (is.null(obj$calculated)) {
    stop("Please run calculate_influence() on the object first.", call. = FALSE)
  }

  # Check if model family is appropriate for binomial plotting
  model_family <- family(obj$model)$family
  if (!grepl("binomial", model_family, ignore.case = TRUE)) {
    warning("This function is designed for binomial models. Current model family: ",
      model_family,
      call. = FALSE
    )
  }

  # --- Data Preparation ---
  focus_var <- obj$focus
  response_var <- as.character(obj$model$formula)[2] # Get response variable name
  obj_data <- obj$data

  # Get observed response data
  observed <- obj_data[[response_var]]
  focus_levels <- obj_data[[focus_var]]

  # Validate response data for binomial context
  if (!all(observed >= 0 & observed <= 1)) {
    stop("Response variable contains values outside [0,1] range. Expected binary (0/1) or proportional data.",
      call. = FALSE
    )
  }

  # Calculate confidence interval multiplier
  alpha <- 1 - confidence_level
  z_value <- qnorm(1 - alpha / 2)

  # --- Panel 1: Raw Proportions ---

  # Calculate raw proportions by focus level
  if (all(observed %in% c(0, 1))) {
    # Binary data - calculate proportions with binomial confidence intervals
    raw_props <- aggregate(observed, list(level = focus_levels), function(x) {
      n <- length(x)
      successes <- sum(x)
      prop <- successes / n

      # Binomial confidence interval (Wilson score interval for better small sample performance)
      if (n > 0) {
        z2 <- z_value^2
        center <- (successes + z2 / 2) / (n + z2)
        margin <- z_value * sqrt((prop * (1 - prop) + z2 / (4 * n)) / (n + z2))
        ci_lower <- pmax(0, center - margin)
        ci_upper <- pmin(1, center + margin)
      } else {
        ci_lower <- ci_upper <- prop
      }

      c(prop = prop, ci_lower = ci_lower, ci_upper = ci_upper, n = n)
    })

    raw_df <- data.frame(
      level = raw_props$level,
      proportion = raw_props$x[, "prop"],
      ci_lower = raw_props$x[, "ci_lower"],
      ci_upper = raw_props$x[, "ci_upper"],
      n = raw_props$x[, "n"]
    )
  } else {
    # Already proportional data - use as is with approximate normal CI
    raw_props <- aggregate(observed, list(level = focus_levels), function(x) {
      n <- length(x)
      prop <- mean(x)
      se <- ifelse(n > 1, sd(x) / sqrt(n), 0)
      ci_lower <- pmax(0, prop - z_value * se)
      ci_upper <- pmin(1, prop + z_value * se)
      c(prop = prop, ci_lower = ci_lower, ci_upper = ci_upper, n = n)
    })

    raw_df <- data.frame(
      level = raw_props$level,
      proportion = raw_props$x[, "prop"],
      ci_lower = raw_props$x[, "ci_lower"],
      ci_upper = raw_props$x[, "ci_upper"],
      n = raw_props$x[, "n"]
    )
  }

  # Convert factor levels to numeric for plotting if needed
  if (is.factor(raw_df$level)) {
    raw_df$level_num <- as.numeric(raw_df$level)
    x_labels <- levels(raw_df$level)
    x_breaks <- seq_along(x_labels)
  } else {
    raw_df$level_num <- as.numeric(raw_df$level)
    x_labels <- sort(unique(raw_df$level))
    x_breaks <- x_labels
  }

  # Prepare data for combined plot
  raw_df$type <- "Observed"

  # --- Panel 2: Model Predictions ---

  # Make predictions directly for each focus level
  new_data <- data.frame(
    unique_level = sort(unique(obj_data[[focus_var]]))
  )
  names(new_data)[1] <- focus_var

  # Add other variables at their median/modal values
  other_vars <- setdiff(names(obj_data), c(focus_var, response_var))
  for (var in other_vars) {
    if (is.numeric(obj_data[[var]])) {
      new_data[[var]] <- median(obj_data[[var]], na.rm = TRUE)
    } else if (is.factor(obj_data[[var]]) || is.character(obj_data[[var]])) {
      # Use most common level
      mode_val <- names(sort(table(obj_data[[var]]), decreasing = TRUE))[1]
      if (is.factor(obj_data[[var]])) {
        new_data[[var]] <- factor(mode_val, levels = levels(obj_data[[var]]))
      } else {
        new_data[[var]] <- mode_val
      }
    }
  }

  # Make predictions
  preds <- predict(obj$model, newdata = new_data, se.fit = TRUE, type = "response")

  pred_df <- data.frame(
    level = new_data[[focus_var]],
    predicted = preds$fit,
    se = preds$se.fit
  )

  # Calculate confidence intervals
  pred_df$ci_lower <- pmax(0, pred_df$predicted - z_value * pred_df$se)
  pred_df$ci_upper <- pmin(1, pred_df$predicted + z_value * pred_df$se)

  # Convert factor levels for plotting
  if (is.factor(pred_df$level)) {
    pred_df$level_num <- as.numeric(pred_df$level)
  } else {
    pred_df$level_num <- as.numeric(pred_df$level)
  }

  # Prepare prediction data for combined plot
  pred_df$type <- "Predicted"

  # Rename columns to match for combining
  names(pred_df)[names(pred_df) == "predicted"] <- "proportion"

  # Combine both datasets
  combined_data <- rbind(
    raw_df[, c("level_num", "proportion", "ci_lower", "ci_upper", "type")],
    pred_df[, c("level_num", "proportion", "ci_lower", "ci_upper", "type")]
  )

  # Create combined plot
  combined_plot <- ggplot2::ggplot(combined_data, ggplot2::aes(
    x = .data$level_num, y = .data$proportion,
    color = .data$type, fill = .data$type
  )) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$ci_lower, ymax = .data$ci_upper),
      alpha = 0.2, color = NA
    ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::scale_color_manual(values = c("Observed" = "black", "Predicted" = "royalblue")) +
    ggplot2::scale_fill_manual(values = c("Observed" = "black", "Predicted" = "royalblue")) +
    ggplot2::labs(
      x = focus_var,
      y = "Proportion",
      color = "Data Type",
      fill = "Data Type"
    )

  return(combined_plot)
}

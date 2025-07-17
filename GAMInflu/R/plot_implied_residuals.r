#' @title Implied Residuals Plot
#' @description Creates a multi-panel plot showing the relationship between the focus term and another variable, including predicted values and observed means/quantiles.
#' @param obj A `gam_influence` object containing calculated indices.
#' @param var The name of the variable for which to compute mean and quantiles (y-axis).
#' @param nbins (Optional) Number of bins if `var` is not a factor.
#' @param islog Logical indicating if the response variable should be log-transformed.
#' @param ylim (Optional) The y-axis limits for the plot.
#' @param n.exclude (Optional) Minimum number of records required to include a group in the plot.
#' @return A ggplot object (multi-panel plot).
#' @export
plot_implied_residuals <- function(obj, var, nbins = 6, islog = NULL, ylim = NULL, n.exclude = 0) {
  # Check object and data
  if (is.null(obj$data)) stop("No data found in gam_influence object.")
  data <- obj$data
  focus_var <- obj$focus
  response_var <- obj$response

  # Determine log transformation
  if (is.null(islog)) islog <- obj$islog

  # Check variable existence
  if (!focus_var %in% names(data)) stop(paste("Focus variable '", focus_var, "' not found in the data."))
  if (!var %in% names(data)) stop(paste("Variable '", var, "' not found in the data."))

  # Bin variable if not a factor
  if (!is.factor(data[[var]])) {
    data[[var]] <- cut_number(data[[var]], n = nbins)
    var_factor <- as.factor(data[[var]])
  } else {
    var_factor <- data[[var]]
  }
  data$var_factor <- var_factor

  # Predict values for the focus variable using the stored model
  predictions <- data.frame(
    focus_var = data[[focus_var]],
    predicted = predict(obj$model, type = "response", newdata = data),
    group = data$var_factor,
    observed = data[[response_var]]
  )

  # Compute mean and quantiles for 'var'
  var_summary <- predictions %>%
    dplyr::group_by(group, focus_var) %>%
    dplyr::summarise(
      predicted = mean(predicted, na.rm = TRUE),
      p_q025 = quantile(predicted, 0.025, na.rm = TRUE),
      p_q975 = quantile(predicted, 0.975, na.rm = TRUE),
      x_mean = mean(observed, na.rm = TRUE),
      x_CI025 = x_mean - 1.96 * sd(observed, na.rm = TRUE) / sqrt(dplyr::n()),
      x_CI975 = x_mean + 1.96 * sd(observed, na.rm = TRUE) / sqrt(dplyr::n()),
      x_q025 = quantile(observed, 0.025, na.rm = TRUE),
      x_q975 = quantile(observed, 0.975, na.rm = TRUE),
      n_records = dplyr::n(),
      .groups = "drop"
    )
  if (!is.null(n.exclude)) {
    var_summary <- var_summary %>% dplyr::filter(n_records > n.exclude)
  }

  # Apply log transformation if needed
  if (islog) {
    var_summary$predicted <- exp(var_summary$predicted)
    var_summary$p_q025 <- exp(var_summary$p_q025)
    var_summary$p_q975 <- exp(var_summary$p_q975)
    var_summary$x_mean <- exp(var_summary$x_mean)
    var_summary$x_CI025 <- exp(var_summary$x_CI025)
    var_summary$x_CI975 <- exp(var_summary$x_CI975)
    var_summary$x_q025 <- exp(var_summary$x_q025)
    var_summary$x_q975 <- exp(var_summary$x_q975)
  }

  # Set y-axis limits
  if (is.null(ylim)) {
    ylim <- c(min(var_summary$x_CI025, na.rm = TRUE), max(var_summary$x_CI975, na.rm = TRUE) * 1.04)
    if (islog) ylim[1] <- 0
  } else {
    var_summary$x_CI025 <- pmax(var_summary$x_CI025, ylim[1])
    var_summary$x_CI975 <- pmin(var_summary$x_CI975, ylim[2])
  }

  # Normalise bar heights for record count
  var_summary$N <- var_summary$n_records / max(var_summary$n_records) * ylim[2] * 0.3

  # Convert focus_var to numeric if possible
  n_focus_var <- suppressWarnings(as.numeric(as.character(var_summary$focus_var)))
  if (!any(is.na(n_focus_var))) var_summary$focus_var <- n_focus_var

  # Create the multi-panel plot
  p <- ggplot(var_summary, aes(x = focus_var, y = x_mean, group = 1)) +
    geom_ribbon(aes(ymin = x_CI025, ymax = x_CI975), fill = "blue", alpha = 0.4) +
    geom_line(colour = "blue") +
    geom_line(aes(y = predicted, colour = "Predicted"), size = 1) +
    geom_col(aes(y = N), fill = "black", alpha = 0.3) +
    facet_wrap(~group, scales = "free_y") +
    coord_cartesian(ylim = ylim) +
    labs(x = focus_var, y = "Predicted Value", colour = "") +
    theme(legend.position = "bottom")
  # Return the plot
  return(p)
}

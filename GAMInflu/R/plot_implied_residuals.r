#' Plot Implied Residuals for GAM Models
#'
#' This function generates a multi-panel plot to visualize the relationship between
#' predicted values of the GAM model's focus variable and the mean and quantiles of
#' another variable in the data.  It is a generic function.
#'
#' @param gam_model A fitted GAM model object (e.g., from `mgcv::gam`).
#' @param data The data frame used to fit the model.
#' @param focus_var The name of the focus variable (x-axis).
#' @param var The name of the variable for which to compute mean and quantiles (y-axis).
#' @param nbins (Optional) If `x_var` is not a factor, the number of bins to create.
#' @param islog Logical indicating if the y variable should be log-transformed.
#' @param ylim (Optional) The y-axis limits for the plot.
#' @param n.exclude (Optional) The minimum number of records required to include a group in the plot.
#' @param ... Additional arguments (currently unused).
#'
#' @return Prints a ggplot object (multi-panel plot), and (optionally) returns a list containing the plot and the data used for plotting.
#'
#' @export
plot_implied_residuals <- function(gam_model, data, focus_var, var, nbins = 6, islog = FALSE, ylim = NULL, n.exclude = 0, ...) {
  UseMethod("plot_implied_residuals")
}

#' @rdname plot_implied_residuals
#' @export
plot_implied_residuals.gam <- function(gam_model, data, focus_var, var, nbins = 6, islog = FALSE, ylim = NULL, n.exclude = 0, ...) {
  # Check if the focus variable exists
  if (!focus_var %in% names(data)) {
    stop(paste("Focus variable '", focus_var, "' not found in the data."))
  }

  # Check if the y variable exists
  if (!var %in% names(data)) {
    stop(paste("Variable '", var, "' not found in the data."))
  }
  data <- data %>% ungroup()

  # Check if variable is a factor.  If not, bin it.
  if (!is.factor(data[[var]])) {
    if (is.null(nbins)) {
      nbins <- 6 # default
    }
    data <- data %>%
      mutate(!!var := cut_number(!!sym(var), n = nbins)) # Use cut_number
    var_factor <- as.factor(data[[var]])
  } else {
    var_factor <- data[[var]]
  }
  data$var_factor <- var_factor

  # 1. Predict values for the focus variable
  new_data <- data
  focus_levels <- unique(new_data[[focus_var]])

  predictions <- data.frame(
    focus_var = new_data[[focus_var]],
    predicted = predict(gam_model, type = "response", newdata = new_data),
    group = new_data$var_factor,
    observed = gam_model$y
  )

  # 2. Compute mean and quantiles for 'variable'
  var_summary <- predictions %>%
    group_by(group, focus_var) %>%
    summarise(
      predicted = mean(predicted, na.rm = TRUE),
      p_q025 = quantile(predicted, 0.025, na.rm = TRUE),
      p_q975 = quantile(predicted, 0.975, na.rm = TRUE),
      x_mean = mean(observed, na.rm = TRUE),
      x_CI025 = x_mean - 1.96 * sd(observed, na.rm = TRUE) / sqrt(n()),
      x_CI975 = x_mean + 1.96 * sd(observed, na.rm = TRUE) / sqrt(n()),
      x_q025 = quantile(observed, 0.025, na.rm = TRUE),
      x_q975 = quantile(observed, 0.975, na.rm = TRUE),
      n_records = n() # Store the number of records.
    )
  if (!is.null(n.exclude)) {
    var_summary <- var_summary %>% filter(n_records > n.exclude)
  } # Exclude groups with fewer than n.exclude records.

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
  if (is.null(ylim)) {
    if (islog) {
      ylim <- c(0, max(var_summary$x_CI975, na.rm = TRUE) * 1.04)
    } else {
      ylim <- c(min(var_summary$x_CI025, na.rm = TRUE), max(var_summary$x_CI975, na.rm = TRUE) * 1.04)
    }
  } else {
    var_summary$x_CI025 <- pmax(var_summary$x_CI025, ylim[1])
    var_summary$x_CI975 <- pmin(var_summary$x_CI975, ylim[2])
  }
  var_summary$N <- var_summary$n_records / max(var_summary$n_records) * ylim[2] * 0.3
  var_summary$n.focus_var <- as.numeric(as.character(var_summary$focus_var))
  if (!any(is.na(var_summary$n.focus_var))) {
    var_summary$focus_var <- var_summary$n.focus_var
  }

  # 3. Create the multi-panel plot
  p <- ggplot(var_summary, aes(x = focus_var, y = x_mean, group = 1)) +
    # geom_point(aes(size = n_records, group = 1), colour = "blue", alpha = 0.3) + # Size by n_records
    geom_ribbon(aes(x = focus_var, ymin = x_CI025, ymax = x_CI975, group = 1), fill = "blue", alpha = 0.4) +
    # geom_ribbon(aes(x = focus_var, ymin = p_q025, ymax = p_q975, group = 1), fill = "black", alpha = 0.4) +
    geom_line(aes(x = focus_var, y = x_mean, group = 1), colour = "blue") +
    geom_line(aes(x = focus_var, y = predicted, group = 1), colour = "black") +
    geom_col(aes(x = focus_var, y = N), fill = "black", alpha = 0.3) +
    # geom_rug(aes(x = focus_var, y = n_records), sides = "b", alpha = 0.6) +
    facet_wrap(~group, scales = "free_y") +
    ylim(ylim) +
    labs(x = focus_var, y = "Predicted Value") +
    theme(legend.position = "bottom")
  print(p)
  invisible(list("plot" = p, "data" = var_summary))
}
# plot_implied_residuals.gam(res1.gam, temp, "myear", "statarea", islog = T, ylim = c(0, 1.5))
# plot_implied_residuals.gam(res1.gam, temp, "myear", "moon.phase", islog = T, ylim = c(0, 1.5))

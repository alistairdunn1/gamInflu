#' @title Implied Residuals Plot
#' @description Creates a multi-panel implied residual plot showing the relationship between the focus term and a non-modelled variable. The implied coefficients are calculated as the normalised focus_term coefficient plus the mean of the standardised residuals for each group of the non-modelled variable.
#' @param obj A `gam_influence` object containing calculated indices and model data.
#' @param var The name of the variable for which to compute implied residual mean and quantiles.
#' @param nbins (Optional) Number of bins if `var` is not a factor.
#' @param islog Logical indicating if the response variable is log-transformed.
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

  # Calculate residuals and standardised residuals
  predictions$residual <- predictions$observed - predictions$predicted
  predictions$std_residual <- predictions$residual / sd(predictions$residual, na.rm = TRUE)
  predictions$focus_coef <- predictions$predicted

  # Calculate mean standardised residual per group
  df <- predictions %>%
    dplyr::group_by(focus_var, group) %>%
    dplyr::summarise(
      mean_std_resid = mean(std_residual, na.rm = TRUE),
      focus_coef = mean(focus_coef, na.rm = TRUE),
      n_records = n()
    )

  # Calculate implied coefficient per group
  df$implied_coef <- df$focus_coef + df$mean_std_resid
  df$q025 <- df$implied_coef - 1.96 * df$mean_std_resid
  df$q975 <- df$implied_coef + 1.96 * df$mean_std_resid

  # Apply log transformation if needed
  if (islog) {
    df$focus_coef <- exp(df$focus_coef)
    focus_mean <- mean(df$focus_coef, na.rm = TRUE)
    df$focus_coef <- df$focus_coef / focus_mean
    df$implied_coef <- exp(df$implied_coef) / focus_mean
    df$q025 <- exp(df$q025) / focus_mean
    df$q975 <- exp(df$q975) / focus_mean
  } else {
    focus_mean <- mean(df$focus_coef, na.rm = TRUE)
    df$focus_coef <- df$focus_coef / focus_mean
    df$implied_coef <- df$implied_coef / focus_mean
    df$q025 <- df$q025 / focus_mean
    df$q975 <- df$q975 / focus_mean
  }

  # Set y-axis limits
  if (is.null(ylim)) {
    ylim <- c(NA, NA)
    if (islog) ylim[1] <- 0
  }

  # Normalise bar heights for record count
  df$N <- df$n_records / max(df$n_records) * ylim[2] # * 0.3

  # Convert focus to numeric if possible
  if (is.factor(df$focus_var) && all(!is.na(as.numeric(as.character(levels(df$focus_var)))))) {
    df$focus_var <- as.numeric(as.character(df$focus_var))
  }

  # Create the multi-panel plot
  # Plot implied_coef by group
  p <- ggplot(df, aes(x = focus_var)) +
    geom_errorbar(aes(ymin = q025, ymax = q975, colour = "Implied"), alpha = 0.6) +
    geom_point(aes(y = implied_coef, colour = "Implied")) +
    geom_line(aes(y = focus_coef, colour = "Expected")) +
    geom_col(aes(y = N, fill = "Observations"), width = 1, alpha = 0.3) +
    scale_colour_manual(values = c("Implied" = "royalblue", "Expected" = "black")) +
    facet_wrap(~group, scales = "fixed") +
    ylab("Implied coefficient") +
    xlab(focus_var) +
    coord_cartesian(ylim = ylim)

  # Return the plot
  return(p)
}

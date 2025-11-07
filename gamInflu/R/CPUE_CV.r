#' Estimate Additional CV for Smooth CPUE Trend
#'
#' @description Estimates the additional coefficient of variation (CV) needed for a GAM-derived
#' standardized index series to follow a smooth trend. Uses a lognormal error model and
#' constrains the maximum rate of change between time points.
#'
#' @param obj A `gam_influence` object containing calculated indices from `calculate_influence()`.
#' @param r Numeric; maximum rate of population increase or change (as a proportion, default 0.1 or 10 percent).
#' @param ylim Numeric vector of length 2 specifying the y-axis limits. If NULL (default), calculated automatically.
#' @param increase.only Logical; if TRUE, only consider positive changes when selecting smoothing parameter.
#' @param add.mean Logical; if TRUE (default), add the smoothed trend line to the plot.
#' @param plot Logical; if TRUE (default), create a plot of the index with uncertainty. If FALSE, just return results without plotting.
#'
#' @return Invisibly returns a list containing:
#'   \itemize{
#'     \item `additional_CV`: The additional CV needed to achieve the smooth trend
#'     \item `total_CV`: Combined model CV and additional smoothing CV
#'     \item `data`: A data frame with the original indices, smoothed values,
#'           model CV, additional CV, and total CV
#'     \item `plot`: The ggplot object (if plot=TRUE)
#'   }
#'
#' @details
#' This function takes standardized indices from a GAM model and estimates the additional
#' process error (CV) needed to achieve a smoothed trend with constrained maximum rate of change.
#' This is useful for:
#'   \enumerate{
#'     \item Identifying excessive inter-annual variability in standardized indices
#'     \item Quantifying additional uncertainty needed for smoother stock assessment inputs
#'     \item Determining if a more constrained smoothing is biologically plausible
#'     \item Presenting smoothed indices with appropriate uncertainty bounds
#'   }
#'
#' The CV is calculated assuming a lognormal error distribution, which is appropriate for
#' most fisheries abundance indices.
#'
#' @examples
#' \dontrun{
#' # First create and fit a gam_influence object
#' gi <- gam_influence(model, focus = "year")
#' gi <- calculate_influence(gi)
#'
#' # Get process error CV with default settings
#' results <- CPUE_CV(gi)
#' cat("Additional CV needed:", results$additional_CV, "\n")
#'
#' # With custom settings for maximum rate of change
#' results <- CPUE_CV(gi, r = 0.2, increase.only = TRUE)
#'
#' # Without plotting
#' results <- CPUE_CV(gi, r = 0.15, plot = FALSE)
#' }
#' @export
CPUE_CV <- function(obj, r = 0.1, ylim = NULL, increase.only = FALSE, add.mean = TRUE, plot = TRUE) {
  # Helper function to calculate lognormal CV
  CV_lognormal <- function(yobs, yfit) {
    lognormal.neg.log.likl <- function(cv, obs, fit) {
      sd <- sqrt(log(1 + cv^2))
      neg.log.likl <- sum(log(sd) + 0.5 * (0.5 * sd + log(obs / fit) / sd)^2)
      return(neg.log.likl)
    }
    est.cv <- optimize(lognormal.neg.log.likl, c(0.01, 0.9), obs = yobs, fit = yfit)$minimum
    return(est.cv)
  }

  # Helper function to calculate CV and create results data frame
  getCV <- function(year, cpue, f) {
    lowess.fit <- lowess(year, cpue, f = f)
    CV <- CV_lognormal(cpue, lowess.fit$y)
    return(data.frame(year = year, cpue = cpue, lowess = lowess.fit$y, cv = rep(CV, length(year))))
  }

  # Validate input
  if (!inherits(obj, "gam_influence")) {
    stop("Object must be a gam_influence object created with calculate_influence()")
  }

  # Extract indices from the gam_influence object
  # Will error if calculate_influence() hasn't been run
  tryCatch(
    {
      indices_df <- extract_indices(obj)
    },
    error = function(e) {
      stop("Could not extract indices. Please run calculate_influence() first")
    }
  )

  # Get focus levels (convert to numeric years), index values and model CV
  year <- as.numeric(as.character(indices_df$level))
  index <- indices_df$index
  model_cv <- indices_df$cv

  # Remove any NA values
  valid <- !is.na(year) & !is.na(index) & !is.na(model_cv)
  year <- year[valid]
  index <- index[valid]
  model_cv <- model_cv[valid]

  # Standardize the index to mean 1
  std <- mean(index)
  index <- index / std

  # Set y-axis limits if not provided
  if (is.null(ylim)) {
    ylim <- c(0, max(ceiling(index), na.rm = TRUE))
  }

  # Test different smoothing parameters
  f <- seq(0.01, 1, length = 20)
  rhat <- rep(NA, length(f))
  CV <- rhat

  # For each smoothing parameter, calculate the max rate of change
  for (i in seq_along(f)) {
    x <- getCV(year, index, f = f[i])
    if (increase.only) {
      # Maximum increase - only consider positive changes
      changes <- diff(x$lowess) / x$lowess[-length(x$lowess)]
      rhat[i] <- max(c(0, changes))
    } else {
      # Maximum rate of change - consider both positive and negative changes
      changes <- diff(x$lowess) / x$lowess[-length(x$lowess)]
      rhat[i] <- max(abs(changes))
    }
    CV[i] <- x$cv[1]
  }

  # Select smoothing parameter that gives desired rate of change
  # Handle case where multiple smoothing parameters give the same rate of change
  # First, remove any NA values
  valid_idx <- !is.na(rhat)
  rhat_valid <- rhat[valid_idx]
  f_valid <- f[valid_idx]

  # Then handle duplicates
  unique_idx <- !duplicated(rhat_valid)
  rhat_unique <- rhat_valid[unique_idx]
  f_unique <- f_valid[unique_idx]

  # Make sure values are sorted for proper interpolation
  sort_idx <- order(rhat_unique)
  fhat <- approx(rhat_unique[sort_idx], f_unique[sort_idx], xout = r, yleft = 1, yright = 0)$y

  # Calculate final results with selected smoothing parameter
  x <- getCV(year, index, f = fhat)

  # Calculate additional process error CV
  # This is the extra CV needed beyond the model CV to fit the smooth curve
  model_var <- log(1 + model_cv^2) # Convert model CV to log variance
  total_var <- log(1 + x$cv[1]^2) # Total variance from smoothing

  # Additional variance = total variance - model variance
  # If model variance > total, no additional CV needed
  additional_var <- pmax(0, total_var - model_var)
  additional_cv <- sqrt(exp(additional_var) - 1)

  # Calculate total CV combining model CV and additional CV
  total_cv <- sqrt((model_cv^2) + (additional_cv^2))

  # Calculate lower and upper bounds using total CV (±2CV)
  lower_bound <- exp(log(x$lowess) - 2 * total_cv)
  upper_bound <- exp(log(x$lowess) + 2 * total_cv)

  # Create complete dataset with all results
  result_data <- data.frame(
    year = x$year,
    index = x$cpue, # Original standardized index
    lowess = x$lowess, # Smoothed index
    model_cv = model_cv, # Original model CV
    additional_cv = additional_cv, # Additional process error CV
    total_cv = total_cv, # Combined CV
    lower = lower_bound, # Lower bound using total CV
    upper = upper_bound # Upper bound using total CV
  )

  # Print additional CV value
  mean_additional_cv <- mean(additional_cv)
  cat("Additional CV needed =", round(mean_additional_cv, 3), "\n")
  cat("Average total CV =", round(mean(total_cv), 3), "\n")

  # Prepare return object
  result <- list(
    additional_CV = mean_additional_cv,
    total_CV = mean(total_cv),
    data = result_data
  )

  # Create plot if requested
  if (plot) {
    # Create ggplot
    p <- ggplot2::ggplot(result$data, ggplot2::aes(x = year, y = index)) +
      # Add confidence interval ribbon
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
        fill = "gray90", alpha = 0.5
      ) +
      # Add points and lines for the index
      ggplot2::geom_line() +
      ggplot2::geom_point(size = 2) +
      # Set theme and labels
      ggplot2::labs(
        x = "Year",
        y = "Index",
        subtitle = paste0(
          "Additional CV = ", round(result$additional_CV, 3),
          ", Total CV = ", round(result$total_CV, 3)
        )
      ) +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) +
      ggplot2::scale_y_continuous(limits = ylim)

    # Add trend line if requested
    if (add.mean) {
      p <- p + ggplot2::geom_line(ggplot2::aes(y = lowess),
        linetype = "dashed",
        color = "royalblue",
        linewidth = 1
      )
    }

    # Add plot to result
    result$plot <- p

    # Print the plot
    print(p)
  }

  # Return data invisibly
  invisible(result)
}

# -----------------------------------------------------------------------------
# Print r2 summary
# -----------------------------------------------------------------------------

#' @title r2
#' @description Calcuate r2 summary statistics
#'
#' @author Alistair Dunn
#' @param obj An object of class `influence_gam`.
#' @param r2_type the type of r2 statistic required.
#' @return The `influence_gam` object with calculation results added.
#' @export
r2 <- function(obj, r2_type = c("r2", "r2Negel", "r2Dev")) {
  UseMethod("r2")
}
#' @rdname r2
#' @export
r2.influence_gam <- function(obj, r2_type = c("r2", "r2Negel", "r2Dev")) {
  if (!inherits(obj, "influence_gam")) {
    stop("Input 'obj' must be of class 'influence_gam'.")
  }
  if (!obj$calculated) {
    stop("Calculations not performed. Run 'calculate_influence()' before summarizing R2 contributions.")
  }
  if (is.null(obj$summary)) {
    stop("Summary table (obj$summary) is missing. Cannot calculate R2 contributions.")
  }

  r2_type <- match.arg(r2_type)

  if (!r2_type %in% names(obj$summary)) {
    available_r2 <- names(obj$summary)[grepl("r2", names(obj$summary))]
    stop("Selected R2 type '", r2_type, "' not found in the summary table. Available: ", paste(available_r2, collapse = ", "))
  }

  summary_df <- obj$summary

  # Check if all values are NA (except possibly intercept)
  if (all(is.na(summary_df[[r2_type]][summary_df$term != "intercept"]))) {
    message("All ", r2_type, " values are NA in the summary table. ")

    # For r2Dev, try to calculate directly from the model
    if (r2_type == "r2Dev") {
      null_dev <- obj$model$null.deviance
      residual_dev <- obj$model$deviance

      if (!is.null(null_dev) && !is.null(residual_dev) &&
        !is.na(null_dev) && !is.na(residual_dev) &&
        null_dev > 0) {
        # Calculate overall R² from deviances
        full_r2dev <- (null_dev - residual_dev) / null_dev
        message("Recalculated full model r2Dev: ", round(full_r2dev, 4))
        message("Computing approximate term contributions...")

        # Option 1: Distribute R² based on degrees of freedom
        terms <- summary_df$term[summary_df$term != "intercept"]
        df <- summary_df$k[summary_df$term != "intercept"]

        # If df is also NA, use equal contribution
        if (all(is.na(df))) {
          n_terms <- length(terms)
          r2_contribution <- rep(full_r2dev / n_terms, n_terms)
        } else {
          # Replace NAs with 1 (conservative assumption)
          df[is.na(df)] <- 1
          # Distribute proportional to degrees of freedom
          r2_contribution <- full_r2dev * (df / sum(df, na.rm = TRUE))
        }

        # Create result data frame
        contributions <- data.frame(
          term = terms,
          r2_contribution = r2_contribution
        )

        # Notify about the approximation
        message("Note: These are approximate contributions based on model degrees of freedom.")
        message("The original summary table did not contain valid r2Dev values.")

        # Add percentage for context
        contributions$r2_percent <- 100 * contributions$r2_contribution / full_r2dev
        attr(contributions, "total_r2") <- full_r2dev

        return(contributions)
      }
    }

    # If we can't calculate for r2Dev or it's a different R² type
    message("Unable to calculate term contributions. Returning zero contributions.")
    contributions <- data.frame(
      term = summary_df$term[summary_df$term != "intercept"],
      r2_contribution = 0
    )
    contributions$r2_percent <- 0
    attr(contributions, "total_r2") <- 0

    return(contributions)
  }

  # Normal case - at least some R² values are not NA
  # Ensure proper ordering
  if (!"intercept" %in% summary_df$term) {
    intercept_row <- summary_df[1, ]
    intercept_row$term <- "intercept"
    intercept_row[[r2_type]] <- 0
    summary_df <- rbind(intercept_row, summary_df)
  }

  # Set intercept R² to zero
  summary_df[[r2_type]][summary_df$term == "intercept"] <- 0

  # Handle NA values more carefully
  r2_values <- summary_df[[r2_type]]
  r2_values_filled <- r2_values
  last_valid_r2 <- 0 # Start with 0 before intercept

  for (i in 1:length(r2_values_filled)) {
    if (is.na(r2_values_filled[i])) {
      r2_values_filled[i] <- last_valid_r2
    } else {
      last_valid_r2 <- r2_values_filled[i]
    }
  }

  # Calculate difference from previous step's filled R2
  r2_diff <- diff(c(0, r2_values_filled))

  # Create result data frame
  contributions <- data.frame(
    term = summary_df$term,
    r2_contribution = r2_diff
  )

  # Handle negative contributions (could happen with NA interpolation)
  if (any(r2_diff < 0, na.rm = TRUE)) {
    contributions$r2_contribution <- pmax(0, contributions$r2_contribution)
  }

  # Remove the intercept row
  contributions <- contributions[contributions$term != "intercept", ]
  rownames(contributions) <- NULL

  # Clean up term names (remove leading '+ ')
  contributions$term <- sub("^\\+ ", "", contributions$term)

  # Calculate percentage
  total_r2 <- sum(contributions$r2_contribution)
  contributions$r2_percent <- ifelse(total_r2 > 0,
    100 * contributions$r2_contribution / total_r2,
    0
  )

  attr(contributions, "total_r2") <- total_r2

  return(contributions)
}

# -----------------------------------------------------------------------------
# R-squared Contribution Summary Function
# -----------------------------------------------------------------------------

#' Summarize R-squared Contribution of Terms
#'
#' Calculates the approximate contribution of each term to the model's R-squared
#' based on the sequential fitting performed during `calculate_influence`.
#'
#' @param obj An object of class `influence_gam` with calculated results.
#' @param r2_type The type of R-squared to use for calculating contributions.
#'   Options are `"r2Dev"` (Deviance Explained, default), `"r2Negel"` (Negelkerke R2),
#'   or `"r2"` (Correlation-based R2).
#'
#' @return A data frame with columns `term` and `r2_contribution`.
#' @export
r2_contribution <- function(obj, r2_type = c("r2Dev", "r2Negel", "r2")) {
  UseMethod("r2_contribution")
}

#' @rdname r2_contribution
#' @export
r2_contribution.influence_gam <- function(obj, r2_type = c("r2Dev", "r2Negel", "r2")) {
  if (!inherits(obj, "influence_gam")) {
    stop("Input 'obj' must be of class 'influence_gam'.")
  }
  if (!obj$calculated) {
    stop("Calculations not performed. Run 'calculate_influence()' before summarizing R2 contributions.")
  }
  if (is.null(obj$summary)){
      stop("Summary table (obj$summary) is missing. Cannot calculate R2 contributions.")
  }

  r2_type <- match.arg(r2_type)

  if (!r2_type %in% names(obj$summary)) {
      available_r2 <- names(obj$summary)[grepl("r2", names(obj$summary))]
      stop("Selected R2 type '", r2_type, "' not found in the summary table. Available: ", paste(available_r2, collapse=", "))
  }

  summary_df <- obj$summary

  # Ensure the intercept row exists and has a valid R2 (should be 0)
  if (!"intercept" %in% summary_df$term) {
      stop("Intercept row missing from summary table.")
  }
  # Ensure intercept R2 is 0, handling potential NA if intercept model failed
  summary_df[[r2_type]][summary_df$term == "intercept"] <- ifelse(is.na(summary_df[[r2_type]][summary_df$term == "intercept"]), NA, 0)


  # Calculate differences
  # Order by the sequence they were added (implicit in summary_df order)
  r2_values <- summary_df[[r2_type]]

  # Handle potential NAs in R2 values before diff calculation
  # Assumption: NA R2 means the step failed or added no improvement over the last valid step.
  r2_values_filled <- r2_values
  last_valid_r2 <- 0 # Start with 0 before intercept
  for(i in 1:length(r2_values_filled)){
      if(is.na(r2_values_filled[i])){
          r2_values_filled[i] <- last_valid_r2 # Carry forward last valid R2 if current is NA
      } else {
          last_valid_r2 <- r2_values_filled[i]
      }
  }

  # Calculate difference from previous step's filled R2
  r2_diff <- diff(c(0, r2_values_filled)) # Prepend 0 for the intercept step

  # Create result data frame
  contributions <- data.frame(
    term = summary_df$term,
    r2_contribution = r2_diff
  )

  # Remove the intercept row as its contribution is implicitly the baseline
  contributions <- contributions[contributions$term != "intercept", ]
  rownames(contributions) <- NULL

  # Clean up term names (remove leading '+ ')
  contributions$term <- sub("^\\+ ", "", contributions$term)

  # Handle cases where contribution might be negative due to NA propagation or model quirks
  contributions$r2_contribution <- pmax(0, contributions$r2_contribution) # Ensure contributions >= 0

  return(contributions)
}
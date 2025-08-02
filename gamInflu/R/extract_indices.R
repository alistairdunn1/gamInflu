#' @title Extract Standardised and Unstandardised Indices
#'
#' @description Extracts the standardised and unstandardised indices for the focus term as a data frame.
#' Works with all supported GLM families (Gaussian, binomial, gamma, Poisson) and returns the underlying
#' data used by \code{plot_standardisation()} for further analysis, export, or custom visualization.
#'
#' @param obj A \code{gam_influence} object containing calculated indices from \code{calculate_influence()}.
#' @return A data frame with columns:
#'   \describe{
#'     \item{focus_term}{Name of the focus term (for reference)}
#'     \item{level}{The levels/values of the focus term}
#'     \item{unstandardised_index}{Raw aggregated values by focus level}
#'     \item{unstandardised_cv}{Coefficient of variation for unstandardised index}
#'     \item{index}{Model-adjusted values accounting for all terms}
#'     \item{se}{Standard error of standardised index (NA for log-transformed data)}
#'     \item{cv}{Coefficient of variation for standardised index}
#'     \item{lower_CI}{Lower confidence bound for standardised index}
#'     \item{upper_CI}{Upper confidence bound for standardised index}
#'   }
#' @details
#' This function provides programmatic access to the standardised and unstandardised indices
#' that are visualized in \code{plot_standardisation()}. The returned data frame contains the same
#' underlying data used for plotting but in a format suitable for further analysis, export,
#' or custom visualization.
#'
#' \strong{Family-specific behaviour:}
#' \itemize{
#'   \item \strong{Gaussian}: Traditional CPUE-style indices with geometric/arithmetic mean aggregation
#'   \item \strong{Binomial}: Probability-based indices for presence/absence data
#'   \item \strong{Gamma}: Positive continuous indices using geometric mean methods
#'   \item \strong{Poisson}: Count-based indices with appropriate statistical methods
#' }
#'
#' \strong{Log-transformed data handling:}
#' For log-transformed data (islog = TRUE), the coefficient of variation (CV) is calculated
#' in log space to provide meaningful measures of relative uncertainty. The standard
#' error (stan_se) is set to NA for log-transformed data as it is not interpretable on the
#' linear scale. Confidence intervals are calculated in log space and then exponentiated.
#'
#' The confidence intervals reflect the uncertainty in the model predictions using the
#' confidence level specified in \code{calculate_influence()} (default 95%).
#' @examples
#' \dontrun{
#' # Basic usage
#' gi <- gam_influence(your_model, focus = "year")
#' gi <- calculate_influence(gi)
#' indices_df <- extract_indices(gi)
#'
#' # View the results including CV information
#' print(indices_df)
#'
#' # Check precision: lower CV indicates more precise estimates
#' summary(indices_df$cv)
#'
#' # Identify years with high uncertainty (high CV)
#' high_cv_years <- indices_df[indices_df$cv > 0.2, ]
#'
#' # Use with different families
#' # Binomial model
#' gi_binom <- calculate_influence(gi_binomial_model)
#' binom_indices <- extract_indices(gi_binom)
#'
#' # Gamma model
#' gi_gamma <- calculate_influence(gi_gamma_model)
#' gamma_indices <- extract_indices(gi_gamma)
#'
#' # Export to CSV including CV information
#' write.csv(indices_df, "focus_indices.csv", row.names = FALSE)
#' }
#' @seealso \code{plot_standardisation}, \code{calculate_influence}, \code{gam_influence}
#' @export
extract_indices <- function(obj) {
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

  # Create clean output data frame with descriptive column names
  result <- data.frame(
    level = df$level,
    unstandardised_index = df$unstan,
    unstandardised_cv = df$unstan_cv,
    index = df$standardised_index,
    se = df$stan_se,
    cv = df$standardised_cv,
    lower_CI = df$stan_lower,
    upper_CI = df$stan_upper,
    focus_term = obj$focus,
    stringsAsFactors = FALSE
  )
  # Rearrange columns in logical order: identification, unstandardised metrics, standardised metrics, confidence bounds
  result <- result[, c(
    "focus_term", "level", "unstandardised_index", "unstandardised_cv",
    "index", "se", "cv", "lower_CI", "upper_CI"
  )]

  # Convert level to numeric if possible (same logic as plot_standardisation)
  if (is.factor(result$level) && all(!is.na(as.numeric(as.character(levels(result$level)))))) {
    result$level <- as.numeric(as.character(result$level))
  }

  # Add class for potential future methods
  class(result) <- c("gam_influence_indices", "data.frame")

  return(result)
}

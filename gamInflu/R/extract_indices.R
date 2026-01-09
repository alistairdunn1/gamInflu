#' @title Extract Standardised and Unstandardised Indices
#'
#' @description Extracts the standardised and unstandardised indices for the focus term as a data frame.
#' This is a generic function that works with both \code{gam_influence} objects (from \code{calculate_influence()})
#' and \code{gam_influence_combined} objects (from \code{combine_indices()}).
#'
#' Works with all supported GLM families (Gaussian, binomial, gamma, Poisson) and returns the underlying
#' data used by \code{plot_standardisation()} for further analysis, export, or custom visualisation.
#'
#' @param obj A \code{gam_influence} or \code{gam_influence_combined} object containing calculated indices.
#' @return A data frame with columns depending on the object type:
#'
#' For \code{gam_influence} objects:
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
#'
#' For \code{gam_influence_combined} objects:
#'   \describe{
#'     \item{focus_term}{Name of the focus term}
#'     \item{level}{The levels/values of the focus term}
#'     \item{binomial_index}{Probability component from binomial model}
#'     \item{binomial_cv}{Coefficient of variation for binomial component}
#'     \item{positive_index}{Conditional catch component from positive model}
#'     \item{positive_cv}{Coefficient of variation for positive component}
#'     \item{combined_index}{Combined delta-GLM index}
#'     \item{combined_se}{Standard error of combined index}
#'     \item{combined_cv}{Coefficient of variation for combined index}
#'     \item{combined_lower_CI}{Lower confidence bound for combined index}
#'     \item{combined_upper_CI}{Upper confidence bound for combined index}
#'     \item{combination_method}{Method used to combine indices}
#'   }
#' @details
#' This function provides programmatic access to the standardised and unstandardised indices
#' that are visualized in \code{plot_standardisation()}. The returned data frame contains the same
#' underlying data used for plotting but in a format suitable for further analysis, export,
#' or custom visualisation.
#'
#' \strong{Family-specific behaviour (gam_influence objects):}
#' \itemize{
#'   \item \strong{Gaussian}: Traditional CPUE-style indices with geometric/arithmetic mean aggregation
#'   \item \strong{Binomial}: Probability-based indices for presence/absence data
#'   \item \strong{Gamma}: Positive continuous indices using geometric mean methods
#'   \item \strong{Poisson}: Count-based indices with appropriate statistical methods
#' }
#'
#' \strong{Combined indices (gam_influence_combined objects):}
#' For objects created by \code{combine_indices()}, this function extracts both the component
#' indices (binomial probability and positive conditional catch) and the combined delta-GLM index.
#' This is particularly useful for:
#' \itemize{
#'   \item Exporting complete delta-GLM results
#'   \item Comparing component contributions
#'   \item Custom visualizations of combined indices
#'   \item Further statistical analysis
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
#' # Basic usage with gam_influence object
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
#' # Use with combined indices (delta-GLM)
#' combined_gi <- combine_indices(gi_binom, gi_pos)
#' combined_df <- extract_indices(combined_gi)
#'
#' # View all components and combined index
#' print(combined_df)
#'
#' # Export to CSV including CV information
#' write.csv(indices_df, "focus_indices.csv", row.names = FALSE)
#' write.csv(combined_df, "combined_indices.csv", row.names = FALSE)
#' }
#' @seealso \code{\link{plot_standardisation}}, \code{\link{calculate_influence}},
#'   \code{\link{gam_influence}}, \code{\link{combine_indices}}
#' @export
extract_indices <- function(obj) {
  UseMethod("extract_indices")
}

#' @title Extract Indices from gam_influence Object
#' @description S3 method for extracting indices from gam_influence objects
#' @param obj A \code{gam_influence} object containing calculated indices from \code{calculate_influence()}.
#' @return A data frame with standardised and unstandardised indices
#' @export
extract_indices.gam_influence <- function(obj) {
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

#' @title Extract Indices from gam_influence_combined Object
#' @description S3 method for extracting component and combined indices from delta-GLM objects
#' @param obj A \code{gam_influence_combined} object from \code{combine_indices()}.
#' @return A data frame with binomial, positive, and combined indices including confidence intervals
#' @export
extract_indices.gam_influence_combined <- function(obj) {
  # Validate input
  if (is.null(obj$combined_indices)) {
    stop("No combined indices found. The gam_influence_combined object appears to be incomplete.",
      call. = FALSE
    )
  }

  df <- obj$combined_indices

  # Check for required columns
  required_cols <- c(
    "level", "standardised_index_binom", "standardised_index_pos",
    "combined_index", "combined_lower", "combined_upper"
  )
  missing_cols <- setdiff(required_cols, names(df))

  if (length(missing_cols) > 0) {
    stop("Combined indices data frame is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (is.null(obj$focus_term)) {
    stop("Focus term is not set in the gam_influence_combined object.", call. = FALSE)
  }

  # Create clean output data frame with descriptive column names
  result <- data.frame(
    level = df$level,
    binomial_index = df$standardised_index_binom,
    binomial_cv = df$standardised_cv_binom,
    positive_index = df$standardised_index_pos,
    positive_cv = df$standardised_cv_pos,
    combined_index = df$combined_index,
    combined_se = df$combined_se,
    combined_cv = df$combined_cv,
    combined_lower_CI = df$combined_lower,
    combined_upper_CI = df$combined_upper,
    focus_term = obj$focus_term,
    combination_method = obj$method,
    stringsAsFactors = FALSE
  )

  # Rearrange columns in logical order: identification, components, combined, metadata
  result <- result[, c(
    "focus_term", "level",
    "binomial_index", "binomial_cv",
    "positive_index", "positive_cv",
    "combined_index", "combined_se", "combined_cv",
    "combined_lower_CI", "combined_upper_CI",
    "combination_method"
  )]

  # Convert level to numeric if possible (consistent with other extract_indices methods)
  if (is.factor(result$level) && all(!is.na(as.numeric(as.character(levels(result$level)))))) {
    result$level <- as.numeric(as.character(result$level))
  } else if (is.character(result$level) && all(!is.na(suppressWarnings(as.numeric(result$level))))) {
    result$level <- as.numeric(result$level)
  }

  # Add class for potential future methods
  class(result) <- c("gam_influence_combined_indices", "data.frame")

  return(result)
}

#' @title Default Method for extract_indices
#' @description Default S3 method that provides an informative error message
#' @param obj An object of unsupported class
#' @return Throws an error
#' @export
extract_indices.default <- function(obj) {
  stop("extract_indices() is not implemented for objects of class '",
    class(obj)[1], "'.\n",
    "Supported classes: 'gam_influence', 'gam_influence_combined'",
    call. = FALSE
  )
}

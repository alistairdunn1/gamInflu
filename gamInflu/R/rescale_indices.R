#' @title Rescale Influence Indices Across Groups
#' @description Rescales influence indices across different groups (e.g., different model types or stocks) using common years as reference points. This function is useful when comparing influence indices from different models or datasets that may have different baseline levels.
#' @param obj A \code{gam_influence} object that has been processed with \code{calculate_influence}, or a data frame containing influence indices.
#' @param group_col Character string specifying the column name that contains the grouping variable (e.g., "stock", "model_type"). Default is "level" for single-group analysis.
#' @param year_col Character string specifying the column name that contains the year/time variable. Default is "level" if the focus variable represents time.
#' @param index_col Character string specifying the column name that contains the index values to rescale. Default is "standardised_index".
#' @param lowerCI_col Character string specifying the column name that contains the lower confidence interval values. Default is "stan_lower".
#' @param upperCI_col Character string specifying the column name that contains the upper confidence interval values. Default is "stan_upper".
#' @param common_years_only Logical. If TRUE (default), only uses years present in all groups for calculating reference means. If FALSE, uses all available years for each group.
#' @return A data frame with rescaled indices. If input was a \code{gam_influence} object, returns the object with rescaled indices in the \code{calculated$indices} data frame. If input was a data frame, returns the rescaled data frame.
#' @details
#' This function rescales influence indices by:
#' 1. Identifying common years across all groups (if \code{common_years_only = TRUE})
#' 2. Calculating the mean index value for each group over the reference period
#' 3. Dividing all index values by their respective group means
#'
#' This approach allows for comparison of relative changes across different groups while accounting for different baseline levels.
#'
#' @examples
#' \dontrun{
#' # Assuming you have multiple models or datasets
#' data_list <- list(
#'   data.frame(year = 2010:2020, index = rnorm(11, 1, 0.1), stock = "StockA"),
#'   data.frame(year = 2010:2020, index = rnorm(11, 1.5, 0.1), stock = "StockB")
#' )
#' combined_data <- do.call(rbind, data_list)
#'
#' # Rescale indices
#' rescaled <- rescale_indices(combined_data,
#'   group_col = "stock",
#'   year_col = "year",
#'   index_col = "index"
#' )
#' }
#' @export
rescale_indices <- function(obj, group_col = NULL, year_col = NULL,
                            index_col = "standardised_index",
                            lowerCI_col = "stan_lower",
                            upperCI_col = "stan_upper",
                            common_years_only = TRUE) {
  # Handle gam_influence object input
  if (inherits(obj, "gam_influence")) {
    if (is.null(obj$calculated) || is.null(obj$calculated$indices)) {
      stop("gam_influence object must be processed with calculate_influence first", call. = FALSE)
    }
    data <- obj$calculated$indices
    is_gam_object <- TRUE

    # Set default column names for gam_influence objects
    if (is.null(group_col)) group_col <- "level"
    if (is.null(year_col)) year_col <- obj$focus
  } else if (is.data.frame(obj)) {
    data <- obj
    is_gam_object <- FALSE

    # Set defaults for data frame input
    if (is.null(group_col)) {
      stop("group_col must be specified when input is a data frame", call. = FALSE)
    }
    if (is.null(year_col)) {
      stop("year_col must be specified when input is a data frame", call. = FALSE)
    }
  } else {
    stop("Input must be a gam_influence object or data frame", call. = FALSE)
  }

  # Check required columns exist
  required_cols <- c(group_col, year_col, index_col)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Required columns not found: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  # Use dplyr for data manipulation
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr package is required for rescale_indices", call. = FALSE)
  }

  # Create symbols for tidy evaluation
  group_sym <- rlang::sym(group_col)
  year_sym <- rlang::sym(year_col)
  index_sym <- rlang::sym(index_col)

  # Step 1: Get total number of unique groups and validate data structure
  group_counts <- data %>%
    dplyr::count(!!group_sym) %>%
    dplyr::pull(n)

  if (length(unique(group_counts)) != 1) {
    stop("All groups must have the same number of observations for rescaling", call. = FALSE)
  }

  n_groups <- length(unique(data[[group_col]]))
  # n_per_group <- unique(group_counts)  # Not used, but kept for potential future use

  # Step 2: Find reference years
  if (common_years_only) {
    # Find years present in all groups
    years_per_group <- data %>%
      dplyr::group_by(!!group_sym, !!year_sym) %>%
      dplyr::summarise(count = dplyr::n(), .groups = "drop")

    reference_years_df <- years_per_group %>%
      dplyr::group_by(!!year_sym) %>%
      dplyr::filter(dplyr::n() == n_groups) %>%
      dplyr::select(!!group_sym, !!year_sym)
  } else {
    # Use all available years
    reference_years_df <- data %>%
      dplyr::select(!!group_sym, !!year_sym) %>%
      dplyr::distinct()
  }

  # Step 3: Calculate mean index over reference period for each group
  mean_reference_df <- data %>%
    dplyr::inner_join(reference_years_df, by = c(group_col, year_col)) %>%
    dplyr::group_by(!!group_sym) %>%
    dplyr::summarise(mean_reference = mean(!!index_sym, na.rm = TRUE), .groups = "drop")

  # Step 4: Join back and rescale all data
  result <- data %>%
    dplyr::left_join(mean_reference_df, by = group_col) %>%
    dplyr::mutate(rescaled_index = !!index_sym / mean_reference) # nolint: object_usage_linter

  # Rescale CI columns if they exist
  if (!is.null(lowerCI_col) && lowerCI_col %in% names(data)) {
    lowerCI_sym <- rlang::sym(lowerCI_col)
    result <- result %>%
      dplyr::mutate(rescaled_lowerCI = !!lowerCI_sym / mean_reference) # nolint: object_usage_linter
  }

  if (!is.null(upperCI_col) && upperCI_col %in% names(data)) {
    upperCI_sym <- rlang::sym(upperCI_col)
    result <- result %>%
      dplyr::mutate(rescaled_upperCI = !!upperCI_sym / mean_reference) # nolint: object_usage_linter
  }

  # Remove temporary column
  result <- result %>%
    dplyr::select(-mean_reference) # nolint: object_usage_linter

  # Return appropriate format
  if (is_gam_object) {
    # Update the gam_influence object
    obj$calculated$indices <- result
    return(obj)
  } else {
    return(result)
  }
}

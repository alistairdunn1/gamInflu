#' @title Analyse Focus Effects by Grouping Variable
#' @description Performs influence analysis on subsets of data grouped by a categorical variable.
#' This allows analysis of how focus term effects (e.g., year trends) vary across different 
#' levels of another variable (e.g., areas, gear types, species).
#' @param model_formula A formula object or character string specifying the model structure.
#' @param focus A character string specifying the focus term (must be a factor).
#' @param grouping_var A character string specifying the grouping variable (must be a factor).
#' @param data A data frame containing all variables referenced in the model formula.
#' @param groups A character vector specifying which groups to analyse. If NULL, analyses all levels.
#' @param family A family specification for the GAM (default: gaussian()).
#' @param islog Logical indicating if response is log-transformed (passed to calculate_influence).
#' @param ... Additional arguments passed to calculate_influence().
#' @return An object of class 'gam_influence_comparison' containing influence analysis results
#'   for each group, with methods for plotting and summarizing results.
#' @details
#' This function addresses the common need to analyse temporal or categorical trends
#' (focus effects) separately within different groups. For example:
#' - Annual CPUE trends by fishing area
#' - Species abundance trends by habitat type  
#' - Treatment effects by experimental site
#'
#' The function:
#' 1. Subsets the data by each level of the grouping variable
#' 2. Fits separate GAM models to each subset
#' 3. Performs influence analysis on each model
#' 4. Returns results in a structured format for comparison
#'
#' **Family Support**: Supports all GLM families (Gaussian, binomial, gamma, Poisson)
#' with automatic family detection in each subset model.
#' @examples
#' \dontrun{
#' library(mgcv)
#' library(gamInflu)
#'
#' # Analyse year effects by area
#' comparison <- analyse_focus_by_group(
#'   model_formula = cpue ~ s(depth) + s(temperature) + year,
#'   focus = "year",
#'   grouping_var = "area", 
#'   data = fisheries_data
#' )
#'
#' # View results
#' summary(comparison)
#' plot(comparison, type = "standardisation")
#' plot(comparison, type = "comparison")
#' }
#' @importFrom mgcv gam
#' @importFrom stats as.formula
#' @export
analyse_focus_by_group <- function(model_formula, focus, grouping_var, data, 
                                   groups = NULL, family = gaussian(), 
                                   islog = NULL, ...) {
  
  # Input validation
  if (!is.data.frame(data)) {
    stop("data must be a data.frame", call. = FALSE)
  }
  
  if (!focus %in% names(data)) {
    stop(paste("Focus variable '", focus, "' not found in data"), call. = FALSE)
  }
  
  if (!grouping_var %in% names(data)) {
    stop(paste("Grouping variable '", grouping_var, "' not found in data"), call. = FALSE)
  }
  
  if (!is.factor(data[[focus]])) {
    stop(paste("Focus variable '", focus, "' must be a factor"), call. = FALSE)
  }
  
  if (!is.factor(data[[grouping_var]])) {
    data[[grouping_var]] <- as.factor(data[[grouping_var]])
    message("Converting grouping variable '", grouping_var, "' to factor")
  }
  
  # Convert formula if needed
  if (is.character(model_formula)) {
    model_formula <- as.formula(model_formula)
  }
  
  # Ensure family is properly specified
  if (is.null(family)) {
    family <- gaussian()
  }
  
  # Determine groups to analyse
  if (is.null(groups)) {
    groups <- levels(data[[grouping_var]])
  } else {
    # Validate requested groups exist
    available_groups <- levels(data[[grouping_var]])
    missing_groups <- setdiff(groups, available_groups)
    if (length(missing_groups) > 0) {
      stop("Groups not found in data: ", paste(missing_groups, collapse = ", "), 
           "\nAvailable groups: ", paste(available_groups, collapse = ", "), call. = FALSE)
    }
  }
  
  message("Analyzing focus effects for ", length(groups), " groups: ", 
          paste(groups, collapse = ", "))
  
  # Analyse each group
  results <- list()
  models <- list()
  
  for (group in groups) {
    message("Processing group: ", group)
    
    # Subset data
    subset_data <- data[data[[grouping_var]] == group, ]
    
    if (nrow(subset_data) == 0) {
      warning("No data for group: ", group, ". Skipping.")
      next
    }
    
    # Check if focus term has multiple levels in subset
    focus_levels_in_subset <- length(levels(droplevels(subset_data[[focus]])))
    if (focus_levels_in_subset < 2) {
      warning("Group '", group, "' has insufficient focus term levels (", 
              focus_levels_in_subset, "). Skipping.")
      next
    }
    
    # Fit model to subset - using explicit gam call
    success <- TRUE
    tryCatch({
      subset_model <- mgcv::gam(formula = model_formula, 
                                data = subset_data, 
                                family = family)
      
      # Ensure the model stores the data for later use
      subset_model$data <- subset_data
      models[[group]] <- subset_model
      
    }, error = function(e) {
      warning("Failed to fit model for group '", group, "': ", e$message)
      success <<- FALSE
    })
    
    if (!success) next
    
    # Perform influence analysis
    tryCatch({
      gi <- gam_influence(subset_model, focus = focus, data = subset_data, islog = islog)
      gi <- calculate_influence(gi, islog = islog, ...)
      
      # Store results with group information
      gi$group <- group
      gi$grouping_var <- grouping_var
      results[[group]] <- gi
      
    }, error = function(e) {
      warning("Failed to perform influence analysis for group '", group, "': ", e$message)
    })
  }
  
  if (length(results) == 0) {
    stop("No successful model fits. Check data and model specification.", call. = FALSE)
  }
  
  # Create comparison object
  comparison <- structure(
    list(
      results = results,
      models = models,
      focus = focus,
      grouping_var = grouping_var,
      groups = names(results),
      formula = model_formula,
      family = family
    ),
    class = c("gam_influence_comparison", "list")
  )
  
  message("Successfully analysed ", length(results), " groups")
  return(comparison)
}

#' @title Compare Focus Effects Across Groups
#' @description A convenience wrapper for analyse_focus_by_group that provides a simpler interface
#' for common comparison scenarios.
#' @param model A fitted GAM model object.
#' @param focus Character string specifying the focus term.
#' @param grouping_var Character string specifying the grouping variable.
#' @param groups Optional character vector of specific groups to analyse.
#' @param ... Additional arguments passed to analyse_focus_by_group.
#' @return A gam_influence_comparison object.
#' @details
#' This function extracts the formula and data from an existing fitted model
#' and performs grouped influence analysis. Useful when you have already fitted
#' a model and want to explore group-specific effects.
#' @examples
#' \dontrun{
#' # Fit full model first
#' mod <- gam(cpue ~ s(depth) + s(temperature) + year * area, data = fisheries_data)
#' 
#' # Compare year effects by area
#' comparison <- compare_focus_by_groups(mod, focus = "year", grouping_var = "area")
#' }
#' @export
compare_focus_by_groups <- function(model, focus, grouping_var, groups = NULL, ...) {
  
  # Extract formula and data from model
  model_formula <- formula(model)
  data <- model$data
  
  if (is.null(data)) {
    stop("Model does not contain data. Please use analyse_focus_by_group with explicit data.", 
         call. = FALSE)
  }
  
  # Extract family - handle different family object structures
  if (is.null(model$family)) {
    family <- gaussian()  # Default fallback
  } else {
    family <- model$family
  }
  
  analyse_focus_by_group(
    model_formula = model_formula,
    focus = focus,
    grouping_var = grouping_var,
    data = data,
    groups = groups,
    family = family,
    ...
  )
}

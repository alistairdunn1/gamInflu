# Updated create_prediction_grid to handle splines
create_prediction_grid <- function(obj, term) {
  # Initialize an empty list to store the grid variables
  grid_vars <- list()
  
  # Always include the focus variable - ensuring factor handling
  focus_var <- obj$data[[obj$focus]]
  if (is.factor(focus_var)) {
    grid_vars[[obj$focus]] <- levels(focus_var)
  } else if (is.numeric(focus_var)) {
    grid_vars[[obj$focus]] <- seq(min(focus_var), max(focus_var), length.out = 20)
  } else {
    grid_vars[[obj$focus]] <- unique(focus_var)
  }
  
  # Check if term is a spline smoother from splines package
  is_spline_smoother <- term %in% names(obj$smoothers) && 
                        obj$smoothers[[term]]$type %in% c("ns", "bs", "poly", "lo", "pspline", "rcs")
  
  if (is_spline_smoother) {
    # For splines package smoothers, include the spline variable in the grid
    var_name <- obj$smoothers[[term]]$vars
    var_data <- obj$data[[var_name]]
    
    if (!(var_name %in% names(grid_vars))) {
      if (is.factor(var_data)) {
        grid_vars[[var_name]] <- levels(var_data)
      } else if (is.numeric(var_data)) {
        # For numeric variables used in splines, we need more points for a smooth curve
        grid_vars[[var_name]] <- seq(min(var_data), max(var_data), length.out = 50)
      } else {
        grid_vars[[var_name]] <- unique(var_data)
      }
    }
  }
  
  # If the term is an interaction, include all variables involved
  if (term %in% names(obj$interactions)) {
    interaction_vars <- obj$interactions[[term]]$vars
    for (var in interaction_vars) {
      if (var != obj$focus && !(var %in% names(grid_vars))) {
        var_data <- obj$data[[var]]
        if (is.factor(var_data)) {
          grid_vars[[var]] <- levels(var_data)
        } else if (is.numeric(var_data)) {
          grid_vars[[var]] <- seq(min(var_data), max(var_data), length.out = 10)
        } else {
          grid_vars[[var]] <- unique(var_data)
        }
      }
    }
  }
  
  # For all other terms involved in interactions with the focus or year
  if (term == obj$focus || term == "year") {
    for (interaction_name in names(obj$interactions)) {
      interaction_info <- obj$interactions[[interaction_name]]
      if (obj$focus %in% interaction_info$vars || "year" %in% interaction_info$vars) {
        for (var in interaction_info$vars) {
          if (var != obj$focus && var != "year" && !(var %in% names(grid_vars))) {
            var_data <- obj$data[[var]]
            if (is.factor(var_data)) {
              grid_vars[[var]] <- levels(var_data)
            } else if (is.numeric(var_data)) {
              grid_vars[[var]] <- seq(min(var_data), max(var_data), length.out = 10)
            } else {
              grid_vars[[var]] <- unique(var_data)
            }
          }
        }
      }
    }
  }
  
  # Handle year variable separately if it's not the focus
  if ("year" %in% names(obj$data) && !"year" %in% names(grid_vars) && "year" != obj$focus) {
    year_data <- obj$data[["year"]]
    if (is.factor(year_data)) {
      grid_vars[["year"]] <- levels(year_data)
    } else if (is.numeric(year_data)) {
      grid_vars[["year"]] <- seq(min(year_data), max(year_data), length.out = 10)
    } else {
      grid_vars[["year"]] <- unique(year_data)
    }
  }
  
  # Add all other variables needed for prediction
  all_vars <- names(obj$data)
  for (var in all_vars) {
    if (!(var %in% names(grid_vars)) && var != obj$response) {
      var_data <- obj$data[[var]]
      # Use median/mode values for variables not directly involved
      if (is.factor(var_data)) {
        # Use the most common level (mode)
        grid_vars[[var]] <- names(sort(table(var_data), decreasing = TRUE)[1])
      } else if (is.numeric(var_data)) {
        # Use median for numeric variables
        grid_vars[[var]] <- median(var_data)
      } else {
        # Use the first value for other types
        grid_vars[[var]] <- var_data[1]
      }
    }
  }
  
  # Create the grid using expand.grid
  pred_grid <- do.call(expand.grid, grid_vars)
  
  # Ensure all factors are preserved as factors with proper levels
  for (var in names(pred_grid)) {
    if (is.factor(obj$data[[var]])) {
      pred_grid[[var]] <- factor(pred_grid[[var]], levels = levels(obj$data[[var]]))
    }
  }
  
  return(pred_grid)
}


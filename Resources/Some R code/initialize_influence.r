# Updated initialize_influence function to handle splines package smoothers
initialize_influence <- function(obj) {
  # If no data was supplied
  if (is.null(obj$data)) {
    # See if the model has a data variable
    obj$data <- obj$model$data
    if (is.null(obj$data)) {
      # See if it is available as a model attribute
      obj$data <- attr(obj$model, "data")
    }
    if (is.null(obj$data)) {
      # For GAMs, use model$model
      if (inherits(obj$model, "gam")) {
        obj$data <- obj$model$model
      } else {
        stop("You need to provide the data that was fitted to for this type of model e.g. influence(model, data=mydata, ...)")
      }
    }
  }
  
  # Get term labels based on model type
  if (inherits(obj$model, "gam")) {
    # For GAMs, handle smooth terms differently
    parameterization <- obj$model$pterms
    obj$terms <- attr(parameterization, "term.labels")
    
    # Add smooth terms from the model
    smooth_terms <- names(obj$model$smooth)
    if (length(smooth_terms) > 0) {
      # Store interaction information
      obj$interactions <- list()
      # Store smoother information (including non-mgcv smoothers)
      obj$smoothers <- list()
      
      for (i in 1:length(obj$model$smooth)) {
        smooth_label <- obj$model$smooth[[i]]$label
        obj$terms <- c(obj$terms, smooth_label)
        
        # Check for interactions with year or focus
        if (inherits(obj$model$smooth[[i]], "tensor.smooth") || 
            inherits(obj$model$smooth[[i]], "tp") || 
            inherits(obj$model$smooth[[i]], "ti")) {
          
          # Get variables involved in this smooth
          vars <- obj$model$smooth[[i]]$term
          
          # Check if this is an interaction involving the focus term or year
          if (length(vars) > 1) {
            # Store interaction information
            obj$interactions[[smooth_label]] <- list(
              type = class(obj$model$smooth[[i]])[1],
              vars = vars,
              by = obj$model$smooth[[i]]$by,
              label = smooth_label
            )
          }
        }
        
        # Store smoother information for all smooth terms
        obj$smoothers[[smooth_label]] <- list(
          type = class(obj$model$smooth[[i]])[1],
          vars = obj$model$smooth[[i]]$term,
          label = smooth_label
        )
      }
    }
  } else {
    # For GLMs and other models
    obj$terms <- attr(obj$model$terms, "term.labels")
    
    # Initialize interaction and smoother tracking
    obj$interactions <- list()
    obj$smoothers <- list()
    
    # Check for non-GAM smoothers (ns, bs, etc.) and interactions in coefficient names
    coef_names <- names(coef(obj$model))
    
    # Function to extract smoother type and variable from coefficient names
    extract_smoother_info <- function(coef_name) {
      # Match smoothers like ns(x, df=4), bs(x, df=4), etc.
      # Both with and without numeric suffix
      smoother_match <- regexpr("^(ns|bs|poly|lo|pspline|rcs)\\(([^,]+)[^)]*\\)($|[0-9]+$)", coef_name)
      
      if (smoother_match > 0) {
        match_text <- regmatches(coef_name, smoother_match)
        
        # Extract smoother type (ns, bs, etc.)
        smoother_type <- gsub("\\(.*$", "", match_text)
        
        # Extract variable name
        var_match <- regexpr("\\(([^,]+)", match_text)
        var_name <- substr(regmatches(match_text, var_match), 2, 1000)
        
        return(list(
          type = smoother_type,
          var = var_name,
          full_name = gsub("([0-9]+$)", "", match_text)  # Remove any numeric suffix
        ))
      }
      
      return(NULL)
    }
    
    # Track detected smoothers
    found_smoothers <- list()
    
    # Check each coefficient name for smoothers
    for (coef_name in coef_names) {
      smoother_info <- extract_smoother_info(coef_name)
      
      if (!is.null(smoother_info)) {
        full_name <- smoother_info$full_name
        
        # Only add each unique smoother once
        if (!(full_name %in% names(found_smoothers))) {
          found_smoothers[[full_name]] <- smoother_info
        }
      }
    }
    
    # Add found smoothers to terms list and smoother info
    for (term_name in names(found_smoothers)) {
      smoother_info <- found_smoothers[[term_name]]
      
      # Add to terms list if not already there
      if (!(term_name %in% obj$terms)) {
        obj$terms <- c(obj$terms, term_name)
      }
      
      # Add smoother information
      obj$smoothers[[term_name]] <- list(
        type = smoother_info$type,
        vars = smoother_info$var,
        label = term_name
      )
    }
    
    # Check for interactions with year or focus in GLM terms
    for (term in obj$terms) {
      if (grepl(":", term, fixed = TRUE)) {
        # This is an interaction term
        parts <- strsplit(term, ":", fixed = TRUE)[[1]]
        obj$interactions[[term]] <- list(
          type = "parametric",
          vars = parts,
          label = term
        )
      }
    }
  }
  
  # Set response and focus if not provided
  if (is.null(obj$response)) {
    obj$response <- as.character(formula(obj$model)[[2]])
  }
  
  if (is.null(obj$focus)) {
    obj$focus <- obj$terms[1]
  }
  
  # Initialize labels and orders
  obj$labels <- setNames(as.list(obj$terms), obj$terms)
  obj$orders <- setNames(as.list(rep("asis", length(obj$terms))), obj$terms)
  
  return(obj)
}


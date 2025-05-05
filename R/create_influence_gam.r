# Ensure necessary packages are available (Ideally handled by DESCRIPTION/NAMESPACE Imports)
# These checks are useful for scripts, but in a package, dependencies are declared.
# We keep them here for script-like usability if files are sourced directly.
if (!requireNamespace("mgcv", quietly = TRUE)) {
  stop("Package 'mgcv' is required but not installed. Please install it.")
}
# Load libraries (Again, Imports are preferred in packages)
library(mgcv)

# -----------------------------------------------------------------------------
# Constructor Function
# -----------------------------------------------------------------------------

#' Create an influence_gam Object
#'
#' Initializes an object for performing influence analysis on a GAM model.
#'
#' @param model A fitted GAM object (from `mgcv::gam`).
#' @param data The data frame used to fit the model. Essential if not stored in the model object.
#'        Ensure variable types (especially factors) are correct *before* fitting the model.
#' @param focus The name (character string) of the focus term (predictor variable)
#'   for which influence is to be calculated. This term must be present in the model.
#'   It can be a factor, numeric, or involved in smooths/interactions.
#' @param response Optional. The name (character string) of the response variable.
#'   If NULL, it will be inferred from the model formula.
#' @param verbose Logical. Print messages during initialization?
#'
#' @return An object of class `influence_gam`.
#' @export
#' @examples
#' \dontrun{
#' # Example usage (requires data and a fitted gam model)
#' set.seed(123)
#' dat <- gamSim(1, n = 200, dist = "normal", scale = 2)
#' dat$fac <- factor(sample(1:3, 200, replace = TRUE))
#' model <- gam(y ~ s(x0) + s(x1) + s(x2, by = fac) + fac, data = dat)
#' influ_obj <- create_influence_gam(model, data = dat, focus = "fac")
#' influ_obj <- calculate_influence(influ_obj)
#' plot(influ_obj, type = "stan")
#' plot(influ_obj, type = "step")
#' plot(influ_obj, type = "influ")
#' plot(influ_obj, type = "cdi", term = "s(x0)") # Note: CDI plot needs adaptation for smooths
#' }
create_influence_gam <- function(model, data = NULL, focus = NULL, response = NULL, verbose = TRUE) {

  # --- Input Validation ---
  if (!inherits(model, "gam")) {
    stop("'model' must be a fitted GAM object from the 'mgcv' package.")
  }

  # Attempt to extract data if not provided
  if (is.null(data)) {
    if (verbose) message("Attempting to extract data from model object.")
    data <- model$model # This usually holds the model frame
    if (is.null(data)) {
        # Try harder for some cases
        data_name <- model$call$data
        if (!is.null(data_name)){
             data <- try(eval(data_name, environment(formula(model))), silent=TRUE)
             if(inherits(data, "try-error")) data <- NULL
        }
    }
     if (is.null(data)) {
         stop("Could not extract data from model. Please provide the 'data' argument explicitly.")
     }
  }
  if (!is.data.frame(data)) {
      stop("'data' must be a data frame.")
  }

  # Infer response if not provided
  if (is.null(response)) {
    response <- as.character(stats::formula(model)[[2]]) # Use stats::formula explicitly
    if (verbose) message("Inferred response variable: ", response)
  }
  if (!response %in% names(data)) {
    stop("Response variable '", response, "' not found in the provided data.")
  }

  # Extract all term labels (including smooths)
  model_terms <- attr(stats::terms(model), "term.labels") # Use stats::terms explicitly

  # Validate focus term
  if (is.null(focus)) {
    focus <- model_terms[1]
    if (verbose) message("No 'focus' term provided. Using the first term: ", focus)
  }

  # Check if focus term exists (can be tricky with by= variables etc.)
  # We check if the focus variable name appears in the data and model terms
  focus_var_name <- gsub("s\\(|\\)|te\\(|ti\\(|by *= *", "", focus) # Basic cleaning
  focus_var_name <- strsplit(focus_var_name, ",")[[1]][1] # Get first part
  focus_var_name <- trimws(focus_var_name)

  if (!focus_var_name %in% names(data)) {
       stop("Focus variable '", focus_var_name, "' (derived from focus term '", focus, "') not found in data.")
  }
  # More robust check: does the focus variable appear in the model terms?
  term_check <- any(grepl(focus_var_name, model_terms, fixed = TRUE))
  if (!term_check) {
      # Check if it's the intercept (though intercept isn't usually a focus)
      if(focus != "intercept"){
          warning("Focus term '", focus, "' or its base variable '", focus_var_name, "' not directly found among model term labels: ", paste(model_terms, collapse=", "), ". Proceeding, but check results carefully.")
      }
  }
  # Warning if focus is not a factor - recommend setting type *before* model fitting
  if (!is.factor(data[[focus_var_name]])) {
      warning("Focus variable '", focus_var_name, "' is not a factor. Influence calculations assume a factor. Please ensure the data used for fitting 'gam' has the correct variable type. Converting to factor internally for calculations.")
      data[[focus_var_name]] <- as.factor(data[[focus_var_name]])
      # Note: This internal conversion happens *after* the model is potentially fitted.
      # Best practice is to ensure data[[focus_var_name]] is a factor *before* calling gam().
  }


  # --- Create Object ---
  obj <- list(
    model = model,
    data = as.data.frame(data), # Ensure it's a standard data frame
    response = response,
    focus = focus,
    focus_var = focus_var_name, # Store the derived variable name
    terms = model_terms,
    # Placeholders for results
    indices = NULL,
    summary = NULL,
    influences = NULL,
    preds = NULL,
    calculated = FALSE # Flag to check if calculations have been run
  )

  class(obj) <- "influence_gam"

  if (verbose) message("influence_gam object created successfully for focus term '", focus, "'. Run calculate_influence() next.")

  return(obj)
}
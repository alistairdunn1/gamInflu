#' @title Create a GAM Influence Object
#' @description Initializes a gam_influence object for a given model and focus term. This is the
#' main constructor for the 'gam_influence' class. It sets up the basic structure,
#' which is then populated by `calculate_influence()`.
#' @param model A fitted model object (supported: `glm`, `gam`).
#' @param focus A character string specifying the name of the focus term. This
#'   term must be a factor in the model's data and is typically the term for
#'   which an annual index is being calculated (e.g., "year").
#' @param data The data frame used to fit the model. This is crucial for models
#'   that do not store the data internally. If NULL, the function will attempt
#'   to extract it from the model object.
#' @param islog A logical value indicating whether the response variable is on a
#'   logarithmic scale.
#' @return An object of S3 class 'gam_influence'. This is a list containing the model,
#'   data, and key term specifications.
#' @importFrom stats terms formula
#' @export
gam_influence <- function(model, focus, data = NULL, islog = NULL) {
  # --- Data Extraction and Validation ---
  if (is.null(data)) {
    data <- model$data
  }
  if (is.null(data)) {
    # Attempt to retrieve data from the model's call if not explicitly provided
    tryCatch(
      {
        data <- eval(model$call$data, environment(formula(model)))
      },
      error = function(e) {
        stop("Could not find model data. Please supply it using the 'data' argument.", call. = FALSE)
      }
    )
  }
  if (is.null(data)) {
    stop("Failed to extract or find the data frame used for the model.", call. = FALSE)
  }

  if (!focus %in% names(data)) {
    stop(paste0("Focus term '", focus, "' not found in the provided data."), call. = FALSE)
  }
  if (!is.factor(data[[focus]])) {
    stop(paste0("The focus term '", focus, "' must be a factor."), call. = FALSE)
  }

  # --- Term and Response Variable Extraction ---
  term_labels <- attr(terms(model), "term.labels")
  response_var <- all.vars(formula(model))[1]

  # --- Structure and Class Assignment ---
  # The main object is a list with the 'gam_influence' class assigned.
  obj <- structure(
    list(
      model = model,
      data = as.data.frame(data),
      focus = focus,
      response = response_var,
      terms = term_labels,
      islog = islog,
      calculated = FALSE # This list will be populated by calculate_influence()
    ),
    class = "gam_influence"
  )

  return(obj)
}

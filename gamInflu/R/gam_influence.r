#' @title Create a GAM Influence Object
#' @description Initializes a gam_influence object for a given model and focus term. This is the
#' main constructor for the 'gam_influence' class that supports multiple GLM families
#' including Gaussian, binomial, gamma, and Poisson distributions. It sets up the basic structure,
#' which is then populated by `calculate_influence()` with family-specific methods.
#' @param model A fitted model object from mgcv (gam) or stats (glm). Supported families:
#'   Gaussian, binomial, Gamma, Poisson, and their quasi variants.
#' @param focus A character string specifying the name of the focus term. This
#'   term must be a factor in the model's data and is typically the term for
#'   which an annual index is being calculated (e.g., "year"). The focus term
#'   should represent the main temporal or categorical variable of interest.
#' @param data The data frame used to fit the model. This is crucial for models
#'   that do not store the data internally. If NULL, the function will attempt
#'   to extract it from the model object.
#' @param islog A logical value indicating whether the response variable is on a
#'   logarithmic scale. If NULL, this will be inferred automatically by
#'   `calculate_influence()` based on the response variable name and model family.
#' @return An object of S3 class 'gam_influence'. This is a list containing the model,
#'   data, focus term, response variable, and key term specifications that will be
#'   used for family-specific influence calculations.
#' @details
#' The gam_influence object stores essential information for subsequent influence analysis:
#' - Model object with family information for automatic family detection
#' - Data validation ensuring focus term is present and properly formatted
#' - Term extraction for stepwise model building
#' - Response variable identification for family-specific calculations
#'
#' **Supported Model Families:**
#' - **Gaussian**: For continuous response data (CPUE, biomass indices)
#' - **Binomial**: For binary or proportion data (presence/absence, catch probability)
#' - **Gamma**: For positive continuous data (biomass, positive CPUE values)
#' - **Poisson**: For count data (fish numbers, catch counts)
#'
#' The focus term must be a factor variable representing the main grouping of interest
#' (typically years for fisheries applications, but could be sites, regions, etc.).
#' @examples
#' \dontrun{
#' library(mgcv)
#'
#' # Gaussian model (traditional CPUE)
#' mod_gauss <- gam(log_cpue ~ s(depth) + s(temp) + year, data = fish_data)
#' gi_gauss <- gam_influence(mod_gauss, focus = "year")
#'
#' # Binomial model (presence/absence)
#' mod_binom <- gam(presence ~ s(depth) + year, family = binomial(), data = fish_data)
#' gi_binom <- gam_influence(mod_binom, focus = "year")
#'
#' # Gamma model (biomass)
#' mod_gamma <- gam(biomass ~ s(effort) + year, family = Gamma(link = "log"), data = survey_data)
#' gi_gamma <- gam_influence(mod_gamma, focus = "year")
#'
#' # Calculate influence with automatic family detection
#' gi_gauss <- calculate_influence(gi_gauss)
#' gi_binom <- calculate_influence(gi_binom)
#' gi_gamma <- calculate_influence(gi_gamma)
#' }
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

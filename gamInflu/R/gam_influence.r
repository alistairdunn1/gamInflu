#' Create a GAM Influence Analysis Object
#'
#' Creates an S3 object for analysing the influence of terms in a GAM model.
#'
#' @param model A fitted GAM model from mgcv (gam, bam, or gamm)
#' @param focus The focus term for influence analysis (smooth or parametric)
#' @param data Optional data frame used to fit the model
#' @param response Optional response variable name
#'
#' @return A 'gam_influence' S3 object
#' @export
#'
#' @examples
#' \dontrun{
#' library(mgcv)
#' data(mcycle)
#' m1 <- gam(accel ~ s(times) + s(times, by = factor), data = mcycle)
#' gi <- gam_influence(m1, focus = "s(times)")
#' }
gam_influence <- function(model, focus = NULL, data = NULL, response = NULL) {
  # Validate input
  if (!inherits(model, c("gam", "bam", "gamm"))) {
    stop("Model must be a GAM fitted with mgcv (gam, bam, or gamm)")
  }

  # Extract model components
  if (inherits(model, "gamm")) {
    gam_model <- model$gam
    data_used <- model$gam$model
  } else {
    gam_model <- model
    data_used <- model$model
  }

  # Get data
  if (is.null(data)) {
    if (is.null(data_used)) {
      stop("Cannot extract data from model. Please provide 'data' argument.")
    }
    data <- data_used
  }

  # Get response variable
  if (is.null(response)) {
    response <- all.vars(formula(gam_model))[1]
  }

  # Extract model terms information
  terms_info <- extract_gam_terms(gam_model)

  # Set focus term
  if (is.null(focus)) {
    focus <- terms_info$smooth_terms[1]
    if (is.na(focus)) focus <- terms_info$parametric_terms[1]
  }

  # Create S3 object
  obj <- list(
    model = gam_model,
    original_model = model,
    data = data,
    response = response,
    focus = focus,
    terms_info = terms_info,
    calculated = FALSE
  )

  obj <- calculate_gam_influence(obj)

  class(obj) <- "gam_influence"
  return(obj)
}

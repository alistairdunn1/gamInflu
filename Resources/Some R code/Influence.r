#' A package for generating step plots, influence plots, CDI plots, and influence metrics for linear models and GAMs
#'
#' The concept of influence in generalised linear models is described in
#' Bentley, N., Kendrick, T. H., Starr, P. J., & Breen, P. A. (2011). Influence plots and metrics: tools for better understanding fisheries catch-per-unit-effort standardisations.
#' This package provides an implementation of the plots and metrics described in that paper. ICES Journal of Marine Science, doi:10.1093/icesjms/fsr174.
#'
#' This package works for \code{glm} and \code{gam} models with log transformed dependent variables. These
#' are the type of models commonly used for one part of the delta-lognormal approach to catch-per-unit-effort (CPUE)
#' standardisation. e.g.\code{model = gam(log(catch)~s(year)+ti(month,vessel)+effort)} or with interactions
#' \code{model = gam(log(catch)~s(year)+ti(month,year)+effort)}
#'
#' @docType package
#' @name influ-package
#' @aliases influ
#' @title Influence in linear models and GAMs with interaction support
#' @author Nokome Bentley (Original), Alistair Dunn (Modified to support GAMs and ggplot2)


#' Create a new Influence object.
#'
#' A new Influence object needs to be created for each combination of a model and focus term.
#' For example, you might compare the influence of the same term in two separate models:
#' \code{
#'    influ1 = influence(model1, 'year')
#'    influ2 = influence(model2, 'year')
#' }
#'
#' @param model The model for which the influence of terms will be determined
#' @param focus The focus term for which influence is to be calculated
#' @param data The data which was used to fit the model (required for some types of models that do not store the data)
#' @param response The response term in the model
#' @param colour The color to use for plots
#' @return A new Influence object (S3)
#' @export
influence <- function(model, focus = NULL, data = NULL, response = NULL, colour = "steelblue") {
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(mgcv)
  library(gridExtra)
  library(patchwork)
  library(stringr)
  
  # Create a new influence object
  obj <- structure(list(
    model = model,
    data = data,
    response = response,
    focus = focus,
    colour = colour,
    terms = NULL,
    labels = NULL,
    orders = NULL,
    indices = NULL,
    summary = NULL,
    preds = NULL,
    influences = NULL,
    interactions = NULL  # New field to store interaction information
  ), class = "influence")
  
  # Initialize the object
  obj <- initialize_influence(obj)
  
  return(obj)
}


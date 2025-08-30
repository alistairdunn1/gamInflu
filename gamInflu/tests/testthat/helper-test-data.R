# Test helper functions for setting up test data and models
library(mgcv)
library(testthat)

#' Create test data for different GLM families
#' @param n Number of observations
#' @param seed Random seed for reproducibility
create_test_data <- function(n = 200, seed = 123) {
  set.seed(seed)

  data <- data.frame(
    year = factor(rep(2015:2019, each = n / 5)),
    depth = runif(n, 10, 100),
    temp = rnorm(n, 15, 3),
    area = factor(sample(c("North", "South"), n, replace = TRUE)),
    vessel = factor(sample(paste0("V", 1:8), n, replace = TRUE))
  )

  # Create base linear predictor
  data$linear_pred <- with(
    data,
    2 + 0.1 * as.numeric(year) + sin(depth / 20) + temp / 10 +
      ifelse(area == "North", 0.5, -0.5)
  )

  # Add noise
  data$noise <- rnorm(n, 0, 0.3)

  # Generate responses for different families
  data$y_gaussian <- data$linear_pred + data$noise
  data$log_catch <- data$linear_pred + data$noise # For lognormal testing
  data$catch <- exp(data$log_catch) # Anti-logged version

  data$y_binomial <- rbinom(n, 1, plogis(data$linear_pred + data$noise))

  # Ensure positive values for Gamma
  gamma_pred <- exp(pmax(data$linear_pred + data$noise, -5))
  data$y_gamma <- rgamma(n, shape = 2, rate = 2 / gamma_pred)

  data$y_poisson <- rpois(n, exp(pmax(data$linear_pred + data$noise - 1, -5)))

  # Weibull response (shape parameter = 2, scale from linear predictor)
  weibull_scale <- exp(pmax(data$linear_pred + data$noise, -5))
  data$y_weibull <- rweibull(n, shape = 2, scale = weibull_scale)

  # Tweedie response (using power = 1.5)
  # For simplicity, approximate with Gamma distribution
  data$y_tweedie <- data$y_gamma

  return(data)
}

#' Fit test models for different families
#' @param data Test data frame
#' @param formula_base Base formula (without response)
fit_test_models <- function(data, formula_base = "~ s(depth) + s(temp) + year + area") {
  models <- list()

  # Store the data globally to ensure gam_influence can find it
  assign("test_data_global", data, envir = .GlobalEnv)

  # Gaussian (normal) model
  models$gaussian <- gam(
    as.formula(paste("y_gaussian", formula_base)),
    data = data,
    family = gaussian()
  )

  # Lognormal models (two versions)
  # Version 1: Pre-logged response with gaussian family
  models$lognormal_prelogged <- gam(
    as.formula(paste("log_catch", formula_base)),
    data = data,
    family = gaussian()
  )

  # Version 2: Raw response with Gamma family and log link
  models$lognormal_gamma <- gam(
    as.formula(paste("catch", formula_base)),
    data = data,
    family = Gamma(link = "log")
  )

  # Binomial model
  models$binomial <- gam(
    as.formula(paste("y_binomial", formula_base)),
    data = data,
    family = binomial()
  )

  # Gamma model
  models$gamma <- gam(
    as.formula(paste("y_gamma", formula_base)),
    data = data,
    family = Gamma(link = "log")
  )

  # Poisson model
  models$poisson <- gam(
    as.formula(paste("y_poisson", formula_base)),
    data = data,
    family = poisson()
  )

  # Weibull model
  models$weibull <- gam(
    as.formula(paste("y_weibull", formula_base)),
    data = data,
    family = weibull_family()
  )

  # Tweedie model (if available)
  if (requireNamespace("tweedie", quietly = TRUE)) {
    models$tweedie <- gam(
      as.formula(paste("y_tweedie", formula_base)),
      data = data,
      family = tw(link = "log")
    )
  }

  # Store data reference in each model for gam_influence
  for (i in seq_along(models)) {
    attr(models[[i]], "test_data") <- data
  }

  return(models)
}

#' Create interaction test models
#' @param data Test data frame
fit_interaction_models <- function(data) {
  models <- list()

  # Store the data globally to ensure gam_influence can find it
  assign("test_data_global", data, envir = .GlobalEnv)

  # Gaussian with interaction
  models$gaussian_interaction <- gam(
    y_gaussian ~ s(depth) + s(temp) + year * area,
    data = data,
    family = gaussian()
  )

  # Gamma with interaction
  models$gamma_interaction <- gam(
    y_gamma ~ s(depth) + s(temp) + year * area,
    data = data,
    family = Gamma(link = "log")
  )

  # Store data reference in each model for gam_influence
  for (i in seq_along(models)) {
    attr(models[[i]], "test_data") <- data
  }

  return(models)
}

#' Create random effects test models
#' @param data Test data frame
fit_random_effects_models <- function(data) {
  models <- list()

  # Store the data globally to ensure gam_influence can find it
  assign("test_data_global", data, envir = .GlobalEnv)

  # Random effects model
  models$random_effects <- gam(
    y_gaussian ~ s(depth) + s(temp) + year + s(vessel, bs = "re"),
    data = data,
    family = gaussian()
  )

  # Random effects with by-variable
  models$random_by_effects <- gam(
    y_gaussian ~ s(depth) + s(temp) + year + s(vessel, by = area, bs = "re"),
    data = data,
    family = gaussian()
  )

  # Store data reference in each model for gam_influence
  for (i in seq_along(models)) {
    attr(models[[i]], "test_data") <- data
  }

  return(models)
}

# Simple comprehensive test for gamInflu package
# This test focuses on the key functionality with proper data handling

library(testthat)
library(mgcv)

test_that("Core gamInflu functionality works with all families and islog settings", {
  # Create test data
  set.seed(123)
  n <- 200
  test_data <- data.frame(
    year = factor(rep(2015:2019, each = 40)),
    depth = runif(n, 10, 100),
    temp = rnorm(n, 15, 3),
    area = factor(sample(c("North", "South"), n, replace = TRUE))
  )

  # Generate different response types
  test_data$linear_pred <- with(test_data, 2 + 0.1 * as.numeric(year) + sin(depth / 20) + temp / 10)
  test_data$log_catch <- test_data$linear_pred + rnorm(n, 0, 0.3)
  test_data$catch <- exp(test_data$log_catch)
  test_data$y_binomial <- rbinom(n, 1, plogis(test_data$linear_pred + rnorm(n, 0, 0.3)))
  test_data$y_gamma <- rgamma(n, shape = 2, rate = 2 / exp(pmax(test_data$linear_pred + rnorm(n, 0, 0.3), -5)))
  test_data$y_gaussian <- test_data$linear_pred + rnorm(n, 0, 0.3)

  # Test 1: Lognormal (pre-logged) with Gaussian family, islog = FALSE
  model_log_gauss <- gam(log_catch ~ s(depth) + s(temp) + year, data = test_data, family = gaussian())
  gi_log_false <- calculate_influence(gam_influence(model_log_gauss, focus = "year", data = test_data), islog = FALSE)

  expect_s3_class(gi_log_false, "gam_influence")
  indices_log_false <- extract_indices(gi_log_false)
  expect_true(all(c("index", "cv", "lower_CI", "upper_CI") %in% names(indices_log_false)))
  expect_equal(nrow(indices_log_false), 5)

  # Test 2: Lognormal (pre-logged) with Gaussian family, islog = TRUE
  gi_log_true <- calculate_influence(gam_influence(model_log_gauss, focus = "year", data = test_data), islog = TRUE)
  indices_log_true <- extract_indices(gi_log_true)

  expect_s3_class(gi_log_true, "gam_influence")
  expect_equal(nrow(indices_log_true), 5)
  # islog=TRUE should anti-log the results
  expect_equal(indices_log_true$index, exp(indices_log_false$index), tolerance = 1e-10)

  # Test 3: Gamma family with log link, islog = FALSE
  model_gamma <- gam(catch ~ s(depth) + s(temp) + year, data = test_data, family = Gamma(link = "log"))
  gi_gamma_false <- calculate_influence(gam_influence(model_gamma, focus = "year", data = test_data), islog = FALSE)
  indices_gamma_false <- extract_indices(gi_gamma_false)

  expect_s3_class(gi_gamma_false, "gam_influence")
  expect_equal(nrow(indices_gamma_false), 5)
  expect_true(all(is.finite(indices_gamma_false$index)))

  # Test 4: Gamma family with log link, islog = TRUE
  gi_gamma_true <- calculate_influence(gam_influence(model_gamma, focus = "year", data = test_data), islog = TRUE)
  indices_gamma_true <- extract_indices(gi_gamma_true)

  expect_s3_class(gi_gamma_true, "gam_influence")
  expect_equal(nrow(indices_gamma_true), 5)
  expect_false(identical(indices_gamma_false$index, indices_gamma_true$index))

  # Test 5: Binomial family, islog = FALSE
  model_binomial <- gam(y_binomial ~ s(depth) + s(temp) + year, data = test_data, family = binomial())
  gi_binomial <- calculate_influence(gam_influence(model_binomial, focus = "year", data = test_data), islog = FALSE)
  indices_binomial <- extract_indices(gi_binomial)

  expect_s3_class(gi_binomial, "gam_influence")
  expect_equal(nrow(indices_binomial), 5)
  expect_true(all(is.finite(indices_binomial$index)))

  # Test 6: Gaussian family, islog = FALSE
  model_gaussian <- gam(y_gaussian ~ s(depth) + s(temp) + year, data = test_data, family = gaussian())
  gi_gaussian <- calculate_influence(gam_influence(model_gaussian, focus = "year", data = test_data), islog = FALSE)
  indices_gaussian <- extract_indices(gi_gaussian)

  expect_s3_class(gi_gaussian, "gam_influence")
  expect_equal(nrow(indices_gaussian), 5)
  expect_true(all(is.finite(indices_gaussian$index)))

  message("✓ Core functionality works with all families and islog settings")
})

test_that("Plotting functions work", {
  # Create simple test data
  set.seed(123)
  n <- 200
  test_data <- data.frame(
    year = factor(rep(2015:2019, each = 40)),
    depth = runif(n, 10, 100),
    temp = rnorm(n, 15, 3)
  )
  test_data$log_catch <- 2 + 0.1 * as.numeric(test_data$year) + sin(test_data$depth / 20) + rnorm(n, 0, 0.3)

  model <- gam(log_catch ~ s(depth) + s(temp) + year, data = test_data, family = gaussian())
  gi <- calculate_influence(gam_influence(model, focus = "year", data = test_data), islog = TRUE)

  # Test main plotting functions
  expect_no_error(plot_standardisation(gi))
  expect_no_error(plot_stepwise_index(gi))
  expect_no_error(plot_step_and_influence(gi))

  # Test term-specific plotting
  terms <- get_terms(gi)
  if (length(terms) > 0) {
    expect_no_error(plot_terms(gi, terms[1]))
    expect_no_error(plot_cdi(gi, terms[1]))
    expect_no_error(plot_term_distribution(gi, terms[1]))
  }

  # Test generic plot method
  expect_no_error(plot(gi, type = "standardisation"))
  expect_no_error(plot(gi, type = "stepwise"))

  # Test residual plots
  expect_no_error(plot_residuals(gi))

  message("✓ Plotting functions work correctly")
})

test_that("Utility functions work", {
  # Test geometric_mean
  x <- c(1, 2, 4, 8)
  expected <- exp(mean(log(x)))
  expect_equal(geometric_mean(x), expected)

  # Create test model for other utilities
  set.seed(123)
  n <- 100
  test_data <- data.frame(
    year = factor(rep(2015:2019, each = 20)),
    x = rnorm(n),
    y = rnorm(n)
  )

  model <- gam(y ~ s(x) + year, data = test_data, family = gaussian())
  gi <- calculate_influence(gam_influence(model, focus = "year", data = test_data))

  # Test get_terms
  terms <- get_terms(gi)
  expect_true(is.character(terms))
  expect_true(length(terms) > 0)

  # Test r2 function
  r2_result <- r2(gi)
  expect_true(is.numeric(r2_result))
  expect_true(all(r2_result >= 0 & r2_result <= 1))

  # Test S3 methods
  expect_no_error(summary(gi))
  expect_no_error(print(gi))

  message("✓ Utility functions work correctly")
})

test_that("Advanced features work", {
  # Create test data
  set.seed(123)
  n <- 200
  test_data <- data.frame(
    year = factor(rep(2015:2019, each = 40)),
    depth = runif(n, 10, 100),
    temp = rnorm(n, 15, 3)
  )
  test_data$y <- 2 + 0.1 * as.numeric(test_data$year) + sin(test_data$depth / 20) + rnorm(n, 0, 0.3)

  model <- gam(y ~ s(depth) + s(temp) + year, data = test_data, family = gaussian())
  gi <- calculate_influence(gam_influence(model, focus = "year", data = test_data))

  # Test residual pattern analysis
  rpa <- analyse_residual_patterns(gi)
  expect_s3_class(rpa, "residual_pattern_analysis")
  expect_true(all(c("linear_results", "recommendations", "analysis_info") %in% names(rpa)))
  expect_no_error(print(rpa))

  message("✓ Advanced features work correctly")
})

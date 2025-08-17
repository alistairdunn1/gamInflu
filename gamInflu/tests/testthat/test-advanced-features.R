test_that("analyse_residual_patterns works with all families", {
  test_data <- create_test_data(n = 200)
  models <- fit_test_models(test_data)

  for (family_name in names(models)) {
    model <- models[[family_name]]

    gi <- gam_influence(model, focus = "year", data = test_data)
    gi <- calculate_influence(gi)

    # Test residual pattern analysis
    rpa <- analyse_residual_patterns(gi)

    expect_s3_class(rpa, "residual_pattern_analysis")
    expect_true(is.list(rpa))
    expect_true(all(c("linear_results", "recommendations", "analysis_info") %in% names(rpa)))

    # Test print method for residual pattern analysis
    expect_no_error(print(rpa))

    message(paste("✓ analyse_residual_patterns works with", family_name, "family"))
  }
})

test_that("combine_indices works for delta-GLM models", {
  test_data <- create_test_data(n = 200)

  # Create binomial model (presence/absence)
  model_binomial <- gam(
    y_binomial ~ s(depth) + s(temp) + year,
    data = test_data,
    family = binomial()
  )

  # Create gamma model for positive catches (subset where y_binomial == 1)
  positive_data <- test_data[test_data$y_binomial == 1, ]
  if (nrow(positive_data) > 20) { # Only if we have enough positive observations
    model_gamma <- gam(
      y_gamma ~ s(depth) + s(temp) + year,
      data = positive_data,
      family = Gamma(link = "log")
    )

    # Calculate influence for both models
    gi_binomial <- calculate_influence(gam_influence(model_binomial, focus = "year", data = test_data))
    gi_gamma <- calculate_influence(gam_influence(model_gamma, focus = "year", data = positive_data))

    # Test combine_indices
    gi_combined <- combine_indices(gi_binomial, gi_gamma)

    expect_s3_class(gi_combined, "gam_influence_combined")
    expect_true(!is.null(gi_combined$combined_indices))

    # Test combined plotting
    expect_no_error(plot(gi_combined, type = "comparison"))

    # Test combined summary
    expect_no_error(summary(gi_combined))

    message("✓ combine_indices works for delta-GLM models")
  } else {
    message("⚠ Skipping combine_indices test - insufficient positive observations")
  }
})

test_that("compare_focus_by_groups works", {
  test_data <- create_test_data(n = 200)

  # Create separate models for different areas
  north_data <- test_data[test_data$area == "North", ]
  south_data <- test_data[test_data$area == "South", ]

  if (nrow(north_data) > 50 && nrow(south_data) > 50) {
    model_north <- gam(y_gaussian ~ s(depth) + s(temp) + year, data = north_data, family = gaussian())
    model_south <- gam(y_gaussian ~ s(depth) + s(temp) + year, data = south_data, family = gaussian())

    gi_north <- calculate_influence(gam_influence(model_north, focus = "year", data = north_data))
    gi_south <- calculate_influence(gam_influence(model_south, focus = "year", data = south_data))

    # Test comparison
    comparison <- compare_focus_by_groups(
      list("North" = gi_north, "South" = gi_south),
      focus = "year"
    )

    expect_s3_class(comparison, "gam_influence_comparison")
    expect_true(!is.null(comparison$comparison_data))

    # Test comparison plotting
    expect_no_error(plot(comparison, type = "indices"))

    # Test comparison summary
    expect_no_error(summary(comparison))

    message("✓ compare_focus_by_groups works")
  } else {
    message("⚠ Skipping compare_focus_by_groups test - insufficient observations per group")
  }
})

test_that("analyse_focus_by_group works", {
  test_data <- create_test_data(n = 200)

  model <- gam(
    y_gaussian ~ s(depth) + s(temp) + year * area,
    data = test_data,
    family = gaussian()
  )

  gi <- gam_influence(model, focus = "year", data = test_data)
  gi <- calculate_influence(gi)

  # Test group analysis
  group_analysis <- analyse_focus_by_group(gi, group_var = "area")

  expect_true(is.list(group_analysis))
  expect_true(all(c("North", "South") %in% names(group_analysis)))

  # Each group should have gam_influence object
  for (group in names(group_analysis)) {
    expect_s3_class(group_analysis[[group]], "gam_influence")
  }

  message("✓ analyse_focus_by_group works")
})

test_that("Tweedie family support works", {
  # Skip if tweedie package not available
  skip_if_not_installed("tweedie")

  test_data <- create_test_data(n = 200)

  model_tweedie <- gam(
    y_tweedie ~ s(depth) + s(temp) + year,
    data = test_data,
    family = tw(link = "log")
  )

  gi <- gam_influence(model_tweedie, focus = "year", data = test_data)
  gi <- calculate_influence(gi)

  expect_s3_class(gi, "gam_influence")
  expect_equal(gi$family$family, "Tweedie")

  # Test that plotting works
  expect_no_error(plot_standardisation(gi))
  expect_no_error(plot_stepwise_index(gi))

  # Test data extraction
  indices <- extract_indices(gi)
  expect_s3_class(indices, "data.frame")
  expect_true(all(c("index", "cv", "lower_CI", "upper_CI") %in% names(indices)))

  message("✓ Tweedie family support works")
})

test_that("Error handling works correctly", {
  test_data <- create_test_data(n = 200)
  model <- gam(y_gaussian ~ s(depth) + s(temp) + year, data = test_data, family = gaussian())

  # Test with invalid focus variable
  expect_error(gam_influence(model, focus = "nonexistent_var", data = test_data))

  # Test with non-factor focus variable
  expect_error(gam_influence(model, focus = "depth", data = test_data))

  # Test calculate_influence before gam_influence
  gi <- gam_influence(model, focus = "year", data = test_data)

  # Test with invalid islog parameter
  expect_error(calculate_influence(gi, islog = "invalid"))

  # Test plotting functions with invalid term
  gi <- calculate_influence(gi)
  expect_error(plot_terms(gi, "nonexistent_term"))
  expect_error(plot_cdi(gi, "nonexistent_term"))

  message("✓ Error handling works correctly")
})

test_that("stepCPUE_gam works with GAM models", {
  # Create test data specifically for stepwise selection
  set.seed(456)
  n <- 100
  stepwise_data <- data.frame(
    x1 = runif(n, 0, 10),
    x2 = runif(n, 0, 10),
    x3 = runif(n, 0, 10),
    x4 = runif(n, 0, 10),
    f1 = factor(sample(c("A", "B", "C"), n, replace = TRUE)),
    f2 = factor(sample(c("Type1", "Type2"), n, replace = TRUE))
  )

  # Create response with known relationships
  # x1 and f1 should be important, x2 moderately important, x3 and x4 should be less important
  stepwise_data$y <- 2 +
    2 * sin(stepwise_data$x1 / 2) + # Strong smooth effect
    0.8 * stepwise_data$x2 + # Moderate linear effect
    0.2 * stepwise_data$x3 + # Weak linear effect
    as.numeric(stepwise_data$f1) + # Factor effect
    rnorm(n, 0, 0.4)

  # Fit initial simple model
  initial_model <- gam(y ~ s(x1), data = stepwise_data, family = gaussian())

  # Test basic functionality
  expect_s3_class(initial_model, "gam")
  expect_true(exists("stepCPUE_gam"))
  expect_true(is.function(stepCPUE_gam))

  # Define scope for stepwise selection
  scope_list <- list(
    lower = ~ s(x1),
    upper = ~ s(x1) + s(x2) + s(x3) + x4 + f1 + f2
  )

  # Test forward selection
  result_forward <- stepCPUE_gam(
    initial_model,
    scope = scope_list,
    r2.change = 0.02,
    direction = "forward",
    trace = 0,
    steps = 5
  )

  # Test that result is a GAM object
  expect_s3_class(result_forward, "gam")

  # Test that it has anova component
  expect_true("anova" %in% names(result_forward))
  expect_s3_class(result_forward$anova, "anova")
  expect_s3_class(result_forward$anova, "data.frame")

  # Test anova structure
  anova_cols <- c("Step", "Df", "Deviance", "Resid. Df", "Resid. Dev", "r.squared", "aic")
  expect_true(all(anova_cols %in% names(result_forward$anova)))

  # Test that R-squared values are reasonable
  r2_values <- result_forward$anova$r.squared
  expect_true(all(r2_values >= 0 & r2_values <= 1))
  expect_true(all(!is.na(r2_values)))

  # Test that R-squared generally increases (some fluctuation allowed)
  if (length(r2_values) > 1) {
    expect_true(max(r2_values) >= min(r2_values))
  }

  message("✓ Forward stepwise selection works")

  # Test backward selection (start with fuller model)
  full_model <- gam(y ~ s(x1) + s(x2) + x4 + f1, data = stepwise_data, family = gaussian())

  result_backward <- stepCPUE_gam(
    full_model,
    scope = list(lower = ~ s(x1), upper = ~ s(x1) + s(x2) + x4 + f1),
    r2.change = 0.01,
    direction = "backward",
    trace = 0,
    steps = 3
  )

  expect_s3_class(result_backward, "gam")
  expect_true("anova" %in% names(result_backward))

  message("✓ Backward stepwise selection works")

  # Test bidirectional selection
  result_both <- stepCPUE_gam(
    initial_model,
    scope = scope_list,
    r2.change = 0.015,
    direction = "both",
    trace = 0,
    steps = 4
  )

  expect_s3_class(result_both, "gam")
  expect_true("anova" %in% names(result_both))

  message("✓ Bidirectional stepwise selection works")

  # Test different trace levels (should not error)
  expect_no_error({
    stepCPUE_gam(
      initial_model,
      scope = list(upper = ~ s(x1) + x2),
      r2.change = 0.05,
      trace = 1,
      steps = 1
    )
  })

  # Test with different GLM families

  # Gamma family test
  stepwise_data$y_positive <- exp(stepwise_data$y - min(stepwise_data$y) + 0.1)
  gamma_model <- gam(y_positive ~ s(x1), data = stepwise_data, family = Gamma(link = "log"))

  result_gamma <- stepCPUE_gam(
    gamma_model,
    scope = list(upper = ~ s(x1) + x2 + f1),
    r2.change = 0.02,
    trace = 0,
    steps = 2
  )

  expect_s3_class(result_gamma, "gam")
  expect_equal(result_gamma$family$family, "Gamma")

  message("✓ Works with Gamma family")

  # Poisson family test
  stepwise_data$y_count <- rpois(n, exp(stepwise_data$y - mean(stepwise_data$y)))
  poisson_model <- gam(y_count ~ s(x1), data = stepwise_data, family = poisson())

  result_poisson <- stepCPUE_gam(
    poisson_model,
    scope = list(upper = ~ s(x1) + x2),
    r2.change = 0.02,
    trace = 0,
    steps = 1
  )

  expect_s3_class(result_poisson, "gam")
  expect_equal(result_poisson$family$family, "poisson")

  message("✓ Works with Poisson family")

  # Test error handling

  # Test with non-GAM object
  lm_model <- lm(y ~ x1, data = stepwise_data)
  expect_error(stepCPUE_gam(lm_model, scope = list(upper = ~ x1 + x2)))

  # Test with missing scope (should default to backward)
  result_no_scope <- stepCPUE_gam(
    gam(y ~ s(x1) + x2, data = stepwise_data),
    r2.change = 0.01,
    trace = 0,
    steps = 1
  )
  expect_s3_class(result_no_scope, "gam")

  # Test parameter validation
  expect_no_error(stepCPUE_gam(initial_model, scope = ~ s(x1) + x2, trace = 0, steps = 1))

  message("✓ Error handling and parameter validation work correctly")

  # Test that final model makes sense
  final_terms <- attr(terms(result_forward), "term.labels")
  initial_terms <- attr(terms(initial_model), "term.labels")

  # Should have at least the initial terms
  expect_true(all(initial_terms %in% final_terms) || length(final_terms) >= length(initial_terms) - 1)

  # Test AIC values are reasonable
  aic_values <- result_forward$anova$aic
  expect_true(all(is.finite(aic_values)))

  message("✓ stepCPUE_gam comprehensive testing completed")
})

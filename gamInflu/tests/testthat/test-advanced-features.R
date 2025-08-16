test_that("analyse_residual_patterns works with all families", {
  test_data <- create_test_data(n = 200)
  models <- fit_test_models(test_data)

  for (family_name in names(models)) {
    model <- models[[family_name]]

    gi <- gam_influence(model, focus = "year")
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
    gi_binomial <- calculate_influence(gam_influence(model_binomial, focus = "year"))
    gi_gamma <- calculate_influence(gam_influence(model_gamma, focus = "year"))

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

    gi_north <- calculate_influence(gam_influence(model_north, focus = "year"))
    gi_south <- calculate_influence(gam_influence(model_south, focus = "year"))

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

  gi <- gam_influence(model, focus = "year")
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

  gi <- gam_influence(model_tweedie, focus = "year")
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
  expect_error(gam_influence(model, focus = "nonexistent_var"))

  # Test with non-factor focus variable
  expect_error(gam_influence(model, focus = "depth"))

  # Test calculate_influence before gam_influence
  gi <- gam_influence(model, focus = "year")

  # Test with invalid islog parameter
  expect_error(calculate_influence(gi, islog = "invalid"))

  # Test plotting functions with invalid term
  gi <- calculate_influence(gi)
  expect_error(plot_terms(gi, "nonexistent_term"))
  expect_error(plot_cdi(gi, "nonexistent_term"))

  message("✓ Error handling works correctly")
})

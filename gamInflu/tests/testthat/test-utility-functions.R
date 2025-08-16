test_that("extract_indices works with all families and islog settings", {
  test_data <- create_test_data(n = 200)
  models <- fit_test_models(test_data)

  for (family_name in names(models)) {
    model <- models[[family_name]]

    # Test with islog = FALSE
    gi_false <- calculate_influence(gam_influence(model, focus = "year"), islog = FALSE)
    indices_false <- extract_indices(gi_false)

    expect_s3_class(indices_false, "data.frame")
    expect_true(all(c("index", "cv", "lower_CI", "upper_CI") %in% names(indices_false)))
    expect_equal(nrow(indices_false), 5) # 5 years
    expect_true(all(is.finite(indices_false$index)))
    expect_true(all(is.finite(indices_false$cv)))
    expect_true(all(indices_false$cv >= 0)) # CV should be non-negative

    # Test with islog = TRUE for applicable families
    if (family_name %in% c("lognormal_prelogged", "lognormal_gamma", "gamma")) {
      gi_true <- calculate_influence(gam_influence(model, focus = "year"), islog = TRUE)
      indices_true <- extract_indices(gi_true)

      expect_s3_class(indices_true, "data.frame")
      expect_true(all(c("index", "cv", "lower_CI", "upper_CI") %in% names(indices_true)))
      expect_equal(nrow(indices_true), 5) # 5 years
      expect_true(all(is.finite(indices_true$index)))

      # For lognormal, islog=TRUE should give different results
      if (family_name %in% c("lognormal_prelogged", "lognormal_gamma")) {
        expect_false(identical(indices_false$index, indices_true$index))
      }
    }

    message(paste("✓ extract_indices works with", family_name, "family"))
  }
})

test_that("get_terms identifies all term types correctly", {
  test_data <- create_test_data(n = 200)

  # Create model with various term types
  model <- gam(
    y_gaussian ~ s(depth) + s(temp) + te(depth, temp) + year + area +
      s(vessel, bs = "re") + s(vessel, by = area, bs = "re"),
    data = test_data,
    family = gaussian()
  )

  gi <- gam_influence(model, focus = "year")
  terms <- get_terms(gi)

  expect_true(is.character(terms))
  expect_true(length(terms) > 0)

  # Should identify various term types
  expect_true(any(grepl("s\\(depth\\)", terms)))
  expect_true(any(grepl("s\\(temp\\)", terms)))
  expect_true(any(grepl("te\\(", terms)))
  expect_true(any(grepl("s\\(vessel", terms)))

  message("✓ get_terms identifies all term types correctly")
})

test_that("geometric_mean function works correctly", {
  # Test with positive values
  x <- c(1, 2, 4, 8)
  expected <- exp(mean(log(x)))
  expect_equal(geometric_mean(x), expected)

  # Test with single value
  expect_equal(geometric_mean(5), 5)

  # Test with identical values
  expect_equal(geometric_mean(c(3, 3, 3, 3)), 3)

  # Test error handling with non-positive values
  expect_error(geometric_mean(c(1, 2, 0, 4)))
  expect_error(geometric_mean(c(1, 2, -1, 4)))

  message("✓ geometric_mean function works correctly")
})

test_that("r2 function works with gam_influence objects", {
  test_data <- create_test_data(n = 200)
  model <- gam(y_gaussian ~ s(depth) + s(temp) + year, data = test_data, family = gaussian())

  gi <- gam_influence(model, focus = "year")
  gi <- calculate_influence(gi)

  r2_result <- r2(gi)

  expect_true(is.numeric(r2_result))
  expect_true(length(r2_result) > 1) # Should return progression of R-squared values
  expect_true(all(r2_result >= 0 & r2_result <= 1)) # R-squared should be between 0 and 1
  expect_true(all(diff(r2_result) >= -1e-10)) # R-squared should be non-decreasing (allowing for small numerical errors)

  message("✓ r2 function works with gam_influence objects")
})

test_that("summary methods work with all families", {
  test_data <- create_test_data(n = 200)
  models <- fit_test_models(test_data)

  for (family_name in names(models)) {
    model <- models[[family_name]]

    gi <- gam_influence(model, focus = "year")
    gi <- calculate_influence(gi)

    # Test summary method
    expect_no_error(summary(gi))
    summary_output <- capture.output(summary(gi))
    expect_true(length(summary_output) > 0)

    message(paste("✓ summary method works with", family_name, "family"))
  }
})

test_that("print methods work with all families", {
  test_data <- create_test_data(n = 200)
  models <- fit_test_models(test_data)

  for (family_name in names(models)) {
    model <- models[[family_name]]

    gi <- gam_influence(model, focus = "year")
    gi <- calculate_influence(gi)

    # Test print method
    expect_no_error(print(gi))
    print_output <- capture.output(print(gi))
    expect_true(length(print_output) > 0)

    message(paste("✓ print method works with", family_name, "family"))
  }
})

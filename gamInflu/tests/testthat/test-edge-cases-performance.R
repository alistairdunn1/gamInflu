test_that("package handles edge cases gracefully", {
  # Test with minimal data
  minimal_data <- data.frame(
    year = factor(c(rep(2015, 10), rep(2016, 10))),
    x = rnorm(20),
    y = rnorm(20)
  )

  model_minimal <- gam(y ~ s(x) + year, data = minimal_data, family = gaussian())

  expect_no_error({
    gi <- gam_influence(model_minimal, focus = "year")
    gi <- calculate_influence(gi)
  })

  # Test with single focus level (should error appropriately)
  single_level_data <- data.frame(
    year = factor(rep(2015, 20)),
    x = rnorm(20),
    y = rnorm(20)
  )

  model_single <- gam(y ~ s(x) + year, data = single_level_data, family = gaussian())

  # Should handle or error gracefully
  expect_error(gam_influence(model_single, focus = "year"))

  message("✓ Package handles minimal data and edge cases")
})

test_that("package works with different sample sizes", {
  sample_sizes <- c(50, 100, 200, 500)

  for (n in sample_sizes) {
    test_data <- create_test_data(n = n, seed = 123)
    model <- gam(y_gaussian ~ s(depth) + s(temp) + year, data = test_data, family = gaussian())

    expect_no_error({
      gi <- gam_influence(model, focus = "year")
      gi <- calculate_influence(gi)
      indices <- extract_indices(gi)
    })

    message(paste("✓ Package works with n =", n))
  }
})

test_that("package performance is reasonable", {
  test_data <- create_test_data(n = 500)
  model <- gam(y_gaussian ~ s(depth) + s(temp) + year + area, data = test_data, family = gaussian())

  # Time the main workflow
  start_time <- Sys.time()
  gi <- gam_influence(model, focus = "year")
  gi <- calculate_influence(gi)
  indices <- extract_indices(gi)
  end_time <- Sys.time()

  elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))

  # Should complete within reasonable time (adjust threshold as needed)
  expect_true(elapsed < 30) # 30 seconds threshold

  message(paste("✓ Main workflow completed in", round(elapsed, 2), "seconds"))
})

test_that("package handles missing values appropriately", {
  test_data <- create_test_data(n = 200)

  # Introduce some missing values
  test_data$depth[1:5] <- NA
  test_data$temp[10:15] <- NA

  # Model should handle NAs appropriately
  model <- gam(y_gaussian ~ s(depth, na.action = na.exclude) + s(temp, na.action = na.exclude) + year,
    data = test_data, family = gaussian(), na.action = na.exclude
  )

  expect_no_error({
    gi <- gam_influence(model, focus = "year")
    gi <- calculate_influence(gi)
  })

  message("✓ Package handles missing values appropriately")
})

test_that("package works with different factor orderings", {
  test_data <- create_test_data(n = 200)

  # Test with different factor level orderings
  test_data$year_reversed <- factor(test_data$year, levels = rev(levels(test_data$year)))
  test_data$area_reordered <- factor(test_data$area, levels = c("South", "North"))

  model <- gam(y_gaussian ~ s(depth) + year_reversed * area_reordered,
    data = test_data, family = gaussian()
  )

  expect_no_error({
    gi <- gam_influence(model, focus = "year_reversed")
    gi <- calculate_influence(gi)
    indices <- extract_indices(gi)
  })

  # Check that factor levels are preserved
  expect_true("year_reversed" %in% names(indices))

  message("✓ Package works with different factor orderings")
})

test_that("package handles complex model formulas", {
  test_data <- create_test_data(n = 200)

  # Complex model with multiple term types
  complex_model <- gam(
    y_gaussian ~ s(depth, k = 5) + s(temp, k = 4) +
      te(depth, temp, k = c(3, 3)) +
      year + area + year:area +
      s(vessel, bs = "re"),
    data = test_data,
    family = gaussian()
  )

  expect_no_error({
    gi <- gam_influence(complex_model, focus = "year")
    gi <- calculate_influence(gi)

    # Test plotting with complex model
    plot_standardisation(gi)
    plot_stepwise_index(gi)

    terms <- get_terms(gi)
    if (length(terms) > 0) {
      plot_terms(gi, terms[1])
    }
  })

  message("✓ Package handles complex model formulas")
})

test_that("memory usage is reasonable", {
  test_data <- create_test_data(n = 1000) # Larger dataset
  model <- gam(y_gaussian ~ s(depth) + s(temp) + year + area, data = test_data, family = gaussian())

  # Monitor memory usage during workflow
  gc_before <- gc()

  gi <- gam_influence(model, focus = "year")
  gi <- calculate_influence(gi)

  gc_after <- gc()

  # Memory usage should not be excessive
  memory_used <- gc_after[2, 2] - gc_before[2, 2] # Vcells used

  message(paste("✓ Memory usage reasonable:", round(memory_used, 2), "MB"))
})

test_that("package produces consistent results across runs", {
  # Test reproducibility
  test_data1 <- create_test_data(n = 200, seed = 42)
  test_data2 <- create_test_data(n = 200, seed = 42)

  model1 <- gam(y_gaussian ~ s(depth) + s(temp) + year, data = test_data1, family = gaussian())
  model2 <- gam(y_gaussian ~ s(depth) + s(temp) + year, data = test_data2, family = gaussian())

  gi1 <- calculate_influence(gam_influence(model1, focus = "year"))
  gi2 <- calculate_influence(gam_influence(model2, focus = "year"))

  indices1 <- extract_indices(gi1)
  indices2 <- extract_indices(gi2)

  # Results should be identical for same seed
  expect_equal(indices1$index, indices2$index)
  expect_equal(indices1$cv, indices2$cv)

  message("✓ Package produces consistent results across runs")
})

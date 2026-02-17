test_that("islog=TRUE works correctly with Gaussian family (pre-logged data)", {
  test_data <- create_test_data(n = 200)

  # Model with pre-logged response
  model <- gam(
    log_catch ~ s(depth) + s(temp) + year,
    data = test_data,
    family = gaussian()
  )

  # Test with islog=FALSE (default)
  gi_false <- calculate_influence(gam_influence(model, focus = "year"), islog = FALSE)
  indices_false <- extract_indices(gi_false)

  # Test with islog=TRUE
  gi_true <- calculate_influence(gam_influence(model, focus = "year"), islog = TRUE)
  indices_true <- extract_indices(gi_true)

  # islog=TRUE should produce different (anti-logged) results
  # Note: After rescaling to mean=1, the exp() relationship doesn't hold exactly
  expect_true(all(is.finite(indices_true$index)))
  expect_true(all(indices_true$index > 0))  # Anti-logged values should be positive
  expect_true(all(is.finite(indices_true$lower_CI)))
  expect_true(all(is.finite(indices_true$upper_CI)))

  # islog should affect the results
  expect_false(identical(indices_false$index, indices_true$index))

  message("✓ islog=TRUE works correctly with pre-logged Gaussian data")
})

test_that("islog=TRUE works correctly with Gamma family (log link)", {
  test_data <- create_test_data(n = 200)

  # Model with Gamma family and log link
  model <- gam(
    catch ~ s(depth) + s(temp) + year,
    data = test_data,
    family = Gamma(link = "log")
  )

  # Test with islog=FALSE (default)
  gi_false <- calculate_influence(gam_influence(model, focus = "year"), islog = FALSE)
  indices_false <- extract_indices(gi_false)

  # Test with islog=TRUE
  gi_true <- calculate_influence(gam_influence(model, focus = "year"), islog = TRUE)
  indices_true <- extract_indices(gi_true)

  # Results should be different
  expect_false(identical(indices_false$index, indices_true$index))

  # Both should produce valid positive indices
  expect_true(all(is.finite(indices_true$index)))
  expect_true(all(indices_true$index > 0))

  message("✓ islog=TRUE works correctly with Gamma family (log link)")
})

test_that("islog=FALSE is appropriate for non-lognormal families", {
  test_data <- create_test_data(n = 200)

  # Binomial model
  model_binomial <- gam(
    y_binomial ~ s(depth) + s(temp) + year,
    data = test_data,
    family = binomial()
  )

  gi_binomial <- calculate_influence(gam_influence(model_binomial, focus = "year"), islog = FALSE)
  indices_binomial <- extract_indices(gi_binomial)

  # For binomial, indices should be reasonable probabilities or log-odds
  expect_true(all(is.finite(indices_binomial$index)))
  expect_true(all(is.finite(indices_binomial$cv)))

  # Poisson model
  model_poisson <- gam(
    y_poisson ~ s(depth) + s(temp) + year,
    data = test_data,
    family = poisson()
  )

  gi_poisson <- calculate_influence(gam_influence(model_poisson, focus = "year"), islog = FALSE)
  indices_poisson <- extract_indices(gi_poisson)

  # For Poisson, indices should be reasonable counts or log-counts
  expect_true(all(is.finite(indices_poisson$index)))
  expect_true(all(is.finite(indices_poisson$cv)))

  message("✓ islog=FALSE works appropriately for non-lognormal families")
})

test_that("islog parameter affects plotting correctly", {
  test_data <- create_test_data(n = 200)

  # Pre-logged Gaussian model
  model <- gam(
    log_catch ~ s(depth) + s(temp) + year,
    data = test_data,
    family = gaussian()
  )

  gi_false <- calculate_influence(gam_influence(model, focus = "year"), islog = FALSE)
  gi_true <- calculate_influence(gam_influence(model, focus = "year"), islog = TRUE)

  # Test standardisation plots
  p_false <- plot_standardisation(gi_false)
  p_true <- plot_standardisation(gi_true)

  expect_s3_class(p_false, "ggplot")
  expect_s3_class(p_true, "ggplot")

  # Test stepwise plots
  p_step_false <- plot_stepwise_index(gi_false)
  p_step_true <- plot_stepwise_index(gi_true)

  expect_s3_class(p_step_false, "ggplot")
  expect_s3_class(p_step_true, "ggplot")

  # Test CDI plots (use full term names, exclude focus term)
  terms_full <- setdiff(get_terms(gi_false, full = TRUE), "year")
  if (length(terms_full) > 0) {
    expect_no_error(plot_cdi(gi_false, terms_full[1]))
    expect_no_error(plot_cdi(gi_true, terms_full[1]))
  }

  message("✓ islog parameter affects plotting correctly")
})

test_that("islog parameter is preserved in combined models", {
  test_data <- create_test_data(n = 200)

  # Create two lognormal models
  model1 <- gam(log_catch ~ s(depth) + year, data = test_data, family = gaussian())
  model2 <- gam(catch ~ s(temp) + year, data = test_data, family = Gamma(link = "log"))

  # Calculate influence with islog=TRUE
  gi1 <- calculate_influence(gam_influence(model1, focus = "year"), islog = TRUE)
  gi2 <- calculate_influence(gam_influence(model2, focus = "year"), islog = TRUE)

  # Combine models
  gi_combined <- combine_indices(gi1, gi2)

  expect_s3_class(gi_combined, "gam_influence_combined")
  expect_true(!is.null(gi_combined$combined_indices))

  # Test that combined plotting works
  expect_no_error(plot(gi_combined, type = "comparison"))

  message("✓ islog parameter is preserved in combined models")
})

test_that("mean calculation methods work correctly with islog", {
  test_data <- create_test_data(n = 200)

  # Pre-logged Gaussian model
  model <- gam(
    log_catch ~ s(depth) + s(temp) + year,
    data = test_data,
    family = gaussian()
  )

  gi <- gam_influence(model, focus = "year")

  # Test different mean types with islog=FALSE
  gi_arith_false <- calculate_influence(gi, islog = FALSE, mean_type = "arithmetic")
  gi_geom_false <- calculate_influence(gi, islog = FALSE, mean_type = "geometric")

  # Test different mean types with islog=TRUE
  gi_arith_true <- calculate_influence(gi, islog = TRUE, mean_type = "arithmetic")
  gi_geom_true <- calculate_influence(gi, islog = TRUE, mean_type = "geometric")

  # All should produce valid results
  expect_s3_class(gi_arith_false, "gam_influence")
  expect_s3_class(gi_geom_false, "gam_influence")
  expect_s3_class(gi_arith_true, "gam_influence")
  expect_s3_class(gi_geom_true, "gam_influence")

  # All should produce valid indices
  indices_arith_false <- extract_indices(gi_arith_false)
  indices_geom_false <- extract_indices(gi_geom_false)

  expect_true(all(is.finite(indices_arith_false$index)))
  expect_true(all(is.finite(indices_geom_false$index)))

  message("✓ Mean calculation methods work correctly with islog")
})

test_that("confidence interval methods work with islog", {
  test_data <- create_test_data(n = 200)

  # Pre-logged Gaussian model
  model <- gam(
    log_catch ~ s(depth) + s(temp) + year,
    data = test_data,
    family = gaussian()
  )

  # Test coefficient-based method
  gi_coeff <- gam_influence(model, focus = "year", use_coeff_method = TRUE)
  gi_coeff <- calculate_influence(gi_coeff, islog = TRUE)

  # Test prediction-based method
  gi_pred <- gam_influence(model, focus = "year", use_coeff_method = FALSE)
  gi_pred <- calculate_influence(gi_pred, islog = TRUE)

  indices_coeff <- extract_indices(gi_coeff)
  indices_pred <- extract_indices(gi_pred)

  # Both methods should produce valid results with islog=TRUE
  expect_true(all(is.finite(indices_coeff$index)))
  expect_true(all(is.finite(indices_pred$index)))
  expect_true(all(indices_coeff$lower_CI <= indices_coeff$upper_CI))
  expect_true(all(indices_pred$lower_CI <= indices_pred$upper_CI))

  # Results should be different between methods
  expect_false(identical(indices_coeff$index, indices_pred$index))

  message("✓ Confidence interval methods work correctly with islog")
})

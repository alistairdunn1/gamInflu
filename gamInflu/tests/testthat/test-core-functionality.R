test_that("gam_influence works with all GLM families", {
  # Setup test data
  test_data <- create_test_data(n = 200)
  models <- fit_test_models(test_data)

  # Test each family
  for (family_name in names(models)) {
    model <- models[[family_name]]

    # Test gam_influence object creation - pass data explicitly
    gi <- gam_influence(model, focus = "year", data = test_data)

    expect_s3_class(gi, "gam_influence")
    expect_equal(gi$focus_var, "year")
    expect_true(is.list(gi$model_data))
    expect_true(!is.null(gi$original_model))

    # Check that family is correctly detected
    expect_true(!is.null(gi$family))

    message(paste("✓ gam_influence works with", family_name, "family"))
  }
})

test_that("calculate_influence works with all GLM families and islog parameter", {
  test_data <- create_test_data(n = 200)
  models <- fit_test_models(test_data)

  for (family_name in names(models)) {
    model <- models[[family_name]]

    # Test with islog = FALSE (default) - pass data explicitly
    gi_false <- gam_influence(model, focus = "year", data = test_data)
    gi_false <- calculate_influence(gi_false, islog = FALSE)

    expect_s3_class(gi_false, "gam_influence")
    expect_true(!is.null(gi_false$indices))
    expect_true(all(c("index", "cv", "lower_CI", "upper_CI") %in% names(gi_false$indices)))

    # Test with islog = TRUE for applicable families
    if (family_name %in% c("lognormal_prelogged", "lognormal_gamma", "gamma", "poisson")) {
      if (family_name != "poisson") { # Skip Poisson for islog=TRUE as it doesn't make sense
        gi_true <- gam_influence(model, focus = "year", data = test_data)
        gi_true <- calculate_influence(gi_true, islog = TRUE)

        expect_s3_class(gi_true, "gam_influence")
        expect_true(!is.null(gi_true$indices))
        expect_true(all(c("index", "cv", "lower_CI", "upper_CI") %in% names(gi_true$indices)))

        # For lognormal models, check that islog=TRUE produces different (anti-logged) results
        if (family_name %in% c("lognormal_prelogged", "lognormal_gamma")) {
          expect_false(identical(gi_false$indices$index, gi_true$indices$index))
          # islog=TRUE results should generally be larger (anti-logged)
          expect_true(mean(gi_true$indices$index) > mean(gi_false$indices$index))
        }
      }
    }

    message(paste("✓ calculate_influence works with", family_name, "family"))
  }
})

test_that("lognormal handling works correctly", {
  test_data <- create_test_data(n = 200)

  # Test pre-logged response with Gaussian family
  model_prelogged <- gam(
    log_catch ~ s(depth) + s(temp) + year,
    data = test_data,
    family = gaussian()
  )

  gi_prelogged_false <- calculate_influence(gam_influence(model_prelogged, focus = "year"), islog = FALSE)
  gi_prelogged_true <- calculate_influence(gam_influence(model_prelogged, focus = "year"), islog = TRUE)

  # Test raw response with Gamma family (log link)
  model_gamma <- gam(
    catch ~ s(depth) + s(temp) + year,
    data = test_data,
    family = Gamma(link = "log")
  )

  gi_gamma_false <- calculate_influence(gam_influence(model_gamma, focus = "year"), islog = FALSE)
  gi_gamma_true <- calculate_influence(gam_influence(model_gamma, focus = "year"), islog = TRUE)

  # Check that islog=TRUE transforms the results appropriately
  expect_true(all(gi_prelogged_true$indices$index > gi_prelogged_false$indices$index))
  expect_true(all(gi_gamma_true$indices$index > gi_gamma_false$indices$index))

  # Check that the transformation is exponential
  expect_equal(gi_prelogged_true$indices$index, exp(gi_prelogged_false$indices$index), tolerance = 1e-10)

  message("✓ Lognormal handling with islog parameter works correctly")
})

test_that("family detection works correctly", {
  test_data <- create_test_data(n = 200)
  models <- fit_test_models(test_data)

  expected_families <- c(
    "gaussian" = "gaussian",
    "lognormal_prelogged" = "gaussian",
    "lognormal_gamma" = "Gamma",
    "binomial" = "binomial",
    "gamma" = "Gamma",
    "poisson" = "poisson"
  )

  for (family_name in names(expected_families)) {
    if (family_name %in% names(models)) {
      model <- models[[family_name]]
      gi <- gam_influence(model, focus = "year")

      expect_equal(gi$family$family, expected_families[family_name])
      message(paste("✓ Family detection correct for", family_name))
    }
  }
})

test_that("interaction terms work with different families", {
  test_data <- create_test_data(n = 200)
  interaction_models <- fit_interaction_models(test_data)

  for (family_name in names(interaction_models)) {
    model <- interaction_models[[family_name]]

    # Test with interaction term as focus
    gi <- gam_influence(model, focus = "year")
    gi <- calculate_influence(gi)

    expect_s3_class(gi, "gam_influence")
    expect_true(!is.null(gi$indices))

    message(paste("✓ Interaction terms work with", family_name))
  }
})

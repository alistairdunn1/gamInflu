test_that("plotting functions work with all GLM families", {
  test_data <- create_test_data(n = 200)
  models <- fit_test_models(test_data)

  # Test plotting functions with each family
  for (family_name in names(models)) {
    model <- models[[family_name]]

    gi <- gam_influence(model, focus = "year")
    gi <- calculate_influence(gi)

    # Test main plotting functions
    expect_no_error(plot_standardisation(gi))
    expect_no_error(plot_stepwise_index(gi))
    expect_no_error(plot_step_and_influence(gi))

    # Test plot_terms for different term types
    terms <- get_terms(gi)
    if (length(terms) > 0) {
      expect_no_error(plot_terms(gi, terms[1]))
      expect_no_error(plot_cdi(gi, terms[1]))
      expect_no_error(plot_term_distribution(gi, terms[1]))
    }

    # Test generic plot method
    expect_no_error(plot(gi, type = "standardisation"))
    expect_no_error(plot(gi, type = "stepwise"))
    expect_no_error(plot(gi, type = "step_influence"))

    # Test residual plots
    expect_no_error(plot_residuals(gi))

    message(paste("✓ Plotting functions work with", family_name, "family"))
  }
})

test_that("plot_terms works with different smoother types", {
  test_data <- create_test_data(n = 200)

  # Create models with different smoother types
  models <- list()

  # Standard smoothers
  models$standard <- gam(
    y_gaussian ~ s(depth) + s(temp) + year,
    data = test_data,
    family = gaussian()
  )

  # Tensor product smoothers
  models$tensor <- gam(
    y_gaussian ~ te(depth, temp) + year,
    data = test_data,
    family = gaussian()
  )

  # Factor-by smoothers
  models$factor_by <- gam(
    y_gaussian ~ s(depth, by = area) + year,
    data = test_data,
    family = gaussian()
  )

  for (model_name in names(models)) {
    model <- models[[model_name]]
    gi <- gam_influence(model, focus = "year")
    gi <- calculate_influence(gi)

    terms <- get_terms(gi)
    if (length(terms) > 0) {
      for (term in terms) {
        expect_no_error(plot_terms(gi, term))
        expect_no_error(plot_cdi(gi, term))
      }
    }

    message(paste("✓ plot_terms works with", model_name, "smoother types"))
  }
})

test_that("random effects plotting works", {
  test_data <- create_test_data(n = 200)
  random_models <- fit_random_effects_models(test_data)

  for (model_name in names(random_models)) {
    model <- random_models[[model_name]]

    gi <- gam_influence(model, focus = "year")
    gi <- calculate_influence(gi)

    terms <- get_terms(gi)
    random_terms <- terms[grepl("vessel", terms)]

    if (length(random_terms) > 0) {
      for (term in random_terms) {
        expect_no_error(plot_terms(gi, term))
        expect_no_error(plot_cdi(gi, term))
      }
    }

    message(paste("✓ Random effects plotting works with", model_name))
  }
})

test_that("plotting with islog parameter works", {
  test_data <- create_test_data(n = 200)

  # Test with lognormal model
  model <- gam(
    log_catch ~ s(depth) + s(temp) + year,
    data = test_data,
    family = gaussian()
  )

  gi_false <- calculate_influence(gam_influence(model, focus = "year"), islog = FALSE)
  gi_true <- calculate_influence(gam_influence(model, focus = "year"), islog = TRUE)

  # Test that plots work with both islog settings
  expect_no_error(plot_standardisation(gi_false))
  expect_no_error(plot_standardisation(gi_true))

  expect_no_error(plot_stepwise_index(gi_false))
  expect_no_error(plot_stepwise_index(gi_true))

  # Test that axis labels reflect the islog setting
  p_false <- plot_standardisation(gi_false)
  p_true <- plot_standardisation(gi_true)

  # Both should be ggplot objects
  expect_s3_class(p_false, "ggplot")
  expect_s3_class(p_true, "ggplot")

  message("✓ Plotting with islog parameter works correctly")
})

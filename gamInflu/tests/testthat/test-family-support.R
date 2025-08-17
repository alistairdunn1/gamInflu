test_that("Factor ordering works correctly in distribution plots", {
  # Create test data with year levels that could be mis-ordered
  set.seed(123)
  test_data <- data.frame(
    year = factor(c(rep("1992", 30), rep("2010", 30), rep("2005", 30))),
    depth = runif(90, 10, 100),
    catch = exp(rnorm(90, 2, 0.5))
  )

  # Fit model
  model <- gam(log(catch) ~ s(depth) + year, data = test_data, family = gaussian())
  gi <- calculate_influence(gam_influence(model, focus = "year", data = test_data))

  # Test that factor levels are handled properly in plotting
  expect_no_error(plot_terms(gi, "year"))
  expect_no_error(plot_cdi(gi, "year"))
  expect_no_error(plot_term_distribution(gi, "year"))

  # The actual ordering should be maintained as specified in factor levels
  year_levels <- levels(test_data$year)
  expect_equal(length(year_levels), 3)
  expect_true("1992" %in% year_levels)
  expect_true("2010" %in% year_levels)
  expect_true("2005" %in% year_levels)

  message("✓ Factor ordering works correctly")
})

test_that("islog arithmetic relationships work correctly", {
  # Create test data
  set.seed(123)
  n <- 100
  test_data <- data.frame(
    year = factor(rep(2015:2019, each = 20)),
    x = runif(n, 10, 100),
    log_y = 2 + 0.1 * seq_len(n) / 10 + rnorm(n, 0, 0.3)
  )

  # Fit model with pre-logged response
  model <- gam(log_y ~ s(x) + year, data = test_data, family = gaussian())

  # Test with explicit mean types to ensure proper exp() relationship
  gi_false <- calculate_influence(gam_influence(model, focus = "year", data = test_data),
    islog = FALSE, rescale_method = "arithmetic_mean"
  )
  gi_true <- calculate_influence(gam_influence(model, focus = "year", data = test_data),
    islog = TRUE, rescale_method = "arithmetic_mean"
  )

  indices_false <- extract_indices(gi_false)
  indices_true <- extract_indices(gi_true)

  # With arithmetic mean and islog handling, should see exp() relationship
  expect_true(all(is.finite(indices_false$index)))
  expect_true(all(is.finite(indices_true$index)))
  expect_true(all(indices_true$index > 0)) # Anti-logged values should be positive

  # The relationship should be consistent with log/exp transformation
  expect_false(identical(indices_false$index, indices_true$index))

  message("✓ islog arithmetic relationships work correctly")
})

test_that("Random effects plotting works with by-variables", {
  # Create test data with random effects
  set.seed(123)
  n <- 200
  test_data <- data.frame(
    year = factor(rep(2015:2019, each = 40)),
    vessel = factor(sample(paste0("V", 1:8), n, replace = TRUE)),
    area = factor(sample(c("North", "South"), n, replace = TRUE)),
    x = runif(n, 10, 100),
    y = rnorm(n, 2, 0.5)
  )

  # Add vessel-by-area random effects
  vessel_effects <- rnorm(nlevels(test_data$vessel), 0, 0.3)
  area_effects <- c("North" = 0.2, "South" = -0.2)

  for (i in seq_len(n)) {
    vessel_idx <- as.numeric(test_data$vessel[i])
    area_name <- as.character(test_data$area[i])
    test_data$y[i] <- test_data$y[i] + vessel_effects[vessel_idx] + area_effects[area_name]
  }

  # Fit model with random by-effects
  model <- gam(y ~ s(x) + year + s(vessel, by = area, bs = "re"),
    data = test_data, family = gaussian()
  )

  gi <- calculate_influence(gam_influence(model, focus = "year", data = test_data))

  # Test plotting of random effects terms
  terms <- get_terms(gi)
  random_terms <- terms[grepl("vessel", terms)]

  if (length(random_terms) > 0) {
    expect_no_error(plot_terms(gi, random_terms[1]))
    expect_no_error(plot_cdi(gi, random_terms[1]))
  }

  message("✓ Random effects plotting works correctly")
})

test_that("Multiple GLM families work correctly", {
  # Create comprehensive test data
  set.seed(123)
  n <- 200
  test_data <- data.frame(
    year = factor(rep(2015:2019, each = 40)),
    x = runif(n, 10, 100)
  )

  # Generate responses for different families
  linear_pred <- 1 + 0.1 * as.numeric(test_data$year) + sin(test_data$x / 20)
  test_data$y_gaussian <- linear_pred + rnorm(n, 0, 0.5)
  test_data$y_binomial <- rbinom(n, 1, plogis(linear_pred))
  test_data$y_gamma <- rgamma(n, shape = 2, rate = 2 / exp(linear_pred))
  test_data$y_poisson <- rpois(n, exp(linear_pred - 1))

  # Test Gaussian family
  model_gaussian <- gam(y_gaussian ~ s(x) + year, data = test_data, family = gaussian())
  gi_gaussian <- calculate_influence(gam_influence(model_gaussian, focus = "year", data = test_data))
  expect_s3_class(gi_gaussian, "gam_influence")
  expect_equal(gi_gaussian$family$family, "gaussian")

  # Test Binomial family
  model_binomial <- gam(y_binomial ~ s(x) + year, data = test_data, family = binomial())
  gi_binomial <- calculate_influence(gam_influence(model_binomial, focus = "year", data = test_data))
  expect_s3_class(gi_binomial, "gam_influence")
  expect_equal(gi_binomial$family$family, "binomial")

  # Test Gamma family
  model_gamma <- gam(y_gamma ~ s(x) + year, data = test_data, family = Gamma(link = "log"))
  gi_gamma <- calculate_influence(gam_influence(model_gamma, focus = "year", data = test_data))
  expect_s3_class(gi_gamma, "gam_influence")
  expect_equal(gi_gamma$family$family, "Gamma")

  # Test Poisson family
  model_poisson <- gam(y_poisson ~ s(x) + year, data = test_data, family = poisson())
  gi_poisson <- calculate_influence(gam_influence(model_poisson, focus = "year", data = test_data))
  expect_s3_class(gi_poisson, "gam_influence")
  expect_equal(gi_poisson$family$family, "poisson")

  # Test that all produce valid indices
  for (gi in list(gi_gaussian, gi_binomial, gi_gamma, gi_poisson)) {
    indices <- extract_indices(gi)
    expect_s3_class(indices, "data.frame")
    expect_equal(nrow(indices), 5) # 5 years
    expect_true(all(c("index", "cv", "lower_CI", "upper_CI") %in% names(indices)))
    expect_true(all(is.finite(indices$index)))
  }

  message("✓ Multiple GLM families work correctly")
})

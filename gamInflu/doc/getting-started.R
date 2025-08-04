## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  warning = FALSE,
  message = FALSE
)

## ----setup--------------------------------------------------------------------
library(gamInflu)
library(mgcv)

## ----eval=FALSE---------------------------------------------------------------
# # Convert focus variable to factor BEFORE fitting the model
# data$year <- factor(data$year)
# data$area <- factor(data$area)

## ----eval=FALSE---------------------------------------------------------------
# # Example: Fit a GAM for fisheries data
# model <- gam(catch ~ s(depth) + s(temperature) + year + area,
#              data = data,
#              family = Gamma(link = "log"))

## ----eval=FALSE---------------------------------------------------------------
# # Create influence analysis object
# gi <- gam_influence(model, focus = "year")

## ----eval=FALSE---------------------------------------------------------------
# # Calculate all influence metrics (auto-detects family)
# gi <- calculate_influence(gi)

## ----eval=FALSE---------------------------------------------------------------
# # Get standardised indices
# indices <- extract_indices(gi)
# print(indices)
# 
# # Export results
# write.csv(indices, "standardised_indices.csv", row.names = FALSE)

## ----eval=FALSE---------------------------------------------------------------
# # Main plots
# plot_standardisation(gi)    # Compare indices
# plot_stepwise_index(gi)     # Model progression
# plot_term_influence(gi)     # Term influences
# plot_residuals(gi, type = "violin")  # Residual diagnostics

## -----------------------------------------------------------------------------
# Set seed for reproducibility
set.seed(123)

# Create simulated fisheries data
n <- 200
data <- data.frame(
  year = factor(rep(2015:2019, each = 40)),
  depth = runif(n, 10, 100),
  temperature = runif(n, 5, 25),
  area = factor(rep(c("North", "South"), length.out = n))
)

# Add year effects (some years have higher catch)
year_effects <- c("2015" = 0.8, "2016" = 1.2, "2017" = 1.0, "2018" = 0.9, "2019" = 1.1)
data$year_effect <- year_effects[as.character(data$year)]

# Simulate catch with depth and temperature effects
data$log_catch <- 2 + 
  log(data$year_effect) +
  0.02 * (data$depth - 50) +
  0.05 * (data$temperature - 15) +
  rnorm(n, 0, 0.3)

data$catch <- exp(data$log_catch)

## -----------------------------------------------------------------------------
# Fit GAM with gamma family (for positive continuous data)
model <- gam(catch ~ s(depth) + s(temperature) + year + area, 
             data = data, 
             family = Gamma(link = "log"))

# Check model summary
summary(model)

## -----------------------------------------------------------------------------
# Create influence object
gi <- gam_influence(model, focus = "year")

# Calculate influence metrics
gi <- calculate_influence(gi)

# Extract indices
indices <- extract_indices(gi)
print(indices)

## -----------------------------------------------------------------------------
# Plot standardised indices
plot_standardisation(gi)

## -----------------------------------------------------------------------------
# Plot model progression
plot_stepwise_index(gi)

## -----------------------------------------------------------------------------
# Plot term influences
plot_term_influence(gi)

## -----------------------------------------------------------------------------
# Residual diagnostics with violin plots
plot_residuals(gi, type = "violin")

## -----------------------------------------------------------------------------
# Create presence/absence data
data$presence <- rbinom(n, 1, plogis(-1 + 0.02 * data$depth - 0.05 * data$temperature))

# Fit binomial model
model_binom <- gam(presence ~ s(depth) + s(temperature) + year, 
                   data = data, 
                   family = binomial())

# Influence analysis
gi_binom <- gam_influence(model_binom, focus = "year")
gi_binom <- calculate_influence(gi_binom)

# Plot results
plot_standardisation(gi_binom)

## -----------------------------------------------------------------------------
# Create count data
lambda <- exp(1 + 0.01 * data$depth + log(data$year_effect))
data$fish_count <- rpois(n, lambda)

# Fit Poisson model
model_pois <- gam(fish_count ~ s(depth) + year, 
                  data = data, 
                  family = poisson())

# Influence analysis
gi_pois <- gam_influence(model_pois, focus = "year")
gi_pois <- calculate_influence(gi_pois)

# Plot results
plot_standardisation(gi_pois)

## -----------------------------------------------------------------------------
# Coefficient-based (default, traditional)
gi_coeff <- gam_influence(model, focus = "year", use_coeff_method = TRUE)
gi_coeff <- calculate_influence(gi_coeff)

# Prediction-based (modern)
gi_pred <- gam_influence(model, focus = "year", use_coeff_method = FALSE)
gi_pred <- calculate_influence(gi_pred)

# Compare results
indices_coeff <- extract_indices(gi_coeff)
indices_pred <- extract_indices(gi_pred)

print("Coefficient-based CVs:")
print(indices_coeff$cv)
print("Prediction-based CVs:")
print(indices_pred$cv)

## -----------------------------------------------------------------------------
# Analyse residual patterns
residual_analysis <- analyse_residual_patterns(gi)
print(residual_analysis)

## ----eval=FALSE---------------------------------------------------------------
# # Export standardised indices
# indices <- extract_indices(gi)
# write.csv(indices, "standardised_indices.csv", row.names = FALSE)
# 
# # Save plots
# ggsave("standardisation_plot.png", plot_standardisation(gi),
#        width = 10, height = 6, dpi = 300)
# ggsave("stepwise_plot.png", plot_stepwise_index(gi),
#        width = 10, height = 6, dpi = 300)


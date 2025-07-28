# gamInflu Family Support Examples
# This example demonstrates the enhanced family support in gamInflu

library(devtools)
load_all("gamInflu")
library(mgcv)

## Example 1: Binomial Model (e.g., presence/absence data)
# Simulate presence/absence data
set.seed(123)
n <- 200
x1 <- runif(n, 0, 10) # Environmental variable
year <- factor(sample(2010:2020, n, replace = TRUE))
# Simulate presence probability declining over time
prob <- plogis(-1 + 0.2 * x1 - 0.1 * as.numeric(as.character(year)) + 2010)
presence <- rbinom(n, 1, prob)

dat_binom <- data.frame(
  presence = presence,
  env_var = x1,
  year = year
)

# Fit binomial GAM
mod_binom <- gam(presence ~ s(env_var) + year,
  data = dat_binom,
  family = binomial()
)

# Create influence object and calculate (family auto-detected)
inf_binom <- gam_influence(mod_binom, focus = "year")
result_binom <- calculate_influence(inf_binom)

cat("Binomial example completed successfully!\n")
cat("Standardized index range:", range(result_binom$standardized_index), "\n")

## Example 2: Gamma Model (e.g., biomass data)
# Simulate positive continuous data (e.g., fish biomass)
y_gamma <- rgamma(n, shape = 2, rate = exp(-2 - 0.3 * x1))
dat_gamma <- data.frame(
  biomass = y_gamma,
  env_var = x1,
  year = year
)

# Fit Gamma GAM with log link
mod_gamma <- gam(biomass ~ s(env_var) + year,
  data = dat_gamma,
  family = Gamma(link = "log")
)

# Calculate influence (family auto-detected)
inf_gamma <- gam_influence(mod_gamma, focus = "year")
result_gamma <- calculate_influence(inf_gamma)

cat("Gamma example completed successfully!\n")
cat("Standardized index range:", range(result_gamma$standardized_index), "\n")

## Example 3: Poisson Model (e.g., count data)
# Simulate count data
count_data <- rpois(n, exp(1 + 0.1 * x1))
dat_pois <- data.frame(
  count = count_data,
  env_var = x1,
  year = year
)

# Fit Poisson GAM
mod_pois <- gam(count ~ s(env_var) + year,
  data = dat_pois,
  family = poisson()
)

# Calculate influence
inf_pois <- gam_influence(mod_pois, focus = "year")
result_pois <- calculate_influence(inf_pois)

cat("Poisson example completed successfully!\n")
cat("Standardized index range:", range(result_pois$standardized_index), "\n")

## Summary
cat("\n=== gamInflu Family Support Summary ===\n")
cat("✓ Gaussian family: Traditional log-normal CPUE standardization\n")
cat("✓ Binomial family: Presence/absence or proportion data\n")
cat("✓ Gamma family: Positive continuous data (biomass, CPUE)\n")
cat("✓ Poisson family: Count data\n")
cat("✓ Automatic family detection from model object\n")
cat("✓ Family-specific index calculation methods\n")
cat("✓ Enhanced validation and error handling\n")

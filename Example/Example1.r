# Example: Using GAMInflu plotting functions with year, depth, month, and vessel

library(mgcv)
library(GAMInflu)

set.seed(123)

# Simulate data
n <- 500
years <- 2010:2020
months <- factor(month.abb, levels = month.abb)
vessels <- paste0("V", 1:10)

dat <- data.frame(
  year = sample(years, n, replace = TRUE),
  depth = runif(n, 10, 200),
  month = factor(sample(months, n, replace = TRUE), levels = months),
  vessel = factor(sample(vessels, n, replace = TRUE), levels = vessels)
)

# Simulate response with effects
dat$y <- 5 +
  0.02 * dat$depth +
  as.numeric(dat$month) * 0.1 +
  rnorm(n, sd = 0.5) +
  rnorm(length(vessels))[dat$vessel] +
  rnorm(length(years))[dat$year - min(years) + 1]

dat$year <- as.factor(dat$year)

# Fit a GAM with depth (continuous), month (factor), vessel (random effect), and year (focus)
m <- gam(y ~ s(depth) + s(month, bs = "re") + s(vessel, bs = "re") + year ,data = dat)

# 1. Plot GAM effects (main effects and random effects)
plots <- plot_gam_effects(
  model = m,
  data = dat,
  terms = c("year", "depth", "month", "vessel", "year"),
  show_rug = TRUE,
  show_ci = TRUE,
  plot_as_grid = TRUE,
  re_style = "panel", # Try "panel", "density", "qqnorm", "caterpillar", or "shrinkage"
  verbose = TRUE
)
print(plots)

# 2. Create influence_gam object
influ_obj <- create_influence_gam(gam_model = m, data = dat,focus = "year")

# 3. Calculate influence
influ_obj <- calculate_influence(influ_obj)

# 4. Plot influence_gam object with all plot types
# Standardization plot
plot(influ_obj, type = "stan")

# Step plot
plot(influ_obj, type = "step")

# Influence plot
plot(influ_obj, type = "influ")

# CDI plot (for a term, e.g., "depth")
plot(influ_obj, type = "cdi", term = "year")

# Step + Influence plot
plot(influ_obj, type = "step_influ")



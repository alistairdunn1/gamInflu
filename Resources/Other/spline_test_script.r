# Test script for updated influence package with ns smoother support
library(mgcv)
library(splines)
library(ggplot2)
library(dplyr)

# Create example data with year as a factor
set.seed(123)
n <- 500
data <- data.frame(
  year = factor(rep(2010:2019, each = n/10)),
  x = runif(n, 0, 10),
  z = rnorm(n),
  group = factor(sample(LETTERS[1:5], n, replace = TRUE)),
  effort = runif(n, 0.5, 10),
  catch = NA
)

# Generate catch data with different effects including natural splines
year_effect <- c(1.0, 1.2, 1.4, 1.3, 0.9, 0.8, 0.7, 0.8, 0.9, 1.0)
for (i in 1:n) {
  year_idx <- as.integer(data$year[i]) - 2009
  
  # Use natural spline effect for x
  x_effect <- 0.8 + 0.5 * sin(pi * data$x[i] / 5)
  
  # Base catch with various effects
  true_catch <- data$effort[i] * year_effect[year_idx] * x_effect * (1 + 0.2 * data$z[i])
  
  # Add noise
  data$catch[i] <- true_catch * exp(rnorm(1, 0, 0.2))
}

# Log-transform catch for modeling
data$log_catch <- log(data$catch)

# Model with natural splines (ns)
ns_model <- glm(log_catch ~ year + ns(x, df=4) + z + group + log(effort), data = data)
summary(ns_model)

# Create influence object and calculate statistics
infl_ns <- influence(ns_model, focus = "year")
print(infl_ns)  # Should show ns as a smoother type
infl_ns <- calc(infl_ns)

# Create plots to verify it works
stan_plot(infl_ns)
step_plot(infl_ns)
influ_plot(infl_ns)
cdi_plot(infl_ns, "ns(x, df=4)")

# Model with B-splines (bs)
bs_model <- glm(log_catch ~ year + bs(x, df=4) + z + group + log(effort), data = data)
summary(bs_model)

# Create influence object and calculate statistics for bs model
infl_bs <- influence(bs_model, focus = "year")
print(infl_bs)  # Should show bs as a smoother type
infl_bs <- calc(infl_bs)

# Create plots to verify bs support works
stan_plot(infl_bs)
step_plot(infl_bs)
influ_plot(infl_bs)
cdi_plot(infl_bs, "bs(x, df=4)")

# Model with both GAM smooth and natural splines to test compatibility
mixed_model <- gam(log_catch ~ s(year) + ns(x, df=4) + z + group + log(effort), data = data)
summary(mixed_model)

# Create influence object and calculate statistics for mixed model
infl_mixed <- influence(mixed_model, focus = "year")
print(infl_mixed)  # Should show both s and ns smoother types
infl_mixed <- calc(infl_mixed)

# Create plots to verify combined support works
stan_plot(infl_mixed)
step_plot(infl_mixed)
influ_plot(infl_mixed)

# Test polynomial smoothers
poly_model <- glm(log_catch ~ year + poly(x, 3) + z + group + log(effort), data = data)
infl_poly <- influence(poly_model, focus = "year")
print(infl_poly)  # Should show poly as a smoother type
infl_poly <- calc(infl_poly)
step_plot(infl_poly)
cdi_plot(infl_poly, "poly(x, 3)")

# Test multiple spline terms in the same model
multi_spline_model <- glm(log_catch ~ year + ns(x, df=4) + bs(z, df=3) + group + log(effort), data = data)
infl_multi <- influence(multi_spline_model, focus = "year")
print(infl_multi)  # Should show both ns and bs smoother types
infl_multi <- calc(infl_multi)
step_plot(infl_multi)

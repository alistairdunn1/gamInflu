# Install required packages if not already installed
# install.packages(c("mgcv", "ggplot2", "patchwork", "tidyr"))

# Load libraries
library(mgcv)
library(ggplot2)
library(patchwork)
library(tidyr) # Required for pivot_longer in plotting functions

# Source the provided R script
#source("Influ-3.r")

# Set seed for reproducibility
set.seed(123)

# Number of observations
n_obs <- 1000

# Simulate data
sim_data <- data.frame(
  year = factor(rep(1990:2019, length.out = n_obs)),
  temperature = rnorm(n_obs, mean = 15, sd = 5),
  region = factor(sample(c("North", "South", "East", "West"), n_obs, replace = TRUE)),
  # Log-transformed response variable
  y = NA
)

# Generate 'y' with some dependencies
# 'year' has a cyclical effect, 'temperature' has a smooth effect, 'region' is a factor effect
sim_data$y <- 5 +
  sin(as.numeric(as.character(sim_data$year)) / 5) * 2 + # Cyclical year effect
  0.1 * sim_data$temperature + # Linear temperature effect
  ifelse(sim_data$region == "North", 0.5,
         ifelse(sim_data$region == "South", -0.5,
                ifelse(sim_data$region == "East", 0.2, -0.2))) + # Region effect
  rnorm(n_obs, mean = 0, sd = 0.5) # Residual error

# Apply log transformation to y for this example, as the script handles log-transformed data
sim_data$log_y <- log(sim_data$y + abs(min(sim_data$y)) + 0.1) # Ensure positive for log

# Fit the GAM model
# s(temperature) for a smooth effect of temperature
# year as a factor (parametric term in this case, but can be smooth if year is numeric)
# region as a factor
gam_model <- gam(log_y ~ year + s(temperature, by = region) + region, data = sim_data)

# Print model summary
summary(gam_model)

# Create the gam_influence object
gi <- gam_influence(model = gam_model, focus = "year", data = sim_data)

# You can inspect the structure of the object
str(gi, max.level = 2)

# Perform calculations
gi <- calculate_influence(gi, islog = TRUE)

# Now, the 'calculated' list within the object should be populated
names(gi$calculated)

# Print the summary of influence analysis
summary(gi)

# look at the pred data


# Get the R-squared and deviance explained progression table
r2_summary <- r2(gi)
print(r2_summary)

# Standardization Plot
plot(gi, type = "stan")

# Stepwise Index Plot
plot(gi, type = "step")

# Combined Step and Influence Plots
plot(gi, type = "all")

# CDI Plot for 's(temperature)'
plot(gi, type = "cdi", term = 2)
plot(gi, type = "cdi", term = 3)

# CDI Plot for 'region' (a factor term)
plot(gi, type = "cdi", term = "region")

#
plot_terms(gi, 3)


plot_term_distribution(gi, 2)



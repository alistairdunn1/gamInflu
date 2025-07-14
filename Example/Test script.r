# Test script for gaminfluence package
# Load required libraries
library(mgcv)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(viridis)

# Source package files
# Load the package files
source("C:/Users/alist/OneDrive/Projects/Software/GAMInflu/GAMInflu/R/calculations.r")
source("C:/Users/alist/OneDrive/Projects/Software/GAMInflu/GAMInflu/R/data_preparation.R") 
source("C:/Users/alist/OneDrive/Projects/Software/GAMInflu/GAMInflu/R/gam_influence_class.R")
source("C:/Users/alist/OneDrive/Projects/Software/GAMInflu/GAMInflu/R/methods.R")
source("C:/Users/alist/OneDrive/Projects/Software/GAMInflu/GAMInflu/R/plot_gam_effects.R")
source("C:/Users/alist/OneDrive/Projects/Software/GAMInflu/GAMInflu/R/plotting.R")
source("C:/Users/alist/OneDrive/Projects/Software/GAMInflu/GAMInflu/R/smooth_parsing.R")
source("C:/Users/alist/OneDrive/Projects/Software/GAMInflu/GAMInflu/R/utils.R")

source("R/utils.R")
source("R/smooth-parsing.R") 
source("R/GAMInfluence-class.R")
source("R/methods.R")
source("R/gam_influence.R")
# Test with mtcars dataset
data(mtcars)

# Convert cyl to factor
mtcars$cyl <- factor(mtcars$cyl)

# Fit a simple GAM
model <- gam(mpg ~ s(wt) + s(hp) + cyl, data = mtcars)

print("Model fitted successfully")
print(summary(model))

# Create influence analysis object
influence <- gam_influence(model, focus = "cyl")

print("GAMInfluence object created successfully")

# Calculate influence metrics
influence$calc()

print("Calculations completed successfully")

# Test plots
print("Creating standardization plot...")
stan_plot <- influence$stan_plot()
print(stan_plot)

print("Creating step plot...")
step_plot <- influence$step_plot(panels = TRUE)
print(step_plot)

print("Creating influence plot...")
influence_plot <- influence$influence_plot()
print(influence_plot)

print("Getting smooth summary...")
smooth_summary <- influence$get_smooth_summary()
print(smooth_summary)

print("All tests completed successfully!")
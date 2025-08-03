#!/usr/bin/env Rscript

# Test the fixed factor level validation for subset analysis

library(mgcv)
set.seed(123)

# Create a dataset where Research_block='486_5' has only one level when subset
data <- data.frame(
  Year = factor(rep(c(2020, 2021, 2022), each = 100)),
  Research_block = factor(c(rep("486_5", 100), rep("487_1", 100), rep("488_2", 100))),
  Line_length = rnorm(300, 50, 10),
  cpue = rnorm(300, 10, 2)
)

# Check data structure
cat("Data structure:\n")
cat("Total observations:", nrow(data), "\n")
cat("Research_block levels in full data:", paste(levels(data$Research_block), collapse = ", "), "\n")
cat("Research_block='486_5' subset size:", sum(data$Research_block == "486_5"), "\n")
cat("Unique Research_block values in subset:", length(unique(data$Research_block[data$Research_block == "486_5"])), "\n\n")

# Fit model with multiple terms
mod <- gam(cpue ~ Year + Research_block + Line_length + Year:Research_block, data = data)

# Load the package functions
source("gamInflu/R/gam_influence.R")
source("gamInflu/R/calculate_influence.R")

# Test subset analysis with Research_block='486_5' (only one level in subset)
cat("Testing subset analysis...\n")
gi <- gam_influence(model = mod, focus = "Year", data = data)
gi_calc <- calculate_influence(gi, subset_var = "Research_block", subset_value = "486_5")

cat("Test completed successfully!\n")
cat("Summary entries:", nrow(gi_calc$calculated$summary), "\n")
print(gi_calc$calculated$summary$term)

# --- Example 4: Interaction term involving focus variable ---
# Load required packages
library(mgcv)     # For gam function and gamSim
library(ggplot2)  # For plotting
library(tidyr)    # For data reshaping
library(patchwork) # For combining plots
library(RColorBrewer) # For color palettes
library(here) # For file path management

DIR <- paste0(here(), "/GAMInflu/R/")

source(paste0(DIR, "calculate_influence.influence_gam.R"))
source(paste0(DIR, "create_influence_gam.R"))
source(paste0(DIR, "plot_gam_effects.r"))
source(paste0(DIR, "plot_influence_gam.R"))
source(paste0(DIR, "print_summary_influence_gam.R"))
source(paste0(DIR, "r2.R"))

print("--- Running Example 4: Interaction term involving focus (year*area) ---")
set.seed(456)

# Check if mgcv is available for gamSim
if (!requireNamespace("mgcv", quietly = TRUE)) {
  stop("The mgcv package is required to run this example.")
}

# Generate example data
dat_ex4 <- gamSim(1, n = 300, dist = "poisson", scale = 0.5)
dat_ex4$year <- factor(sample(2001:2005, 300, replace = TRUE))
dat_ex4$area <- factor(sample(c("North", "South"), 300, replace = TRUE))
dat_ex4$fac <- factor(sample(LETTERS[1:3], 300, replace = TRUE))

# Fit the model
model4 <- gam(y ~ year * area + s(x0) + s(x1, by = fac) + fac, data = dat_ex4, family = poisson)

# Create and calculate influence
influ4 <- create_influence_gam(model4, data = dat_ex4, focus = "year", verbose = FALSE)
influ4 <- calculate_influence(influ4, verbose = TRUE)

summary(influ4)

# Get R2 contribution summary with error handling
print("R2 Contribution Summary (Example 4):")
#r2_contribution(influ4, r2_type = "r2Dev")
r2_contribution.influence_gam(influ4, r2_type = "r2Dev")

print("Influence terms calculated (should include year:area interaction):")
print(names(influ4$influences))

# Check for interaction terms with robust pattern matching
interaction_term_name <- grep("year:area|area:year", names(influ4$influences), value = TRUE)
print(head(influ4$influences[, c("level", interaction_term_name)]))

# Plotting with error handling
print("Generating plots for Example 4...")

# Function to safely generate and print plots
plot(influ4, type = "stan", main = "Ex4: Stan Plot (focus=year, year*area model)")

plot(influ4, type = "influ", panels = TRUE, main = "Ex4: Influence Plot (focus=year, year*area model)")

plot(influ4, type = "step_influ", main = "Ex4: Step & Influence (focus=year, year*area model)")

plot(influ4, type = "interaction", term = interaction_term_name[1], 
           main = paste("Ex4: Interaction Plot (", interaction_term_name[1], ")"))

# Add CDI plot example for Example 4 (e.g., for term 'area')
print("Influence terms for CDI plot (Ex4):")
print(names(influ4$influences))

if ("area" %in% names(influ4$influences)) {
    print("Generating CDI plot for term 'area' in Example 4...")
    plot(influ4, type = "cdi", term = "area", main = "Ex4: CDI Plot (area)")
} else {
    print("Term 'area' not found in influences for CDI plot in Example 4.")
    # Try to find an alternative term for CDI plot
    other_terms <- setdiff(names(influ4$influences), c("level", interaction_term_name))
    if (length(other_terms) > 0) {
        print(paste("Attempting CDI plot with alternative term:", other_terms[1]))
        plot(influ4, type = "cdi", term = other_terms[1])
    }
}

plot_gam_effects(model4, dat_ex4)


print("--- Example 4 Finished ---")

print("--- Validation Script Completed ---")

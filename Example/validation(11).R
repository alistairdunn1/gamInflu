# Validation script for influence_gam.R

# Source the rewritten code
source("/home/ubuntu/influence_gam.R")

# Load necessary library
library(mgcv)

# --- Example 1: Poisson GAM with factor interaction ---
print("--- Running Example 1: Poisson GAM with factor interaction (focus = fac) ---")
set.seed(123)
dat <- gamSim(1, n = 200, dist = "poisson", scale = 0.5)
dat$fac <- factor(sample(LETTERS[1:4], 200, replace = TRUE))
dat$fac2 <- factor(sample(c("low", "high"), 200, replace = TRUE))

# Model with smooths and factor interaction
model1 <- gam(y ~ s(x0) + s(x1) + s(x2, by = fac) + fac + fac2, data = dat, family = poisson)

# Create and calculate influence object
influ1 <- create_influence_gam(model1, data = dat, focus = "fac", verbose = FALSE)
print(influ1)
influ1 <- calculate_influence(influ1, verbose = FALSE)
print(influ1)
summary(influ1)

# Get R2 contribution summary
print("R2 Contribution Summary (Example 1):")
print(r2_contribution(influ1, r2_type = "r2Dev"))

# Generate plots
print("Generating plots for Example 1...")
print(plot(influ1, type = "stan", main = "Ex1: Stan Plot (focus=fac)"))
print(plot(influ1, type = "step", panels = TRUE, main = "Ex1: Step Plot (focus=fac, panelled)"))
print(plot(influ1, type = "influ", panels = TRUE, main = "Ex1: Influence Plot (focus=fac, panelled)"))
print(plot(influ1, type = "step_influ", main = "Ex1: Step & Influence (focus=fac)"))

# CDI plot (experimental)
print("Influence terms for CDI plot:")
print(names(influ1$influences))
if ("s(x0)" %in% names(influ1$influences)) {
    print(plot(influ1, type = "cdi", term = "s(x0)", main = "Ex1: CDI Plot (s(x0))"))
} else {
    print("Term 's(x0)' not found in influences for CDI plot.")
}

print("--- Example 1 Finished ---")

# --- Example 2: Focus on a different factor --- 
print("--- Running Example 2: Focus on fac2 ---")
influ1b <- create_influence_gam(model1, data = dat, focus = "fac2", verbose = FALSE)
influ1b <- calculate_influence(influ1b, verbose = FALSE)
summary(influ1b)

# Get R2 contribution summary
print("R2 Contribution Summary (Example 2):")
print(r2_contribution(influ1b, r2_type = "r2Dev"))

print("Generating plots for Example 2...")
print(plot(influ1b, type = "stan", main = "Ex2: Stan Plot (focus=fac2)"))
print(plot(influ1b, type = "influ", panels=FALSE, main = "Ex2: Influence Plot (focus=fac2, overlaid)"))

print("--- Example 2 Finished ---")


# --- Example 3: Continuous focus (requires cutting) & Smooth Interaction ---
print("--- Running Example 3: Continuous focus (x1_cut) & Smooth Interaction --- ")
dat$x1_cut <- cut_number(dat$x1, 4) # Create 4 groups based on quantiles
model2b <- gam(y ~ s(x0) + s(x1) + s(x2) + ti(x1, x2) + x1_cut, data = dat, family = poisson)

influ2 <- create_influence_gam(model2b, data = dat, focus = "x1_cut", verbose = FALSE)
influ2 <- calculate_influence(influ2, verbose = FALSE)
summary(influ2)

# Get R2 contribution summary
print("R2 Contribution Summary (Example 3):")
print(r2_contribution(influ2, r2_type = "r2Dev"))

print("Generating plots for Example 3...")
print(plot(influ2, type = "stan", main = "Ex3: Stan Plot (focus=x1_cut)"))
print(plot(influ2, type = "influ", panels=FALSE, main = "Ex3: Influence Plot (focus=x1_cut, overlaid)"))

# CDI plot for interaction term
print("Influence terms for CDI plot (Ex3):")
print(names(influ2$influences))
if ("ti(x1,x2)" %in% names(influ2$influences)) {
    tryCatch({
        print(plot(influ2, type = "cdi", term = "ti(x1,x2)", main = "Ex3: CDI Plot (ti(x1,x2))"))
    }, error = function(e) { print(paste("CDI plot failed for ti(x1,x2):", e$message)) })
} else {
    print("Term 'ti(x1,x2)' not found in influences for CDI plot.")
}

print("--- Example 3 Finished ---")


# --- Example 4: Interaction term involving focus variable ---
print("--- Running Example 4: Interaction term involving focus (year*area) ---")
set.seed(456)
dat_ex4 <- gamSim(1, n = 300, dist = "poisson", scale = 0.5)
dat_ex4$year <- factor(sample(2001:2005, 300, replace = TRUE))
dat_ex4$area <- factor(sample(c("North", "South"), 300, replace = TRUE))
dat_ex4$fac <- factor(sample(LETTERS[1:3], 300, replace = TRUE))

model4 <- gam(y ~ year * area + s(x0) + s(x1, by = fac) + fac, data = dat_ex4, family = poisson)

influ4 <- create_influence_gam(model4, data = dat_ex4, focus = "year", verbose = FALSE)
influ4 <- calculate_influence(influ4, verbose = TRUE)

print("Summary for Example 4:")
summary(influ4)

# Get R2 contribution summary
print("R2 Contribution Summary (Example 4):")
print(r2_contribution.influence_gam(influ4, r2_type = "r2Dev"))

print("Influence terms calculated (should include year:area interaction):")
print(names(influ4$influences))

interaction_term_name <- grep("year:area|area:year", names(influ4$influences), value = TRUE)
if (length(interaction_term_name) > 0) {
    print(paste("Interaction term found:", interaction_term_name))
    print("Head of influences including interaction:")
    print(head(influ4$influences[, c("level", interaction_term_name)]))
} else {
    print("WARNING: Interaction term year:area not found in influences!")
}

print("Generating plots for Example 4...")
print(plot(influ4, type = "stan", main = "Ex4: Stan Plot (focus=year, year*area model)"))
print(plot(influ4, type = "influ", panels = TRUE, main = "Ex4: Influence Plot (focus=year, year*area model)"))
print(plot(influ4, type = "step_influ", main = "Ex4: Step & Influence (focus=year, year*area model)"))

print("--- Example 4 Plotting Finished ---")

# Add plot for the interaction term using the new plot type
if (length(interaction_term_name) > 0) {
    print(paste("Generating interaction plot for:", interaction_term_name[1]))
    print(plot(influ4, type = "interaction", term = interaction_term_name[1], 
         main = paste("Ex4: Interaction Plot (", interaction_term_name[1], ")")))
} else {
    print("Skipping interaction plot as term not found.")
}

# Add CDI plot example for Example 4 (e.g., for term 'area')
print("Influence terms for CDI plot (Ex4):")
print(names(influ4$influences))
if ("area" %in% names(influ4$influences)) {
    print("Generating CDI plot for term 'area' in Example 4...")
    print(plot(influ4, type = "cdi", term = "area", main = "Ex4: CDI Plot (area)"))
} else {
    print("Term 'area' not found in influences for CDI plot in Example 4.")
}

# Note on CDI for specific interaction levels:
# The standard predict(type="terms") provides the effect for the whole interaction term (e.g., year:area).
# Generating a separate CDI plot for each specific level combination (e.g., year2001:areaNorth)
# is not directly supported by this structure. Use the interaction plot above to visualize these effects.

print("--- Example 4 Finished ---")

print("--- Validation Script Completed ---")


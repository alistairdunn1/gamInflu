source("C:/Users/alist/OneDrive/Projects/Software/GAMInflu/GAMInflu/R/calculate_influence.influence_gam.R")
source("C:/Users/alist/OneDrive/Projects/Software/GAMInflu/GAMInflu/R/create_influence_gam.R")
source("C:/Users/alist/OneDrive/Projects/Software/GAMInflu/GAMInflu/R/plot_influence_gam.R")
source("C:/Users/alist/OneDrive/Projects/Software/GAMInflu/GAMInflu/R/print_summary_influence_gam.R")
source("C:/Users/alist/OneDrive/Projects/Software/GAMInflu/GAMInflu/R/print_summary_influence_gam.R")




# --- Example 4: Interaction term involving focus variable ---
# Load required packages
library(mgcv)     # For gam function and gamSim
library(ggplot2)  # For plotting
library(tidyr)    # For data reshaping
library(patchwork) # For combining plots
library(RColorBrewer) # For color palettes

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
model4 <- tryCatch({
  gam(y ~ year * area + s(x0) + s(x1, by = fac) + fac, data = dat_ex4, family = poisson)
}, error = function(e) {
  stop("Error fitting GAM model: ", e$message)
})

# Create and calculate influence
tryCatch({
  influ4 <- create_influence_gam(model4, data = dat_ex4, focus = "year", verbose = FALSE)
  influ4 <- calculate_influence(influ4, verbose = TRUE)
}, error = function(e) {
  stop("Error in influence calculation: ", e$message)
})

print("Summary for Example 4:")
summary(influ4)

# Get R2 contribution summary with error handling
print("R2 Contribution Summary (Example 4):")
r2_result <- tryCatch({
  r2_contribution(influ4, r2_type = "r2Dev")
}, error = function(e) {
  message("Error calculating R2 contribution: ", e$message)
  return(NULL)
})

if (!is.null(r2_result)) {
  print(r2_result)
}

print("Influence terms calculated (should include year:area interaction):")
print(names(influ4$influences))

# Check for interaction terms with robust pattern matching
interaction_term_name <- grep("year:area|area:year", names(influ4$influences), value = TRUE)
if (length(interaction_term_name) > 0) {
    print(paste("Interaction term found:", interaction_term_name))
    print("Head of influences including interaction:")
    print(head(influ4$influences[, c("level", interaction_term_name)]))
} else {
    print("WARNING: Interaction term year:area not found in influences!")
}

# Plotting with error handling
print("Generating plots for Example 4...")

# Function to safely generate and print plots
safe_plot <- function(expr, plot_desc) {
  result <- tryCatch({
    p <- eval(expr)
    print(p)
    TRUE
  }, error = function(e) {
    message(paste("Error generating", plot_desc, ":", e$message))
    FALSE
  })
  return(result)
}

safe_plot(
  expr = quote(plot(influ4, type = "stan", main = "Ex4: Stan Plot (focus=year, year*area model)")),
  plot_desc = "standardization plot"
)

safe_plot(
  expr = quote(plot(influ4, type = "influ", panels = TRUE, main = "Ex4: Influence Plot (focus=year, year*area model)")),
  plot_desc = "influence plot"
)

safe_plot(
  expr = quote(plot(influ4, type = "step_influ", main = "Ex4: Step & Influence (focus=year, year*area model)")),
  plot_desc = "step & influence plot"
)

print("--- Example 4 Plotting Finished ---")

# Add plot for the interaction term using the new plot type
if (length(interaction_term_name) > 0) {
    print(paste("Generating interaction plot for:", interaction_term_name[1]))
    safe_plot(
      expr = quote(plot(influ4, type = "interaction", term = interaction_term_name[1], 
           main = paste("Ex4: Interaction Plot (", interaction_term_name[1], ")"))),
      plot_desc = paste("interaction plot for", interaction_term_name[1])
    )
} else {
    print("Skipping interaction plot as term not found.")
}

# Add CDI plot example for Example 4 (e.g., for term 'area')
print("Influence terms for CDI plot (Ex4):")
print(names(influ4$influences))

if ("area" %in% names(influ4$influences)) {
    print("Generating CDI plot for term 'area' in Example 4...")
    safe_plot(
      expr = quote(plot(influ4, type = "cdi", term = "area", main = "Ex4: CDI Plot (area)")),
      plot_desc = "CDI plot for term 'area'"
    )
} else {
    print("Term 'area' not found in influences for CDI plot in Example 4.")
    # Try to find an alternative term for CDI plot
    other_terms <- setdiff(names(influ4$influences), c("level", interaction_term_name))
    if (length(other_terms) > 0) {
        print(paste("Attempting CDI plot with alternative term:", other_terms[1]))
        safe_plot(
          expr = quote(plot(influ4, type = "cdi", term = other_terms[1], 
               main = paste("Ex4: CDI Plot (", other_terms[1], ")"))),
          plot_desc = paste("CDI plot for term", other_terms[1])
        )
    }
}

plot_gam_effects(model4, dat_ex4)


print("--- Example 4 Finished ---")

print("--- Validation Script Completed ---")
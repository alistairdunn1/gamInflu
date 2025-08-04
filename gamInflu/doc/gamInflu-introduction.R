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
# # Prepare data: focus term must be a factor
# data$year <- factor(data$year)
# 
# # Create influence object (default: coefficient-based CIs)
# gi <- gam_influence(model, focus = "year")
# 
# # Alternative: use prediction-based CIs
# gi <- gam_influence(model, focus = "year", use_coeff_method = FALSE)

## ----eval=FALSE---------------------------------------------------------------
# # Calculate influence metrics (automatic family detection)
# gi <- calculate_influence(gi)
# 
# # Or with prediction-based CIs
# gi <- calculate_influence(gi, use_coeff_method = FALSE)

## ----eval=FALSE---------------------------------------------------------------
# # Extract results with standardised column names
# indices <- extract_indices(gi)
# print(indices)
# 
# # Column names follow standardised convention:
# # - focus_term: The name of the focus variable
# # - level: The levels of the focus variable
# # - index: Standardised index values (main result)
# # - cv: Coefficient of variation
# # - lower_CI: Lower confidence interval bound
# # - upper_CI: Upper confidence interval bound

## ----eval=FALSE---------------------------------------------------------------
# # Main plotting functions
# plot_standardisation(gi)    # Index comparison
# plot_stepwise_index(gi)     # Model progression
# plot_term_influence(gi)     # Term influences
# plot_residuals(gi)          # Residual diagnostics
# 
# # Alternative: Use generic plot method
# plot(gi, type = "stan")     # Same as plot_standardisation(gi)
# plot(gi, type = "step")     # Same as plot_stepwise_index(gi)
# plot(gi, type = "influ")    # Same as plot_term_influence(gi)

## ----eval=FALSE---------------------------------------------------------------
# # Default method - coefficient-based
# gi <- gam_influence(model, focus = "year")  # use_coeff_method = TRUE by default
# gi <- calculate_influence(gi)

## ----eval=FALSE---------------------------------------------------------------
# # Prediction-based method
# gi <- gam_influence(model, focus = "year", use_coeff_method = FALSE)
# gi <- calculate_influence(gi)

## ----eval=FALSE---------------------------------------------------------------
# # Prepare data
# fish_data$year <- factor(fish_data$year)
# 
# # Fit binomial model
# mod_binom <- gam(presence ~ s(depth) + s(temp) + year,
#                  data = fish_data, family = binomial())
# 
# # Influence analysis
# gi_binom <- gam_influence(mod_binom, focus = "year")
# gi_binom <- calculate_influence(gi_binom)  # Auto-detects binomial
# 
# # Visualise results
# plot_standardisation(gi_binom)

## ----eval=FALSE---------------------------------------------------------------
# # Prepare data
# survey_data$year <- factor(survey_data$year)
# 
# # Fit gamma model
# mod_gamma <- gam(biomass ~ s(effort) + s(sst) + year,
#                  data = survey_data, family = Gamma(link="log"))
# 
# # Influence analysis
# gi_gamma <- gam_influence(mod_gamma, focus = "year")
# gi_gamma <- calculate_influence(gi_gamma)  # Auto-detects gamma
# 
# # Extract indices
# indices_gamma <- extract_indices(gi_gamma)

## ----eval=FALSE---------------------------------------------------------------
# # Prepare data
# catch_data$year <- factor(catch_data$year)
# 
# # Fit Poisson model
# mod_pois <- gam(fish_count ~ s(depth) + s(longitude) + year,
#                 data = catch_data, family = poisson())
# 
# # Influence analysis
# gi_pois <- gam_influence(mod_pois, focus = "year")
# gi_pois <- calculate_influence(gi_pois)   # Auto-detects Poisson

## ----eval=FALSE---------------------------------------------------------------
# # Prepare data
# fisheries_data$year <- factor(fisheries_data$year)
# 
# # Fit Tweedie model
# mod_tweedie <- gam(catch_kg ~ s(depth) + s(vessel_power) + year,
#                    data = fisheries_data, family = tw())
# 
# # Influence analysis
# gi_tweedie <- gam_influence(mod_tweedie, focus = "year")
# gi_tweedie <- calculate_influence(gi_tweedie)  # Auto-detects Tweedie

## ----eval=FALSE---------------------------------------------------------------
# # Analyse residual patterns for model improvement
# residual_analysis <- analyse_residual_patterns(gi)
# print(residual_analysis)
# 
# # This shows variables that may improve model fit
# # if significant patterns are found in residuals

## ----eval=FALSE---------------------------------------------------------------
# # Subset analysis for specific area
# gi_subset <- calculate_influence(gi, subset_var = "area", subset_value = "North")

## ----eval=FALSE---------------------------------------------------------------
# # Fit separate models for catch probability and positive catch
# mod_binom <- gam(presence ~ s(depth) + year, data = data, family = binomial())
# mod_gamma <- gam(catch ~ s(depth) + year, data = data[data$catch > 0,], family = Gamma())
# 
# # Create influence objects
# gi_binom <- gam_influence(mod_binom, focus = "year")
# gi_gamma <- gam_influence(mod_gamma, focus = "year")
# 
# # Calculate influences
# gi_binom <- calculate_influence(gi_binom)
# gi_gamma <- calculate_influence(gi_gamma)
# 
# # Combine for delta-GLM analysis
# gi_combined <- combine_indices(gi_binom, gi_gamma)
# plot_standardisation(gi_combined)

## ----eval=FALSE---------------------------------------------------------------
# # Index comparison with confidence intervals
# plot_standardisation(gi)
# 
# # Show only standardised indices (cleaner)
# plot_standardisation(gi, show_unstandardised = FALSE)
# 
# # Model progression showing term additions
# plot_stepwise_index(gi)
# 
# # Term influence on final standardised index
# plot_term_influence(gi)

## ----eval=FALSE---------------------------------------------------------------
# # Standard GAM diagnostic plots
# plot_residuals(gi, type = "standard")
# 
# # Violin plots by focus levels (shows distribution)
# plot_residuals(gi, type = "violin")
# 
# # Faceted violin plots by another variable
# plot_residuals(gi, type = "violin", by = "area")

## ----eval=FALSE---------------------------------------------------------------
# # CDI plot for specific smooth term
# plot(gi, type = "cdi", term = "s(temp)")
# 
# # Data distribution for specific term
# plot(gi, type = "distribution", term = "s(temp)")

## ----eval=FALSE---------------------------------------------------------------
# # Always convert focus term to factor BEFORE fitting model
# data$year <- factor(data$year)
# data$area <- factor(data$area)
# 
# # Fit model with factor variables
# model <- gam(response ~ s(continuous_var) + factor_var, data = data)

## ----eval=FALSE---------------------------------------------------------------
# # Let the package auto-detect family (recommended)
# gi <- calculate_influence(gi)
# 
# # Only override if needed
# gi <- calculate_influence(gi, family_method = "binomial")

## ----eval=FALSE---------------------------------------------------------------
# # Use coefficient-based for traditional fisheries work
# gi <- gam_influence(model, focus = "year")  # Default
# 
# # Use prediction-based for modern GAM applications
# gi <- gam_influence(model, focus = "year", use_coeff_method = FALSE)

## ----eval=FALSE---------------------------------------------------------------
# # Always check residual patterns
# residual_analysis <- analyse_residual_patterns(gi)
# print(residual_analysis)
# 
# # Use diagnostic plots
# plot_residuals(gi, type = "violin")


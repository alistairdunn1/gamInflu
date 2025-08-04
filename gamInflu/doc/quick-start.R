## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(gamInflu)

## ----eval=FALSE---------------------------------------------------------------
# # Basic workflow
# library(gamInflu)
# library(mgcv)
# 
# # 1. Prepare data (focus must be factor)
# data$year <- factor(data$year)
# 
# # 2. Fit model
# model <- gam(response ~ s(x) + year, data = data)
# 
# # 3. Influence analysis
# gi <- gam_influence(model, focus = "year")
# gi <- calculate_influence(gi)
# 
# # 4. Results
# indices <- extract_indices(gi)
# plot_standardisation(gi)


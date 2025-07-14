# Install dependencies
install.packages(c("mgcv", "ggplot2", "dplyr", "tidyr", "gridExtra", "viridis", "R6", "rlang"))

# Load the package files
source("C:/Users/alist/OneDrive/Projects/Software/GAMInflu/GAMInflu/R/imports.R")
source("C:/Users/alist/OneDrive/Projects/Software/GAMInflu/GAMInflu/R/gam_influence_private.R")
source("C:/Users/alist/OneDrive/Projects/Software/GAMInflu/GAMInflu/R/gam_influence_calculations.R") 
source("C:/Users/alist/OneDrive/Projects/Software/GAMInflu/GAMInflu/R/gam_influence_cdi_plots.R")
source("C:/Users/alist/OneDrive/Projects/Software/GAMInflu/GAMInflu/R/gam_influence_smooth_plots.R")
source("C:/Users/alist/OneDrive/Projects/Software/GAMInflu/GAMInflu/R/gam_influence_class.R")
source("C:/Users/alist/OneDrive/Projects/Software/GAMInflu/GAMInflu/R/gam_influence.R")

# Example usage
library(mgcv)
model <- gam(mpg ~ s(wt) + te(wt, hp) + factor(cyl), data = mtcars)
influence <- gam_influence(model, focus = "factor(cyl)")
influence$calc()
influence$stan_plot() 
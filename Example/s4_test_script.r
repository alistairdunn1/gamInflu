# Test Script for S4 GAM Influence Implementation
# This script tests the S4 implementation of the GAM influence analysis package

# Load required libraries
library(mgcv)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(viridis)
library(rlang)
library(methods)

# Source all package files in correct order
source("C:/Users/alist/OneDrive/Projects/Software/GAMInflu/GAMInflu/R/utils.R")
source("C:/Users/alist/OneDrive/Projects/Software/GAMInflu/GAMInflu/R/smooth-parsing.R") 
source("C:/Users/alist/OneDrive/Projects/Software/GAMInflu/GAMInflu/R/calculations.R")
source("C:/Users/alist/OneDrive/Projects/Software/GAMInflu/GAMInflu/R/data-preparation.R")
source("C:/Users/alist/OneDrive/Projects/Software/GAMInflu/GAMInflu/R/plotting.R")
source("C:/Users/alist/OneDrive/Projects/Software/GAMInflu/GAMInflu/R/GAMInfluence-class.R")
source("C:/Users/alist/OneDrive/Projects/Software/GAMInflu/GAMInflu/R/methods.R")
source("C:/Users/alist/OneDrive/Projects/Software/GAMInflu/GAMInflu/R/gam_influence.R")

cat("=== Testing S4 GAM Influence Implementation ===\n\n")

# Test 1: Basic object creation
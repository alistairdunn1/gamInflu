# gamInflu

**gamInflu** provides comprehensive influence analysis tools for Generalised Additive Models (GAMs) fitted with the `mgcv` package in R. Enhanced with **multi-family support** for Gaussian, binomial, gamma, and Poisson distributions. It supports all smoother types (`s()`, `te()`, `ti()`, `t2()`, and `by=` terms), and generates stepwise index plots, term influence plots, coefficient-distribution-influence (CDI) plots, diagnostics for random effects, and family-specific standardised indices to help you understand model structure and the influence of each term.

[![R Package](https://img.shields.io/badge/R-package-blue.svg)](https://www.r-project.org/)
[![Version](https://img.shields.io/badge/version-0.1-orange.svg)](https://github.com/alistairdunn1/gamInflu)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

---

## Family Support

**gamInflu**  supports multiple GLM families with automatic detection and family-specific statistical methods:

- **Gaussian** 📈 Traditional log-normal CPUE standardisation
- **Binomial** 🎯 Presence/absence, catch probability, proportion data
- **Gamma** 📊 Positive continuous data (biomass, CPUE without zeros)
- **Poisson** 🔢 Count data (fish numbers, abundance indices)

**Key Features:**

- ✅ Automatic family detection from model objects
- ✅ Family-appropriate index calculation methods
- ✅ Enhanced statistical rigour with geometric mean support
- ✅ Comprehensive validation and error handling
- ✅ Backward compatibility with existing workflows

---

## Installation

```r
# If you have the source directory:
devtools::install_github("alistairdunn1/gamInflu", subdir = "gamInflu")
# Or, if you have a tar.gz:
install.packages("gamInflu_0.1.0.tar.gz", repos = NULL, type = "source")
```

---

## Quick Start

### Basic Workflow (All Families)

```r
library(gamInflu)

# Prepare data: focus term must be a factor
your_data$year <- factor(your_data$year)

# 1. Create influence object (works with any GLM family)
gi <- gam_influence(your_gam_model, focus = "year")

# 2. Calculate influence metrics (automatic family detection)
gi <- calculate_influence(gi)

# 3. Visualize results
plot_standardisation(gi)    # Index comparison
plot_stepwise_index(gi)     # Model progression
plot_term_influence(gi)     # Term influences
```

### Family-Specific Examples

```r
library(mgcv)

# Important: Focus term must be a factor!
# Best practice: Convert to factor in data beforehand

# Binomial model (presence/absence)
fish_data$year <- factor(fish_data$year)
mod_binom <- gam(presence ~ s(depth) + s(temp) + year, 
                 data = fish_data, family = binomial())
gi_binom <- gam_influence(mod_binom, focus = "year")
gi_binom <- calculate_influence(gi_binom)  # Auto-detects binomial

# Gamma model (biomass)
survey_data$year <- factor(survey_data$year)
mod_gamma <- gam(biomass ~ s(effort) + s(sst) + year, 
                 data = survey_data, family = Gamma(link="log"))
gi_gamma <- gam_influence(mod_gamma, focus = "year")
gi_gamma <- calculate_influence(gi_gamma)  # Auto-detects gamma

# Poisson model (count data)
catch_data$year <- factor(catch_data$year)
mod_pois <- gam(fish_count ~ s(depth) + s(longitude) + year, 
                data = catch_data, family = poisson())
gi_pois <- gam_influence(mod_pois, focus = "year")
gi_pois <- calculate_influence(gi_pois)   # Auto-detects poisson
```

---

## Family Support Details

### Automatic Family Detection

By default, `calculate_influence()` automatically detects the GLM family from your model:

```r
# No need to specify family - automatically detected
gi <- calculate_influence(gi)
```

### Manual Family Specification

Override automatic detection when needed:

```r
# Force specific family method
gi <- calculate_influence(gi, family_method = "binomial")
```

### Family-Specific Methods

| Family       | Data Type           | Index Calculation         | Use Cases                           |
| ------------ | ------------------- | ------------------------- | ----------------------------------- |
| **Gaussian** | Continuous          | Geometric/arithmetic mean | Log-normal CPUE, linear models      |
| **Binomial** | Binary/proportions  | Proportion-based          | Presence/absence, catch probability |
| **Gamma**    | Positive continuous | Geometric mean            | Biomass, positive CPUE              |
| **Poisson**  | Non-negative counts | Count-appropriate         | Fish numbers, abundance counts      |

---

## Main Functions and Examples

### Create a gam_influence object

Initializes the analysis object for your fitted GAM.

```r
# Assuming your data has been prepared with year as factor
gi <- gam_influence(mod, focus = "year")
# or explicitly provide data if needed:
# gi <- gam_influence(mod, focus = "year", data = mydata)
```

### Calculate influence metrics

Computes all indices, predictions, and influence statistics.

```r
gi <- calculate_influence(gi)
```

### Plot standardisation

Shows the unstandardised and standardised indices for the focus term.

```r
plot_standardisation(gi)
```

### Plot stepwise index

Shows how the index changes as each term is added to the model.

```r
plot_stepwise_index(gi)
```

### Plot term influence

Shows the influence of each non-focus term on the focus index.

```r
plot_term_influence(gi)
```

### Plot combined step and influence

Displays both the stepwise index and term influence plots side-by-side.

```r
plot_step_and_influence(gi)
```

### Plot coefficient-distribution-influence (CDI) for a term

Multi-panel plot showing the effect, data distribution, and influence for a specified term.

```r
plot_cdi(gi, term = "s(temp)")
```

or use an integer index to refer to the term,

```r
plot_cdi(gi, term = 2)
```

### Plot predicted effects for all terms

Visualises the predicted effects for each term in the model.

```r
plot_terms(gi)
```

### Plot predicted effects for a single term (including random effects and by-variable panels)

```r
plot_terms(gi, term = "s(site, bs='re')", re_type = "points")
plot_terms(gi, term = "s(temp, by=season)")
plot_terms(gi, term = "te(lon, lat)")
```

or

```r
plot_terms(gi, term = 3, re_type = "points")
plot_terms(gi, term = 4)
plot_terms(gi, term = 5)
```

### Plot data distribution for a term

Shows the distribution of the data for a specified model term and the focus variable.

```r
plot_term_distribution(gi, term = "s(temp)")
```

### Plot random effect diagnostics (QQ plot, histogram, caterpillar plot, points)

```r
plot_re(gi, term = "site", re_type = "qq")
plot_re(gi, term = "site", re_type = "hist")
plot_re(gi, term = "site", re_type = "caterpillar")
plot_re(gi, term = "site", re_type = "points")
```

### Extract model terms

Get the variable names or full term expressions used in the model.

```r
get_terms(gi)            # Variable names
get_terms(gi, full=TRUE) # Full term expressions (e.g., s(temp), te(lon,lat))
```

### Model progression statistics

Extracts a summary table of model progression statistics (AIC, R-squared, deviance explained).

```r
r2(gi)
```

### Print a summary of term contributions

Prints a formatted summary table of term contributions and influence metrics.

```r
summary(gi)
```

---

## Advanced Usage

#### Plotting with by-variable panels

```r
plot_terms(gi, term = "s(temp, by=season)")
```

#### Plotting random effects

```r
plot_terms(gi, term = "s(site, bs='re')", re_type = "points")
```

#### Plotting 2D tensor smooths

```r
plot_terms(gi, term = "te(lon, lat)")
```

---

## Requirements

- R >= 4.0.0
- mgcv (for GAM fitting)
- ggplot2 (for plotting)
- patchwork (plot composition)
- dplyr, tidyr (data manipulation)
- rlang (non-standard evaluation)

---

## Family-Specific Examples

### Binomial Models (Presence/Absence)

```r
library(mgcv)
library(gamInflu)

# Fit binomial GAM for presence/absence data
# Prepare data with factor
survey_data$year <- factor(survey_data$year)
mod_presence <- gam(fish_present ~ s(depth) + s(temperature) + year, 
                    data = survey_data, 
                    family = binomial())

# Analyse presence probability trends
gi_presence <- gam_influence(mod_presence, focus = "year")
gi_presence <- calculate_influence(gi_presence)  # Auto-detects binomial

# Visualize results
plot_standardisation(gi_presence)     # Presence probability index
plot_stepwise_index(gi_presence)      # Model building effects  
plot_term_influence(gi_presence)      # Environmental influences
```

### Gamma Models (Biomass Data)

```r
# Fit Gamma GAM for positive biomass data
# Prepare data with factor
trawl_data$year <- factor(trawl_data$year)
mod_biomass <- gam(biomass ~ s(effort) + s(sst) + te(lon, lat) + year,
                   data = trawl_data,
                   family = Gamma(link = "log"))

# Analyse biomass index trends
gi_biomass <- gam_influence(mod_biomass, focus = "year") 
gi_biomass <- calculate_influence(gi_biomass)  # Auto-detects gamma

# Comprehensive analysis
plot_step_and_influence(gi_biomass)   # Combined stepwise + influence
plot_cdi(gi_biomass, term = "s(sst)") # Sea surface temperature effects
plot_terms(gi_biomass)                # All model terms
```

### Poisson Models (Count Data)

```r
# Fit Poisson GAM for fish count data
# Prepare data with factor
acoustic_data$year <- factor(acoustic_data$year)
mod_counts <- gam(fish_count ~ s(depth) + s(current_speed) + area + year,
                  data = acoustic_data,
                  family = poisson())

# Analyse abundance index
gi_counts <- gam_influence(mod_counts, focus = "year")
gi_counts <- calculate_influence(gi_counts)  # Auto-detects poisson

# Random effects diagnostics  
plot_re(gi_counts, term = "area", re_type = "caterpillar")
plot_re(gi_counts, term = "area", re_type = "qq")
```

### Advanced Family Options

```r
# Manual family specification (override auto-detection)
gi_manual <- calculate_influence(gi_object, 
                                family_method = "gamma",
                                rescale_method = "geometric_mean",
                                confidence_level = 0.90)

# Custom rescaling for specific applications
gi_custom <- calculate_influence(gi_object,
                                rescale_method = "custom", 
                                custom_rescale_value = 100)  # Scale to 100

# Enhanced confidence intervals
gi_wide <- calculate_influence(gi_object, confidence_level = 0.99)
```

---

## Citation

To cite the **gamInflu** package in publications, please use:

> Dunn, A. (2025). *gamInflu: Influence Analysis Tools for Generalized Additive Models in R*. R package version 0.1.0. Enhanced with multi-family support for Gaussian, binomial, gamma, and Poisson distributions. Available at: https://github.com/alistairdunn1/gamInflu

Please also cite the foundational methodology:

> Bentley, N., Kendrick, T. H., Starr, P. J., & Breen, P. A. (2012). Influence plots and metrics: tools for better understanding fisheries catch-per-unit-effort standardisations. *ICES Journal of Marine Science*, 69(1), 84–88. https://doi.org/10.1093/icesjms/fsr174

You can generate a citation in R with:

```r
citation("gamInflu")
```

---

## Support

For questions, bug reports, or feature requests, please contact the maintainer

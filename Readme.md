# gamInflu

**gamInflu** provides influence analysis tools for Generalised Additive Models (GAMs) fitted with the `mgcv` package in R. The package supports Gaussian, binomial, gamma, Poisson, and Tweedie distributions with automatic family detection. It offers both traditional coefficient-based confidence intervals (compatible with influ.r) and modern prediction-based methods. The package handles all smoother types (`s()`, `te()`, `ti()`, `t2()`, and `by=` terms) and generates stepwise index plots, term influence plots, coefficient-distribution-influence (CDI) plots, diagnostics for random effects, and family-specific standardised indices to understand model structure and the influence of each term.

[![R Package](https://img.shields.io/badge/R-package-blue.svg)](https://www.r-project.org/)
[![Version](https://img.shields.io/badge/version-0.1-orange.svg)](https://github.com/alistairdunn1/gamInflu)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

---

## Family Support

**gamInflu** supports multiple GLM families with automatic detection and family-specific statistical methods:

- **Gaussian** 📈 Traditional log-normal CPUE standardisation
- **Binomial** 🎯 Presence/absence, catch probability, proportion data
- **Gamma** 📊 Positive continuous data (biomass, CPUE without zeros)
- **Poisson** 🔢 Count data (fish numbers, abundance indices)
- **Tweedie** 🎲 Semi-continuous data (fisheries catch with exact zeros)

---

## Installation

```r
# Install from source directory:
devtools::install_github("alistairdunn1/gamInflu", subdir = "gamInflu")
# Or install from tar.gz file:
install.packages("gamInflu_0.1.0.tar.gz", repos = NULL, type = "source")
```

---

## Quick Start

### Basic Workflow (All Families)

```r
library(gamInflu)

# Prepare data: focus term must be a factor
data$year <- factor(data$year)

# 1. Create influence object (works with any GLM family)
gi <- gam_influence(model, focus = "year")  # Default: coefficient-based CIs
# or: gi <- gam_influence(model, focus = "year", use_coeff_method = FALSE)  # prediction-based CIs

# 2. Calculate influence metrics (automatic family detection)
gi <- calculate_influence(gi)  # Default: coefficient-based CIs (influ.r approach)
# or: gi <- calculate_influence(gi, use_coeff_method = FALSE)  # prediction-based CIs

# 3. Extract results with standardised column names
indices <- extract_indices(gi)  # Includes: index, cv, lower_CI, upper_CI, etc.

# 4. Visualise results
plot_standardisation(gi)    # Index comparison
plot_stepwise_index(gi)     # Model progression
plot_term_influence(gi)     # Term influences

# 5. Subset analysis for interaction models
gi_subset <- calculate_influence(gi, subset_var = "area", subset_value = "North")
```

### Family-Specific Examples

```r
library(mgcv)

# Important: Focus term must be a factor
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

# Tweedie model (semi-continuous data with exact zeros)
fisheries_data$year <- factor(fisheries_data$year)
mod_tweedie <- gam(catch_kg ~ s(depth) + s(vessel_power) + year,
                   data = fisheries_data, family = tw())
gi_tweedie <- gam_influence(mod_tweedie, focus = "year")
gi_tweedie <- calculate_influence(gi_tweedie)  # Auto-detects Tweedie

# CI method examples (works with any family):
# Default coefficient-based approach (influ.r style)
gi_coeff <- gam_influence(mod_gamma, focus = "year")  # use_coeff_method = TRUE default
gi_coeff <- calculate_influence(gi_coeff)

# Modern prediction-based approach  
gi_pred <- gam_influence(mod_gamma, focus = "year", use_coeff_method = FALSE)
gi_pred <- calculate_influence(gi_pred)
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

| Family       | Data Type           | Index Calculation         | Use Cases                             |
| ------------ | ------------------- | ------------------------- | ------------------------------------- |
| **Gaussian** | Continuous          | Geometric/arithmetic mean | Log-normal CPUE, linear models        |
| **Binomial** | Binary/proportions  | Proportion-based          | Presence/absence, catch probability   |
| **Gamma**    | Positive continuous | Geometric mean            | Biomass, positive CPUE                |
| **Poisson**  | Non-negative counts | Count-appropriate         | Fish numbers, abundance counts        |
| **Tweedie**  | Semi-continuous     | Gamma-like treatment      | Fisheries catch data with exact zeros |

---

## Confidence Interval Calculation Methods

**gamInflu** provides two methods for calculating confidence intervals, allowing users to choose between traditional and modern approaches:

### Coefficient-Based Method (Default)

The **coefficient-based approach** follows the traditional influ.r methodology, providing backwards compatibility and familiar results for users transitioning from influ.r.

```r
# Default method - coefficient-based (influ.r approach)
gi <- gam_influence(model, focus = "year")  # use_coeff_method = TRUE by default
gi <- calculate_influence(gi)

# Or explicitly specify
gi <- gam_influence(model, focus = "year", use_coeff_method = TRUE)
gi <- calculate_influence(gi)
```

**Characteristics:**
- Uses model coefficients directly for relative effect calculations
- Follows established influ.r methodology using the **Francis method**
- Often provides more pronounced between-group differences
- Ideal for users familiar with traditional fisheries standardisation approaches
- Confidence intervals calculated as: `exp(coefficients ± CI_multiplier × SEs)`

**Francis Method Implementation:**
The coefficient-based approach implements the Francis method with proper variance-covariance matrix transformation to ensure **non-zero confidence intervals for all levels**, including reference levels. This is achieved through:
- Q matrix transformation: `Q %*% vcov %*% t(Q)` where Q transforms contrasts to relative-to-mean effects
- Guaranteed non-zero CI ranges for all factor levels (no zero-width intervals)
- Proper uncertainty propagation from model coefficients to standardised indices
- Consistent with established fisheries assessment methodology

### Prediction-Based Method (Modern)

The **prediction-based approach** uses modern mgcv prediction methods with full uncertainty propagation.

```r
# Prediction-based method (modern approach)
gi <- gam_influence(model, focus = "year", use_coeff_method = FALSE)
gi <- calculate_influence(gi)
```

**Characteristics:**
- Uses model predictions with complete uncertainty propagation  
- Incorporates all sources of model uncertainty
- More conservative confidence intervals
- Applies delta method for log-link models automatically
- Modern statistical approach following mgcv best practices

### Method Comparison

| Aspect                      | Coefficient-Based (Default)  | Prediction-Based (Modern)   |
| --------------------------- | ---------------------------- | --------------------------- |
| **Backwards Compatibility** | ✅ influ.r equivalent         | ❌ Different from influ.r    |
| **Index Differences**       | More pronounced              | More conservative           |
| **Uncertainty Handling**    | Francis method with Q matrix | Full prediction uncertainty |
| **Statistical Approach**    | Traditional (coefficient)    | Modern mgcv methods         |
| **CI Coverage**             | Non-zero CIs for all levels  | Non-zero CIs for all levels |
| **Use Case**                | Familiar to influ.r users    | Modern GAM best practices   |

### Choosing the Right Method

- **Use coefficient-based (default)** when:
  - Transitioning from influ.r
  - Wanting familiar, established results with Francis method guarantees
  - Working with traditional fisheries applications
  - Need more pronounced index differences
  - Require non-zero confidence intervals for all levels (including reference)

- **Use prediction-based** when:
  - Following modern GAM best practices  
  - Want complete uncertainty propagation
  - Prefer conservative confidence intervals
  - Working with complex model structures

---

## Main Functions and Examples

### Create a gam_influence object

Initializes the analysis object for your fitted GAM.

```r
# Assuming your data has been prepared with year as factor
gi <- gam_influence(mod, focus = "year")  # Default: coefficient-based CIs
# or explicitly provide data and method:
# gi <- gam_influence(mod, focus = "year", data = mydata, use_coeff_method = FALSE)
```

**Parameters:**
- `model`: A fitted GAM or GLM model object
- `focus`: Character string specifying the focus term (must be a factor)
- `data`: Optional data frame (extracted from model if not provided)
- `use_coeff_method`: Logical, use coefficient-based CIs (TRUE, default) or prediction-based CIs (FALSE)

### Calculate influence metrics

Computes all indices, predictions, and influence statistics.

```r
gi <- calculate_influence(gi)
```

### Get standardised and unstandardised indices

Extract the indices as a data frame for further analysis or export. Includes coefficient of variation (CV) for assessing precision. The package uses a standardised column naming convention for consistency:

```r
indices_df <- extract_indices(gi)
print(indices_df)

# Column names follow standardised convention:
# - focus_term: The name of the focus variable
# - level: The levels of the focus variable
# - unstandardised_index: Raw index values
# - unstandardised_cv: Coefficient of variation for raw index
# - index: Standardised index values (main result)
# - se: Standard error of standardised index
# - cv: Coefficient of variation of standardised index
# - lower_CI: Lower confidence interval bound
# - upper_CI: Upper confidence interval bound

# Check precision: lower CV indicates more precise estimates
summary(indices_df$cv)

# Export including all uncertainty information
write.csv(indices_df, "focus_indices.csv", row.names = FALSE)
```

The package uses `type="response"` predictions for calculating confidence intervals, providing accurate uncertainty estimates that include all model uncertainty, not just partial effects.

### Plot standardisation

Shows the unstandardised and standardised indices for the focus term. Option to show only standardised for cleaner presentation.

```r
plot_standardisation(gi)                          # Default: both indices with legend
plot_standardisation(gi, show_unstandardised = FALSE)  # Clean: standardised only, no legend
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

### Subset Analysis

Perform influence analysis on a subset of the data whilst using the full model for predictions. This is particularly useful for models with interaction terms (e.g., `year:area`) where focus term effects within specific subsets are of interest:

```r
# For models with interactions like year + area + year:area
# Analyse year effects specifically in the North area
gi_subset <- calculate_influence(gi, subset_var = "area", subset_value = "North")

# The subset analysis:
# - Uses only North area data for calculations
# - Applies the full model (including interactions) for predictions
# - Preserves the original data in the returned object
# - Provides focus term effects specific to that subset

# Extract indices for the subset
subset_indices <- extract_indices(gi_subset)

# Compare with full analysis
full_indices <- extract_indices(gi)
```

### Plotting with by-variable panels

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

### Tweedie Models (Semi-Continuous Data)

```r
# Fit Tweedie GAM for fisheries catch data with exact zeros
# Prepare data with factor
fisheries_data$year <- factor(fisheries_data$year)
mod_catch <- gam(catch_kg ~ s(depth) + s(vessel_power) + area + year,
                 data = fisheries_data,
                 family = tw())

# Analyse catch index trends
gi_catch <- gam_influence(mod_catch, focus = "year")
gi_catch <- calculate_influence(gi_catch)  # Auto-detects Tweedie

# Tweedie-specific analysis
plot_standardisation(gi_catch)        # Catch probability × abundance index
plot_term_influence(gi_catch)         # Environmental and operational influences
extract_indices(gi_catch)             # Export with proper zero handling
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

### Statistical Methodology

### Coefficient of Variation (CV) Calculation

`gamInflu` uses mathematically appropriate methods for calculating coefficients of variation based on the model family:

**For Log-Link Models** (Gamma, Poisson, Tweedie with log link):

- Uses the **delta method**: `CV = √(exp(σ²) - 1)`
- Where σ is the standard error on the log scale
- Properly accounts for the non-linear transformation from log to response scale
- `islog = FALSE` (log link ≠ pre-logged response)

**For Pre-Logged Response Models** (response names starting with `log_`, `ln_`, or `log(`):

- Uses standard definition: `CV = SE / μ` on the log scale
- `islog = TRUE` (automatically detected or user-specified)
- Plotting functions will exponentiate effects for display

**For Linear Models** (Gaussian, identity link):

- Uses the standard definition: `CV = SE / μ`
- Where SE is the standard error and μ is the mean on the response scale
- `islog = FALSE`

**For Other Link Functions**:

- Uses response-scale calculations: `CV = SE_response / μ_response`
- Ensures consistency across different GLM families

This approach ensures that CVs are always positive, mathematically correct, and comparable across different model types and confidence levels.

---

## Citation

To cite the **gamInflu** package in publications:

> Dunn, A. (2025). *gamInflu: Influence Analysis Tools for Generalised Additive Models in R*. R package version 0.1.0. Available at: https://github.com/alistairdunn1/gamInflu

Please also cite the foundational methodology:

> Bentley, N., Kendrick, T. H., Starr, P. J., & Breen, P. A. (2012). Influence plots and metrics: tools for better understanding fisheries catch-per-unit-effort standardisations. *ICES Journal of Marine Science*, 69(1), 84–88. https://doi.org/10.1093/icesjms/fsr174

Generate a citation in R:

```r
citation("gamInflu")
```

---

## Support

For questions, bug reports, or feature requests, please contact the maintainer.

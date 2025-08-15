# gamInflu

**gamInflu** provides influence analysis tools for Generalised Additive Models (GAMs) fitted with the `mgcv` package in R. The package supports Gaussian, binomial, gamma, Poisson, and Tweedie distributions with automatic family detection. It offers both traditional coefficient-based confidence intervals and modern prediction-based methods for the model terms. The package handles smoother types (`s()`, `te()`, `ti()`, `t2()`, and `by=` terms) and generates stepwise index plots, term influence plots, coefficient-distribution-influence (CDI) plots, residual diagnostics, residual pattern analysis for model adequacy assessment, delta-GLM analysis (combined indices) for fisheries data, diagnostics for random effects, and family-specific standardised indices to understand model structure and the influence of each term.

[![R Package](https://img.shields.io/badge/R-package-blue.svg)](https://www.r-project.org/)
[![Version](https://img.shields.io/badge/version-0.2-orange.svg)](https://github.com/alistairdunn1/gamInflu)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

---

## Installation

```r
# Install from source directory:
devtools::install_github("alistairdunn1/gamInflu", subdir = "gamInflu")
# Or install from tar.gz file:
install.packages("gamInflu_0.2.0.tar.gz", repos = NULL, type = "source")
```

---

## Quick Start

### Basic Workflow

```r
library(gamInflu)

# Prepare data: focus term (e.g., year) must be a factor
data$year <- factor(data$year)

# 1. Create influence object (works with any GLM family)
gi <- gam_influence(model, focus = "year")  # Default: coefficient-based CIs

# 2. Calculate influence metrics (automatic family detection)
gi <- calculate_influence(gi)

# 3. Extract results with standardised column names
indices <- extract_indices(gi)  # Includes: index, cv, lower_CI, upper_CI, etc.

# 4. Visualise results
plot_standardisation(gi)    # Index comparison
plot_stepwise_index(gi)     # Model progression
plot_term_influence(gi)     # Term influences
plot_residuals(gi)          # Residual diagnostics

# 5. Check model adequacy through residual pattern analysis
residual_analysis <- analyse_residual_patterns(gi)
print(residual_analysis)    # Shows variables that may improve model fit
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
gi_pois <- calculate_influence(gi_pois)   # Auto-detects Poisson

# Tweedie model (semi-continuous data with exact zeros)
fisheries_data$year <- factor(fisheries_data$year)
mod_tweedie <- gam(catch_kg ~ s(depth) + s(vessel_power) + year,
                   data = fisheries_data, family = tw())
gi_tweedie <- gam_influence(mod_tweedie, focus = "year")
gi_tweedie <- calculate_influence(gi_tweedie)  # Auto-detects Tweedie
```

---

## Family Support

**gamInflu** supports multiple GLM families with automatic detection and family-specific statistical methods:

- **Gaussian** 📈 Traditional log-normal CPUE standardisation  
- **Binomial** 🎯 Presence/absence, catch probability, proportion data
- **Gamma** 📊 Positive continuous data (biomass, CPUE without zeros)
- **Poisson** 🔢 Count data (fish numbers, abundance indices)
- **Tweedie** 🎲 Semi-continuous data (fisheries catch with exact zeros)

### Automatic Family Detection

By default, `calculate_influence()` automatically detects the GLM family from your model, and if the response is denoted as 'log_' or 'log(...)' will interpret the model as being fitted in log space with an identity link:

```r
# No need to specify family - automatically detected
gi <- calculate_influence(gi)
```

### Manual Family Specification

Override automatic detection when needed:

```r
# Force specific family method
gi <- calculate_influence(gi, family = "binomial")
# or for Gaussian with log transformation:
gi <- calculate_influence(gi, family = "gaussian", islog = TRUE)
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

The **coefficient-based approach** follows traditional fisheries methodology, providing familiar results for users working with established standardisation approaches.

```r
# Default method - coefficient-based
gi <- gam_influence(model, focus = "year")  # use_coeff_method = TRUE by default
gi <- calculate_influence(gi)

# Or explicitly specify
gi <- gam_influence(model, focus = "year", use_coeff_method = TRUE)
gi <- calculate_influence(gi)
```

**Characteristics:**
- Uses model coefficients directly for relative effect calculations
- Follows established fisheries methodology using the **Francis method**
- Often provides more pronounced between-group differences
- Confidence intervals calculated as: `exp(coefficients ± CI_multiplier × SEs)`

**Francis Method Implementation:**
The coefficient-based approach implements the Francis method with the variance-covariance matrix transformation to ensure **non-zero confidence intervals for all levels**, including reference levels. This is achieved through:
- Q matrix transformation: `Q %*% vcov %*% t(Q)` where Q transforms contrasts to relative-to-mean effects
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
- Incorporates full model uncertainty
- More conservative confidence intervals
- Applies delta method for log-link models automatically
- Modern statistical approach following mgcv best practices

### Method Comparison

| Aspect                      | Coefficient-Based (Default)  | Prediction-Based (Modern)   |
| --------------------------- | ---------------------------- | --------------------------- |
| **Backwards Compatibility** | ✅ Traditional equivalent     | ❌ Modern approach           |
| **Index Differences**       | More pronounced              | More conservative           |
| **Uncertainty Handling**    | Francis method with Q matrix | Full prediction uncertainty |
| **Statistical Approach**    | Traditional (coefficient)    | Modern mgcv methods         |
| **CI Coverage**             | Non-zero CIs for all levels  | Non-zero CIs for all levels |
| **Use Case**                | Traditional fisheries users  | Modern GAM best practices   |

### Choosing the Right Method

- **Use coefficient-based (default)** when:
  - Working with traditional fisheries applications
  - Wanting established results with Francis method guarantees
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

### Plot standardisation

Shows the unstandardised and standardised indices for the focus term. Option to show only standardised for cleaner presentation.

```r
plot_standardisation(gi)                          # Default: both indices with legend
plot_standardisation(gi, show_unstandardised = FALSE)  # Standardised index only
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

### Plot residual diagnostics

Generate comprehensive residual diagnostic plots for model validation. Supports both traditional 4-panel diagnostics and violin plots showing residual distributions by focus term levels. Includes automatic numeric axis conversion for sequential focus terms (e.g., years) and optional faceting by grouping variables.

```r
plot_residuals(gi)                                    # Default: standard 4-panel diagnostics
plot_residuals(gi, type = "violin")                  # Violin plots by focus levels
plot_residuals(gi, type = "both")                    # Both standard and violin plots
plot_residuals(gi, residual_type = "deviance")       # Specify residual type
plot_residuals(gi, type = "violin", by = "area")     # Faceted violin plots by grouping variable
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

# Or using the generic plot method:
plot(gi, type = "distribution", term = "s(temp)")
```

### Plot random effect diagnostics (QQ plot, histogram, caterpillar plot, points)

```r
plot_re(gi, term = "site", re_type = "qq")
plot_re(gi, term = "site", re_type = "hist")
plot_re(gi, term = "site", re_type = "caterpillar")
plot_re(gi, term = "site", re_type = "points")
```

### Analyse residual patterns for model improvement

Identify potential covariates that should be included in the model by analysing patterns in residuals.

```r
# Basic residual pattern analysis
rpa <- analyse_residual_patterns(gi)

# View results and recommendations
print(rpa)
summary(rpa)

# Customised analysis
rpa_custom <- analyse_residual_patterns(gi,
  candidate_vars = c("bottom_temp", "salinity", "moon_phase", "vessel_length"),
  exclude_vars = c("trip_id", "date"),
  residual_type = "deviance",
  smooth_terms = TRUE,
  significance_level = 0.01,
  min_r_squared = 0.05,
  plot_results = TRUE
)

# Access detailed results
linear_results <- rpa_custom$linear_results
smooth_results <- rpa_custom$smooth_results
significant_vars <- rpa_custom$significant_vars

# View diagnostic plots for significant variables
if (length(rpa_custom$plots) > 0) {
  # Display first plot
  print(rpa_custom$plots[[1]])
  
  # Access specific variable plot
  if ("bottom_temp" %in% names(rpa_custom$plots)) {
    print(rpa_custom$plots$bottom_temp)
  }
}

# Example of acting on recommendations:
# If the analysis identifies bottom_temp as significant with non-linear pattern:
# Re-fit model including the recommended term
mod_improved <- gam(response ~ s(depth) + s(temp) + s(bottom_temp) + year, 
                    family = gamma(link = "log"), data = data)

# Re-run analysis to check improvement
gi_improved <- calculate_influence(gam_influence(mod_improved, focus = "year"))
rpa_check <- analyse_residual_patterns(gi_improved)
```

### Delta-GLM Analysis for Fisheries Data

Combine binomial (catch probability) and positive catch GAMs for comprehensive fisheries indices.

```r
# Fit separate models for catch probability and positive catch amounts
mod_binomial <- gam(presence ~ year + s(depth), family = binomial(), data = data)
mod_positive <- gam(positive_catch ~ year + s(depth), family = Gamma(link="log"), 
                    data = data[data$presence == 1, ])

# Create influence objects
gi_binom <- calculate_influence(gam_influence(mod_binomial, focus = "year"))
gi_positive <- calculate_influence(gam_influence(mod_positive, focus = "year"))

# Combine indices using different methods
gi_combined <- combine_indices(gi_binom, gi_positive, method = "multiplicative")  # Default
gi_combined_geo <- combine_indices(gi_binom, gi_positive, method = "geometric")
gi_combined_arith <- combine_indices(gi_binom, gi_positive, method = "arithmetic")

# Different confidence interval methods
gi_combined_boot <- combine_indices(gi_binom, gi_positive, 
                                   confidence_method = "bootstrap", bootstrap_n = 1000)
gi_combined_indep <- combine_indices(gi_binom, gi_positive, 
                                    confidence_method = "independent")

# Visualise combined results
plot(gi_combined, type = "combined")    # Combined index only
plot(gi_combined, type = "components")  # Individual components
plot(gi_combined, type = "comparison")  # All indices together

# Extract combined indices
combined_indices <- extract_indices(gi_combined)
summary(gi_combined)
```

### Compare focus effects across multiple groups

Compare the focus term effects across different groups or datasets.

```r
# Create influence objects for different areas
gi_north <- calculate_influence(gi, subset_var = "area", subset_value = "North")
gi_south <- calculate_influence(gi, subset_var = "area", subset_value = "South")
gi_east <- calculate_influence(gi, subset_var = "area", subset_value = "East")

# Compare focus effects across areas
comparison <- compare_focus_by_groups(
  list(gi_north, gi_south, gi_east), 
  group_names = c("North", "South", "East")
)

# Visualise and summarise comparisons
plot(comparison)
summary(comparison)
```

### Analyse focus term by group levels

Perform influence analysis separately for each level of a grouping variable.

```r
# Analyse year effects separately by area
area_analysis <- analyse_focus_by_group(gi, group_var = "area")

# Extract results for each area
area_results <- extract_indices(area_analysis)
summary(area_analysis)
```

### Advanced residual analysis

Generate implied residual plots for model diagnostic purposes.

```r
# Advanced residual diagnostics
plot_implied_residuals(gi, var)

# Can be combined with other diagnostic plots
plot_residuals(gi, residual_type = "combined")
```

### Utility functions

Calculate geometric mean for rescaling or other statistical purposes.

```r
# Calculate geometric mean (handles zeros and negative values appropriately)
geom_mean_value <- geometric_mean(c(1, 2, 4, 8), na.rm = TRUE)

# Use in custom rescaling
gi_custom <- calculate_influence(gi, rescale_method = "custom",
                                custom_rescale_value = geom_mean_value)
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

### Combine indices for delta-GLM analysis

Combine binomial and positive catch GAMs for comprehensive fisheries CPUE indices.

```r
# Combine indices from binomial and positive models
gi_combined <- combine_indices(gi_binomial, gi_positive)
gi_combined <- combine_indices(gi_binomial, gi_positive, method = "geometric")
gi_combined <- combine_indices(gi_binomial, gi_positive, confidence_method = "bootstrap")

# Extract and visualise combined results
plot(gi_combined)
summary(gi_combined)
combined_indices <- extract_indices(gi_combined)
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

## Complete Function Reference

### Core Analysis Functions

| Function                      | Purpose                                    | Usage Example                                |
| ----------------------------- | ------------------------------------------ | -------------------------------------------- |
| `gam_influence()`             | Initialize influence analysis object       | `gi <- gam_influence(model, focus = "year")` |
| `calculate_influence()`       | Compute all indices and influence metrics  | `gi <- calculate_influence(gi)`              |
| `analyse_residual_patterns()` | Identify missing covariates from residuals | `analyse_residual_patterns(gi)`              |
| `extract_indices()`           | Extract standardised results as data frame | `indices <- extract_indices(gi)`             |
| `get_terms()`                 | Get model term names or full expressions   | `get_terms(gi, full = TRUE)`                 |
| `r2()`                        | Extract model progression statistics       | `r2(gi)`                                     |

### Plotting Functions

The package provides both specific plotting functions and a generic `plot()` method. The generic method supports the following types:
- `"stan"`: Standardisation plot (`plot_standardisation()`)
- `"step"`: Stepwise plot (`plot_stepwise_index()`) 
- `"influ"`: Influence plot (`plot_term_influence()`)
- `"cdi"`: CDI plot (`plot_cdi()`, requires `term` argument)
- `"distribution"`: Data distribution plot (`plot_term_distribution()`, requires `term` argument)
- `"all"`: Combined step and influence plots (`plot_step_and_influence()`)

| Function                    | Purpose                                            | Usage Example                                  |
| --------------------------- | -------------------------------------------------- | ---------------------------------------------- |
| `plot_standardisation()`    | Compare unstandardised vs standardised indices     | `plot_standardisation(gi)`                     |
| `plot_stepwise_index()`     | Show index changes as terms are added              | `plot_stepwise_index(gi)`                      |
| `plot_term_influence()`     | Display influence of each term on focus            | `plot_term_influence(gi)`                      |
| `plot_step_and_influence()` | Combined stepwise and influence plots              | `plot_step_and_influence(gi)`                  |
| `plot_cdi()`                | Coefficient-Distribution-Influence plot for a term | `plot_cdi(gi, term = "s(temp)")`               |
| `plot_terms()`              | Plot predicted effects for model terms             | `plot_terms(gi, term = "s(depth)")`            |
| `plot_residuals()`          | Comprehensive residual diagnostic plots            | `plot_residuals(gi, type = "violin")`          |
| `plot_re()`                 | Random effects diagnostics                         | `plot_re(gi, term = "site", re_type = "qq")`   |
| `plot_term_distribution()`  | Data distribution for a specific term              | `plot_term_distribution(gi, term = "s(temp)")` |
| `plot_implied_residuals()`  | Advanced residual analysis                         | `plot_implied_residuals(gi)`                   |

### Delta-GLM Functions

| Function                    | Purpose                                    | Usage Example                                    |
| --------------------------- | ------------------------------------------ | ------------------------------------------------ |
| `combine_indices()`         | Combine binomial and positive catch models | `combine_indices(gi_binom, gi_positive)`         |
| `compare_focus_by_groups()` | Compare focus effects across groups        | `compare_focus_by_groups(gi_list, group_names)`  |
| `analyse_focus_by_group()`  | Analyse focus term by group levels         | `analyse_focus_by_group(gi, group_var = "area")` |

### Utility Functions

| Function           | Purpose                                | Usage Example                          |
| ------------------ | -------------------------------------- | -------------------------------------- |
| `geometric_mean()` | Calculate geometric mean for rescaling | `geometric_mean(values, na.rm = TRUE)` |

### S3 Methods

| Method                             | Purpose                                    | Usage Example                              |
| ---------------------------------- | ------------------------------------------ | ------------------------------------------ |
| `plot.gam_influence()`             | Generic plot method with type selection    | `plot(gi, type = "cdi", term = "s(temp)")` |
| `plot.gam_influence_combined()`    | Plot method for combined delta-GLM objects | `plot(gi_combined, type = "comparison")`   |
| `summary.gam_influence()`          | Formatted summary of analysis results      | `summary(gi)`                              |
| `summary.gam_influence_combined()` | Summary for combined objects               | `summary(gi_combined)`                     |

---

## Testing and Validation

### Comprehensive Function Testing

The package includes comprehensive testing for all user-facing functions:

**Core Analysis Functions:**
- ✅ `gam_influence()` - Object creation with all GLM families
- ✅ `calculate_influence()` - Influence calculations across family types
- ✅ `analyse_residual_patterns()` - Residual pattern analysis for model improvement
- ✅ `extract_indices()` - Data extraction and formatting
- ✅ `get_terms()` - Term identification and parsing
- ✅ `r2()` - Model progression statistics

**Plotting Functions:**
- ✅ `plot_standardisation()` - Index comparison plots
- ✅ `plot_stepwise_index()` - Stepwise progression visualization
- ✅ `plot_term_influence()` - Term influence plots  
- ✅ `plot_step_and_influence()` - Combined plotting
- ✅ `plot_cdi()` - Coefficient-Distribution-Influence plots
- ✅ `plot_terms()` - Term effect visualization (all smoother types)
- ✅ `plot_residuals()` - Comprehensive residual diagnostics
- ✅ `plot_re()` - Random effects diagnostics
- ✅ `plot_term_distribution()` - Data distribution plots
- ✅ `plot_implied_residuals()` - Advanced residual analysis
- ✅ `plot.gam_influence()` - Generic plot method with type selection

**Delta-GLM Functions:**
- ✅ `combine_indices()` - Multiple combination methods
- ✅ `compare_focus_by_groups()` - Multi-group comparisons
- ✅ `analyse_focus_by_group()` - Group-wise analysis

**Utility Functions:**
- ✅ `geometric_mean()` - Statistical calculations

**S3 Methods:**
- ✅ All `plot.*()` methods with type specifications
- ✅ All `summary.*()` methods for formatted output

### Test Coverage Areas

**Family Support Testing:**
- Gaussian, binomial, gamma, Poisson, and Tweedie distributions
- Automatic family detection and appropriate statistical methods
- Log-link and identity-link model handling
- Pre-logged response variable detection

**Smoother Type Testing:**
- `s()` - Standard smoothers (splines, random effects)
- `te()` - Tensor product smoothers  
- `ti()` - Tensor product interactions
- `t2()` - Scaled tensor products
- `by=` terms - Factor-smoother interactions
- Random effects with `bs="re"`

**Advanced Feature Testing:**
- Subset analysis with interaction terms
- Coefficient-based vs prediction-based confidence intervals
- Bootstrap confidence intervals for combined models
- Custom rescaling methods and confidence levels
- Multi-variable interaction handling in CDI plots

### Validation Examples

```r
# Test basic functionality across families
library(mgcv)
library(gamInflu)

# Create test data
set.seed(123)
n <- 200
test_data <- data.frame(
  year = factor(rep(2015:2019, each = 40)),
  depth = runif(n, 10, 100),
  temp = rnorm(n, 15, 3),
  area = factor(sample(c("North", "South"), n, replace = TRUE))
)

# Generate responses for different families
test_data$y_gaussian <- with(test_data, 2 + 0.1*as.numeric(year) + sin(depth/20) + rnorm(n, 0, 0.5))
test_data$y_binomial <- rbinom(n, 1, plogis(test_data$y_gaussian - 2))
test_data$y_gamma <- rgamma(n, shape = 2, rate = 2/exp(test_data$y_gaussian))
test_data$y_poisson <- rpois(n, exp(test_data$y_gaussian - 1))

# Test all families
families <- list(
  gaussian = gaussian(),
  binomial = binomial(), 
  gamma = Gamma(link = "log"),
  poisson = poisson()
)

# Validate each family
for(fam_name in names(families)) {
  cat("Testing", fam_name, "family...\n")
  
  # Fit model
  response_var <- paste0("y_", fam_name)
  model <- gam(
    formula = as.formula(paste(response_var, "~ s(depth) + s(temp) + year")),
    data = test_data,
    family = families[[fam_name]]
  )
  
  # Test complete workflow
  gi <- gam_influence(model, focus = "year")
  gi <- calculate_influence(gi)
  
  # Test all plot functions
  plot_standardisation(gi)
  plot_stepwise_index(gi)  
  plot_term_influence(gi)
  plot_cdi(gi, term = "s(depth)")
  plot_term_distribution(gi, term = "s(temp)")
  
  # Test generic plot method with all types
  plot(gi, type = "stan")
  plot(gi, type = "step") 
  plot(gi, type = "influ")
  plot(gi, type = "cdi", term = "s(depth)")
  plot(gi, type = "distribution", term = "s(temp)")
  plot(gi, type = "all")
  plot_residuals(gi)
  
  # Test residual pattern analysis
  rpa <- analyse_residual_patterns(gi)
  stopifnot(class(rpa) == "residual_pattern_analysis")
  stopifnot(all(c("linear_results", "recommendations", "analysis_info") %in% names(rpa)))
  
  # Test data extraction
  indices <- extract_indices(gi)
  stopifnot(nrow(indices) == 5)  # 5 years
  stopifnot(all(c("index", "cv", "lower_CI", "upper_CI") %in% names(indices)))
  
  cat("✅", fam_name, "family validation complete\n\n")
}

# Test delta-GLM workflow
mod_binom <- gam(y_binomial ~ s(depth) + year, data = test_data, family = binomial())
mod_gamma <- gam(y_gamma ~ s(depth) + year, data = test_data[test_data$y_binomial == 1,], family = Gamma(link="log"))

gi_binom <- calculate_influence(gam_influence(mod_binom, focus = "year"))
gi_gamma <- calculate_influence(gam_influence(mod_gamma, focus = "year"))
gi_combined <- combine_indices(gi_binom, gi_gamma)

plot(gi_combined, type = "comparison")
summary(gi_combined)

cat("✅ All function validation tests passed\n")
```

---

## Citation

To cite the **gamInflu** package in publications:

> Dunn, A. (2025). *gamInflu: Influence Analysis Tools for Generalised Additive Models in R*. R package version 0.2.0. Available at: https://github.com/alistairdunn1/gamInflu

Please also cite the foundational methodology:

> Bentley, N., Kendrick, T. H., Starr, P. J., & Breen, P. A. (2012). Influence plots and metrics: tools for better understanding fisheries catch-per-unit-effort standardisations. *ICES Journal of Marine Science*, 69(1), 84–88. https://doi.org/10.1093/icesjms/fsr174

Generate a citation in R:

```r
citation("gamInflu")
```

---

## Support

For questions, bug reports, or feature requests, please contact the maintainer.

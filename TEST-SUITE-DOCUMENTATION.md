# gamInflu Package - Comprehensive Unit Test Suite

## Overview

This document describes the comprehensive unit test suite created for the gamInflu package, which provides extensive testing for GAM (Generalized Additive Model) influence analysis across different GLM families and the `islog` parameter.

## Test Coverage

### Core Functionality Tests (`test-simple-comprehensive.R`)

The main comprehensive test covers all key functionality:

1. **Lognormal Models with Pre-logged Data**
   - Tests `islog = FALSE` (default handling of log-scaled data)
   - Tests `islog = TRUE` (anti-logged results using geometric mean)
   - Verifies proper exp() transformation behavior

2. **Gamma Family with Log Link**
   - Tests `islog = FALSE` and `islog = TRUE`
   - Both use geometric mean internally (family-appropriate)
   - Verifies positive, finite results

3. **Binomial Family**
   - Tests `islog = FALSE` with probability scale preservation
   - Uses raw rescaling method for binomial family
   - Verifies finite indices on appropriate scale

4. **Gaussian Family**
   - Tests standard Gaussian family with `islog = FALSE`
   - Uses arithmetic mean rescaling
   - Verifies basic functionality

### GLM Family Support

The package correctly handles all major GLM families:

- **Gaussian**: `family = gaussian()`
  - Auto-selects arithmetic mean when `islog = FALSE`
  - Auto-selects geometric mean when `islog = TRUE`
  
- **Gamma**: `family = Gamma(link = "log")`
  - Always uses geometric mean (appropriate for positive data)
  - Handles both `islog = TRUE/FALSE` properly
  
- **Binomial**: `family = binomial()`
  - Uses raw rescaling to preserve probability scale (0-1)
  - Appropriate for proportion/presence-absence data
  
- **Poisson**: `family = poisson()` (tested in extended suite)
  - Uses geometric mean for count data
  
- **Tweedie**: `family = tw(link = "log")` (if package available)
  - Compound Poisson-Gamma distribution support

### islog Parameter Testing

The `islog` parameter is comprehensively tested:

1. **When `islog = FALSE`** (default):
   - Data treated as on original scale
   - Uses family-appropriate mean type
   - No exp() transformation applied

2. **When `islog = TRUE`**:
   - Data treated as log-transformed
   - Applies exp() transformation during calculation
   - Results returned on anti-logged scale
   - Particularly important for pre-logged catch data

### Mean Type and Rescaling Method Testing

The package supports multiple rescaling approaches:

- **`arithmetic_mean`**: Standard arithmetic averaging
- **`geometric_mean`**: Appropriate for log-normal data and positive values
- **`raw`**: No rescaling, preserves original scale (important for binomial)
- **`auto`**: Intelligent selection based on family and islog parameter

### Plotting Function Tests

All major plotting functions are tested:

- `plot_standardisation()`: Index comparison plots
- `plot_stepwise_index()`: Model progression visualization  
- `plot_step_and_influence()`: Combined plotting
- `plot_terms()`: Term effect visualization
- `plot_cdi()`: Coefficient-Distribution-Influence plots
- `plot_residuals()`: Comprehensive residual diagnostics

### Utility Function Tests

Core utility functions are verified:

- `geometric_mean()`: Geometric mean calculation
- `get_terms()`: Term identification from GAM models
- `extract_indices()`: Data extraction with proper formatting
- `r2()`: Model progression R-squared values (returns data.frame)
- S3 methods: `summary()`, `print()` for all object types

### Advanced Feature Tests

Advanced functionality is tested where applicable:

- **Residual Pattern Analysis**: `analyse_residual_patterns()` 
  - Identifies potential missing variables
  - Provides model improvement recommendations
  - Handles cases with no suitable candidates gracefully

## Running the Tests

### Quick Test Suite

```r
# Run comprehensive functionality test
source("run-comprehensive-unit-tests.R")
```

### Individual Test Components  

```r
# Load required packages
library(testthat)
library(gamInflu) 
library(mgcv)

# Run specific test file
test_file("gamInflu/tests/testthat/test-simple-comprehensive.R")
```

### Using testthat Framework

```r
# Run all tests in testthat format
library(testthat)
test_dir("gamInflu/tests/testthat")
```

## Test Results Summary

✅ **All Core Tests Pass**

The comprehensive test suite validates:

- ✅ Lognormal data handling with `islog = TRUE/FALSE`
- ✅ All GLM family support (Gaussian, Gamma, Binomial, Poisson, Tweedie*)
- ✅ Proper `islog` parameter transformation behavior
- ✅ Family-appropriate mean type selection
- ✅ All plotting functions work correctly
- ✅ Data extraction and utility functions
- ✅ S3 methods (summary, print) for all object types
- ✅ Advanced features (residual pattern analysis)
- ✅ Error handling and edge cases

*Tweedie requires the `tweedie` package

## Key Testing Insights

1. **islog Parameter Behavior**: 
   - The relationship between `islog=FALSE` and `islog=TRUE` results depends on the mean type used
   - When the same mean type is used, `islog=TRUE` applies exp() transformation
   - Auto-selection may choose different mean types, leading to different relationships

2. **Family-Specific Behavior**:
   - Binomial models preserve probability scale with raw rescaling
   - Gamma models always use geometric mean (appropriate for positive data)
   - Gaussian models choose mean type based on islog parameter

3. **Robust Error Handling**:
   - Package handles edge cases gracefully
   - Provides informative messages about family detection and method selection
   - Fails safely when required conditions aren't met

## Files in Test Suite

- `gamInflu/tests/testthat/test-simple-comprehensive.R`: Main comprehensive tests
- `gamInflu/tests/testthat/helper-test-data.R`: Test data generation functions
- `run-comprehensive-unit-tests.R`: Standalone test runner
- `gamInflu/tests/testthat.R`: testthat integration file

This comprehensive test suite ensures the gamInflu package works reliably across all supported GLM families and parameter combinations, providing confidence for users working with diverse ecological and fisheries datasets.

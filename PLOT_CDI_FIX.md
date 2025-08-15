# Fix for plot_cdi Error: Influence Data Matching Issue

## Problem Description

The `plot_cdi()` function was throwing warnings:
```
Warning messages:
1: In min(influ_data$influence, na.rm = TRUE) :
  no non-missing arguments to min; returning Inf
2: In max(influ_data$influence, na.rm = TRUE) :
  no non-missing arguments to max; returning -Inf
```

When called with:
```r
p1 <- plot_cdi(r1, 4, re_type = "points")
```

## Root Cause

The issue was in the `subplot_influence()` function where influence data was being filtered incorrectly:

**BEFORE (broken):**
```r
# Extract variable names from term (e.g., "vessel" from "s(vessel, bs = 're')")
term_vars_in_model <- term_vars[term_vars %in% model_terms]

# Try to match variable names against full term names - THIS FAILS!
influ_data <- subset(obj$calculated$influences, term %in% term_vars_in_model)
```

**Problem:** 
- `term_vars_in_model` contains variable names like `"vessel"`
- `obj$calculated$influences$term` contains full term expressions like `"s(vessel, bs = \"re\")"`
- The `%in%` comparison fails, returning no rows
- Empty dataframe leads to NA values after transformations
- `min()` and `max()` get no valid values, causing warnings

## Solution

**AFTER (fixed):**
```r
# Match against the full term name directly
influ_data <- obj$calculated$influences[obj$calculated$influences$term == term, ]

# Apply transformations
influ_data$influence <- exp(influ_data$influence)
influ_data$influence <- influ_data$influence / mean(influ_data$influence)

# Safety checks
if (nrow(influ_data) == 0) {
  stop("No influence data found for term: ", term, ". Available terms: ", 
       paste(unique(obj$calculated$influences$term), collapse = ", "))
}

if (all(is.na(influ_data$influence))) {
  warning("All influence values are NA for term: ", term)
  influ_data$influence <- rep(1, nrow(influ_data))  # Default to no influence
}
```

## Key Changes

1. **Direct term matching:** Instead of extracting variable names and trying to match them, we now match the full term name directly
2. **Explicit transformations:** Removed dplyr pipe operations to avoid potential namespace issues  
3. **Safety checks:** Added validation to catch cases where no data is found or all values are NA
4. **Better error messages:** Provide helpful information about available terms when matching fails

## Files Modified

- `gamInflu/R/subplot_influence.R`: Fixed the term matching logic and added safety checks

This fix ensures that `plot_cdi()` can correctly find and plot influence data for any valid model term, including random effects like `s(vessel, bs = "re")`.

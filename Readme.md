# gamInflu

**gamInflu** provides influence analysis tools for Generalized Additive Models (GAMs) fitted with the `mgcv` package in R. It supports all smoother types (`s()`, `te()`, `ti()`, `t2()`, and `by=` terms), and generates stepwise index plots, term influence plots, coefficient-distribution-influence (CDI) plots, and more to help you understand model structure and the influence of each term.

---

## Installation

```r
# If you have the source directory:
devtools::install_local("path/to/gamInflu")
# Or, if you have a tar.gz:
install.packages("gamInflu_0.1.0.tar.gz", repos = NULL, type = "source")
```

---

## Quick Start

```r
library(mgcv)
library(gamInflu)

# Fit a GAM model
mod <- gam(y ~ year + s(temp) + s(site, bs="re"), data = mydata, family = poisson())

# Create a gam_influence object
gi <- gam_influence(mod, focus = "year", data = mydata)

# Calculate influence metrics
gi <- calculate_influence(gi)

# Plot standardisation
plot_standardisation(gi)

# Plot stepwise index
plot_stepwise_index(gi)

# Plot term influence
plot_term_influence(gi)

# Plot combined step and influence
plot_step_and_influence(gi)

# Plot coefficient-distribution-influence for a term
plot_cdi(gi, term = "s(temp)")
# or reference the term by its order in the GAM
plot_cdi(gi, term = 2)

# Plot predicted effects for all terms
plot_terms(gi)

# Plot predicted effects for a single term (with random effect options)
plot_terms(gi, term = "s(site, bs='re')", type = "bar")
# or reference the term by its order in the GAM
plot_terms(gi, term = 3, type = "bar")

# Plot data distribution for a term
plot_term_distribution(gi, term = "s(temp)")
```

---

## Main Functions

- `gam_influence(model, focus, data = NULL, islog = NULL)`: Create a gam_influence object.
- `calculate_influence(obj)`: Calculate all influence metrics and indices.
- `plot_standardisation(obj)`: Plot unstandardised and standardised indices for the focus term.
- `plot_stepwise_index(obj, show_previous = FALSE)`: Plot how the index changes as each term is added.
- `plot_term_influence(obj)`: Plot the influence of each non-focus term.
- `plot_step_and_influence(obj)`: Combined stepwise and influence plot.
- `plot_cdi(obj, term)`: Coefficient-Distribution-Influence plot for a specified term.
- `plot_terms(obj, term = NULL, type = "point")`: Plot predicted effects for each term, including random effects and by-variable panels.
- `plot_term_distribution(obj, term)`: Plot the data distribution for a specified model term.
- `get_terms(obj, full = FALSE)`: Return the terms used in the model.

---

## Advanced Usage

### Plotting with by-variable panels

```r
plot_terms(gi, term = "s(temp, by=season)")
```

### Plotting random effects

```r
plot_terms(gi, term = "s(site, bs='re')", type = "violin")
```

### Extracting model terms

```r
get_terms(gi)            # Variable names
get_terms(gi, full=TRUE) # Full term expressions (e.g., s(temp), te(lon,lat))
```

---

## Requirements

- R >= 4.0.0
- mgcv
- ggplot2
- patchwork
- dplyr
- tidyr
- rlang

---

## Citation

To cite the **gamInflu** package in publications, please use:

> Dunn, A. (2025). *gamInflu: Influence Analysis Tools for Generalized Additive Models in R*. R package version 0.1.0. Available at: https://github.com/alistairdunn1/gamInflu

Please also cite:

> Bentley, N., Kendrick, T. H., Starr, P. J., & Breen, P. A. (2012). Influence plots and metrics: tools for better understanding fisheries catch-per-unit-effort standardizations. *ICES Journal of Marine Science*, 69(1), 84–88. https://doi.org/10.1093/icesjms/fsr176

You can generate a citation in R with:

```r
citation("gamInflu")
```

---

## Support

For questions, bug reports, or feature requests, please contact the maintainer


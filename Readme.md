# gamInflu

**gamInflu** provides influence analysis tools for Generalized Additive Models (GAMs) fitted with the `mgcv` package in R. It supports all smoother types (`s()`, `te()`, `ti()`, `t2()`, and `by=` terms), and generates stepwise index plots, term influence plots, coefficient-distribution-influence (CDI) plots, diagnostics for random effects, and more to help you understand model structure and the influence of each term.

---

## Installation

```r
# If you have the source directory:
devtools::install_local("path/to/gamInflu")
# Or, if you have a tar.gz:
install.packages("gamInflu_0.1.0.tar.gz", repos = NULL, type = "source")
```

---

## Main Functions and Examples

### Create a gam_influence object

Initializes the analysis object for your fitted GAM.

```r
gi <- gam_influence(mod, focus = "year", data = mydata)
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

Visualizes the predicted effects for each term in the model.

```r
plot_terms(gi)
```

### Plot predicted effects for a single term (including random effects and by-variable panels)

```r
plot_terms(gi, term = "s(site, bs='re')", type = "bar")
plot_terms(gi, term = "s(temp, by=season)")
plot_terms(gi, term = "te(lon, lat)")
```
or
```r
plot_terms(gi, term = 3, type = "bar")
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
plot_re(gi, term = "site", type = "qq")
plot_re(gi, term = "site", type = "hist")
plot_re(gi, term = "site", type = "caterpillar")
plot_re(gi, term = "site", type = "points")
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
plot_terms(gi, term = "s(site, bs='re')", type = "violin")
```

#### Plotting 2D tensor smooths

```r
plot_terms(gi, term = "te(lon, lat)")
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

> Bentley, N., Kendrick, T. H., Starr, P. J., & Breen, P. A. (2012). Influence plots and metrics: tools for better understanding fisheries catch-per-unit-effort standardizations. *ICES Journal of Marine Science*, 69(1), 84–88. https://doi.org/10.1093/icesjms/fsr174

You can generate a citation in R with:

```r
citation("gamInflu")
```

---

## Support

For questions, bug reports, or feature requests, please contact the maintainer

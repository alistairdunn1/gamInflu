# gamInflu: Influence Analysis for Generalized Additive Models

**gamInflu** is an R package for comprehensive influence analysis and visualization of Generalized Additive Models (GAMs) fitted with the `mgcv` package. It supports all smoother types (including `s()`, `te()`, `ti()`, `t2()`, and `by=` terms), and provides a suite of diagnostic and influence plots to help you understand model structure and the influence of each term.

---

## Features

- **Step plots**: Visualize how the focus term's fitted values change as each model term is added.
- **Influence plots**: Quantify and visualize the influence of each term on the focus variable.
- **CDI plots**: Coefficient-Distribution-Influence plots for deep dives into term effects.
- **Effect plots**: Visualize smooth and parametric effects for all model terms.
- **Support for random effects**: Specialized plots for random effect terms.
- **Easy export of influence summaries and per-level influences.**
- **Quick testing and validation functions.**

---

## Installation

```r
# Install from GitHub (if not on CRAN)
# install.packages("devtools")
devtools::install_github("yourusername/gamInflu")
```

---

## Quick Start Example

```r
library(mgcv)
library(gamInflu)

# Simulate data
set.seed(123)
n <- 300
dat <- data.frame(
  year = factor(sample(2010:2015, n, replace = TRUE)),
  depth = runif(n, 10, 200),
  month = factor(sample(month.abb, n, replace = TRUE), levels = month.abb),
  vessel = factor(sample(paste0("V", 1:8), n, replace = TRUE))
)
dat$y <- 5 + 0.02 * dat$depth + as.numeric(dat$month) * 0.1 +
  rnorm(n, sd = 0.5) + rnorm(length(levels(dat$vessel)))[dat$vessel] +
  rnorm(length(levels(dat$year)))[dat$year]

# Fit a GAM
m <- mgcv::gam(y ~ year + s(depth) + s(month, bs = "re") + s(vessel, bs = "re"), data = dat)
# Create a gam_influence object, focusing on "year"
gi <- gam_influence(m, focus = "year", data = dat)

# Print summary
print(gi)
summary(gi)
```

---

## Plotting Options

### 1. **Step Plot**
Shows how the focus term's fitted values change as each model term is added.

```r
plot_step(gi, panels = TRUE)         # Faceted by step
plot_step(gi, panels = FALSE)        # All steps in one panel
```

### 2. **Influence Plot**
Visualizes the influence of each term on the focus variable.

```r
plot_influence(gi, panels = TRUE)    # Faceted by term
plot_influence(gi, panels = FALSE)   # All terms in one panel
```

### 3. **Step + Influence Plot**
Side-by-side step and influence plots.

```r
plot_step_influence(gi)
```

### 4. **CDI Plot**
Coefficient-Distribution-Influence plot for a specific term.

```r
plot_cdi(gi, term = "s(depth)")
plot_cdi(gi, term = "s(month)")
```

### 5. **Effect Plots**
Visualize smooth and parametric effects for all model terms.

```r
plot_effects(gi)                     # All terms
plot_effects(gi, terms = c("s(depth)", "s(month)")) # Selected terms
```

---

## Exporting Results

Export a summary data frame of influence results:

```r
results <- export_results(gi)
head(results)
```

---

## Utility Functions

- **Get model summary:**  
  ```r
  get_model_summary(gi)
  ```
- **Get per-level influences:**  
  ```r
  get_per_level_influences(gi)
  ```
- **Reset calculations:**  
  ```r
  gi <- reset_calculations(gi)
  ```
- **Quick test example:**  
  ```r
  quick_test_example()
  ```
- **Test package installation:**  
  ```r
  test_package()
  ```

---

## Advanced Options

- **Random effect visualization:**  
  Use the `re_style` argument in plotting functions to choose between `"panel"`, `"qqnorm"`, `"caterpillar"`, or `"shrinkage"` for random effects.

- **Customizing plots:**  
  All plotting functions return `ggplot2` objects, so you can further customize them with additional `ggplot2` layers.

---

## Documentation

- All functions are documented with examples. See the package help pages, e.g.:
  ```r
  ?plot_step
  ?plot_influence
  ?plot_cdi
  ?plot_effects
  ```

---

## Contributing

Contributions, bug reports, and suggestions are welcome! Please open an issue or pull request on GitHub.

---

## License

GPL-3

---

## Acknowledgements

This package builds on the excellent `mgcv` and `ggplot2` packages.

---

**Enjoy exploring your GAMs with gamInflu!**
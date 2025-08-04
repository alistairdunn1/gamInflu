# gamInflu Vignettes

This directory contains vignettes for the **gamInflu** package:

## Available Vignettes

1. **quick-start.Rmd** - Basic introduction and simple usage examples
2. **getting-started.Rmd** - Comprehensive tutorial with working examples  
3. **gamInflu-introduction.Rmd** - Complete package overview based on README

## Building Vignettes

Vignettes are automatically built when the package is installed with:

```r
devtools::install("gamInflu", build_vignettes = TRUE)
```

Or when building from source:

```bash
R CMD build gamInflu
```

## Accessing Vignettes

After installation with vignettes:

```r
# List available vignettes
vignette(package = "gamInflu")

# Open specific vignettes
vignette("quick-start", package = "gamInflu")
vignette("getting-started", package = "gamInflu")
vignette("gamInflu-introduction", package = "gamInflu")
```

## Development Notes

- All vignettes use `eval=FALSE` for most code chunks to avoid requiring example data
- The `getting-started` vignette includes working examples with simulated data
- Vignettes are written in R Markdown format with knitr engine

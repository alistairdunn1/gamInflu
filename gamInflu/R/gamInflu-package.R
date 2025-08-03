#' @title gamInflu: Influence Analysis for Generalized Additive Models
#' @description Provides comprehensive influence analysis tools for Generalized Additive Models (GAMs)
#' fitted with the mgcv package. Supports all smoother types (s(), te(), ti(), t2(), by= terms) and
#' multiple GLM families including Gaussian, binomial, gamma, and Poisson distributions. Generates
#' stepwise index plots, term influence plots, coefficient-distribution-influence (CDI) plots,
#' diagnostics for random effects, and standardised indices to help understand model structure
#' and the influence of each term.
#'
#' @details
#' **Key Features:**
#' - **Multi-family Support**: Automatic detection and family-specific calculations for Gaussian,
#'   binomial, gamma, and Poisson GLM families
#' - **Comprehensive Plotting**: Stepwise indices, term influences, CDI plots, random effect diagnostics
#' - **Statistical Rigour**: Proper confidence intervals, geometric mean calculations, family-appropriate aggregation
#' - **Flexible Analysis**: Supports all mgcv smoother types and by-variable interactions
#'
#' **Main Functions:**
#' - `gam_influence()`: Create influence analysis object
#' - `calculate_influence()`: Compute indices and influence metrics with family support
#' - `plot_standardisation()`: Compare unstandardised vs standardised indices
#' - `plot_stepwise_index()`: Show model building progression
#' - `plot_term_influence()`: Visualize term-specific influences
#' - `plot_cdi()`: Coefficient-distribution-influence plots
#' - `plot_terms()`: Predicted effects visualisation
#' - `plot_re()`: Random effects diagnostics
#'
#' **GLM Family Support:**
#' - **Gaussian**: Traditional log-normal CPUE standardisation
#' - **Binomial**: Presence/absence, catch probability, proportion data
#' - **Gamma**: Positive continuous data (biomass, CPUE without zeros)
#' - **Poisson**: Count data (fish numbers, abundance indices)
#'
#' The package builds upon influence analysis concepts from Bentley et al. (2012)
#' "Influence plots and metrics: tools for better understanding fisheries catch-per-unit-effort
#' standardisations" (ICES Journal of Marine Science).
#'
#' @references
#' Bentley, N., Kendrick, T. H., Starr, P. J., & Breen, P. A. (2012).
#' Influence plots and metrics: tools for better understanding fisheries catch-per-unit-effort standardisations.
#' ICES Journal of Marine Science, 69(1), 84-88.
#'
#' @author Alistair Dunn
#' @keywords package GAM influence analysis mgcv
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL

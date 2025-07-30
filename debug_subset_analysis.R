# Simple debug test for subset analysis
library(devtools)
load_all("gamInflu")
library(mgcv)

cat("Debug test for subset analysis...\n")

# Create very simple test data
set.seed(123)
n <- 120
test_data <- data.frame(
  year = factor(rep(2018:2020, each = 40)),
  area = factor(rep(c("A", "B", "C"), length.out = n)),
  depth = runif(n, 10, 100),
  cpue = exp(rnorm(n, 2, 0.5))
)

cat("Test data structure:\n")
str(test_data)

cat("\nTesting single subset manually...\n")
subset_a <- test_data[test_data$area == "A", ]
cat("Subset A has", nrow(subset_a), "rows\n")

# Test fitting model to subset
tryCatch({
  mod_a <- mgcv::gam(cpue ~ s(depth) + year, data = subset_a, family = Gamma(link = "log"))
  cat("✓ Manual model fitting: SUCCESS\n")
  cat("Model family:", mod_a$family$family, "\n")
  
  # Test gamInflu on subset
  gi_a <- gam_influence(mod_a, focus = "year", data = subset_a)
  gi_a <- calculate_influence(gi_a)
  cat("✓ Manual gamInflu analysis: SUCCESS\n")
  
}, error = function(e) {
  cat("✗ Manual test FAILED:", e$message, "\n")
})

cat("\nTesting analyze_focus_by_group with simple data...\n")
tryCatch({
  result <- analyze_focus_by_group(
    model_formula = cpue ~ s(depth) + year,
    focus = "year",
    grouping_var = "area",
    data = test_data,
    family = Gamma(link = "log")
  )
  cat("✓ analyze_focus_by_group: SUCCESS\n")
  print(result)
}, error = function(e) {
  cat("✗ analyze_focus_by_group FAILED:", e$message, "\n")
})

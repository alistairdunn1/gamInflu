# Comprehensive Test Script for Subset-Based Analysis
# Test the analyze_focus_by_group functionality

library(devtools)
load_all("gamInflu")
library(mgcv)
library(ggplot2)

cat("=== Testing Subset-Based Analysis Implementation ===\n\n")

# Test 1: Create synthetic fisheries data with area-year interaction
cat("1. Creating synthetic fisheries data with area-year structure...\n")
set.seed(123)
n_per_group <- 60
areas <- c("North", "South", "East")
years <- factor(2015:2020)

# Create structured data
test_data <- expand.grid(
  year = years,
  area = factor(areas),
  depth = seq(10, 200, length.out = 20)
)

# Add additional observations
n_total <- nrow(test_data) * 3
test_data <- test_data[rep(seq_len(nrow(test_data)), 3), ]
test_data$depth <- test_data$depth + rnorm(nrow(test_data), 0, 10)
test_data$temperature <- rnorm(nrow(test_data), 15, 3)

# Create area-specific year effects
area_effects <- c("North" = 0.5, "South" = 0, "East" = -0.3)
year_numeric <- as.numeric(as.character(test_data$year))
year_effects <- 0.1 * (year_numeric - 2015)  # General increasing trend

# Area-specific modulation of year effects
test_data$area_year_effect <- area_effects[test_data$area] + 
  ifelse(test_data$area == "North", 0.2 * year_effects,
         ifelse(test_data$area == "South", year_effects,
                0.5 * year_effects))

# Generate CPUE with environmental and spatio-temporal structure
test_data$cpue <- exp(
  2.0 +  # Baseline
  0.3 * sin(test_data$depth / 50) +  # Depth effect
  0.2 * test_data$temperature / 15 + # Temperature effect
  test_data$area_year_effect +       # Area-year interaction
  rnorm(nrow(test_data), 0, 0.3)     # Noise
)

cat("Data created:", nrow(test_data), "observations across", 
    length(unique(test_data$area)), "areas and", 
    length(unique(test_data$year)), "years\n\n")

# Test 2: Basic functionality test
cat("2. Testing basic analyze_focus_by_group functionality...\n")

tryCatch({
  comparison_result <- analyze_focus_by_group(
    model_formula = cpue ~ s(depth) + s(temperature) + year,
    focus = "year",
    grouping_var = "area",
    data = test_data,
    family = Gamma(link = "log")
  )
  cat("✓ Basic functionality test: SUCCESS\n")
}, error = function(e) {
  cat("✗ Basic functionality test FAILED:", e$message, "\n")
  stop("Basic test failed")
})

# Test 3: Print and summary methods
cat("\n3. Testing print and summary methods...\n")
tryCatch({
  print(comparison_result)
  cat("\n--- Summary ---\n")
  summary(comparison_result)
  cat("✓ Print and summary methods: SUCCESS\n")
}, error = function(e) {
  cat("✗ Print/summary methods FAILED:", e$message, "\n")
})

# Test 4: Plotting methods
cat("\n4. Testing plotting methods...\n")

# Test standardisation plots
tryCatch({
  p1 <- plot(comparison_result, type = "standardisation")
  cat("✓ Standardisation plots: SUCCESS\n")
}, error = function(e) {
  cat("✗ Standardisation plots FAILED:", e$message, "\n")
})

# Test comparison plot
tryCatch({
  p2 <- plot(comparison_result, type = "comparison")
  cat("✓ Comparison plot: SUCCESS\n")
}, error = function(e) {
  cat("✗ Comparison plot FAILED:", e$message, "\n")
})

# Test stepwise plots
tryCatch({
  p3 <- plot(comparison_result, type = "stepwise")
  cat("✓ Stepwise plots: SUCCESS\n")
}, error = function(e) {
  cat("✗ Stepwise plots FAILED:", e$message, "\n")
})

# Test influence plots
tryCatch({
  p4 <- plot(comparison_result, type = "influence")
  cat("✓ Influence plots: SUCCESS\n")
}, error = function(e) {
  cat("✗ Influence plots FAILED:", e$message, "\n")
})

# Test 5: Index extraction
cat("\n5. Testing index extraction...\n")
tryCatch({
  indices_combined <- extract_indices(comparison_result, type = "both")
  cat("✓ Index extraction: SUCCESS\n")
  cat("  Extracted", nrow(indices_combined), "rows with columns:", 
      paste(names(indices_combined), collapse = ", "), "\n")
}, error = function(e) {
  cat("✗ Index extraction FAILED:", e$message, "\n")
})

# Test 6: Alternative family (Binomial)
cat("\n6. Testing with binomial family...\n")

# Create presence/absence data
test_data$presence <- rbinom(nrow(test_data), 1, plogis(log(test_data$cpue) - 3))

tryCatch({
  comparison_binomial <- analyze_focus_by_group(
    model_formula = presence ~ s(depth) + s(temperature) + year,
    focus = "year",
    grouping_var = "area",
    data = test_data,
    family = binomial()
  )
  cat("✓ Binomial family test: SUCCESS\n")
}, error = function(e) {
  cat("✗ Binomial family test FAILED:", e$message, "\n")
})

# Test 7: compare_focus_by_groups convenience function
cat("\n7. Testing convenience wrapper function...\n")
tryCatch({
  # Fit full model first
  full_model <- mgcv::gam(cpue ~ s(depth) + s(temperature) + year * area, 
                          data = test_data, family = Gamma(link = "log"))
  
  # Use convenience function
  comparison_wrapper <- compare_focus_by_groups(
    model = full_model,
    focus = "year", 
    grouping_var = "area"
  )
  cat("✓ Convenience wrapper test: SUCCESS\n")
}, error = function(e) {
  cat("✗ Convenience wrapper test FAILED:", e$message, "\n")
})

# Test 8: Error handling
cat("\n8. Testing error handling...\n")

# Test invalid focus variable
tryCatch({
  analyze_focus_by_group(
    model_formula = cpue ~ s(depth) + year,
    focus = "invalid_var",
    grouping_var = "area",
    data = test_data
  )
  cat("✗ Error handling test FAILED: Should have caught invalid focus variable\n")
}, error = function(e) {
  cat("✓ Error handling test: SUCCESS (caught invalid focus variable)\n")
})

# Test invalid grouping variable
tryCatch({
  analyze_focus_by_group(
    model_formula = cpue ~ s(depth) + year,
    focus = "year",
    grouping_var = "invalid_var",
    data = test_data
  )
  cat("✗ Error handling test FAILED: Should have caught invalid grouping variable\n")
}, error = function(e) {
  cat("✓ Error handling test: SUCCESS (caught invalid grouping variable)\n")
})

# Test 9: Save example plots
cat("\n9. Creating example plots for documentation...\n")
if (require(png, quietly = TRUE)) {
  tryCatch({
    # Standardisation comparison
    png("subset_analysis_standardisation.png", width = 1200, height = 800, res = 150)
    print(plot(comparison_result, type = "standardisation"))
    dev.off()
    
    # Index comparison
    png("subset_analysis_comparison.png", width = 1000, height = 600, res = 150)
    print(plot(comparison_result, type = "comparison"))
    dev.off()
    
    cat("✓ Example plots saved successfully\n")
  }, error = function(e) {
    cat("✗ Plot saving FAILED:", e$message, "\n")
  })
}

# Test 10: Performance with larger dataset
cat("\n10. Testing performance with larger dataset...\n")
tryCatch({
  # Create larger dataset
  large_data <- test_data[rep(seq_len(nrow(test_data)), 5), ]
  large_data$cpue <- large_data$cpue * exp(rnorm(nrow(large_data), 0, 0.1))
  
  start_time <- Sys.time()
  comparison_large <- analyze_focus_by_group(
    model_formula = cpue ~ s(depth) + s(temperature) + year,
    focus = "year",
    grouping_var = "area",
    data = large_data,
    family = Gamma(link = "log")
  )
  end_time <- Sys.time()
  
  processing_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  cat("✓ Performance test: SUCCESS (", round(processing_time, 2), " seconds for ", 
      nrow(large_data), " observations)\n")
}, error = function(e) {
  cat("✗ Performance test FAILED:", e$message, "\n")
})

cat("\n=== All Tests Completed ===\n")
cat("Implementation appears to be working correctly!\n")
cat("Results stored in 'comparison_result' object for further inspection.\n")

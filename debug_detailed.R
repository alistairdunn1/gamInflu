# More detailed debugging
library(devtools)
load_all("gamInflu")
library(mgcv)

# Create test data
set.seed(123)
test_data <- data.frame(
  year = factor(rep(2018:2020, each = 40)),
  area = factor(rep(c("A", "B", "C"), length.out = 120)),
  depth = runif(120, 10, 100),
  cpue = exp(rnorm(120, 2, 0.5))
)

cat("Testing step by step...\n")

# Test 1: Check if basic components work
subset_a <- test_data[test_data$area == "A", ]
cat("Subset A:", nrow(subset_a), "rows\n")

# Test 2: Check family object
gamma_family <- Gamma(link = "log")
cat("Gamma family object:\n")
print(str(gamma_family))

# Test 3: Try fitting model with explicit family
cat("\nFitting model with Gamma family...\n")
tryCatch({
  mod_test <- mgcv::gam(cpue ~ s(depth) + year, data = subset_a, family = gamma_family)
  cat("✓ Model fit successful\n")
  cat("Model family:", class(mod_test$family), "\n")
  
  # Store data explicitly
  mod_test$data <- subset_a
  
  # Test gamInflu
  gi_test <- gam_influence(mod_test, focus = "year", data = subset_a)
  cat("✓ gamInflu object created\n")
  
}, error = function(e) {
  cat("✗ Error:", e$message, "\n")
})

# Test 4: Step through analyze_focus_by_group logic manually
cat("\nTesting analyze_focus_by_group logic manually...\n")

model_formula <- cpue ~ s(depth) + year
focus <- "year"
grouping_var <- "area"
data <- test_data
groups <- levels(data[[grouping_var]])
family <- Gamma(link = "log")

cat("Groups:", paste(groups, collapse = ", "), "\n")
cat("Family:", family$family, "\n")

results <- list()
for (group in groups) {
  cat("Processing group:", group, "\n")
  subset_data <- data[data[[grouping_var]] == group, ]
  cat("  Subset size:", nrow(subset_data), "\n")
  
  tryCatch({
    cat("  Fitting model...\n")
    subset_model <- mgcv::gam(model_formula, data = subset_data, family = family)
    subset_model$data <- subset_data
    cat("  ✓ Model fitted\n")
    
    cat("  Creating gam_influence object...\n") 
    gi <- gam_influence(subset_model, focus = focus, data = subset_data)
    cat("  ✓ gam_influence created\n")
    
    cat("  Calculating influence...\n")
    gi <- calculate_influence(gi)
    cat("  ✓ Influence calculated\n")
    
    results[[group]] <- gi
    
  }, error = function(e) {
    cat("  ✗ Error:", e$message, "\n")
  })
}

cat("Results:", length(results), "successful fits\n")

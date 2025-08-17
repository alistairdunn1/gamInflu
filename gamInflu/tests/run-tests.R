# Comprehensive Test Runner for gamInflu Package
# Run this during development to test all major functionality
# For formal package testing, use devtools::test() or R CMD check

library(testthat)

# Try to load the package - during R CMD check, the package is already available
if (requireNamespace("devtools", quietly = TRUE) && interactive()) {
  cat("Loading package with devtools...\n")
  devtools::load_all(quiet = TRUE)
} else {
  # During R CMD check, just load the package normally
  cat("Loading package...\n")
  library(gamInflu, quietly = TRUE)
}

cat("=======================================================\n")
cat("         Running gamInflu Package Tests\n")
cat("=======================================================\n\n")

# Option 1: Run all testthat tests
cat("Running all testthat tests...\n")
cat("=============================\n")

# Find the correct path to testthat directory
testthat_path <- if (file.exists("testthat")) {
  "testthat" # Running from tests/ directory
} else if (file.exists("tests/testthat")) {
  "tests/testthat" # Running from package root
} else {
  # Last resort - find it relative to this script
  file.path(dirname(sys.frame(1)$ofile), "testthat")
}

if (!dir.exists(testthat_path)) {
  cat("Error: Could not find testthat directory\n")
  cat("Current working directory:", getwd(), "\n")
  cat("Looking for:", testthat_path, "\n")
  stop("testthat directory not found")
}

test_results <- test_dir(testthat_path, reporter = "summary")

cat("\n=======================================================\n")
cat("              Test Results Summary\n")
cat("=======================================================\n")

if (length(test_results) > 0) {
  cat("✅ All tests completed!\n")
  cat("Check output above for any warnings or failures.\n")
} else {
  cat("⚠️  No test results returned.\n")
}

cat("\nTo run individual test categories:\n")
cat("- Core functionality: test_active_file('tests/testthat/test-core-functionality.R')\n")
cat("- Plotting functions: test_active_file('tests/testthat/test-plotting-functions.R')\n")
cat("- islog parameter:    test_active_file('tests/testthat/test-islog-parameter.R')\n")
cat("- Utility functions:  test_active_file('tests/testthat/test-utility-functions.R')\n")
cat("- Family support:     test_active_file('tests/testthat/test-family-support.R')\n")
cat("- Edge cases:         test_active_file('tests/testthat/test-edge-cases.R')\n")

cat("\n=======================================================\n")
cat("Tests complete. Use devtools::test() for standard testing.\n")
cat("=======================================================\n")

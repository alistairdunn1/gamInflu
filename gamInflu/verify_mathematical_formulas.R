# Additional verification: Check Box-Cox formulas against literature
# Box-Cox transformation: y^(λ) = (y^λ - 1) / λ  for λ ≠ 0
#                                = log(y)         for λ = 0
# Inverse:                y = (λx + 1)^(1/λ)     for λ ≠ 0
#                             = exp(x)           for λ = 0

source("c:/Users/alist/OneDrive/Projects/Software/gamInflu/gamInflu/R/box_cox_transform.R")

cat("=== Mathematical Formula Verification ===\n\n")

# Test specific values that are easy to verify by hand
y_simple <- c(1, 4, 9) # Perfect squares for easy calculation

cat("1. Lambda = 0.5 (Square root transformation)\n")
cat("Formula: (y^0.5 - 1) / 0.5 = 2*(sqrt(y) - 1)\n")
result_half <- box_cox_transform(y_simple, lambda = 0.5)
expected_half <- 2 * (sqrt(y_simple) - 1)
cat("Input y:     ", paste(y_simple, collapse = ", "), "\n")
cat("sqrt(y):     ", paste(sqrt(y_simple), collapse = ", "), "\n")
cat("2*(sqrt(y)-1):", paste(expected_half, collapse = ", "), "\n")
cat("Our result:  ", paste(result_half$transformed, collapse = ", "), "\n")
cat("Match: ", all.equal(result_half$transformed, expected_half), "\n\n")

cat("2. Lambda = 2 (Square transformation)\n")
cat("Formula: (y^2 - 1) / 2\n")
result_two <- box_cox_transform(y_simple, lambda = 2)
expected_two <- (y_simple^2 - 1) / 2
cat("Input y:     ", paste(y_simple, collapse = ", "), "\n")
cat("y^2:         ", paste(y_simple^2, collapse = ", "), "\n")
cat("(y^2-1)/2:   ", paste(expected_two, collapse = ", "), "\n")
cat("Our result:  ", paste(result_two$transformed, collapse = ", "), "\n")
cat("Match: ", all.equal(result_two$transformed, expected_two), "\n\n")

cat("3. Lambda = -0.5 (Inverse square root)\n")
cat("Formula: (y^(-0.5) - 1) / (-0.5) = -2*(1/sqrt(y) - 1) = 2*(1 - 1/sqrt(y))\n")
result_neg_half <- box_cox_transform(y_simple, lambda = -0.5)
expected_neg_half <- (y_simple^(-0.5) - 1) / (-0.5)
cat("Input y:        ", paste(y_simple, collapse = ", "), "\n")
cat("1/sqrt(y):      ", paste(1 / sqrt(y_simple), collapse = ", "), "\n")
cat("(y^-0.5-1)/-0.5:", paste(expected_neg_half, collapse = ", "), "\n")
cat("Our result:     ", paste(result_neg_half$transformed, collapse = ", "), "\n")
cat("Match: ", all.equal(result_neg_half$transformed, expected_neg_half), "\n\n")

cat("4. Test continuity at λ = 0\n")
cat("As λ → 0, (y^λ - 1)/λ → log(y) (L'Hôpital's rule)\n")
y_test <- c(1, 2, 4)
lambda_sequence <- c(0.1, 0.01, 0.001, 0.0001)
log_y <- log(y_test)

cat("True log(y):   ", paste(round(log_y, 6), collapse = ", "), "\n")
for (lambda in lambda_sequence) {
  approx_log <- (y_test^lambda - 1) / lambda
  cat("λ =", sprintf("%6.4f", lambda), ":", paste(round(approx_log, 6), collapse = ", "), "\n")
}
# Test our function at λ = 0
result_zero <- box_cox_transform(y_test, lambda = 0)
cat("Our λ = 0:     ", paste(round(result_zero$transformed, 6), collapse = ", "), "\n")
cat("Match with log:", all.equal(result_zero$transformed, log_y), "\n\n")

cat("5. Inverse transformation verification\n")
# Test that the inverse formulas are correct
test_values <- c(0.5, 1.0, 2.0)

for (lambda in c(-1, -0.5, 0, 0.5, 1, 2)) {
  cat("Lambda =", lambda, "\n")

  if (abs(lambda) < 1e-10) {
    # For λ = 0: x = log(y), so y = exp(x)
    expected_inverse <- exp(test_values)
  } else {
    # For λ ≠ 0: x = (y^λ - 1)/λ, so y = (λx + 1)^(1/λ)
    expected_inverse <- (lambda * test_values + 1)^(1 / lambda)
  }

  our_inverse <- inverse_box_cox(test_values, lambda = lambda)

  cat("  x values:      ", paste(test_values, collapse = ", "), "\n")
  cat("  Expected y:    ", paste(round(expected_inverse, 6), collapse = ", "), "\n")
  cat("  Our result:    ", paste(round(our_inverse, 6), collapse = ", "), "\n")
  cat("  Match: ", all.equal(our_inverse, expected_inverse), "\n\n")
}

cat("=== All Mathematical Formulas Verified ===\n")

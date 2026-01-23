#!/usr/bin/env Rscript

# Simple test after package reload
library(RLDM)

cat("Checking pfilter...\n")
cat("Exists:", exists("pfilter"), "\n")
cat("Is generic:", isGeneric("pfilter"), "\n")

# Create minimal model
tmpl <- tmpl_stsp_full(1, 1, 1, sigma_L = "chol")
model <- r_model(tmpl, bpoles = 1, sd = 0.5)
data <- sim(model, n.obs = 10)
y <- data$y

cat("\nTrying pfilter...\n")
result <- try(pfilter(model, y, N_particles = 100))
if (inherits(result, "try-error")) {
  cat("Error:", result[1], "\n")
} else {
  cat("Success! Log-likelihood:", result$log_likelihood, "\n")
}
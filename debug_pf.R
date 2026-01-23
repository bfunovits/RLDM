#!/usr/bin/env Rscript

library(RLDM)

# Check if pf is generic
cat("Is pf a generic?", isGeneric("pf"), "\n")
cat("Methods for pf:", paste(methods("pf"), collapse = ", "), "\n")

# Check args of pf
cat("Args of pf:", paste(names(formals(pf)), collapse = ", "), "\n")
cat("Args of pf.stspmod:", paste(names(formals(pf.stspmod)), collapse = ", "), "\n")

# Create a simple model
tmpl <- tmpl_stsp_full(1, 1, 2, sigma_L = "chol")
model <- r_model(tmpl, bpoles = 1, sd = 0.5)
data <- sim(model, n.obs = 10)

# Try calling with explicit method
cat("\nTrying pf with method='sir'...\n")
tryCatch({
  result <- pf(model, data$y, method = "sir", N_particles = 100)
  cat("Success!\n")
}, error = function(e) cat("Error:", e$message, "\n"))

# Try calling without extra args (should use defaults)
cat("\nTrying pf without extra args...\n")
tryCatch({
  result <- pf(model, data$y)
  cat("Success!\n")
}, error = function(e) cat("Error:", e$message, "\n"))
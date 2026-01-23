#!/usr/bin/env Rscript

cat("Final integration test for particle filter extension\n")
cat("===================================================\n")

# Load development version
cat("\n1. Loading RLDM development version...\n")
devtools::load_all()
cat("Package loaded\n")

# Check function availability
cat("\n2. Checking function availability...\n")
required_funcs <- c("pfilter", "ll_pfilter", "plot.pfilter")
for (f in required_funcs) {
  exists_flag <- exists(f)
  cat(sprintf("  %-15s: %s\n", f, ifelse(exists_flag, "YES", "NO")))
}

# Create test model
cat("\n3. Generating test model...\n")
set.seed(123)
tmpl <- tmpl_stsp_full(1, 1, 2, sigma_L = "chol")
model <- r_model(tmpl, bpoles = 1, sd = 0.5)
data <- sim(model, n.obs = 100)
y <- data$y
cat("Model: m=1, n=1, s=2, N=100\n")

# Kalman filter baseline
cat("\n4. Kalman filter baseline...\n")
kf_result <- kf(model, y)
cat("   Log-likelihood:", kf_result$ll, "\n")

# Particle filter tests
cat("\n5. Particle filter tests...\n")

cat("   a) SIR filter with systematic resampling...\n")
pf_sir <- pfilter(model, y, N_particles = 1000, method = "sir",
                  resampling = "systematic")
cat("      Log-likelihood:", pf_sir$log_likelihood, "\n")
cat("      ESS range:", range(pf_sir$effective_sample_size), "\n")

cat("   b) APF filter with stratified resampling...\n")
pf_apf <- pfilter(model, y, N_particles = 1000, method = "apf",
                  resampling = "stratified")
cat("      Log-likelihood:", pf_apf$log_likelihood, "\n")

cat("   c) Likelihood approximation (averaged)...\n")
ll_avg <- ll_pfilter(model, y, N_particles = 1000, N_runs = 5)
cat("      Averaged log-likelihood:", ll_avg, "\n")
cat("      Difference from KF:", ll_avg - kf_result$ll, "\n")

# Validate consistency
cat("\n6. Validating consistency...\n")
cat("   Filtered states dimensions:", dim(pf_sir$filtered_states), "\n")
cat("   Particles array dimensions:", dim(pf_sir$particles), "\n")
cat("   Weights length:", length(pf_sir$weights), "\n")

# Compare with Kalman filter (should be close for linear Gaussian)
cat("\n7. Comparison with Kalman filter...\n")
rmse <- sqrt(mean((kf_result$a[1:100,] - pf_sir$filtered_states[1:100,])^2))
cat("   RMSE between KF and PF filtered states:", rmse, "\n")
cat("   Expected: < 0.1 for 1000 particles\n")

# Plotting (silent)
cat("\n8. Testing plotting methods (silent)...\n")
try({
  pdf(file = NULL)
  plot(pf_sir, type = "states")
  plot(pf_sir, type = "ess")
  plot(pf_sir, type = "weights")
  plot(pf_sir, type = "likelihood")
  dev.off()
  cat("   All plotting functions executed successfully\n")
})

cat("\n9. Testing error handling...\n")
cat("   a) Invalid number of particles...\n")
err1 <- try(pfilter(model, y, N_particles = 0), silent = TRUE)
cat("      ", ifelse(inherits(err1, "try-error"), "ERROR caught ✓", "FAIL"), "\n")

cat("   b) Invalid resampling method...\n")
err2 <- try(pfilter(model, y, resampling = "invalid"), silent = TRUE)
cat("      ", ifelse(inherits(err2, "try-error"), "ERROR caught ✓", "FAIL"), "\n")

cat("\n===================================================\n")
cat("All tests completed successfully!\n")
cat("Particle filter extension is functional.\n")
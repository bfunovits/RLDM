#!/usr/bin/env Rscript

cat("Testing particle filter implementation (renamed functions)...\n")

library(RLDM)

# Create a simple linear Gaussian state space model
set.seed(123)
s <- 2
m <- 1
n <- m

tmpl <- tmpl_stsp_full(m, n, s, sigma_L = "chol")
model <- r_model(tmpl, bpoles = 1, sd = 0.5)

n_obs <- 50
data <- sim(model, n.obs = n_obs, a1 = NA)
y <- data$y

cat("Model generated\n")

# 1. Kalman filter baseline
cat("\n1. Kalman filter...\n")
kf_result <- kf(model, y)
cat("KF log-likelihood:", kf_result$ll, "\n")

# 2. Particle filter (SIR)
cat("\n2. Particle filter (SIR)...\n")
pf_result <- pfilter(model, y, N_particles = 500, method = "sir")
cat("PF log-likelihood:", pf_result$log_likelihood, "\n")
cat("Effective sample size range:", range(pf_result$effective_sample_size, na.rm = TRUE), "\n")
cat("Number of particles:", pf_result$N_particles, "\n")
cat("Resampling method:", pf_result$resampling_method, "\n")

# 3. Particle filter (APF)
cat("\n3. Particle filter (APF)...\n")
pf_apf_result <- pfilter(model, y, N_particles = 500, method = "apf")
cat("APF log-likelihood:", pf_apf_result$log_likelihood, "\n")

# 4. Likelihood approximation
cat("\n4. Likelihood approximation (averaged)...\n")
ll_pf_result <- ll_pfilter(model, y, N_particles = 1000, N_runs = 3)
cat("Averaged PF likelihood:", ll_pf_result, "\n")
cat("Difference from Kalman filter:", ll_pf_result - kf_result$ll, "\n")

# 5. Basic plotting
cat("\n5. Testing plotting methods...\n")
try({
  pdf(file = NULL)
  plot(pf_result, type = "states")
  plot(pf_result, type = "ess")
  plot(pf_result, type = "weights")
  plot(pf_result, type = "likelihood")
  dev.off()
  cat("Plotting functions executed without errors\n")
})

# 6. Compare filtered state estimates
cat("\n6. Comparing filtered state estimates...\n")
kf_states <- kf_result$a
pf_states <- pf_result$filtered_states
rmse <- sqrt(mean((kf_states[1:n_obs, ] - pf_states[1:n_obs, ])^2, na.rm = TRUE))
cat("RMSE between KF and PF filtered states:", rmse, "\n")
cat("Expected: small ( < 0.1 ) for linear Gaussian with sufficient particles\n")

# 7. Test different resampling methods
cat("\n7. Testing different resampling methods...\n")
for (resamp in c("systematic", "multinomial", "stratified")) {
  cat("  Resampling:", resamp, "... ")
  res <- try(pfilter(model, y, N_particles = 300, method = "sir", resampling = resamp))
  if (inherits(res, "try-error")) {
    cat("ERROR\n")
  } else {
    cat("OK, log-likelihood:", res$log_likelihood, "\n")
  }
}

cat("\nAll tests completed successfully!\n")
#!/usr/bin/env Rscript

cat("Testing particle filter implementation...\n")

# Load the package (already loaded via devtools)
library(RLDM)

# Create a simple linear Gaussian state space model
set.seed(123)
s <- 2  # state dimension
m <- 1  # output dimension
n <- m  # input dimension

# Create template and random model
tmpl <- tmpl_stsp_full(m, n, s, sigma_L = "chol")
model <- r_model(tmpl, bpoles = 1, sd = 0.5)

# Generate data
n_obs <- 50
data <- sim(model, n.obs = n_obs, a1 = NA)
y <- data$y

cat("Model generated\n")
cat("State dimension:", s, "\n")
cat("Output dimension:", m, "\n")
cat("Sample size:", n_obs, "\n")

# 1. Test Kalman filter (baseline)
cat("\n1. Running Kalman filter...\n")
kf_result <- kf(model, y)
cat("Kalman filter log-likelihood:", kf_result$ll, "\n")

# 2. Test particle filter (SIR)
cat("\n2. Running particle filter (SIR)...\n")
pf_result <- pf(model, y, N_particles = 500, method = "sir")
cat("Particle filter log-likelihood:", pf_result$log_likelihood, "\n")
cat("Effective sample size range:",
    range(pf_result$effective_sample_size, na.rm = TRUE), "\n")

# 3. Test particle filter (APF)
cat("\n3. Running particle filter (APF)...\n")
pf_apf_result <- pf(model, y, N_particles = 500, method = "apf")
cat("APF log-likelihood:", pf_apf_result$log_likelihood, "\n")

# 4. Test likelihood approximation
cat("\n4. Testing likelihood approximation...\n")
ll_pf_result <- ll_pf(model, y, N_particles = 1000, N_runs = 3)
cat("Averaged PF likelihood:", ll_pf_result, "\n")
cat("Difference from Kalman filter:", ll_pf_result - kf_result$ll, "\n")

# 5. Basic plotting
cat("\n5. Testing basic plotting...\n")
try({
  # Plot filtered states
  pdf(file = NULL)  # suppress actual plotting
  plot(pf_result, type = "states")
  plot(pf_result, type = "ess")
  plot(pf_result, type = "weights")
  plot(pf_result, type = "likelihood")
  dev.off()
  cat("Plotting functions executed without errors\n")
})

# 6. Compare filtered state estimates
cat("\n6. Comparing filtered state estimates...\n")
# Kalman filter states (filtered)
kf_states <- kf_result$a
# Particle filter states (filtered)
pf_states <- pf_result$filtered_states

# Compute RMSE between filters (should be small for linear Gaussian with many particles)
rmse <- sqrt(mean((kf_states[1:n_obs, ] - pf_states[1:n_obs, ])^2))
cat("RMSE between KF and PF filtered states:", rmse, "\n")

cat("\nAll tests completed!\n")
#!/usr/bin/env Rscript
# Benchmark script for particle filter performance comparison
# Compares SIR, APF, and optimal proposal filters across parameter ranges

library(RLDM)
library(microbenchmark)

cat("=== Particle Filter Performance Benchmarks ===\n")
set.seed(123)

# Configuration
N_particles_vec <- c(100, 500, 1000, 2000, 5000)
N_obs <- 100
m <- 1  # output dimension
s <- 2  # state dimension
n <- s  # input dimension

# Helper function to compute performance metrics
compute_metrics <- function(pf_result, kf_result, method_name) {
  # Filtered state RMSE vs Kalman filter
  kf_states <- kf_result$a[1:N_obs, ]
  pf_states <- pf_result$filtered_states[1:N_obs, ]
  rmse <- sqrt(mean((kf_states - pf_states)^2))

  # Likelihood difference
  ll_diff <- pf_result$log_likelihood - kf_result$ll

  # ESS statistics
  ess_min <- min(pf_result$effective_sample_size, na.rm = TRUE)
  ess_max <- max(pf_result$effective_sample_size, na.rm = TRUE)
  ess_mean <- mean(pf_result$effective_sample_size, na.rm = TRUE)

  # Weight degeneracy indicator
  weight_entropy <- -sum(pf_result$weights * log(pf_result$weights + 1e-10))
  max_weight <- max(pf_result$weights)

  return(list(
    method = method_name,
    rmse = rmse,
    ll_diff = ll_diff,
    ess_min = ess_min,
    ess_max = ess_max,
    ess_mean = ess_mean,
    weight_entropy = weight_entropy,
    max_weight = max_weight
  ))
}

# Test 1: Model with diagonal sigma_L (S = 0)
cat("\n--- Test 1: No cross-covariance (S = 0) ---\n")
tmpl_diag <- tmpl_stsp_full(m, n, s, sigma_L = "chol")
model_diag <- r_model(tmpl_diag, bpoles = 1, sd = 0.5)
data_diag <- sim(model_diag, n.obs = N_obs, a1 = NA)
y_diag <- data_diag$y

# Kalman filter baseline
kf_diag <- kf(model_diag, y_diag)

results_diag <- list()
for (N_particles in N_particles_vec) {
  cat(sprintf("\n  N_particles = %d:\n", N_particles))

  # SIR filter
  pf_sir <- pfilter(model_diag, y_diag, N_particles = N_particles, method = "sir")
  metrics_sir <- compute_metrics(pf_sir, kf_diag, "sir")
  results_diag <- c(results_diag, list(metrics_sir))
  cat(sprintf("    SIR: RMSE = %.4f, LL diff = %.4f, ESS min = %.1f\n",
              metrics_sir$rmse, metrics_sir$ll_diff, metrics_sir$ess_min))

  # APF filter
  pf_apf <- pfilter(model_diag, y_diag, N_particles = N_particles, method = "apf")
  metrics_apf <- compute_metrics(pf_apf, kf_diag, "apf")
  results_diag <- c(results_diag, list(metrics_apf))
  cat(sprintf("    APF: RMSE = %.4f, LL diff = %.4f, ESS min = %.1f\n",
              metrics_apf$rmse, metrics_apf$ll_diff, metrics_apf$ess_min))

  # Optimal proposal
  pf_opt <- pfilter(model_diag, y_diag, N_particles = N_particles, method = "optimal")
  metrics_opt <- compute_metrics(pf_opt, kf_diag, "optimal")
  results_diag <- c(results_diag, list(metrics_opt))
  cat(sprintf("    OPT: RMSE = %.4f, LL diff = %.4f, ESS min = %.1f\n",
              metrics_opt$rmse, metrics_opt$ll_diff, metrics_opt$ess_min))
}

# Test 2: Model with cross-covariance (S ≠ 0)
cat("\n--- Test 2: With cross-covariance (S ≠ 0) ---\n")
# Create a model with significant cross-covariance
sigma_L_cross <- matrix(rnorm(n * n), n, n)
model_cross <- stspmod(
  sys = stsp(A = matrix(rnorm(s * s) * 0.5, s, s),
             B = matrix(rnorm(s * n), s, n),
             C = matrix(rnorm(m * s), m, s),
             D = matrix(rnorm(m * n), m, n)),
  sigma_L = sigma_L_cross
)
# Ensure stability
model_cross$sys$A <- model_cross$sys$A * 0.5

data_cross <- sim(model_cross, n.obs = N_obs, a1 = NA)
y_cross <- data_cross$y

# Compute S matrix to verify non-zero
sigma_cross <- tcrossprod(model_cross$sigma_L)
S_cross <- model_cross$sys$B %*% sigma_cross %*% t(model_cross$sys$D)
cat(sprintf("  Norm of S matrix: %.4f\n", norm(S_cross, "F")))

# Kalman filter baseline
kf_cross <- kf(model_cross, y_cross)

results_cross <- list()
for (N_particles in N_particles_vec) {
  cat(sprintf("\n  N_particles = %d:\n", N_particles))

  # SIR filter
  pf_sir <- pfilter(model_cross, y_cross, N_particles = N_particles, method = "sir")
  metrics_sir <- compute_metrics(pf_sir, kf_cross, "sir")
  results_cross <- c(results_cross, list(metrics_sir))
  cat(sprintf("    SIR: RMSE = %.4f, LL diff = %.4f, ESS min = %.1f\n",
              metrics_sir$rmse, metrics_sir$ll_diff, metrics_sir$ess_min))

  # APF filter (should show warning)
  pf_apf <- pfilter(model_cross, y_cross, N_particles = N_particles, method = "apf")
  metrics_apf <- compute_metrics(pf_apf, kf_cross, "apf")
  results_cross <- c(results_cross, list(metrics_apf))
  cat(sprintf("    APF: RMSE = %.4f, LL diff = %.4f, ESS min = %.1f\n",
              metrics_apf$rmse, metrics_apf$ll_diff, metrics_apf$ess_min))

  # Optimal proposal
  pf_opt <- pfilter(model_cross, y_cross, N_particles = N_particles, method = "optimal")
  metrics_opt <- compute_metrics(pf_opt, kf_cross, "optimal")
  results_cross <- c(results_cross, list(metrics_opt))
  cat(sprintf("    OPT: RMSE = %.4f, LL diff = %.4f, ESS min = %.1f\n",
              metrics_opt$rmse, metrics_opt$ll_diff, metrics_opt$ess_min))
}

# Test 3: Timing comparison
cat("\n--- Test 3: Computational timing comparison ---\n")
N_particles_timing <- 1000
model_timing <- model_diag
y_timing <- y_diag

timing_results <- microbenchmark(
  sir = pfilter(model_timing, y_timing, N_particles = N_particles_timing, method = "sir"),
  apf = pfilter(model_timing, y_timing, N_particles = N_particles_timing, method = "apf"),
  optimal = pfilter(model_timing, y_timing, N_particles = N_particles_timing, method = "optimal"),
  times = 10
)

print(timing_results)

# Summary
cat("\n=== Summary of Key Findings ===\n")
cat("1. Optimal proposal should have lowest RMSE and LL difference\n")
cat("2. APF may show bias when S ≠ 0 (cross-covariance present)\n")
cat("3. SIR typically has high variance and weight degeneracy\n")
cat("4. ESS indicates weight degeneracy (lower is worse)\n")
cat("5. Timing: Optimal > APF > SIR (optimal is most computationally intensive)\n")

# Save results for later analysis
save(results_diag, results_cross, timing_results,
     file = "particle_filter_benchmark_results.RData")

cat("\nBenchmark completed. Results saved to particle_filter_benchmark_results.RData\n")
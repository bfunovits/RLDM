#!/usr/bin/env Rscript
devtools::load_all()

cat("=== Optimal Proposal Validation ===\n")

set.seed(123)
m <- 1
s <- 2
n <- s
tmpl <- tmpl_stsp_full(m, n, s, sigma_L = "chol")
model <- r_model(tmpl, bpoles = 1, sd = 0.5)
data <- sim(model, n.obs = 50, a1 = NA)
y <- data$y

cat("Model: m =", m, "s =", s, "n.obs =", nrow(y), "\n")

# Kalman filter (optimal)
kf_result <- kf(model, y)
ll_kf <- ll_kf(model, y)
cat("Kalman filter log-likelihood:", ll_kf, "\n")

# Particle filters
N_particles <- 1000
cat("\nRunning particle filters with N =", N_particles, "particles...\n")

# SIR filter
pf_sir <- pfilter(model, y, N_particles = N_particles, method = "sir")
ll_sir <- ll_pfilter(model, y, N_particles = N_particles, filter_type = "sir", N_runs = 5)
ess_sir <- min(pf_sir$effective_sample_size)
cat("SIR filter:\n")
cat("  Log-likelihood:", ll_sir, " (diff from KF:", ll_sir - ll_kf, ")\n")
cat("  Min ESS:", ess_sir, "\n")

# APF filter
pf_apf <- pfilter(model, y, N_particles = N_particles, method = "apf")
ll_apf <- ll_pfilter(model, y, N_particles = N_particles, filter_type = "apf", N_runs = 5)
ess_apf <- min(pf_apf$effective_sample_size)
cat("APF filter:\n")
cat("  Log-likelihood:", ll_apf, " (diff from KF:", ll_apf - ll_kf, ")\n")
cat("  Min ESS:", ess_apf, "\n")

# Optimal proposal filter
pf_opt <- pfilter(model, y, N_particles = N_particles, method = "optimal")
ll_opt <- ll_pfilter(model, y, N_particles = N_particles, filter_type = "optimal", N_runs = 5)
ess_opt <- min(pf_opt$effective_sample_size)
cat("Optimal proposal filter:\n")
cat("  Log-likelihood:", ll_opt, " (diff from KF:", ll_opt - ll_kf, ")\n")
cat("  Min ESS:", ess_opt, "\n")

# Compare filtered states RMSE
rmse_sir <- sqrt(mean((pf_sir$filtered_states - kf_result$a)^2))
rmse_apf <- sqrt(mean((pf_apf$filtered_states - kf_result$a)^2))
rmse_opt <- sqrt(mean((pf_opt$filtered_states - kf_result$a)^2))
cat("\nRMSE of filtered states vs Kalman filter:\n")
cat("  SIR:", rmse_sir, "\n")
cat("  APF:", rmse_apf, "\n")
cat("  Optimal:", rmse_opt, "\n")

# Check weight degeneracy
cat("\nWeight degeneracy indicators (ESS < N_particles/2):\n")
cat("  SIR: ESS < threshold at", sum(pf_sir$effective_sample_size < N_particles/2), "time points\n")
cat("  APF: ESS < threshold at", sum(pf_apf$effective_sample_size < N_particles/2), "time points\n")
cat("  Optimal: ESS < threshold at", sum(pf_opt$effective_sample_size < N_particles/2), "time points\n")

cat("\n=== Conclusion ===\n")
cat("Optimal proposal should have:\n")
cat("1. Log-likelihood close to Kalman filter (small difference)\n")
cat("2. High ESS (little weight degeneracy)\n")
cat("3. Small RMSE for filtered states\n")
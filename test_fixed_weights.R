devtools::load_all()
set.seed(123)
m <- 1
s <- 2
n <- s
tmpl <- tmpl_stsp_full(m, n, s, sigma_L = "chol")
model <- r_model(tmpl, bpoles = 1, sd = 0.5)
data <- sim(model, n.obs = 10, a1 = NA)
y <- data$y
cat("Kalman filter log-likelihood:", ll_kf(model, y), "\n")

# Test optimal filter with 100 particles
pf_opt <- pfilter(model, y, N_particles = 100, method = "optimal")
cat("Optimal filter log-likelihood:", pf_opt$log_likelihood, "\n")
cat("ESS min:", min(pf_opt$effective_sample_size), "\n")

# Test SIR filter with 100 particles
pf_sir <- pfilter(model, y, N_particles = 100, method = "sir")
cat("SIR filter log-likelihood:", pf_sir$log_likelihood, "\n")
cat("ESS min:", min(pf_sir$effective_sample_size), "\n")

# Test APF filter with 100 particles
pf_apf <- pfilter(model, y, N_particles = 100, method = "apf")
cat("APF filter log-likelihood:", pf_apf$log_likelihood, "\n")
cat("ESS min:", min(pf_apf$effective_sample_size), "\n")

# Run ll_pfilter with multiple runs
cat("\nll_pfilter averages (N_runs=5):\n")
cat("Optimal:", ll_pfilter(model, y, N_particles = 100, filter_type = "optimal", N_runs = 5), "\n")
cat("SIR:", ll_pfilter(model, y, N_particles = 100, filter_type = "sir", N_runs = 5), "\n")
cat("APF:", ll_pfilter(model, y, N_particles = 100, filter_type = "apf", N_runs = 5), "\n")
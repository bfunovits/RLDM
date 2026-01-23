devtools::load_all()
set.seed(123)
m <- 1
s <- 2
n <- s
tmpl <- tmpl_stsp_full(m, n, s, sigma_L = "chol")
model <- r_model(tmpl, bpoles = 1, sd = 0.5)
data <- sim(model, n.obs = 5, a1 = NA)
y <- data$y
cat("Running optimal filter with 5 particles\n")
pf <- pfilter(model, y, N_particles = 5, method = "optimal")
cat("\nFiltered states:\n")
print(pf$filtered_states)
cat("\nWeight trajectories:\n")
print(pf$weight_trajectories)
cat("\nESS:\n")
print(pf$effective_sample_size)
cat("\nLog-likelihood contributions:\n")
print(pf$log_likelihood_contributions)
cat("\nTotal log-likelihood:", pf$log_likelihood, "\n")
# Kalman filter for comparison
kf_result <- kf(model, y)
cat("\nKalman filtered states:\n")
print(kf_result$a)
cat("\nKalman log-likelihood:", ll_kf(model, y), "\n")
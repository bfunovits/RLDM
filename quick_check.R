devtools::load_all(quiet=TRUE)
set.seed(123)
m <- 1; s <- 2; n <- s
tmpl <- tmpl_stsp_full(m, n, s, sigma_L = "chol")
model <- r_model(tmpl, bpoles = 1, sd = 0.5)
data <- sim(model, n.obs = 3, a1 = NA)
y <- data$y

# Kalman filter
ll_kf <- ll_kf(model, y)
cat("Kalman LL:", ll_kf, "\n")

# Optimal PF with 1000 particles
pf <- pfilter(model, y, N_particles = 1000, method = "optimal")
cat("Optimal PF LL:", pf$log_likelihood, "\n")
cat("ESS min:", min(pf$effective_sample_size), "\n")

# Compare per-step contributions
kf_result <- kf(model, y)
cat("\nKalman innovations and F:\n")
for (t in 1:3) {
  v <- kf_result$v[t,1]
  F_val <- kf_result$F[,,t][1,1]
  ll_contrib <- -0.5 * (log(F_val) + v^2/F_val + log(2*pi))
  cat(sprintf("t=%d: v=%.4f, F=%.4f, LL_contrib=%.4f\n", t, v, F_val, ll_contrib))
}
cat("Sum of Kalman contributions:", sum(sapply(1:3, function(t) {
  v <- kf_result$v[t,1]; F_val <- kf_result$F[,,t][1,1]
  -0.5 * (log(F_val) + v^2/F_val + log(2*pi))
})), "\n")

cat("\nPF contributions:", pf$log_likelihood_contributions, "\n")
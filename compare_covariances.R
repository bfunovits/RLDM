devtools::load_all(quiet=TRUE)
set.seed(123)
m <- 1; s <- 2; n <- s
tmpl <- tmpl_stsp_full(m, n, s, sigma_L = "chol")
model <- r_model(tmpl, bpoles = 1, sd = 0.5)
data <- sim(model, n.obs = 3, a1 = NA)
y <- data$y

# Extract matrices
A <- model$sys$A
C <- model$sys$C
sigma <- tcrossprod(model$sigma_L)
R <- model$sys$D %*% sigma %*% t(model$sys$D)
S <- model$sys$B %*% sigma %*% t(model$sys$D)
Q <- model$sys$B %*% sigma %*% t(model$sys$B)

# Initial covariance P1 (from model)
P1 <- lyapunov(A, tcrossprod(model$sys$B %*% model$sigma_L), non_stable = 'stop')
cat("P1 (initial state covariance):\n"); print(P1)

# Kalman filter step 1: F_1 = C P1 C' + R
F1 <- C %*% P1 %*% t(C) + R
cat("\nF1 = C P1 C' + R:\n"); print(F1)

# Our predictive covariance S_cov = C Q C' + R + C S + S' C'
S_cov <- C %*% Q %*% t(C) + R + C %*% S + t(S) %*% t(C)
cat("\nS_cov = C Q C' + R + C S + S' C':\n"); print(S_cov)

cat("\nRatio F1 / S_cov:", F1[1,1] / S_cov[1,1], "\n")
cat("F1 is", ifelse(F1[1,1] > S_cov[1,1], "larger", "smaller"), "than S_cov\n")

# Compute Kalman gain for first step: K_1 = (A P1 C' + S) (C P1 C' + R)^{-1}
K1 <- (A %*% P1 %*% t(C) + S) %*% solve(F1)
cat("\nK1 (first step Kalman gain):\n"); print(K1)

# Our K = (Q C' + S) (C Q C' + R + C S + S' C')^{-1}
K_our <- (Q %*% t(C) + S) %*% solve(S_cov)
cat("\nK_our (our proposal gain):\n"); print(K_our)

# Compute steady-state Kalman gain
P_steady <- P1
for (iter in 1:1000) {
  F_iter <- C %*% P_steady %*% t(C) + R
  K_iter <- (A %*% P_steady %*% t(C) + S) %*% solve(F_iter)
  P_next <- A %*% P_steady %*% t(A) + Q - K_iter %*% F_iter %*% t(K_iter)
  if (max(abs(P_next - P_steady)) < 1e-12) break
  P_steady <- P_next
}
F_steady <- C %*% P_steady %*% t(C) + R
K_steady <- (A %*% P_steady %*% t(C) + S) %*% solve(F_steady)
cat("\nSteady-state solution:\n")
cat("P_steady:\n"); print(P_steady)
cat("F_steady:\n"); print(F_steady)
cat("K_steady:\n"); print(K_steady)

# Compare with Q
cat("\nQ:\n"); print(Q)
cat("P_steady - Q:\n"); print(P_steady - Q)

# The optimal proposal uses conditional distribution p(x_t | x_{t-1}, y_t)
# This uses the model noise covariance Q, not the filtering covariance P
# So our formulas using Q are correct for the optimal proposal.

# But the weight p(y_t | x_{t-1}) uses S_cov = C Q C' + R + C S + S' C'
# While the Kalman filter uses F_t = C P_{t|t-1} C' + R
# These are different! p(y_t | x_{t-1}) is the conditional distribution given a specific x_{t-1}
# p(y_t | y_{1:t-1}) is the predictive distribution marginalizing over x_{t-1}
# The particle filter approximates the latter by averaging over particles.

# So maybe the particle filter likelihood should indeed be different from Kalman filter?
# No, as N → ∞, the average should converge to the marginal.

# Let's compute the Kalman filter likelihood contributions exactly
kf_result <- kf(model, y)
cat("\nKalman innovations e:\n"); print(kf_result$e[,1])

# Compute likelihood contributions from innovations
ll_contrib_kf <- numeric(3)
for (t in 1:3) {
  # Need F_t, not stored in kf_result. Let's run Kalman filter manually.
}
cat("\nKalman total ll:", kf_result$ll, "\n")

# Run optimal PF with many particles
pf <- pfilter(model, y, N_particles = 10000, method = "optimal")
cat("\nPF with 10000 particles:\n")
cat("LL:", pf$log_likelihood, "\n")
cat("ESS min:", min(pf$effective_sample_size), "\n")

# Compute variance of estimator
N_runs <- 20
ll_values <- numeric(N_runs)
for (r in 1:N_runs) {
  pf_run <- pfilter(model, y, N_particles = 1000, method = "optimal")
  ll_values[r] <- pf_run$log_likelihood
}
cat("\nPF LL mean (20 runs, N=1000):", mean(ll_values), "\n")
cat("PF LL sd:", sd(ll_values), "\n")
cat("Bias from Kalman:", mean(ll_values) - kf_result$ll, "\n")
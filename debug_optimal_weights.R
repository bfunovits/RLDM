devtools::load_all()
set.seed(123)
m <- 1
s <- 2
n <- s
tmpl <- tmpl_stsp_full(m, n, s, sigma_L = "chol")
model <- r_model(tmpl, bpoles = 1, sd = 0.5)
data <- sim(model, n.obs = 3, a1 = NA)  # Very short for debugging
y <- data$y

# Extract matrices
A <- model$sys$A
C <- model$sys$C
sigma <- tcrossprod(model$sigma_L)
R <- model$sys$D %*% sigma %*% t(model$sys$D)
S <- model$sys$B %*% sigma %*% t(model$sys$D)
Q <- model$sys$B %*% sigma %*% t(model$sys$B)

cat("=== Model Matrices ===\n")
cat("A:\n"); print(A)
cat("C:\n"); print(C)
cat("Q:\n"); print(Q)
cat("R:\n"); print(R)
cat("S:\n"); print(S)

# Compute predictive covariance
S_cov <- C %*% Q %*% t(C) + R + C %*% S + t(S) %*% t(C)
cat("\n=== Predictive Covariance S_cov ===\n")
cat("CQC':\n"); print(C %*% Q %*% t(C))
cat("R:\n"); print(R)
cat("CS:\n"); print(C %*% S)
cat("S'C':\n"); print(t(S) %*% t(C))
cat("S_cov total:\n"); print(S_cov)

# Compute Kalman gain
K <- (Q %*% t(C) + S) %*% solve(S_cov)
cat("\n=== Kalman Gain K ===\n"); print(K)

# Compute proposal covariance Sigma
Sigma <- Q - K %*% (C %*% Q + t(S))
cat("\n=== Proposal Covariance Sigma ===\n"); print(Sigma)

# Run Kalman filter for comparison
kf_result <- kf(model, y)
cat("\n=== Kalman Filter Results ===\n")
cat("Log-likelihood:", ll_kf(model, y), "\n")
cat("Innovation covariances F_t:\n")
for (t in 1:nrow(y)) {
  cat(sprintf("t=%d: F_t = %f\n", t, kf_result$F[,,t]))
}
cat("Kalman gains K_t (first few):\n")
for (t in 1:min(3, nrow(y))) {
  cat(sprintf("t=%d: K_t = [%f, %f]\n", t, kf_result$K[1,1,t], kf_result$K[2,1,t]))
}

# Run optimal particle filter with small N
N_particles <- 5
cat(sprintf("\n=== Optimal Particle Filter (N=%d) ===\n", N_particles))
pf <- pfilter(model, y, N_particles = N_particles, method = "optimal")
cat("PF log-likelihood:", pf$log_likelihood, "\n")
cat("Log-likelihood contributions:\n"); print(pf$log_likelihood_contributions)
cat("ESS:\n"); print(pf$effective_sample_size)

# Let's manually compute weights for first time step
cat("\n=== Manual weight computation for t=1 ===\n")
# Initial particles and weights
particles_t0 <- pf$particles[,,1]  # t=0
weights_t0 <- pf$weight_trajectories[1,]
cat("Initial particles (columns):\n"); print(particles_t0)
cat("Initial weights:\n"); print(weights_t0)

# Compute predictive means for each particle: μ_i = A x_{t-1} + K (y_t - C A x_{t-1})
y_obs <- y[1,]
mu <- matrix(0, s, N_particles)
for (i in 1:N_particles) {
  x_pred <- A %*% particles_t0[,i]
  mu[,i] <- x_pred + K %*% (y_obs - C %*% x_pred)
}
cat("Predictive means mu_i (columns):\n"); print(mu)

# Compute predictive likelihood p(y_t | x_{t-1}^i)
log_pdf <- numeric(N_particles)
for (i in 1:N_particles) {
  x_pred <- A %*% particles_t0[,i]
  y_pred <- C %*% x_pred
  residual <- y_obs - y_pred
  # Using S_cov (should match Kalman filter innovation variance at steady state)
  log_pdf[i] <- -0.5 * (m * log(2*pi) + log(det(S_cov)) + t(residual) %*% solve(S_cov) %*% residual)
}
cat("Log p(y_t | x_{t-1}^i):\n"); print(log_pdf)

# Unnormalized weights: w_t^i ∝ w_{t-1}^i * p(y_t | x_{t-1}^i)
log_weights <- log(weights_t0) + log_pdf
cat("Log weights (log w_{t-1} + log p(y|x)):\n"); print(log_weights)

# Normalize
max_log <- max(log_weights)
shifted <- log_weights - max_log
unnorm <- exp(shifted)
sum_unnorm <- sum(unnorm)
norm_weights <- unnorm / sum_unnorm
cat("Normalized weights:\n"); print(norm_weights)

# Likelihood contribution: log( (1/N) Σ exp(log_weights) )
ll_contribution <- max_log + log(sum_unnorm) - log(N_particles)
cat("Likelihood contribution (manual):", ll_contribution, "\n")
cat("PF likelihood contribution:", pf$log_likelihood_contributions[1], "\n")

# Compare with Kalman filter likelihood contribution
# For Kalman filter, log-likelihood contribution at t is:
# -0.5 * (log|F_t| + v_t' F_t^{-1} v_t + m*log(2π))
# where v_t = y_t - C a_{t|t-1} is innovation
cat("\n=== Kalman filter comparison ===\n")
cat("Kalman innovation v_1:", kf_result$v[1,1], "\n")
cat("Kalman F_1:", kf_result$F[,,1], "\n")
kf_ll_contribution <- -0.5 * (log(kf_result$F[,,1]) + kf_result$v[1,1]^2 / kf_result$F[,,1] + m*log(2*pi))
cat("Kalman LL contribution:", kf_ll_contribution, "\n")

# Check if S_cov matches steady-state innovation variance
cat("\n=== Steady-state check ===\n")
# Compute steady-state Kalman filter covariance via DARE
library(control)
# For time-invariant system, innovation covariance F converges to steady state
# Solve DARE: P = A P A' + Q - (A P C' + S) (C P C' + R)^{-1} (C P A' + S')
# Then F = C P C' + R
# Let's try to compute steady-state P
P_steady <- Q  # initial guess
for (iter in 1:100) {
  F_iter <- C %*% P_steady %*% t(C) + R
  K_iter <- (A %*% P_steady %*% t(C) + S) %*% solve(F_iter)
  P_next <- A %*% P_steady %*% t(A) + Q - K_iter %*% F_iter %*% t(K_iter)
  if (max(abs(P_next - P_steady)) < 1e-10) break
  P_steady <- P_next
}
F_steady <- C %*% P_steady %*% t(C) + R
cat("Steady-state innovation covariance F_steady:\n"); print(F_steady)
cat("S_cov (our predictive covariance):\n"); print(S_cov)
cat("Difference:", F_steady - S_cov, "\n")

# Note: S_cov should be C Q C' + R + C S + S' C'
# But steady-state F is C P C' + R where P is solution to DARE
# They are not equal unless P = Q (which is not true in general)
# The optimal proposal uses p(y_t | x_{t-1}) ~ N(C A x_{t-1}, C Q C' + R + C S + S' C')
# This is the one-step predictive distribution given x_{t-1}, not the steady-state Kalman innovation variance!
# The Kalman filter uses p(y_t | y_{1:t-1}) ~ N(C a_{t|t-1}, F_t) where a_{t|t-1} = A a_{t-1|t-1}
# And F_t = C P_{t|t-1} C' + R where P_{t|t-1} evolves by Riccati equation

# So maybe the issue is that we're using the wrong predictive covariance!
# The optimal proposal weight should use p(y_t | x_{t-1}) which is indeed N(C A x_{t-1}, C Q C' + R + C S + S' C')
# But the Kalman filter likelihood uses p(y_t | y_{1:t-1}) which is different.
# Wait, but for particle filter with optimal proposal, the likelihood approximation should converge to the true likelihood as N→∞.
# The particle filter likelihood is (1/N) Σ w_t^i where w_t^i = w_{t-1}^i * p(y_t | x_{t-1}^i)
# This is an unbiased estimate of p(y_t | y_{1:t-1})? Actually, it's an approximation of the marginal likelihood.
# Hmm, need to think more.
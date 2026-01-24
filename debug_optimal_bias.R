devtools::load_all(quiet=TRUE)
set.seed(123)

# Simple scalar linear Gaussian model
m <- 1
s <- 1
n <- m
tmpl <- tmpl_stsp_full(m, n, s, sigma_L = "chol")
model <- r_model(tmpl, bpoles = 0.5, sd = 0.5)  # stable A
data <- sim(model, n.obs = 5, a1 = NA)
y <- data$y

# Extract matrices
A <- model$sys$A
C <- model$sys$C
sigma <- tcrossprod(model$sigma_L)
R <- model$sys$D %*% sigma %*% t(model$sys$D)
S <- model$sys$B %*% sigma %*% t(model$sys$D)
Q <- model$sys$B %*% sigma %*% t(model$sys$B)

# Initial state covariance
P1 <- lyapunov(A, tcrossprod(model$sys$B %*% model$sigma_L), non_stable = 'stop')
a1 <- double(s)

cat("=== Model ===\n")
cat("A:", A, "\n")
cat("C:", C, "\n")
cat("Q:", Q, "\n")
cat("R:", R, "\n")
cat("S:", S, "\n")
cat("P1:", P1, "\n")

# Kalman filter to get true likelihood contributions
kf_result <- kf(model, y)
kf_ll_total <- kf_result$ll * nrow(y)
cat("\n=== Kalman Filter ===\n")
cat("Total log-likelihood:", kf_ll_total, "\n")

# Compute per-step contributions using Kalman recursion
N <- nrow(y)
v <- numeric(N)  # innovations
F <- numeric(N)  # innovation variances
Pt <- P1
for (t in 1:N) {
  # Prediction
  a_pred <- A %*% a1
  y_pred <- C %*% a_pred
  v[t] <- y[t] - y_pred
  F[t] <- C %*% Pt %*% t(C) + R
  # Update
  K <- (A %*% Pt %*% t(C) + S) %*% solve(F[t])
  a1 <- a_pred + K %*% v[t]
  Pt <- A %*% Pt %*% t(A) + Q - K %*% F[t] %*% t(K)
}
ll_contrib_kf <- -0.5 * (log(F) + v^2/F + m*log(2*pi))
cat("Kalman contributions:", ll_contrib_kf, "\n")
cat("Sum:", sum(ll_contrib_kf), "\n")

# Optimal particle filter with many particles
N_particles <- 10000
pf <- pfilter(model, y, N_particles = N_particles, method = "optimal")
cat("\n=== Particle Filter (optimal) ===\n")
cat("PF total ll:", pf$log_likelihood, "\n")
cat("PF contributions:", pf$log_likelihood_contributions, "\n")
cat("Difference total:", pf$log_likelihood - kf_ll_total, "\n")
cat("Difference per step:", pf$log_likelihood_contributions - ll_contrib_kf, "\n")

# Analyze weight degeneracy
cat("\n=== ESS ===\n")
cat("ESS min:", min(pf$effective_sample_size), "\n")
cat("ESS mean:", mean(pf$effective_sample_size), "\n")

# Compute theoretical optimal proposal parameters
S_cov <- C %*% Q %*% t(C) + R + C %*% S + t(S) %*% t(C)
K_opt <- (Q %*% t(C) + S) %*% solve(S_cov)
Sigma_opt <- Q - K_opt %*% (C %*% Q + t(S))
cat("\n=== Optimal Proposal Parameters ===\n")
cat("S_cov:", S_cov, "\n")
cat("K_opt:", K_opt, "\n")
cat("Sigma_opt:", Sigma_opt, "\n")

# For each time step, compute average p(y_t | x_{t-1}) using particles
# and compare with Kalman contribution
cat("\n=== Per-step analysis ===\n")
for (t in 1:N) {
  # Particles at time t-1
  if (t == 1) {
    particles_tm1 <- pf$particles[,,1]
    weights_tm1 <- pf$weight_trajectories[1,]
  } else {
    particles_tm1 <- pf$particles[,,t]
    weights_tm1 <- pf$weight_trajectories[t,]
  }
  # Compute p(y_t | x_{t-1}) for each particle
  log_pdf <- numeric(N_particles)
  for (i in 1:N_particles) {
    x_pred <- A %*% particles_tm1[,i]
    y_pred <- C %*% x_pred
    residual <- y[t] - y_pred
    log_pdf[i] <- -0.5 * (m*log(2*pi) + log(S_cov) + t(residual) %*% solve(S_cov) %*% residual)
  }
  avg <- sum(weights_tm1 * exp(log_pdf))
  log_avg <- log(avg)
  cat(sprintf("t=%d: PF contrib=%.6f, Kalman contrib=%.6f, Avg p(y|x)=%.6f, log avg=%.6f\n",
              t, pf$log_likelihood_contributions[t], ll_contrib_kf[t], avg, log_avg))
}

# Check if weight update matches theoretical
cat("\n=== Weight update check ===\n")
t <- 1
particles_t0 <- pf$particles[,,1]
weights_t0 <- pf$weight_trajectories[1,]
y_obs <- y[1,]
log_pdf <- numeric(N_particles)
for (i in 1:N_particles) {
  x_pred <- A %*% particles_t0[,i]
  y_pred <- C %*% x_pred
  residual <- y_obs - y_pred
  log_pdf[i] <- -0.5 * (m*log(2*pi) + log(S_cov) + t(residual) %*% solve(S_cov) %*% residual)
}
unnorm_weights <- weights_t0 * exp(log_pdf)
sum_unnorm <- sum(unnorm_weights)
log_sum <- log(sum_unnorm)
cat("Sum unnorm weights:", sum_unnorm, "\n")
cat("Log sum:", log_sum, "\n")
cat("PF contribution:", pf$log_likelihood_contributions[1], "\n")
cat("Expected contribution (log sum):", log_sum, "\n")
cat("Difference:", pf$log_likelihood_contributions[1] - log_sum, "\n")
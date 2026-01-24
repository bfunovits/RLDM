# Example: Particle Filter for Nonlinear State Transition
# Model: x_t = sin(x_{t-1}) + η_t, η_t ~ N(0, σ_η^2)
#        y_t = x_t + ε_t, ε_t ~ N(0, σ_ε^2)
# Particle filter handles nonlinearity naturally, while Kalman filter requires linearization.

library(RLDM)
source("inst/examples/helper_pf.R")  # Load helper functions

set.seed(123)

# Parameters
sigma_eta <- 0.5   # state noise
sigma_eps <- 0.2   # observation noise
N_time <- 100      # time steps
N_particles <- 1000

# Simulate true states and observations
x_true <- numeric(N_time + 1)
y_obs <- numeric(N_time)

x_true[1] <- rnorm(1, 0, sigma_eta)  # initial state
for (t in 1:N_time) {
  x_true[t+1] <- sin(x_true[t]) + sigma_eta * rnorm(1)
  y_obs[t] <- x_true[t] + sigma_eps * rnorm(1)
}
# Remove initial state for plotting
x_true <- x_true[-1]
time_idx <- 1:N_time

# Define model functions for particle filter
transition <- function(x_prev) {
  sin(x_prev) + sigma_eta * rnorm(1)
}

observation_density <- function(x, y) {
  # log p(y | x) = log N(y; x, σ_ε^2)
  dnorm(y, mean = x, sd = sigma_eps, log = TRUE)
}

# Prior for initial state: N(0, σ_η^2)
x0_prior <- function(N) {
  rnorm(N, 0, sigma_eta)
}

log_prior_density <- function(x) {
  dnorm(x, mean = 0, sd = sigma_eta, log = TRUE)
}

# Run particle filter
cat("Running particle filter with", N_particles, "particles...\n")
pf_result <- particle_filter_r(matrix(y_obs, ncol = 1),
                               transition = transition,
                               observation_density = observation_density,
                               N_particles = N_particles,
                               resampling_method = "systematic",
                               ess_threshold = 0.5,
                               x0_prior = x0_prior,
                               log_prior_density = log_prior_density)

cat("Particle filter log-likelihood:", pf_result$log_likelihood, "\n")
cat("Minimum ESS:", min(pf_result$effective_sample_size), "\n")

# Extended Kalman Filter (EKF) for comparison
# We linearize the state transition around the current estimate
ekf_filter <- function(y, sigma_eta, sigma_eps) {
  N <- length(y)
  x_est <- numeric(N + 1)
  P <- numeric(N + 1)

  # Initialization
  x_est[1] <- 0
  P[1] <- sigma_eta^2

  for (t in 1:N) {
    # Prediction step (linearize f(x) = sin(x) ≈ sin(x_est) + cos(x_est)*(x - x_est))
    # EKF uses Jacobian F = cos(x_est)
    F <- cos(x_est[t])
    x_pred <- sin(x_est[t])
    P_pred <- F^2 * P[t] + sigma_eta^2

    # Update step (linear observation: H = 1)
    H <- 1
    y_pred <- x_pred
    innovation <- y[t] - y_pred
    S <- H^2 * P_pred + sigma_eps^2
    K <- P_pred * H / S
    x_est[t+1] <- x_pred + K * innovation
    P[t+1] <- (1 - K * H) * P_pred
  }
  list(filtered = x_est[-1], P = P[-1])
}

cat("\nRunning Extended Kalman Filter...\n")
ekf_result <- ekf_filter(y_obs, sigma_eta, sigma_eps)

# Compute RMSE
pf_rmse <- sqrt(mean((pf_result$filtered_states[-1] - x_true)^2))
ekf_rmse <- sqrt(mean((ekf_result$filtered - x_true)^2))

cat("\nRoot Mean Square Error:\n")
cat("Particle Filter:", pf_rmse, "\n")
cat("Extended Kalman Filter:", ekf_rmse, "\n")

# Plot results
par(mfrow = c(2, 2))

# 1. True state vs filtered estimates
plot(time_idx, x_true, type = "l", lwd = 2, col = "black",
     xlab = "Time", ylab = "State",
     main = "True State and Filtered Estimates",
     ylim = range(c(x_true, pf_result$filtered_states[-1], ekf_result$filtered)))
lines(time_idx, pf_result$filtered_states[-1], col = "blue", lwd = 1.5)
lines(time_idx, ekf_result$filtered, col = "red", lwd = 1.5, lty = 2)
legend("topright", legend = c("True", "Particle Filter", "EKF"),
       col = c("black", "blue", "red"), lty = c(1, 1, 2), lwd = c(2, 1.5, 1.5))

# 2. Effective sample size
plot(0:N_time, pf_result$effective_sample_size, type = "l",
     xlab = "Time", ylab = "ESS",
     main = "Effective Sample Size")
abline(h = 0.5 * N_particles, col = "red", lty = 2)
legend("topright", legend = c("ESS", "Threshold"),
       col = c("black", "red"), lty = c(1, 2))

# 3. Particle weights histogram at final time
hist(pf_result$weights[, N_time + 1], breaks = 30,
     main = "Particle Weights (Final Time)",
     xlab = "Weight", col = "lightblue")
abline(v = 1/N_particles, col = "red", lty = 2)

# 4. Log-likelihood contributions
plot(time_idx, pf_result$log_likelihood_contributions, type = "l",
     xlab = "Time", ylab = "Log-likelihood contribution",
     main = "Log-likelihood Contributions")

par(mfrow = c(1, 1))

cat("\nExample completed. Particle filter handles nonlinear transition,\n")
cat("while EKF relies on linearization. For strong nonlinearities,\n")
cat("particle filter often outperforms EKF.\n")
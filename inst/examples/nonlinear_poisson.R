# Example: Particle Filter for Non-Gaussian Observation
# Model: x_t = A x_{t-1} + η_t, η_t ~ N(0, σ_η^2)  (linear Gaussian state)
#        y_t ~ Poisson(λ = exp(C x_t))            (non-Gaussian observation)
# Particle filter handles non-Gaussian observations naturally.

library(RLDM)
source("inst/examples/helper_pf.R")

set.seed(456)

# Parameters
A <- 0.8          # state transition
C <- 1.0          # observation coefficient
sigma_eta <- 0.3  # state noise
N_time <- 100
N_particles <- 2000

# Simulate true states and Poisson observations
x_true <- numeric(N_time + 1)
y_obs <- numeric(N_time)

x_true[1] <- rnorm(1, 0, sigma_eta / sqrt(1 - A^2))  # stationary distribution
for (t in 1:N_time) {
  x_true[t+1] <- A * x_true[t] + sigma_eta * rnorm(1)
  lambda <- exp(C * x_true[t])
  y_obs[t] <- rpois(1, lambda)
}
x_true <- x_true[-1]
time_idx <- 1:N_time

# Particle filter functions
transition <- function(x_prev) {
  A * x_prev + sigma_eta * rnorm(1)
}

observation_density <- function(x, y) {
  # log p(y | x) = log Poisson(y; λ = exp(C*x))
  lambda <- exp(C * x)
  dpois(y, lambda = lambda, log = TRUE)
}

# Prior for initial state: stationary distribution N(0, σ_η^2/(1-A^2))
stationary_var <- sigma_eta^2 / (1 - A^2)
x0_prior <- function(N) {
  rnorm(N, 0, sqrt(stationary_var))
}

log_prior_density <- function(x) {
  dnorm(x, mean = 0, sd = sqrt(stationary_var), log = TRUE)
}

cat("Running particle filter for Poisson observations...\n")
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

# Gaussian approximation (Kalman filter) for comparison
# Approximate Poisson as Gaussian with variance equal to mean (variance-stabilizing)
# This is suboptimal, especially for small counts.
kf_filter <- function(y, A, C, sigma_eta) {
  N <- length(y)
  x_est <- numeric(N + 1)
  P <- numeric(N + 1)

  stationary_var <- sigma_eta^2 / (1 - A^2)
  x_est[1] <- 0
  P[1] <- stationary_var

  for (t in 1:N) {
    # Prediction
    x_pred <- A * x_est[t]
    P_pred <- A^2 * P[t] + sigma_eta^2

    # Update using Gaussian approximation
    # Observation: y_t ~ N(exp(C x_t), exp(C x_t)) (mean = variance)
    # Linearize around x_pred: h(x) ≈ exp(C x_pred) + C exp(C x_pred) (x - x_pred)
    lambda_pred <- exp(C * x_pred)
    H <- C * lambda_pred  # Jacobian
    y_pred <- lambda_pred
    innovation <- y[t] - y_pred
    # Variance of observation approximation: use predicted lambda
    R <- lambda_pred  # variance = mean for Poisson
    S <- H^2 * P_pred + R
    K <- P_pred * H / S
    x_est[t+1] <- x_pred + K * innovation
    P[t+1] <- (1 - K * H) * P_pred
  }
  list(filtered = x_est[-1], P = P[-1])
}

cat("\nRunning Gaussian approximation (EKF for Poisson)...\n")
kf_result <- kf_filter(y_obs, A, C, sigma_eta)

# Compute RMSE (state estimation error)
pf_rmse <- sqrt(mean((pf_result$filtered_states[-1] - x_true)^2))
kf_rmse <- sqrt(mean((kf_result$filtered - x_true)^2))

cat("\nRoot Mean Square Error (state estimation):\n")
cat("Particle Filter:", pf_rmse, "\n")
cat("Gaussian Approximation:", kf_rmse, "\n")

# Plot results
par(mfrow = c(2, 2))

# 1. True state vs filtered estimates
plot(time_idx, x_true, type = "l", lwd = 2, col = "black",
     xlab = "Time", ylab = "State",
     main = "True State and Filtered Estimates",
     ylim = range(c(x_true, pf_result$filtered_states[-1], kf_result$filtered)))
lines(time_idx, pf_result$filtered_states[-1], col = "blue", lwd = 1.5)
lines(time_idx, kf_result$filtered, col = "red", lwd = 1.5, lty = 2)
legend("topright", legend = c("True", "Particle Filter", "Gaussian Approx"),
       col = c("black", "blue", "red"), lty = c(1, 1, 2), lwd = c(2, 1.5, 1.5))

# 2. Observations vs predicted means
lambda_true <- exp(C * x_true)
lambda_pf <- exp(C * pf_result$filtered_states[-1])
lambda_kf <- exp(C * kf_result$filtered)
plot(time_idx, y_obs, pch = 16, cex = 0.5, col = "darkgray",
     xlab = "Time", ylab = "Count / Predicted λ",
     main = "Observations and Predicted Means")
lines(time_idx, lambda_true, col = "black", lwd = 2)
lines(time_idx, lambda_pf, col = "blue", lwd = 1.5)
lines(time_idx, lambda_kf, col = "red", lwd = 1.5, lty = 2)
legend("topright", legend = c("Observations", "True λ", "PF λ", "KF λ"),
       col = c("darkgray", "black", "blue", "red"),
       pch = c(16, NA, NA, NA), lty = c(NA, 1, 1, 2), lwd = c(NA, 2, 1.5, 1.5))

# 3. Effective sample size
plot(0:N_time, pf_result$effective_sample_size, type = "l",
     xlab = "Time", ylab = "ESS",
     main = "Effective Sample Size")
abline(h = 0.5 * N_particles, col = "red", lty = 2)

# 4. Particle weights histogram at final time
hist(pf_result$weights[, N_time + 1], breaks = 30,
     main = "Particle Weights (Final Time)",
     xlab = "Weight", col = "lightblue")
abline(v = 1/N_particles, col = "red", lty = 2)

par(mfrow = c(1, 1))

cat("\nExample completed. Particle filter handles Poisson observations exactly,\n")
cat("while Gaussian approximation is biased, especially for small counts.\n")
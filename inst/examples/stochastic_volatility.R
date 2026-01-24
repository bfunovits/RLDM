# Example: Particle Filter for Stochastic Volatility Model
# Model: x_t = α + β x_{t-1} + σ_v v_t, v_t ~ N(0,1)  (log volatility)
#        y_t = exp(x_t/2) ε_t, ε_t ~ N(0,1)           (returns)
# This is a classic nonlinear/non-Gaussian model in finance.

library(RLDM)
source("inst/examples/helper_pf.R")

set.seed(789)

# Parameters (typical values from literature)
alpha <- -0.1   # intercept in log volatility
beta <- 0.95    # persistence
sigma_v <- 0.2  # volatility of log volatility
N_time <- 200
N_particles <- 3000

# Simulate true log-volatility and returns
x_true <- numeric(N_time + 1)
y_obs <- numeric(N_time)

# Initial state from stationary distribution
stationary_var <- sigma_v^2 / (1 - beta^2)
x_true[1] <- rnorm(1, alpha / (1 - beta), sqrt(stationary_var))

for (t in 1:N_time) {
  x_true[t+1] <- alpha + beta * x_true[t] + sigma_v * rnorm(1)
  y_obs[t] <- exp(x_true[t] / 2) * rnorm(1)
}
x_true <- x_true[-1]
time_idx <- 1:N_time

# Particle filter functions
transition <- function(x_prev) {
  alpha + beta * x_prev + sigma_v * rnorm(1)
}

observation_density <- function(x, y) {
  # log p(y | x) = log N(y; 0, exp(x))
  # y_t ~ N(0, exp(x_t))
  variance <- exp(x)
  dnorm(y, mean = 0, sd = sqrt(variance), log = TRUE)
}

# Prior for initial state: stationary distribution
x0_prior <- function(N) {
  rnorm(N, alpha / (1 - beta), sqrt(stationary_var))
}

log_prior_density <- function(x) {
  dnorm(x, mean = alpha / (1 - beta), sd = sqrt(stationary_var), log = TRUE)
}

cat("Running particle filter for stochastic volatility model...\n")
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

# Common approximation: log-squared transformation
# Let z_t = log(y_t^2) ≈ x_t + log(ε_t^2)
# The noise log(ε_t^2) follows log-chi-square(1) distribution, often approximated as Gaussian.
# This enables a linear Gaussian state space model (approximate).
# We'll implement this approximate Kalman filter for comparison.

log_sq_approx <- function(y) {
  # Add small constant to avoid log(0)
  eps <- 1e-4
  log(y^2 + eps)
}

# Kalman filter for the approximate model: z_t = x_t + w_t, w_t ~ N(μ_w, σ_w^2)
# where μ_w and σ_w^2 are mean and variance of log(chi-square(1))
mu_w <- digamma(0.5) - log(0.5)  # ≈ -1.27
sigma_w_sq <- trigamma(0.5)      # ≈ 4.93

kf_approx <- function(y, alpha, beta, sigma_v) {
  z <- log_sq_approx(y)
  N <- length(z)
  x_est <- numeric(N + 1)
  P <- numeric(N + 1)

  stationary_var <- sigma_v^2 / (1 - beta^2)
  x_est[1] <- alpha / (1 - beta)
  P[1] <- stationary_var

  for (t in 1:N) {
    # Prediction
    x_pred <- alpha + beta * x_est[t]
    P_pred <- beta^2 * P[t] + sigma_v^2

    # Update (observation equation: z_t = x_t + w_t)
    H <- 1
    y_pred <- x_pred + mu_w
    innovation <- z[t] - y_pred
    S <- H^2 * P_pred + sigma_w_sq
    K <- P_pred * H / S
    x_est[t+1] <- x_pred + K * innovation
    P[t+1] <- (1 - K * H) * P_pred
  }
  list(filtered = x_est[-1], P = P[-1])
}

cat("\nRunning approximate Kalman filter (log-squared transformation)...\n")
kf_result <- kf_approx(y_obs, alpha, beta, sigma_v)

# Compute RMSE for log-volatility estimation
pf_rmse <- sqrt(mean((pf_result$filtered_states[-1] - x_true)^2))
kf_rmse <- sqrt(mean((kf_result$filtered - x_true)^2))

cat("\nRoot Mean Square Error (log-volatility):\n")
cat("Particle Filter:", pf_rmse, "\n")
cat("Approximate KF:", kf_rmse, "\n")

# Plot results
par(mfrow = c(2, 2))

# 1. True log-volatility vs filtered estimates
plot(time_idx, x_true, type = "l", lwd = 2, col = "black",
     xlab = "Time", ylab = "Log-Volatility",
     main = "True Log-Volatility and Estimates",
     ylim = range(c(x_true, pf_result$filtered_states[-1], kf_result$filtered)))
lines(time_idx, pf_result$filtered_states[-1], col = "blue", lwd = 1.5)
lines(time_idx, kf_result$filtered, col = "red", lwd = 1.5, lty = 2)
legend("topright", legend = c("True", "Particle Filter", "Approx KF"),
       col = c("black", "blue", "red"), lty = c(1, 1, 2), lwd = c(2, 1.5, 1.5))

# 2. Observations (returns)
plot(time_idx, y_obs, type = "h", col = "gray",
     xlab = "Time", ylab = "Returns",
     main = "Observed Returns")
abline(h = 0, col = "black", lty = 2)

# 3. Effective sample size
plot(0:N_time, pf_result$effective_sample_size, type = "l",
     xlab = "Time", ylab = "ESS",
     main = "Effective Sample Size")
abline(h = 0.5 * N_particles, col = "red", lty = 2)

# 4. Estimated volatility (exp(x/2))
vol_true <- exp(x_true / 2)
vol_pf <- exp(pf_result$filtered_states[-1] / 2)
vol_kf <- exp(kf_result$filtered / 2)
plot(time_idx, vol_true, type = "l", lwd = 2, col = "black",
     xlab = "Time", ylab = "Volatility",
     main = "Estimated Volatility",
     ylim = range(c(vol_true, vol_pf, vol_kf)))
lines(time_idx, vol_pf, col = "blue", lwd = 1.5)
lines(time_idx, vol_kf, col = "red", lwd = 1.5, lty = 2)
legend("topright", legend = c("True", "Particle Filter", "Approx KF"),
       col = c("black", "blue", "red"), lty = c(1, 1, 2), lwd = c(2, 1.5, 1.5))

par(mfrow = c(1, 1))

cat("\nExample completed. Particle filter handles stochastic volatility model directly,\n")
cat("while approximate KF relies on log-squared transformation with Gaussian approximation.\n")
cat("Particle filter is more accurate but computationally more expensive.\n")
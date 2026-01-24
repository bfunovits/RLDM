# Helper functions for particle filter examples
# These are simplified R implementations for educational purposes

#' Systematic resampling
#'
#' @param weights Normalized importance weights (vector of length N)
#' @param N Number of particles (default: length(weights))
#' @return Indices of selected particles (integer vector of length N)
systematic_resampling <- function(weights, N = length(weights)) {
  u <- runif(1, 0, 1/N)
  cum_weights <- cumsum(weights)
  indices <- numeric(N)
  j <- 1
  for (i in 1:N) {
    while (cum_weights[j] < u + (i-1)/N) {
      j <- j + 1
    }
    indices[i] <- j
  }
  return(indices)
}

#' Multinomial resampling
multinomial_resampling <- function(weights, N = length(weights)) {
  sample.int(length(weights), size = N, replace = TRUE, prob = weights)
}

#' Stratified resampling
stratified_resampling <- function(weights, N = length(weights)) {
  u <- (runif(N) + (0:(N-1))) / N
  cum_weights <- cumsum(weights)
  findInterval(u, cum_weights) + 1
}

#' Effective sample size
effective_sample_size <- function(log_weights) {
  # log_weights are unnormalized log-weights
  max_log <- max(log_weights)
  scaled <- exp(log_weights - max_log)
  sum_scaled <- sum(scaled)
  weights <- scaled / sum_scaled
  1 / sum(weights^2)
}

#' Normalize log-weights to get normalized weights
normalize_log_weights <- function(log_weights) {
  max_log <- max(log_weights)
  scaled <- exp(log_weights - max_log)
  scaled / sum(scaled)
}

#' Generic bootstrap particle filter (SIR) in R
#'
#' @param y Observations matrix (N_time x m)
#' @param transition Function(x_prev) returning a new state sample
#' @param observation_density Function(x, y) returning log p(y | x)
#' @param N_particles Number of particles
#' @param resampling_method One of "systematic", "multinomial", "stratified"
#' @param ess_threshold Resample when ESS < ess_threshold * N_particles
#' @param x0_prior Function() returning initial state samples
#' @param log_prior_density Function(x) returning log prior density p(x_0)
#' @return List with filtered states, particles, weights, ESS, log-likelihood
particle_filter_r <- function(y, transition, observation_density,
                              N_particles = 1000,
                              resampling_method = c("systematic", "multinomial", "stratified"),
                              ess_threshold = 0.5,
                              x0_prior = NULL,
                              log_prior_density = NULL) {

  resampling_method <- match.arg(resampling_method)
  N_time <- nrow(y)
  m <- ncol(y)

  # Storage
  particles <- array(NA, dim = c(N_particles, N_time + 1))
  weights <- matrix(NA, nrow = N_particles, ncol = N_time + 1)
  log_weights <- matrix(NA, nrow = N_particles, ncol = N_time + 1)
  ess <- numeric(N_time + 1)
  filtered_states <- numeric(N_time + 1)
  log_likelihood_contributions <- numeric(N_time)

  # Initialization
  if (is.null(x0_prior)) {
    stop("Initial prior sampling function x0_prior must be provided")
  }
  particles[, 1] <- x0_prior(N_particles)

  if (is.null(log_prior_density)) {
    # Assume uniform prior (log density constant)
    log_weights[, 1] <- 0
  } else {
    log_weights[, 1] <- log_prior_density(particles[, 1])
  }

  # Normalize initial weights
  weights[, 1] <- normalize_log_weights(log_weights[, 1])
  ess[1] <- effective_sample_size(log_weights[, 1])
  filtered_states[1] <- sum(particles[, 1] * weights[, 1])

  # Main loop
  for (t in 1:N_time) {
    # Resample if needed
    if (ess[t] < ess_threshold * N_particles) {
      if (resampling_method == "systematic") {
        idx <- systematic_resampling(weights[, t], N_particles)
      } else if (resampling_method == "multinomial") {
        idx <- multinomial_resampling(weights[, t], N_particles)
      } else {
        idx <- stratified_resampling(weights[, t], N_particles)
      }
      particles[, t] <- particles[idx, t]
      # Reset weights to uniform
      log_weights[, t] <- 0
      weights[, t] <- 1/N_particles
    }

    # Propagate particles
    for (i in 1:N_particles) {
      particles[i, t+1] <- transition(particles[i, t])
    }

    # Update weights with observation y[t, ]
    for (i in 1:N_particles) {
      log_weights[i, t+1] <- observation_density(particles[i, t+1], y[t, ])
    }

    # Compute log-likelihood contribution
    max_log <- max(log_weights[, t+1])
    scaled <- exp(log_weights[, t+1] - max_log)
    sum_scaled <- sum(scaled)
    log_likelihood_contributions[t] <- max_log + log(sum_scaled) - log(N_particles)

    # Normalize weights
    weights[, t+1] <- scaled / sum_scaled
    ess[t+1] <- effective_sample_size(log_weights[, t+1])
    filtered_states[t+1] <- sum(particles[, t+1] * weights[, t+1])
  }

  # Total log-likelihood (average per time step)
  log_likelihood <- sum(log_likelihood_contributions)

  list(
    filtered_states = filtered_states,
    particles = particles,
    weights = weights,
    log_weights = log_weights,
    effective_sample_size = ess,
    log_likelihood = log_likelihood,
    log_likelihood_contributions = log_likelihood_contributions,
    N_particles = N_particles,
    resampling_method = resampling_method
  )
}
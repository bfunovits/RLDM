# 5) Particle Filter Estimation ----

#' Sequential Monte Carlo (Particle Filter) Methods
#'
#' Implementation of Sequential Monte Carlo (SMC) methods, also known as particle filters,
#' for state space models. These methods extend Kalman filtering to nonlinear and/or
#' non-Gaussian state space models.
#'
#' The particle filter approximates the filtering distribution \eqn{p(x_t | y_{1:t})}
#' using a set of weighted particles (samples). The basic algorithm is the
#' Sampling Importance Resampling (SIR) filter, also known as the bootstrap filter.
#'
#' The model considered is
#' \deqn{x_{t+1} = f(x_t, u_t)}{x[t+1] = f(x[t], u[t])}
#' \deqn{y_t = h(x_t, v_t)}{y[t] = h(x[t], v[t])}
#' where \eqn{f} is the state transition function, \eqn{h} is the observation function,
#' and \eqn{u_t}, \eqn{v_t} are noise processes.
#'
#' For linear Gaussian models (the default), these reduce to
#' \deqn{x_{t+1} = A x_t + B u_t}{x[t+1] = A x[t] + B u[t]}
#' \deqn{y_t = C x_t + D v_t}{y[t] = C x[t] + D v[t]}
#' with \eqn{u_t \sim N(0, Q)}{u[t] ~ N(0, Q)} and \eqn{v_t \sim N(0, R)}{v[t] ~ N(0, R)}.
#'
#' @param model [stspmod()] object, which represents the state space model.
#'   For nonlinear models, additional parameters may be required.
#' @param y sample, i.e. an \eqn{(N,m)} dimensional matrix,
#'   or a "time series" object (i.e. `as.matrix(y)` should return an
#'   \eqn{(N,m)}-dimensional numeric matrix). Missing values (`NA`, `NaN` and
#'   `Inf`) are **not** supported.
#' @param method Character string specifying the particle filter algorithm.
#'   Options: `"sir"` (Sampling Importance Resampling, default),
#'   `"apf"` (Auxiliary Particle Filter), `"optimal"` (Optimal proposal for linear Gaussian).
#'   Note: The APF method may produce biased likelihood estimates for models with
#'   cross-covariance between state and observation noise (S != 0). For linear Gaussian
#'   models with cross-correlation, the optimal proposal is recommended.
#' @param N_particles Number of particles to use (default: 1000).
#' @param resampling Resampling method: `"multinomial"`, `"systematic"` (default),
#'   or `"stratified"`.
#' @param ess_threshold Effective sample size threshold for triggering resampling
#'   (default: 0.5). Resampling occurs when ESS < ess_threshold * N_particles.
#' @param P1 \eqn{(s,s)} dimensional covariance matrix of the error of the initial state estimate,
#'   i.e. \eqn{\Pi_{1|0}}{P[1|0]}. If `NULL`, then the state covariance
#'   \eqn{P = APA'+B\Sigma B'} is used.
#' @param a1 \eqn{s} dimensional vector, which holds the initial estimate \eqn{a_{1|0}}{a[1|0]}
#'   for the state at time \eqn{t=1}. If `a1=NULL`, then a zero vector is used.
#' @param ... Additional arguments passed to filter implementations.
#'
#' @return List with components
#' \item{filtered_states}{\eqn{(N+1,s)} dimensional matrix with the filtered state estimates.
#'   The \eqn{t}-th row corresponds to \eqn{E[x_t | y_{1:t}]}.}
#' \item{predicted_states}{\eqn{(N+1,s)} dimensional matrix with the predicted state estimates.
#'   The \eqn{t}-th row corresponds to \eqn{E[x_t | y_{1:t-1}]}.}
#' \item{particles}{\eqn{(s, N\_particles, N+1)} dimensional array containing the particle trajectories.}
#' \item{weight_trajectories}{\eqn{(N+1, N\_particles)} dimensional matrix of particle weights over time.
#'   The \eqn{t}-th row corresponds to weights at time \eqn{t}.}
#' \item{weights}{\eqn{(N\_particles)} dimensional vector of normalized particle weights at final time.}
#' \item{log_likelihood}{Particle approximation of the log-likelihood.}
#' \item{log_likelihood_contributions}{\eqn{(N)} dimensional vector of log-likelihood contributions per time step.}
#' \item{effective_sample_size}{\eqn{(N+1)} dimensional vector of effective sample sizes.}
#' \item{N_particles}{Number of particles used.}
#' \item{resampling_method}{Resampling method used.}
#' \item{filter_type}{Type of particle filter used.}
#'
#' @seealso [kf()] for Kalman filter (optimal for linear Gaussian models),
#'   \code{\link[=ll_pfilter]{ll_pfilter()}} for particle filter approximation of log-likelihood.
#'
#' @name pfilter
#' @rdname pfilter
#' @importFrom graphics lines legend abline hist
#' @export
#'
#' @examples
#' # Linear Gaussian example: compare particle filter with Kalman filter
#' set.seed(123)
#' s = 2  # state dimension
#' m = 1  # number of outputs
#' n = m  # number of inputs
#' n.obs = 100 # sample size
#'
#' # Generate a stable state space model
#' tmpl = tmpl_stsp_full(m, n, s, sigma_L = "chol")
#' model = r_model(tmpl, bpoles = 1, sd = 0.5)
#' # Generate a sample
#' data = sim(model, n.obs = n.obs, a1 = NA)
#'
#' # Run particle filter
#' pf_result = pfilter(model, data$y, N_particles = 500)
#'
#' # Compare with Kalman filter
#' kf_result = kf(model, data$y)
#'
#' # Plot filtered states comparison
#' plot(pf_result$filtered_states[,1], type = "l", col = "blue",
#'      main = "Filtered State Estimates")
#' lines(kf_result$a[,1], col = "red", lty = 2)
#' legend("topright", legend = c("Particle Filter", "Kalman Filter"),
#'        col = c("blue", "red"), lty = 1:2)
pfilter = function(model, y, method = c('sir', 'apf', 'optimal'),
              N_particles = 1000,
              resampling = c('systematic', 'multinomial', 'stratified'),
              ess_threshold = 0.5,
              P1 = NULL, a1 = NULL, ...) {
  UseMethod("pfilter")
}


#' @name pfilter
#' @export
pfilter.stspmod = function(model, y, method = c('sir', 'apf', 'optimal'),
                      N_particles = 1000,
                      resampling = c('systematic', 'multinomial', 'stratified'),
                      ess_threshold = 0.5,
                      P1 = NULL, a1 = NULL, ...) {

  method = match.arg(method)
  resampling = match.arg(resampling)

  if (!inherits(model, 'stspmod')) stop('model must be an "stspmod" object!')

  d = dim(model$sys)
  m = d[1]
  n = d[2]
  s = d[3]

  if ((m == 0) || (n < m)) stop('only wide, non empty state space systems (n>=m>0) are supported.')

  # Check 'y'
  y = try(as.matrix(y))
  if ( inherits(y, 'try-error') || (!is.numeric(y)) || (!is.matrix(y)) ) {
    stop('input "y" must be a data object which may be coerced to a matrix with "as.matrix(y)')
  }
  if (length(y) == 0) stop('"y" contains no data')
  n.obs = nrow(y)  # sample size
  if (ncol(y) != m) stop('data "y" is not compatible with the model!')
  if (any(!is.finite(y))) stop('"y" contains NAs/NaNs/Infs!')

  # Check 'P1'
  if (is.null(P1)) {
    if (s > 0) {
      P1 = lyapunov(model$sys$A, tcrossprod(model$sys$B %*% model$sigma_L),
                    non_stable = 'stop')
    } else {
      P1 = matrix(0, nrow = s, ncol = s)
    }
  }
  if (!is.numeric(P1))  stop('parameter "P1" is not numeric')
  if ( is.vector(P1) ) {
    if (length(P1) == s) P1 = diag(P1, nrow = s, ncol = s)
    if (length(P1) == (s^2)) P1 = matrix(P1, nrow = s, ncol = s)
  }
  if ( (!is.matrix(P1)) || any(dim(P1) != s) ) {
    stop('parameter "P1" is not compatible with the model')
  }

  # Check 'a1'
  if (is.null(a1)) a1 = double(s)
  a1 = as.vector(a1)
  if (length(a1) != s) stop('parameter "a1" is not compatible with the model')

  sigma = tcrossprod(model$sigma_L)
  R = model$sys$D %*% sigma %*% t(model$sys$D)
  S = model$sys$B %*% sigma %*% t(model$sys$D)
  Q = model$sys$B %*% sigma %*% t(model$sys$B)

  # Warning for APF with cross-covariance S != 0
  if (method == 'apf' && norm(S, "F") > 1e-10) {
    warning("APF method may produce biased likelihood estimates with cross-covariance S != 0. Consider using method = 'optimal' for linear Gaussian models with cross-correlation.")
  }

  if (method == 'sir') {
    out = .Call(`_RLDM_pf_sir_cpp`, model$sys$A, model$sys$C, Q, R, S,
                t(y), P1, a1, N_particles, resampling, ess_threshold)
  } else if (method == 'apf') {
    out = .Call(`_RLDM_pf_apf_cpp`, model$sys$A, model$sys$C, Q, R, S,
                t(y), P1, a1, N_particles, resampling, ess_threshold)
  } else if (method == 'optimal') {
    out = .Call(`_RLDM_pf_optimal_cpp`, model$sys$A, model$sys$C, Q, R, S,
                t(y), P1, a1, N_particles, resampling, ess_threshold)
  } else {
    stop('Unknown particle filter method')
  }

  # Convert array dimensions to be consistent with R conventions
  # C++ returns particles as (s, N_particles, N+1) which is fine for R
  # but we might want to transpose for easier indexing
  # For now, leave as is

  class(out) <- "pfilter"
  return(out)
}


#' Particle Filter Approximation of Log-Likelihood
#'
#' Approximates the log-likelihood of a state space model using particle filters.
#' This function is useful for nonlinear/non-Gaussian models where the exact
#' Kalman filter likelihood is not available.
#'
#' @inheritParams pfilter
#' @param N_runs Number of independent particle filter runs to average over
#'   (default: 10). Averaging reduces the variance of the likelihood estimator.
#' @param filter_type Type of particle filter to use for likelihood approximation.
#'   Options: `"sir"` (default), `"apf"`.
#'
#' @return Particle filter approximation of the log-likelihood (scaled by \eqn{1/N}).
#'
#' @seealso [ll_kf()] for exact Kalman filter likelihood (linear Gaussian models),
#'   [pf()] for particle filtering.
#'
#' @name ll_pfilter
#' @rdname ll_pfilter
#' @export
#'
#' @examples
#' # Linear Gaussian example: compare particle filter likelihood with Kalman filter
#' set.seed(123)
#' s = 2
#' m = 1
#' n = m
#' n.obs = 100
#'
#' tmpl = tmpl_stsp_full(m, n, s, sigma_L = "chol")
#' model = r_model(tmpl, bpoles = 1, sd = 0.5)
#' data = sim(model, n.obs = n.obs)
#'
#' # Kalman filter likelihood (exact)
#' ll_exact = ll_kf(model, data$y)
#'
#' # Particle filter likelihood approximation
#' ll_approx = ll_pfilter(model, data$y, N_particles = 1000, N_runs = 5)
#'
#' # Compare (should be close for linear Gaussian model with enough particles)
#' cat("Exact (Kalman):", ll_exact, "\n")
#' cat("Approx (PF):", ll_approx, "\n")
#' cat("Difference:", ll_approx - ll_exact, "\n")
ll_pfilter = function(model, y, N_particles = 1000,
                 filter_type = c('sir', 'apf', 'optimal'),
                 resampling = c('systematic', 'multinomial', 'stratified'),
                 ess_threshold = 0.5,
                 N_runs = 10, ...) {
  UseMethod("ll_pfilter")
}


#' @name ll_pfilter
#' @export
ll_pfilter.stspmod = function(model, y, N_particles = 1000,
                         filter_type = c('sir', 'apf', 'optimal'),
                         resampling = c('systematic', 'multinomial', 'stratified'),
                         ess_threshold = 0.5,
                         N_runs = 10, ...) {

  filter_type = match.arg(filter_type)
  resampling = match.arg(resampling)

  if (!inherits(model, 'stspmod')) stop('model must be an "stspmod" object!')

  d = dim(model$sys)
  m = d[1]
  n = d[2]
  s = d[3]

  if ((m == 0) || (n < m)) stop('only wide, non empty state space systems (n>=m>0) are supported.')

  # Check 'y'
  y = try(as.matrix(y))
  if ( inherits(y, 'try-error') || (!is.numeric(y)) || (!is.matrix(y)) ) {
    stop('input "y" must be a data object which may be coerced to a matrix with "as.matrix(y)')
  }
  if (length(y) == 0) stop('"y" contains no data')
  n.obs = nrow(y)  # sample size
  if (ncol(y) != m) stop('data "y" is not compatible with the model!')
  if (any(!is.finite(y))) stop('"y" contains NAs/NaNs/Infs!')

  # Compute P1 if not provided
  P1 = matrix(0, nrow = s, ncol = s)
  if (s > 0) {
    P1 = lyapunov(model$sys$A, tcrossprod(model$sys$B %*% model$sigma_L),
                  non_stable = 'stop')
  }

  a1 = double(s)

  sigma = tcrossprod(model$sigma_L)
  R = model$sys$D %*% sigma %*% t(model$sys$D)
  S = model$sys$B %*% sigma %*% t(model$sys$D)
  Q = model$sys$B %*% sigma %*% t(model$sys$B)

  ll = .Call(`_RLDM_ll_pf_cpp`, model$sys$A, model$sys$C, Q, R, S,
             t(y), P1, a1, N_particles, filter_type, resampling,
             ess_threshold, N_runs)

  return(ll)
}


#' Diagnostics for Particle Filter Results
#'
#' Plotting and diagnostic functions for particle filter output.
#'
#' @param x Particle filter result object from [pf()].
#' @param type Type of diagnostic plot:
#'   `"states"` (filtered states with credibility intervals),
#'   `"ess"` (effective sample size over time),
#'   `"weights"` (histogram of final weights),
#'   `"likelihood"` (log-likelihood contributions over time).
#' @param ... Additional arguments passed to plotting functions.
#'
#' @return Plot (invisibly returns `x`).
#'
#' @name plot.pfilter
#' @examples
#' # See examples in pf()
#' @export
#'
plot.pfilter = function(x, type = c("states", "ess", "weights", "likelihood"), ...) {
  # This is a placeholder - implementation would be in 06_visualization.R
  # For now, provide basic plotting
  type = match.arg(type)

  if (type == "states") {
    n_states = ncol(x$filtered_states)
    n_time = nrow(x$filtered_states) - 1
    time_idx = 0:n_time

    plot(NA, xlim = range(time_idx), ylim = range(x$filtered_states),
         xlab = "Time", ylab = "State", main = "Filtered State Estimates")
    for (i in 1:n_states) {
      lines(time_idx, x$filtered_states[, i], col = i, lwd = 2)
    }
    legend("topright", legend = paste("State", 1:n_states),
           col = 1:n_states, lwd = 2)

  } else if (type == "ess") {
    n_time = length(x$effective_sample_size) - 1
    time_idx = 0:n_time
    plot(time_idx, x$effective_sample_size, type = "l",
         xlab = "Time", ylab = "Effective Sample Size",
         main = "Effective Sample Size Over Time")
    abline(h = x$ess_threshold * x$N_particles, col = "red", lty = 2)
    legend("topright", legend = c("ESS", "Threshold"),
           col = c("black", "red"), lty = 1:2)

  } else if (type == "weights") {
    hist(x$weights, main = "Particle Weights (Final Time)",
         xlab = "Weight", col = "lightblue",
         breaks = min(30, length(x$weights) %/% 10))
    abline(v = 1/length(x$weights), col = "red", lty = 2)

  } else if (type == "likelihood") {
    n_time = length(x$log_likelihood_contributions)
    time_idx = 1:n_time
    plot(time_idx, x$log_likelihood_contributions, type = "l",
         xlab = "Time", ylab = "Log-Likelihood Contribution",
         main = "Log-Likelihood Contributions Over Time")
    abline(h = 0, col = "gray", lty = 2)
  }

  invisible(x)
}

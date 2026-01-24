// pf.cpp
// Sequential Monte Carlo (Particle Filter) implementations for RLDM package
// Extends Kalman filtering to nonlinear/non-Gaussian state space models

#include <RcppArmadillo.h>
using namespace Rcpp;

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// we use rationalmatrices::lyapunov_impl from the header-only implementation
// [[Rcpp::depends(rationalmatrices)]]
#include <rationalmatrices_lyapunov.h>

// To make your C++ code callable from C++ code in other packages.
// This will generate a header file, inst/include/mypackage.h that
// can be included by other packages
// [[Rcpp::interfaces(r, cpp)]]


//' @name pf
//' @rdname pf
//' @export
// [[Rcpp::export]]
Rcpp::List pf_sir_cpp(const arma::mat& A, const arma::mat& C,
                      const arma::mat& Q, const arma::mat& R, const arma::mat& S,
                      const arma::mat& y_t, const arma::mat& P1, const arma::colvec& a1,
                      int N_particles = 1000,
                      const std::string& resampling = "systematic",
                      double ess_threshold = 0.5) {

  // note y_t is (m X N)!
  unsigned long int m = y_t.n_rows;
  unsigned long int N = y_t.n_cols;
  unsigned long int s = A.n_rows;

  // Validate inputs
  if (N_particles <= 0) {
    stop("Number of particles must be positive");
  }

  // Initialize particles
  arma::mat particles = arma::zeros<arma::mat>(s, N_particles);
  arma::colvec weights = arma::ones<arma::colvec>(N_particles) / N_particles;
  arma::colvec log_weights = arma::zeros<arma::colvec>(N_particles);

  // Initialize with prior distribution: a1 ~ N(a1, P1)
  // Generate initial particles using Cholesky decomposition of P1
  arma::mat P1_chol;
  if (s > 0) {
    P1_chol = chol(P1).t(); // Left Cholesky factor
    for (int i = 0; i < N_particles; i++) {
      particles.col(i) = a1 + P1_chol * arma::randn<arma::colvec>(s);
    }
  } else {
    // Zero-dimensional state space
    particles = arma::zeros<arma::mat>(0, N_particles);
  }

  // Storage for outputs
  arma::mat filtered_states = arma::zeros<arma::mat>(s, N+1);
  arma::mat predicted_states = arma::zeros<arma::mat>(s, N+1);
  arma::cube particle_trajectories = arma::zeros<arma::cube>(s, N_particles, N+1);
  arma::mat effective_sample_size = arma::zeros<arma::mat>(N+1, 1);
  arma::mat log_likelihood_contributions = arma::zeros<arma::mat>(N, 1);
  arma::mat weight_trajectories = arma::zeros<arma::mat>(N_particles, N+1);

  // Store initial particles
  particle_trajectories.slice(0) = particles;
  weight_trajectories.col(0) = weights;

  // Compute initial filtered state estimate (weighted average)
  filtered_states.col(0) = particles * weights;

  // Cholesky decompositions for noise covariances with regularization for singular matrices
  arma::mat Q_chol, R_chol;
  if (s > 0) {
    arma::mat Q_reg = Q;
    double eps = 1e-10 * arma::norm(Q, "fro");
    if (eps == 0) eps = 1e-10;
    Q_reg.diag() += eps;
    Q_chol = chol(Q_reg).t();
  }
  arma::mat R_reg = R;
  double eps_r = 1e-10 * arma::norm(R, "fro");
  if (eps_r == 0) eps_r = 1e-10;
  R_reg.diag() += eps_r;
  R_chol = chol(R_reg).t();

  // Cross-covariance S is not directly used in particle filter propagation
  // but may be needed for certain proposal distributions
  // For bootstrap filter, we ignore S and use transition prior as proposal

  // Main particle filter loop
  for (unsigned long int t = 0; t < N; t++) {

    // 1. Propagate particles: x_t ~ p(x_t | x_{t-1}) = N(A * x_{t-1}, Q)
    arma::mat new_particles = arma::zeros<arma::mat>(s, N_particles);
    for (int i = 0; i < N_particles; i++) {
      arma::colvec noise = Q_chol * arma::randn<arma::colvec>(s);
      new_particles.col(i) = A * particles.col(i) + noise;
    }

    // 2. Compute weights: w_t ∝ p(y_t | x_t) = N(C * x_t, R)
    arma::colvec y_obs = y_t.col(t);
    arma::colvec new_log_weights = arma::zeros<arma::colvec>(N_particles);

    for (int i = 0; i < N_particles; i++) {
      arma::colvec y_pred = C * new_particles.col(i);
      arma::colvec residual = y_obs - y_pred;

      // Compute log Gaussian density: -0.5 * (m*log(2π) + log|R| + residual' * R^{-1} * residual)
      // Using Cholesky decomposition for efficiency
      arma::colvec scaled_residual = solve(trimatl(R_chol), residual);
      double log_pdf = -0.5 * (m * std::log(2 * arma::datum::pi) +
                               2 * sum(log(R_chol.diag())) +
                               dot(scaled_residual, scaled_residual));

      new_log_weights(i) = std::log(weights(i)) + log_pdf;
    }

    // 3. Normalize weights with log-sum-exp trick for numerical stability
    double max_log_weight = new_log_weights.max();
    arma::colvec shifted_log_weights = new_log_weights - max_log_weight;
    arma::colvec new_weights = exp(shifted_log_weights);
    double sum_weights = sum(new_weights);

    if (sum_weights <= 0 || !std::isfinite(sum_weights)) {
      // Numerical issues - reset to uniform weights
      new_weights = arma::ones<arma::colvec>(N_particles) / N_particles;
      log_likelihood_contributions(t) = -std::numeric_limits<double>::infinity();
    } else {
      new_weights = new_weights / sum_weights;
      // Store log-likelihood contribution for this time step
      log_likelihood_contributions(t) = max_log_weight + std::log(sum_weights);
    }

    // 4. Compute effective sample size
    double ess = 1.0 / sum(square(new_weights));
    effective_sample_size(t) = ess;

    // 5. Resampling if ESS is below threshold
    arma::uvec indices = arma::regspace<arma::uvec>(0, N_particles-1);

    if (ess < ess_threshold * N_particles) {
      if (resampling == "multinomial") {
        // Multinomial resampling
        arma::vec cum_weights = cumsum(new_weights);
        indices = arma::zeros<arma::uvec>(N_particles);
        for (int i = 0; i < N_particles; i++) {
          double u = R::runif(0, 1);
          indices(i) = as_scalar(find(cum_weights >= u, 1));
        }
      } else if (resampling == "systematic") {
        // Systematic resampling
        indices = arma::zeros<arma::uvec>(N_particles);
        double u = R::runif(0, 1.0 / N_particles);
        arma::vec cum_weights = cumsum(new_weights);
        unsigned long int j = 0;
        for (int i = 0; i < N_particles; i++) {
          while (cum_weights(j) < u + static_cast<double>(i) / N_particles) {
            j++;
          }
          indices(i) = j;
        }
      } else if (resampling == "stratified") {
        // Stratified resampling
        indices = arma::zeros<arma::uvec>(N_particles);
        arma::vec cum_weights = cumsum(new_weights);
        for (int i = 0; i < N_particles; i++) {
          double u = (R::runif(0, 1) + i) / N_particles;
          indices(i) = as_scalar(find(cum_weights >= u, 1));
        }
      } else {
        stop("Unknown resampling method");
      }

      // Resample particles
      new_particles = new_particles.cols(indices);
      new_weights = arma::ones<arma::colvec>(N_particles) / N_particles;
    }

    // 6. Update for next iteration
    particles = new_particles;
    weights = new_weights;

    // 7. Store results
    particle_trajectories.slice(t+1) = particles;
    weight_trajectories.col(t+1) = weights;
    filtered_states.col(t+1) = particles * weights;

    // Predicted state (before seeing observation) is mean of propagated particles
    predicted_states.col(t+1) = arma::mean(new_particles, 1);
  }

  // Compute overall log-likelihood approximation
  double log_likelihood = accu(log_likelihood_contributions);

  // Return results in similar format to kf_cpp for consistency
  return Rcpp::List::create(
    Rcpp::Named("filtered_states") = filtered_states.t(),
    Rcpp::Named("predicted_states") = predicted_states.t(),
    Rcpp::Named("particles") = particle_trajectories,
    Rcpp::Named("weight_trajectories") = weight_trajectories.t(),
    Rcpp::Named("weights") = weights,
    Rcpp::Named("log_likelihood") = log_likelihood,
    Rcpp::Named("log_likelihood_contributions") = log_likelihood_contributions,
    Rcpp::Named("effective_sample_size") = effective_sample_size,
    Rcpp::Named("N_particles") = N_particles,
    Rcpp::Named("resampling_method") = resampling
  );
}

//' @name pf
//' @rdname pf
//' @export
// [[Rcpp::export]]
Rcpp::List pf_apf_cpp(const arma::mat& A, const arma::mat& C,
                      const arma::mat& Q, const arma::mat& R, const arma::mat& S,
                      const arma::mat& y_t, const arma::mat& P1, const arma::colvec& a1,
                      int N_particles = 1000,
                      const std::string& resampling = "systematic",
                      double ess_threshold = 0.5) {

  // note y_t is (m X N)!
  unsigned long int m = y_t.n_rows;
  unsigned long int N = y_t.n_cols;
  unsigned long int s = A.n_rows;

  // Validate inputs
  if (N_particles <= 0) {
    stop("Number of particles must be positive");
  }

  // Initialize particles
  arma::mat particles = arma::zeros<arma::mat>(s, N_particles);
  arma::colvec weights = arma::ones<arma::colvec>(N_particles) / N_particles;

  // Initialize with prior distribution
  arma::mat P1_chol;
  if (s > 0) {
    P1_chol = chol(P1).t();
    for (int i = 0; i < N_particles; i++) {
      particles.col(i) = a1 + P1_chol * arma::randn<arma::colvec>(s);
    }
  }

  // Storage for outputs
  arma::mat filtered_states = arma::zeros<arma::mat>(s, N+1);
  arma::cube particle_trajectories = arma::zeros<arma::cube>(s, N_particles, N+1);
  arma::mat effective_sample_size = arma::zeros<arma::mat>(N+1, 1);
  arma::mat log_likelihood_contributions = arma::zeros<arma::mat>(N, 1);
  arma::mat weight_trajectories = arma::zeros<arma::mat>(N_particles, N+1);

  particle_trajectories.slice(0) = particles;
  weight_trajectories.col(0) = weights;
  filtered_states.col(0) = particles * weights;

  // Cholesky decompositions with regularization for singular matrices
  arma::mat Q_chol, R_chol;
  arma::mat Q_reg = Q;
  if (s > 0) {
    double eps = 1e-10 * arma::norm(Q, "fro");
    if (eps == 0) eps = 1e-10;
    Q_reg.diag() += eps;
    Q_chol = chol(Q_reg).t();
  }
  arma::mat R_reg = R;
  double eps_r = 1e-10 * arma::norm(R, "fro");
  if (eps_r == 0) eps_r = 1e-10;
  R_reg.diag() += eps_r;
  R_chol = chol(R_reg).t();

  // Covariance for p(y_t | μ_t) where μ_t = A x_{t-1}
  arma::mat S_cov = C * Q * C.t() + R;
  arma::mat S_cov_reg = S_cov;
  double eps_s = 1e-10 * arma::norm(S_cov, "fro");
  if (eps_s == 0) eps_s = 1e-10;
  S_cov_reg.diag() += eps_s;
  arma::mat S_chol = chol(S_cov_reg).t();


  // Auxiliary Particle Filter (APF) - two-stage resampling
  for (unsigned long int t = 0; t < N; t++) {
    arma::colvec y_obs = y_t.col(t);

    // First stage: compute auxiliary weights and store log densities
    arma::colvec log_pdf_mu = arma::zeros<arma::colvec>(N_particles);
    arma::colvec auxiliary_weights = arma::zeros<arma::colvec>(N_particles);
    for (int i = 0; i < N_particles; i++) {
      arma::colvec mu = A * particles.col(i);
      arma::colvec y_pred = C * mu;
      arma::colvec residual = y_obs - y_pred;
      arma::colvec scaled_residual = solve(trimatl(S_chol), residual);
      double log_pdf = -0.5 * (m * std::log(2 * arma::datum::pi) +
                               2 * sum(log(S_chol.diag())) +
                               dot(scaled_residual, scaled_residual));
      log_pdf_mu(i) = log_pdf;
      auxiliary_weights(i) = weights(i) * exp(log_pdf);
    }

    // Normalize auxiliary weights
    auxiliary_weights = auxiliary_weights / sum(auxiliary_weights);

    // Resample indices based on auxiliary weights
    arma::uvec indices = arma::zeros<arma::uvec>(N_particles);
    if (resampling == "systematic") {
      // Systematic resampling
      double u = R::runif(0, 1.0 / N_particles);
      arma::vec cum_weights = cumsum(auxiliary_weights);
      unsigned long int j = 0;
      for (int i = 0; i < N_particles; i++) {
        while (cum_weights(j) < u + static_cast<double>(i) / N_particles) {
          j++;
        }
        indices(i) = j;
      }
    } else {
      // Default to multinomial
      arma::vec cum_weights = cumsum(auxiliary_weights);
      for (int i = 0; i < N_particles; i++) {
        double u = R::runif(0, 1);
        indices(i) = as_scalar(find(cum_weights >= u, 1));
      }
    }

    // Resample particles and corresponding log_pdf_mu
    particles = particles.cols(indices);
    arma::colvec log_pdf_mu_resampled = log_pdf_mu.elem(indices);
    weights = arma::ones<arma::colvec>(N_particles) / N_particles;

    // Second stage: propagate resampled particles
    arma::mat new_particles = arma::zeros<arma::mat>(s, N_particles);
    for (int i = 0; i < N_particles; i++) {
      arma::colvec noise = Q_chol * arma::randn<arma::colvec>(s);
      new_particles.col(i) = A * particles.col(i) + noise;
    }

    // Compute weights: w_t ∝ p(y_t | x_t) / p(y_t | μ_t)
    arma::colvec log_pdf_x = arma::zeros<arma::colvec>(N_particles);
    for (int i = 0; i < N_particles; i++) {
      arma::colvec y_pred = C * new_particles.col(i);
      arma::colvec residual = y_obs - y_pred;
      arma::colvec scaled_residual = solve(trimatl(R_chol), residual);
      double log_pdf = -0.5 * (m * std::log(2 * arma::datum::pi) +
                               2 * sum(log(R_chol.diag())) +
                               dot(scaled_residual, scaled_residual));
      log_pdf_x(i) = log_pdf;
    }

    // Compute unnormalized weights: exp(log_pdf_x - log_pdf_mu_resampled)
    arma::colvec log_weights = log_pdf_x - log_pdf_mu_resampled;
    double max_log_weight = log_weights.max();
    arma::colvec shifted_log_weights = log_weights - max_log_weight;
    arma::colvec unnormalized_weights = exp(shifted_log_weights);
    double sum_unnormalized = sum(unnormalized_weights);

    arma::colvec new_weights;
    if (sum_unnormalized <= 0 || !std::isfinite(sum_unnormalized)) {
      // Numerical issues - reset to uniform weights
      new_weights = arma::ones<arma::colvec>(N_particles) / N_particles;
      log_likelihood_contributions(t) = -std::numeric_limits<double>::infinity();
    } else {
      new_weights = unnormalized_weights / sum_unnormalized;
      // Store log-likelihood contribution: log( (1/N) Σ exp(log_weights) )
      log_likelihood_contributions(t) = max_log_weight + std::log(sum_unnormalized) - std::log(N_particles);
    }

    // Compute effective sample size
    double ess = 1.0 / sum(square(new_weights));
    effective_sample_size(t) = ess;

    // Update particles and weights
    particles = new_particles;
    weights = new_weights;

    // Store results
    particle_trajectories.slice(t+1) = particles;
    weight_trajectories.col(t+1) = weights;
    filtered_states.col(t+1) = particles * weights;
  }

  // Compute overall log-likelihood
  double log_likelihood = accu(log_likelihood_contributions);

  return Rcpp::List::create(
    Rcpp::Named("filtered_states") = filtered_states.t(),
    Rcpp::Named("particles") = particle_trajectories,
    Rcpp::Named("weight_trajectories") = weight_trajectories.t(),
    Rcpp::Named("weights") = weights,
    Rcpp::Named("log_likelihood") = log_likelihood,
    Rcpp::Named("log_likelihood_contributions") = log_likelihood_contributions,
    Rcpp::Named("effective_sample_size") = effective_sample_size,
    Rcpp::Named("N_particles") = N_particles,
    Rcpp::Named("resampling_method") = resampling,
    Rcpp::Named("filter_type") = "auxiliary"
  );
}

//' @name pf
//' @rdname pf
//' @export
// [[Rcpp::export]]
Rcpp::List pf_optimal_cpp(const arma::mat& A, const arma::mat& C,
                      const arma::mat& Q, const arma::mat& R, const arma::mat& S,
                      const arma::mat& y_t, const arma::mat& P1, const arma::colvec& a1,
                      int N_particles = 1000,
                      const std::string& resampling = "systematic",
                      double ess_threshold = 0.5) {

  // note y_t is (m X N)!
  unsigned long int m = y_t.n_rows;
  unsigned long int N = y_t.n_cols;
  unsigned long int s = A.n_rows;

  // Validate inputs
  if (N_particles <= 0) {
    stop("Number of particles must be positive");
  }

  // Initialize particles
  arma::mat particles = arma::zeros<arma::mat>(s, N_particles);
  arma::colvec weights = arma::ones<arma::colvec>(N_particles) / N_particles;
  arma::colvec log_weights = arma::zeros<arma::colvec>(N_particles);

  // Initialize with prior distribution: a1 ~ N(a1, P1)
  arma::mat P1_chol;
  if (s > 0) {
    P1_chol = chol(P1).t(); // Left Cholesky factor
    for (int i = 0; i < N_particles; i++) {
      particles.col(i) = a1 + P1_chol * arma::randn<arma::colvec>(s);
    }
  } else {
    // Zero-dimensional state space
    particles = arma::zeros<arma::mat>(0, N_particles);
  }

  // Storage for outputs
  arma::mat filtered_states = arma::zeros<arma::mat>(s, N+1);
  arma::mat predicted_states = arma::zeros<arma::mat>(s, N+1);
  arma::cube particle_trajectories = arma::zeros<arma::cube>(s, N_particles, N+1);
  arma::mat effective_sample_size = arma::zeros<arma::mat>(N+1, 1);
  arma::mat log_likelihood_contributions = arma::zeros<arma::mat>(N, 1);
  arma::mat weight_trajectories = arma::zeros<arma::mat>(N_particles, N+1);

  // Store initial particles and weights
  particle_trajectories.slice(0) = particles;
  weight_trajectories.col(0) = weights;
  filtered_states.col(0) = particles * weights;

  // Cholesky decompositions for noise covariances with regularization for singular matrices
  arma::mat Q_chol, R_chol;
  if (s > 0) {
    arma::mat Q_reg = Q;
    double eps = 1e-10 * arma::norm(Q, "fro");
    if (eps == 0) eps = 1e-10;
    Q_reg.diag() += eps;
    Q_chol = chol(Q_reg).t();
  }
  arma::mat R_reg = R;
  double eps_r = 1e-10 * arma::norm(R, "fro");
  if (eps_r == 0) eps_r = 1e-10;
  R_reg.diag() += eps_r;
  R_chol = chol(R_reg).t();

  // Covariance for predictive likelihood: p(y_t | x_{t-1}) ~ N(C A x_{t-1}, C Q C' + R + C S + S' C')
  // Note: S is cross-covariance between state and observation noise
  arma::mat S_cov = C * Q * C.t() + R + C * S + S.t() * C.t();
  arma::mat S_cov_reg = S_cov;
  double eps_s = 1e-10 * arma::norm(S_cov, "fro");
  if (eps_s == 0) eps_s = 1e-10;
  S_cov_reg.diag() += eps_s;
  arma::mat S_chol = chol(S_cov_reg).t();

  // Optimal proposal parameters: Kalman gain and covariance (with cross-covariance S)
  arma::mat K, Sigma, Sigma_chol;
  if (s > 0) {
    // Kalman gain: K = (Q C' + S) (C Q C' + R + C S + S' C')^{-1}
    K = (Q * C.t() + S) * inv(S_cov_reg); // using regularized covariance for stability
    // Proposal covariance: Σ = Q - K (C Q + S')
    Sigma = Q - K * (C * Q + S.t());
    // Regularize Sigma
    arma::mat Sigma_reg = Sigma;
    double eps_sigma = 1e-10 * arma::norm(Sigma, "fro");
    if (eps_sigma == 0) eps_sigma = 1e-10;
    Sigma_reg.diag() += eps_sigma;
    Sigma_chol = chol(Sigma_reg).t();
  }

  // Main particle filter loop
  for (unsigned long int t = 0; t < N; t++) {
    arma::colvec y_obs = y_t.col(t);

    // 1. Compute predictive means for each particle: μ_i = A x_{t-1} + K (y_t - C A x_{t-1})
    arma::mat mu = arma::zeros<arma::mat>(s, N_particles);
    for (int i = 0; i < N_particles; i++) {
      arma::colvec x_pred = A * particles.col(i);
      mu.col(i) = x_pred + K * (y_obs - C * x_pred);
    }

    // 2. Propagate particles using optimal proposal: x_t ~ N(μ_i, Σ)
    arma::mat new_particles = arma::zeros<arma::mat>(s, N_particles);
    for (int i = 0; i < N_particles; i++) {
      arma::colvec noise = Sigma_chol * arma::randn<arma::colvec>(s);
      new_particles.col(i) = mu.col(i) + noise;
    }

    // 3. Compute weights using predictive likelihood: w_t ∝ p(y_t | x_{t-1})
    arma::colvec new_log_weights = arma::zeros<arma::colvec>(N_particles);
    for (int i = 0; i < N_particles; i++) {
      arma::colvec x_pred = A * particles.col(i);
      arma::colvec y_pred = C * x_pred;
      arma::colvec residual = y_obs - y_pred;
      arma::colvec scaled_residual = solve(trimatl(S_chol), residual);
      double log_pdf = -0.5 * (m * std::log(2 * arma::datum::pi) +
                               2 * sum(log(S_chol.diag())) +
                               dot(scaled_residual, scaled_residual));
      new_log_weights(i) = std::log(weights(i)) + log_pdf;
    }

    // 4. Normalize weights with log-sum-exp trick
    double max_log_weight = new_log_weights.max();
    arma::colvec shifted_log_weights = new_log_weights - max_log_weight;
    arma::colvec new_weights = exp(shifted_log_weights);
    double sum_weights = sum(new_weights);

    if (sum_weights <= 0 || !std::isfinite(sum_weights)) {
      // Numerical issues - reset to uniform weights
      new_weights = arma::ones<arma::colvec>(N_particles) / N_particles;
      log_likelihood_contributions(t) = -std::numeric_limits<double>::infinity();
    } else {
      new_weights = new_weights / sum_weights;
      // Store log-likelihood contribution for this time step
      log_likelihood_contributions(t) = max_log_weight + std::log(sum_weights);
    }

    // 5. Compute effective sample size
    double ess = 1.0 / sum(square(new_weights));
    effective_sample_size(t) = ess;

    // 6. Resampling if ESS is below threshold
    arma::uvec indices = arma::regspace<arma::uvec>(0, N_particles-1);

    if (ess < ess_threshold * N_particles) {
      if (resampling == "multinomial") {
        // Multinomial resampling
        arma::vec cum_weights = cumsum(new_weights);
        indices = arma::zeros<arma::uvec>(N_particles);
        for (int i = 0; i < N_particles; i++) {
          double u = R::runif(0, 1);
          indices(i) = as_scalar(find(cum_weights >= u, 1));
        }
      } else if (resampling == "systematic") {
        // Systematic resampling
        indices = arma::zeros<arma::uvec>(N_particles);
        double u = R::runif(0, 1.0 / N_particles);
        arma::vec cum_weights = cumsum(new_weights);
        unsigned long int j = 0;
        for (int i = 0; i < N_particles; i++) {
          while (cum_weights(j) < u + static_cast<double>(i) / N_particles) {
            j++;
          }
          indices(i) = j;
        }
      } else if (resampling == "stratified") {
        // Stratified resampling
        indices = arma::zeros<arma::uvec>(N_particles);
        arma::vec cum_weights = cumsum(new_weights);
        for (int i = 0; i < N_particles; i++) {
          double u = (R::runif(0, 1) + i) / N_particles;
          indices(i) = as_scalar(find(cum_weights >= u, 1));
        }
      } else {
        stop("Unknown resampling method");
      }

      // Resample particles
      new_particles = new_particles.cols(indices);
      new_weights = arma::ones<arma::colvec>(N_particles) / N_particles;
    }

    // 7. Update for next iteration
    particles = new_particles;
    weights = new_weights;

    // 8. Store results
    particle_trajectories.slice(t+1) = particles;
    weight_trajectories.col(t+1) = weights;
    filtered_states.col(t+1) = particles * weights;

    // Predicted state (before seeing observation) is mean of propagated particles
    predicted_states.col(t+1) = arma::mean(mu, 1);
  }

  // Compute overall log-likelihood approximation
  double log_likelihood = accu(log_likelihood_contributions);

  // Return results in similar format to other filters
  return Rcpp::List::create(
    Rcpp::Named("filtered_states") = filtered_states.t(),
    Rcpp::Named("predicted_states") = predicted_states.t(),
    Rcpp::Named("particles") = particle_trajectories,
    Rcpp::Named("weight_trajectories") = weight_trajectories.t(),
    Rcpp::Named("weights") = weights,
    Rcpp::Named("log_likelihood") = log_likelihood,
    Rcpp::Named("log_likelihood_contributions") = log_likelihood_contributions,
    Rcpp::Named("effective_sample_size") = effective_sample_size,
    Rcpp::Named("N_particles") = N_particles,
    Rcpp::Named("resampling_method") = resampling,
    Rcpp::Named("filter_type") = "optimal"
  );
}


//' @name ll_pf
//' @rdname ll_pf
//' @export
// [[Rcpp::export]]
double ll_pf_cpp(const arma::mat& A, const arma::mat& C,
                 const arma::mat& Q, const arma::mat& R, const arma::mat& S,
                 const arma::mat& y_t, const arma::mat& P1, const arma::colvec& a1,
                 int N_particles = 1000,
                 const std::string& filter_type = "sir",
                 const std::string& resampling = "systematic",
                 double ess_threshold = 0.5,
                 unsigned int N_runs = 10) {

  // Run particle filter multiple times and average log-likelihood estimates
  // to reduce variance of the estimator

  double total_log_likelihood = 0.0;

  for (unsigned int run = 0; run < N_runs; run++) {
    Rcpp::List pf_result;
    if (filter_type == "sir") {
      pf_result = pf_sir_cpp(A, C, Q, R, S, y_t, P1, a1,
                             N_particles, resampling, ess_threshold);
    } else if (filter_type == "apf") {
      pf_result = pf_apf_cpp(A, C, Q, R, S, y_t, P1, a1,
                             N_particles, resampling, ess_threshold);
    } else if (filter_type == "optimal") {
      pf_result = pf_optimal_cpp(A, C, Q, R, S, y_t, P1, a1,
                             N_particles, resampling, ess_threshold);
    } else {
      stop("Unknown filter type");
    }

    total_log_likelihood += Rcpp::as<double>(pf_result["log_likelihood"]);
  }

  return total_log_likelihood / N_runs;
}
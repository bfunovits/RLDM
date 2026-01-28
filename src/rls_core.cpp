// rls_core.cpp
// Recursive Least Squares implementations for RLDM package

#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

/**
 * @brief Debug utility to display matrix contents
 *
 * @param X Matrix to display
 * @param name Name to print before matrix
 */
// [[Rcpp::export]]
void showMatrix(arma::mat X, const char* name) {
  Rcout << name << std::endl << X << std::endl;
}

/**
 * @brief Multivariate Recursive Least Squares with exponential forgetting
 *
 * Implements the Recursive Least Squares (RLS) algorithm with exponentially
 * down-weighted past or fixed window size. Allows for multivariate regressions
 * with multiple left-hand-side variables y. The algorithm assumes the same
 * regressors for each component of the LHS variable:
 *   Y = XB + U
 * where:
 *   - Y and U are (T x n_y) dimensional
 *   - X is (T x n_x) dimensional
 *   - B is (n_x x n_y) dimensional
 *
 * In a VAR(p) context where the endogenous variable y is n-dimensional,
 * this corresponds to dimensions (T x n), (T x n p), (n p x n) for Y, X, B
 * respectively. Using the same regressors for each LHS component simplifies
 * the updating formula.
 *
 * @param Y Left-hand-side variables in multivariate regression (T x n_y)
 * @param X Right-hand-side variables (T x n_x)
 * @param r Forgetting factors matrix (n_r x 1). Each observation has the same
 *        weight (different components have same forgetting factor)
 * @param n_init Number of observations used for initial estimate of regression
 *        coefficients. Default is twice the number of regressors in X.
 *        Forecasts and forecast errors are only produced AFTER initial estimate.
 * @param allow_neg If false, negative forecasts are set to zero (default: true)
 * @param debug_flag If true, output is printed to console (default: false)
 * @return Rcpp::List containing:
 *   - y_pred: Cube (n_obs x n_y x n_r) with one-step-ahead predictions for all
 *             forgetting factors. For windowed RLS, third dimension is discarded.
 *   - fe_honest: Cube (n_obs x n_y x n_r) with one-step-ahead forecast errors.
 *                For windowed RLS, third dimension is discarded.
 *   - forgetting: Input forgetting factors r
 *
 * @note Inputs are not checked - validation should be done in R before calling
 * @note Main reference: Chapter 4 in Young (2012)
 */
// [[Rcpp::export]]
Rcpp::List rls_exp_cpp(const arma::mat Y, const arma::mat X,
                       const arma::mat r, int n_init,
                       bool allow_neg = true, bool debug_flag = false) {

  // For debugging: Rprintf("\n %d \n", i); showMatrix(r_mat, "x: ");

  // Dimensions
  int n_x = X.n_cols;
  int n_obs = Y.n_rows;
  int n_y = Y.n_cols;
  int n_r = r.n_rows;  // Number of forgetting factors

  // Initial estimator
  arma::mat y_init = arma::zeros(n_init, n_y);
  arma::mat X_init = arma::zeros(n_init, n_x);
  arma::mat P = arma::zeros(n_x, n_x);
  arma::mat B = arma::zeros(n_x, n_y);

  X_init = X(arma::span(0, n_init - 1), arma::span(0, n_x - 1));
  y_init = Y(arma::span(0, n_init - 1), arma::span(0, n_y - 1));
  P = pinv(X_init.t() * X_init);
  B = P * X_init.t() * y_init;

  // Initialization for RLS
  arma::cube y_pred = arma::zeros(n_obs, n_y, n_r);
  arma::cube fe_honest = arma::zeros(n_obs, n_y, n_r);

  // Initializing stuff for for-loop
  int idx_r;
  int idx_t;
  int idx_y;

  arma::mat x = arma::zeros(1, n_x);
  arma::mat g = arma::zeros(n_x, 1);
  arma::mat Px = arma::zeros(n_x, 1);

  arma::mat tmp_neg_val = arma::zeros(1, n_y);

  for (idx_r = 0; idx_r < n_r; idx_r++) {
    for (idx_t = n_init; idx_t < n_obs; idx_t++) {
      // Define quantities as in Young (2012) Chapter 4
      // If the RHS variables were not the same for all LHS variables,
      // e.g., G below would be a matrix (instead of a vector)
      x = X(idx_t, arma::span());
      Px = (P * x.t());
      g = Px.each_row() / (r(idx_r, 0) + x * P * x.t());

      // Prediction
      y_pred.slice(idx_r)(idx_t, arma::span()) = x * B;

      // Predictions if negative values not allowed
      if (allow_neg == false) {
        tmp_neg_val = x * B;
        for (idx_y = 0; idx_y < n_y; idx_y++) {
          if (tmp_neg_val(0, idx_y) < 0) {
            y_pred.slice(idx_r)(idx_t, idx_y) = 0;  // comment_bf: types ok?
          }
        }
      }

      // Forecast errors
      fe_honest.slice(idx_r)(idx_t, arma::span()) =
          Y(idx_t, arma::span()) - y_pred.slice(idx_r)(idx_t, arma::span());

      // Update regression coefficients and (X'X)^{-1} equivalent
      B = B + g * fe_honest.slice(idx_r)(idx_t, arma::span());
      P = 1 / r(idx_r, 0) * (P - g * x * P);
    }
  }

  return Rcpp::List::create(Rcpp::Named("y_pred") = y_pred,
                            Rcpp::Named("fe_honest") = fe_honest,
                            Rcpp::Named("forgetting") = r);
}

/**
 * @brief Windowed Recursive Least Squares
 *
 * Implements RLS with fixed window size. Uses sliding window of observations
 * to compute estimates, dropping oldest observation when adding new one.
 *
 * @param Y Left-hand-side variables (T x n_y)
 * @param X Right-hand-side variables (T x n_x)
 * @param ws Window size (number of observations in estimation window)
 * @param allow_neg If false, negative forecasts are set to zero (default: true)
 * @param debug_flag If true, output is printed to console (default: false)
 * @return Rcpp::List containing:
 *   - y_pred: Matrix (n_obs x n_y) with one-step-ahead predictions
 *   - fe_honest: Matrix (n_obs x n_y) with one-step-ahead forecast errors
 */
// [[Rcpp::export]]
Rcpp::List rls_window_cpp(const arma::mat Y, const arma::mat X, int ws,
                          bool allow_neg = true, bool debug_flag = false) {

  // Rprintf("\n %d \n", i); showMatrix(X, "x: ");

  // Dimensions
  int n_obs = Y.n_rows;
  int n_y = Y.n_cols;
  int n_x = X.n_cols;

  // Initial estimators for residuals and regression coefficients
  arma::mat Y_init = arma::zeros(ws, n_y);
  arma::mat X_init = arma::zeros(ws, n_x);
  arma::mat P = arma::zeros(n_x, n_x);
  arma::mat A = arma::zeros(n_x, n_y);

  X_init = X(arma::span(0, ws - 1), arma::span(0, n_x - 1));
  Y_init = Y(arma::span(0, ws - 1), arma::span(0, n_y - 1));
  P = pinv(X_init.t() * X_init);
  A = P * X_init.t() * Y_init;

  // Initialization for RLS
  arma::mat y_pred = arma::zeros(n_obs, n_y);
  arma::mat fe_honest = arma::zeros(n_obs, n_y);

  // Initializing stuff for for-loop
  int idx_t;
  int idx_y;

  arma::mat xf = arma::zeros(1, n_x);
  arma::mat xb = arma::zeros(1, n_x);
  arma::mat g = arma::zeros(n_x, 1);
  arma::mat Px = arma::zeros(n_x, 1);  // Corresponds to P(k-1) x(k) in Young p52

  arma::mat tmp_neg_val = arma::zeros(1, n_y);

  for (idx_t = ws; idx_t < n_obs; idx_t++) {
    // 1) Update with new observation

    // Simplifying notation
    xf = X(idx_t, arma::span());
    Px = (P * xf.t());

    // "Kalman-Gain"
    g = Px.each_row() / (1 + xf * Px);

    // Update regression coefficients and state covariance
    A = A + g * (Y(idx_t, arma::span(0, n_y - 1)) - xf * A);
    P = P - g * xf * P;

    // 2) Subtracting last observation in (old) window

    // Simplifying notation
    xb = X(idx_t - ws, arma::span());
    Px = P * xb.t();

    // Updates
    g = Px.each_row() / (-1 + xb * Px);
    A = A + g * (Y(idx_t - ws, arma::span(0, n_y - 1)) - xb * A);
    P = P - g * xb * P;

    // Prediction
    y_pred(idx_t, arma::span()) = xf * A;

    // Predictions if negative values not allowed
    if (allow_neg == false) {
      tmp_neg_val = y_pred(idx_t, arma::span());
      for (idx_y = 0; idx_y < n_y; idx_y++) {
        if (tmp_neg_val(0, idx_y) < 0) {
          y_pred(idx_t, idx_y) = 0;  // comment_bf: types ok?
        }
      }
    }

    // Forecast errors
    fe_honest(idx_t, arma::span()) =
        Y(idx_t, arma::span()) - y_pred(idx_t, arma::span());
  }

  return Rcpp::List::create(Rcpp::Named("y_pred") = y_pred,
                            Rcpp::Named("fe_honest") = fe_honest);
}
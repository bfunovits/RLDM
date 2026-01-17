#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
void showMatrix(arma::mat X, const char* name) {
  Rcout << name << std::endl << X << std::endl;
}


//' (Multivariate) Recursive Least Squares
//'
//' This function implements the Recursive Least Squares (RLS) algorithm with exponentially down-weighted past of fixed window size.
//' It allows for multivariate regressions, i.e. for multiple left-hand-side variables \code{y}.
//' However, the algorithm is formulated such that the regressors are the same for each individual component of the LHS variable,
//' i.e. \deqn{Y = XB + U}
//' where
//' \itemize{
//' \item \code{Y} and \code{U} are \eqn{(T \times n_{y})} dimensional,
//' \item \code{X} is \eqn{(T \times n_{x})} dimensional,
//' \item \code{B} is \eqn{(n_{x} \times n_{y})} dimensional,
//' }
//' In a VAR(p) context where the endogenous variable \code{y} is \code{n}-dimensional,
//' this corresponds to dimensions \eqn{(T \times n)}, \eqn{(T \times n p)}, \eqn{( n p \times n)} for \code{Y}, \code{X}, \code{B} respectively.
//' The fact that we use the same regressors for each component of the LHS simplifies the updating formula in RLS a bit.
//'
//' Inputs are not checked, so this should be done within R before calling this function.
//'
//' The main reference is Chapter 4 in Young (2012).
//'
//' @param Y,X Matrices of doubles.
//'   Left-hand-side and right-hand-side variables in possibly multivariate regression.
//'   The number of rows corresponds to the number of observations
//' @param r Matrix of doubles of dimension \eqn{(n_{r} \times 1)}.
//'   Contains forgetting factors. Note that each observation has the same weight (i.e. different components )
//' @param n_init Integer. Number of observations used for initial estimate of regression coefficients.
//'   Default is set to twice the number of regressors in \code{X}.
//'   Note that forecasts and forecast errors are only produced AFTER the initial estimate.
//' @param ws Integer. Fixed number of observation to be used to calculate estimates.
//'   Number of observations used for initial estimate of regression coefficients.
//' @param allow_neg Boolean. Default set to true. If false, negative forecasts are not allowed and set to zero.
//' @param debug_flag Boolean. Default set to false. If true, output is printed to the console.
//' @return List containing
//'   \itemize{
//'     \item \code{y_pred}: Cube of doubles of dimension \eqn{(n_{obs} \times n_y \times n_{r})}.
//'       Contains (one-step-ahead) predictions for all forgetting factors.
//'       In case of windowed RLS, the third dimension (containing forgetting factors) is discarded (such that the output is a matrix).
//'     \item \code{fe_honest:} Cube of doubles of dimension \eqn{(n_{obs} \times n_y \times n_{r})}.
//'       Contains (one-step-ahead) forecast errors for all forgetting factors.
//'       In case of windowed RLS, the third dimension (containing forgetting factors) is discarded (such that the output is a matrix).
//'   }
//' @export
// [[Rcpp::export]]
Rcpp::List rls_exp_cpp(const arma::mat Y,
                        const arma::mat X,
                        const arma::mat r,
                        int n_init,
                        bool allow_neg = true,
                        bool debug_flag = false){

  // For debugging: Rprintf("\n %d \n", i); showMatrix(r_mat, "x: ");

  // dimensions
  int n_x = X.n_cols;
  int n_obs = Y.n_rows;
  int n_y = Y.n_cols;
  int n_r = r.n_rows; // number of forgetting factors

  // initial estimator
  mat y_init = zeros(n_init, n_y);
  mat X_init = zeros(n_init, n_x);
  mat P = zeros(n_x, n_x);
  mat B = zeros(n_x, n_y);

  X_init = X(span(0, n_init-1), span(0,n_x-1));
  y_init = Y(span(0, n_init-1), span(0,n_y-1));
  P = pinv(X_init.t() * X_init );
  B = P * X_init.t() * y_init;

  // Initialization for RLS
  cube y_pred = zeros(n_obs, n_y, n_r);
  cube fe_honest = zeros(n_obs, n_y, n_r);


  // Initializing stuff for for-loop
  int idx_r;
  int idx_t;
  int idx_y;

  mat x = zeros(1, n_x);
  mat g = zeros(n_x, 1);
  mat Px = zeros(n_x, 1);

  mat tmp_neg_val = zeros(1, n_y);

  for (idx_r = 0; idx_r < n_r; idx_r++){
    for (idx_t = n_init; idx_t < n_obs; idx_t++){
      // Define quantities as in Young (2012) Chapter 4
      // If the RHS variables were not the same for all LHS variables, e.g., G below would be a matrix (instead of a vector)
      x = X(idx_t, span());
      Px = (P * x.t());
      g = Px.each_row() / ( r(idx_r,0) + x * P * x.t());

      // Prediction
      y_pred.slice(idx_r)(idx_t, span()) =  x * B;

      // Predictions if negative values not allowed
      if (allow_neg == false){
        tmp_neg_val = x * B;
        for (idx_y = 0; idx_y < n_y; idx_y++){
          if (tmp_neg_val(0, idx_y) < 0) {
            y_pred.slice(idx_r)(idx_t, idx_y) =  0; // comment_bf: types ok?
          }

        }
      }

      // Forecast errors
      fe_honest.slice(idx_r)(idx_t, span()) = Y(idx_t, span()) - y_pred.slice(idx_r)(idx_t, span());

      // Update regression coefficients and (X'X)^{-1} equivalent
      B = B + g * fe_honest.slice(idx_r)(idx_t, span());
      P = 1/r(idx_r,0) * ( P - g * x * P );
    }
  }

  return Rcpp::List::create(Rcpp::Named("y_pred") = y_pred,
                            Rcpp::Named("fe_honest") = fe_honest,
                            Rcpp::Named("forgetting") = r);
}




//' @rdname rls_exp_cpp
// [[Rcpp::export]]
Rcpp::List rls_window_cpp(const arma::mat Y,
                     const arma::mat X,
                     int ws,
                     bool allow_neg = true,
                     bool debug_flag = false){

  // Rprintf("\n %d \n", i); showMatrix(X, "x: ");

  // dimensions
  int n_obs = Y.n_rows;
  int n_y = Y.n_cols;
  int n_x = X.n_cols;

  // initial estimators for residuals and regression coefficients
  mat Y_init = zeros(ws, n_y);
  mat X_init = zeros(ws, n_x);
  mat P = zeros(n_x, n_x);
  mat A = zeros(n_x, n_y);

  X_init = X(span(0, ws - 1), span(0, n_x - 1));
  Y_init = Y(span(0, ws - 1), span(0, n_y - 1));
  P = pinv(X_init.t() * X_init );
  A = P * X_init.t() * Y_init;

  // initialization for RLS
  mat y_pred = zeros(n_obs, n_y);
  mat fe_honest = zeros(n_obs, n_y);

  // Initializing stuff for for-loop
  int idx_t;
  int idx_y;

  mat xf = zeros(1, n_x);
  mat xb = zeros(1, n_x);
  mat g = zeros(n_x, 1);
  mat Px = zeros(n_x, 1); // corresponds to P(k-1) x(k) in Young p52 (RWP-2)

  mat tmp_neg_val = zeros(1, n_y);

  for (idx_t = ws; idx_t < n_obs; idx_t++){
    // 1) update with new observation

      // Simplifying notation
      xf = X(idx_t, span());
      Px = (P * xf.t());

      // "Kalman-Gain"
      g = Px.each_row() / ( 1 + xf * Px);

      // Update regression coefficients and state covariance
      A = A + g * (Y(idx_t, span(0, n_y - 1)) - xf * A);
      P = P - g * xf * P;

    // 2) Subtracting last observation in (old) window

      // Simplifying notation
      xb = X(idx_t - ws, span());
      Px = P * xb.t();

      // Updates
      g = Px.each_row() / ( -1 + xb * Px);
      A = A + g * (Y(idx_t - ws, span(0, n_y - 1)) - xb * A);
      P = P - g * xb * P;

    // Prediction
    y_pred(idx_t, span()) =  xf * A;

      // Predictions if negative values not allowed
      if (allow_neg == false){
        tmp_neg_val = y_pred(idx_t, span());
        for (idx_y = 0; idx_y < n_y; idx_y++){
          if (tmp_neg_val(0, idx_y) < 0) {
            y_pred(idx_t, idx_y) =  0; // comment_bf: types ok?
          }
        }
      }

    // Forecast errors
    fe_honest(idx_t, span()) = Y(idx_t, span()) - y_pred(idx_t, span());

  }

  return Rcpp::List::create(Rcpp::Named("y_pred") = y_pred,
                            Rcpp::Named("fe_honest") = fe_honest);
}



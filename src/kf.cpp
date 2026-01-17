// kf.cpp
// basiert auf "parametrization/src/kf.cpp".
// ist das die letzte version?

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
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

// via the exports attribute we tell Rcpp to make this function
// available from R



//' @name kf
//' @rdname kf
//' @export
// [[Rcpp::export]]
Rcpp::List kf_cpp(const arma::mat& A, const arma::mat& C,
                  const arma::mat& Q, const arma::mat& R, const arma::mat& S,
                  const arma::mat& y_t, const arma::mat& P1, const arma::colvec& a1) {

  // note y_t is (m X N)!
  unsigned long int m = y_t.n_rows;
  unsigned long int N = y_t.n_cols;
  unsigned long int s = A.n_rows;

  arma::mat K = arma::zeros<arma::mat>(s,m);
  arma::mat tK = K.t();
  arma::mat tC = C.t();
  arma::mat tA = A.t();
  arma::mat tS = S.t();
  arma::mat sigma = C * P1 * tC + R;
  // left Cholesky factor of Sigma[t|t-1]
  arma::mat sigma_L = chol(sigma).t();//' @name solve_rmfd_cpp


  double ll = N * m * log(2 * arma::datum::pi);

  // filtered states
  arma::mat a = arma::zeros<arma::mat>(s,(N+1));
  // scaled residuals, initialize with outputs y!
  arma::mat e = y_t;

  // P[t|t-1]
  arma::mat Pt = P1;

  for (unsigned long int t = 0; t < N; t++) {
    // fprintf(stdout, "t=%3d m=%2d N=%2d\n", t, e.n_rows, e.n_cols);

    // e[t] = Sigma[t|t-1]^(-1/2) ( y[t] - C a[t|t-1] )
    // note e is created with e = y_t = y.t();
    e.col(t) = solve( trimatl(sigma_L), e.col(t) - C * a.col(t));

    // add log det Sigma[t|t-1] to log Likelihood
    ll += 2*log(prod(abs( sigma_L.diag() )));

    // K Sigma[t|t-1]^(-1/2)  is the Kalman gain
    tK = solve( trimatl(sigma_L), C * Pt * tA + tS );
    K = tK.t();

    a.col(t+1) = A * a.col(t) + K * e.col(t);

    Pt = A * Pt * tA + Q - K * tK;
    // make sure that Pt is symetric / is this necessary?
    Pt = (Pt + Pt.t()) / 2;

    sigma = C * Pt * tC + R;
    sigma_L = chol(sigma).t();
  }

  // add residual sum of squares to log Likelihood
  ll += accu(square(e));

  return Rcpp::List::create(Rcpp::Named("e") = e.t(),
                            Rcpp::Named("a") = a.t(),
                            Rcpp::Named("ll") = -ll/(2*N),
                            Rcpp::Named("P1") = Pt);

}

//' @name kf
//' @rdname kf
//' @export
// [[Rcpp::export]]
Rcpp::List kf2_cpp(const arma::mat& A, const arma::mat& C,
                   const arma::mat& H_t,
                   const arma::mat& y_t, arma::mat& P1_R, const arma::colvec& a1) {

  // note y_t is (m x N)!
  unsigned long int m = y_t.n_rows;
  unsigned long int N = y_t.n_cols;
  unsigned long int s = A.n_rows;

  double ll = N * m * std::log(2 * arma::datum::pi);

  // tCA = cbind(t(C), t(A)) ( <=> [C', A'] )
  arma::mat tCA = join_rows(C.t(), A.t());

  // M = rbind( P1_R %*% t(CA), H_t )
  arma::mat M = join_cols( P1_R * tCA, H_t );

  arma::mat L = arma::zeros<arma::mat>(s, s);
  arma::mat sigma_L = arma::zeros<arma::mat>(m, m);

  arma::mat Q, R = arma::zeros<arma::mat>(m+s, m+s);

  // filtered states
  arma::mat a = arma::zeros<arma::mat>(s,(N+1));
  // scaled residuals, initialize with outputs y!
  arma::mat e = y_t;

  for (unsigned long int t = 0; t < N; t++) {
    // QR decomposition of M
    qr_econ(Q, R, M);

    // make sure that diagonal entries are positive
    for (unsigned long int i = 0; i < (m+s); i++) {
      if (R(i,i) < 0) {
        R.row(i) = -R.row(i);
      }
    }

    // (right) Cholesky factor of P[t+1|t]
    P1_R = R.submat(m, m, m+s-1, m+s-1);
    // (left) Cholesky factor of Sigma[t|t-1]
    sigma_L = R.submat(0, 0, m-1, m-1).t();

    // K = L Sigma[t|t-1]^(-1/2) is the Kalman gain
    L = R.submat(0, m, m-1, m+s-1).t();

    // e[t] = Sigma[t|t-1]^(-1/2) ( y[t] - C a[t|t-1] )
    // note e is created with e = y.t()
    e.col(t) = solve( trimatl(sigma_L), e.col(t) - C * a.col(t));
    // state estimate a[t+1|t] = A a[t|t-1] + K ( y[t] - C a[t|t-1] ) = A a[t|t-1] + L e[t]
    a.col(t+1) = A * a.col(t) + L * e.col(t);

    // add log det Sigma[t|t-1] to log Likelihood
    ll += 2*log(prod(abs( sigma_L.diag() )));

    M.rows(0, s-1) = P1_R * tCA;

    // fprintf(stdout, "t=%3d ll=%1.9f\n ", t, ll);
  }

  // add residual sum of squares to log Likelihood
  ll += accu(square(e));

  return Rcpp::List::create(Rcpp::Named("e") = e.t(),
                            Rcpp::Named("a") = a.t(),
                            Rcpp::Named("ll") = -ll/(2*N),
                            Rcpp::Named("P1") = P1_R.t() * P1_R);

}


//' @name ll_kf
//' @rdname ll_kf
//' @export
// [[Rcpp::export]]
double ll_kf_cpp(const arma::mat& A, const arma::mat& C,
                 const arma::mat& Q, const arma::mat& R, const arma::mat& S,
                 const arma::mat& y_t, const arma::mat& P1, const arma::colvec& a1,
                 double tol) {

  // note y_t is (m X N)!
  unsigned long int m = y_t.n_rows;
  unsigned long int N = y_t.n_cols;
  unsigned long int s = A.n_rows;
  // unsigned long int t0;

  arma::mat K = arma::zeros<arma::mat>(s, m);
  arma::mat tK = K.t();
  arma::mat tC = C.t();
  arma::mat tA = A.t();
  arma::mat tS = S.t();
  arma::mat sigma = C * P1 * tC + R;
  // left Cholesky factor of Sigma[t|t-1]
  arma::mat sigma_L = chol(sigma).t();

  double ll;
  // the case s=0 is simple, no need to run the Kalman filter
  if (s == 0) {
    ll = m * log(2 * arma::datum::pi);
    ll = ll + 2*log( prod( abs( sigma_L.diag() ) ) );
    ll = ll + accu( square( solve( trimatl(sigma_L), y_t) ) ) / N;
    ll = -0.5*ll;
    return(ll);
  }

  ll = N * m * log(2 * arma::datum::pi);

  // filtered states
  arma::colvec a = a1;
  // scaled residuals
  arma::colvec e;

  // P[t|t-1]
  arma::mat Pt = P1;

  unsigned long int t = 0;
  bool converged = false;
  arma::mat Pt0 = arma::zeros<arma::mat>(s, s);
  arma::mat K0 = arma::zeros<arma::mat>(s, m);

  while (!converged) {
    Pt0 = Pt;
    K0 = K;

    // e[t] = Sigma[t|t-1]^(-1/2) ( y[t] - C a[t|t-1] )
    e = solve( trimatl(sigma_L), y_t.col(t) - C * a);

    // add log det Sigma[t|t-1] and sum of squared, scaled resiuals to log Likelihood
    ll += 2*log(prod(abs( sigma_L.diag() ))) + accu(square(e));

    // K Sigma[t|t-1]^(-1/2)  is the Kalman gain
    tK = solve( trimatl(sigma_L), C * Pt * tA + tS );
    K = tK.t();

    a = A * a + K * e;

    Pt = A * Pt * tA + Q - K * tK;
    // make sure that Pt is symetric / is this necessary?
    Pt = (Pt + Pt.t()) / 2;

    sigma = C * Pt * tC + R;
    sigma_L = chol(sigma).t();

    t++;
    // fprintf(stdout, "t=%3d ll=%1.9f diff=%1.9f | %1.9f\\n ", t, ll,
    //        (abs(K - K0)).max(), (abs(Pt - Pt0)).max());

    // to be sure check P and K!
    converged = ((t == N) or ( ((abs(Pt - Pt0)).max() < tol) and ((abs(K - K0)).max() < tol) ));
  }

  // unsigned long int t0 = t;
  while (t < N) {
    // finish with fixed Kalman gain K

    // e[t] = Sigma[t|t-1]^(-1/2) ( y[t] - C a[t|t-1] )
    e = solve( trimatl(sigma_L), y_t.col(t) - C * a);

    ll += 2*log(prod(abs( sigma_L.diag() ))) + accu(square(e));

    a = A * a + K * e;

    t++;
  }


  return (-ll/(2*N));
}


//' @name ll_kf
//' @rdname ll_kf
//' @export
// [[Rcpp::export]]
double ll_kf2_cpp(arma::mat& A, arma::mat& C, arma::mat& H_t,
                  arma::mat& y_t, arma::mat& P1_R, arma::colvec& a1,
                  double tol) {

  // ACHTUNG: y_t is (m x N)!
  unsigned long int m = y_t.n_rows;
  unsigned long int N = y_t.n_cols;
  unsigned long int s = A.n_rows;
  // unsigned long int t0;

  double ll;

  // the case s=0 is simple, no need to run the Kalman filter
  if (s == 0) {
    // left Cholesky factor of Sigma = var(y[t])
    arma::mat sigma_L = H_t.t() * H_t;
    sigma_L = chol(sigma_L).t();

    ll = m * log(2 * arma::datum::pi);
    ll = ll + 2*log( prod( abs( sigma_L.diag() ) ) );
    ll = ll + accu( square( solve( trimatl(sigma_L), y_t) ) ) / N;
    ll = -0.5*ll;
    return(ll);
  }

  ll = N * m * log(2 * arma::datum::pi);

  // tCA = cbind(t(C), t(A)) ( <=> [C', A'] )
  arma::mat tCA = join_rows(C.t(), A.t());

  // M = rbind( P1_R %*% t(CA), H_t )
  arma::mat M = join_cols( P1_R * tCA, H_t );

  arma::mat L = arma::zeros<arma::mat>(s, s);
  arma::mat sigma_L = arma::zeros<arma::mat>(m, m);

  arma::mat Q, R = arma::zeros<arma::mat>(m+s, m+s);

  // filtered states
  arma::colvec a = a1;
  // scaled residuals
  arma::colvec e;


  unsigned long int t = 0;
  bool converged = false;
  arma::mat P1_R0 = arma::zeros<arma::mat>(s, s);

  while (!converged) {
    P1_R0 = P1_R;

    // QR decomposition of M
    qr_econ(Q, R, M);

    // (right) Cholesky factor of P[t+1|t]
    P1_R = R.submat(m, m, m+s-1, m+s-1);
    // (left) Cholesky factor of Sigma[t|t-1]
    sigma_L = R.submat(0, 0, m-1, m-1).t();

    // K = L Sigma[t|t-1]^(-1/2) is the Kalman gain
    L = R.submat(0, m, m-1, m+s-1).t();

    // e[t] = Sigma[t|t-1]^(-1/2) ( y[t] - C a[t|t-1] )
    e = solve( trimatl(sigma_L), y_t.col(t) - C * a);
    // state estimate a[t+1|t] = A a[t|t-1] + K ( y[t] - C a[t|t-1] ) = A a[t|t-1] + L e[t]
    a = A * a + L * e;

    // add log det Sigma[t|t-1] to log Likelihood
    ll += 2*log(prod(abs( sigma_L.diag() ))) + accu(square(e));

    M.rows(0, s-1) = P1_R * tCA;

    t++;
    // fprintf(stdout, "t=%3d ll=%1.9f diff=%1.9f\n ", t, ll, (abs(P1_R - P1_R0)).max());

    converged = ((t == N) or ((abs(P1_R - P1_R0)).max() < tol));
  }

  // unsigned long int t0 = t;
  while (t < N) {
    // finish with fixed Kalman gain K

    // e[t] = Sigma[t|t-1]^(-1/2) ( y[t] - C a[t|t-1] )
    e = solve( trimatl(sigma_L), y_t.col(t) - C * a);
    // state estimate a[t+1|t] = A a[t|t-1] + K ( y[t] - C a[t|t-1] ) = A a[t|t-1] + L e[t]
    a = A * a + L * e;

    ll += 2*log(prod(abs( sigma_L.diag() ))) + accu(square(e));

    t++;
  }

  return (-ll/(2*N));
}

//' Compute the log likelihood for a statespace system
//' described by a model template.
//'
//' This is an internal helper function, used by the function factory \code{\link{ll_FUN}}. For a more detailed
//' documentation of the log Likelihood, see \code{\link{ll_kf}}.
//'
//' @param theta \eqn{(K)} dimensional vector of "deep" parameters.
//' @param y \eqn{(m,N)} matrix with the observed outputs:
//'        \eqn{(y_1,y_2,\ldots,y_N)}{(y[1],y[2],...,y[N])}.
//' @param SYS \eqn{(m+s,m+s)} matrix, is overwritten with the system matrix
//'        \eqn{[A,B | C,D]}.
//' @param H_SYS \eqn{(m+s)^2, K)} matrix.
//' @param h_SYS \eqn{((m+s)^2)}-dimensional vector. Note that \code{vec(SYS) = H_SYS*theta + h_SYS}.
//' @param sigma_L \eqn{(m,m)} matrix, is overwritten with the left square root of the
//'        noise covariance matrix.
//' @param H_sigma_L \eqn{(m^2, K)} matrix.
//' @param h_sigma_L \eqn{(m^2)}-dimensional vector. Note that
//'        \code{vec(sigma_L) = H_sigma_L*theta + h_sigma_L}.
//' @param VAR \eqn{(m+s,m+s)} matrix, is overwritten with the
//'        covariance matrix \eqn{[Q,S | S',R] = [B | C] sigma_L sigma_L' [B', C']}
//' @param P1 \eqn{(s,s)} matrix, is overwritten with the
//'        initial state covariance matrix (computed via a Lyapunov equation).
//' @param tol (double) tolerance used by ll_kf_cpp.
//' @param err (double) return err, if the computation of P1 fails.
//' @export
// [[Rcpp::export]]
double ll_kf_theta_cpp(const arma::vec& theta, const arma::mat& y,
                       arma::mat& SYS, const arma::mat& H_SYS, const arma::vec& h_SYS,
                       arma::mat& sigma_L, const arma::mat& H_sigma_L, const arma::vec& h_sigma_L,
                       arma::mat& VAR, arma::mat& P1,
                       double tol, double err) {

  unsigned long int s = P1.n_rows;
  unsigned long int m = y.n_rows;

  arma::vec a1(s);
  a1.zeros();
  arma::vec lambda_r(s);
  arma::vec lambda_i(s);

  // http://arma.sourceforge.net/docs.html#memptr
  // vec(ptr_aux_mem, number_of_elements, copy_aux_mem = true, strict = false)
  // SYS_vec is the vectorised version of the matrix SYS (shares the memory with SYS)
  arma::vec SYS_vec = arma::vec(SYS.memptr(), (m+s)*(m+s), false, true);
  // sigma_L_vec is the vectorised version of the matrix sigma_L (shares the memory with sigma_L)
  arma::vec sigma_L_vec = arma::vec(sigma_L.memptr(), (m)*(m), false, true);

  SYS_vec = h_SYS + H_SYS * theta;
  sigma_L_vec = h_sigma_L + H_sigma_L * theta;

  // A = SYS.submat(0,0,s-1,s-1);
  // B = SYS.submat(0,s,s-1,s+m-1);
  // C = SYS.submat(s,0,s+m-1,s-1);
  // D = SYS.submat(s,s,s+m-1,s+m-1);

  // VAR is the covariance matrix [Q  S]
  //                              [S' R]
  // R = D * sigma * D'
  // S = B * sigma * D'
  // Q = B * sigma * B'
  SYS.cols(s, s+m-1) = SYS.cols(s, s+m-1) * sigma_L;
  VAR = SYS.cols(s, s+m-1) * SYS.cols(s, s+m-1).t();

  if (s > 0) {
    // Using header-only implementation from rationalmatrices for cross-package linking
    bool ok = rationalmatrices::lyapunov_impl(SYS.submat(0,0,s-1,s-1), VAR.submat(0,0,s-1,s-1),
                                              P1, lambda_r, lambda_i, true);
    if (!ok) {
      // computation of P1 fails if the system is not stable!!!
      return(err);
    }
  }

  // Rcout << SYS.submat(0,0,s-1,s-1) << VAR.submat(0,0,s-1,s-1) << P1 << std::endl;


  double ll;
  if (s > 0) {
    //             A,                       C
    ll = ll_kf_cpp(SYS.submat(0,0,s-1,s-1), SYS.submat(s,0,s+m-1,s-1),
    //             Q                        R                            S
                   VAR.submat(0,0,s-1,s-1), VAR.submat(s,s,s+m-1,s+m-1), VAR.submat(0,s,s-1,s+m-1),
                   y, P1, a1, tol);
  } else {
    //             A,              C
    ll = ll_kf_cpp(arma::mat(0,0), arma::mat(m,0),
    //             Q               R    S
                   arma::mat(0,0), VAR, arma::mat(0,m),
                   y, P1, a1, tol);
  }

  return(ll);
}

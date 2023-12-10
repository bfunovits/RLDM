// #include "RcppArmadillo.h"
// using namespace Rcpp;
//
// //' @rdname logLik_RMFD_cpp
// //' @export
// //' @keywords internal
// // [[Rcpp::export]]
// double logLik_RMFD_cpp(const arma::vec& th, int m, int p, int q, int skip, // th = deep parameters / m = dimension of outputs = dim of inputs // (p,q) = degrees of a(z), b(z)
//                        const arma::vec& hb0, const arma::mat& Hb0, arma::mat& b0,
//                        const arma::vec& hb, const arma::mat& Hb, arma::mat& b,
//                        const arma::vec& ha, const arma::mat& Ha, arma::mat& a0, arma::mat& a,
//                        const arma::vec& hSigmaL, const arma::mat& HSigmaL, arma::mat& Sigma, // these parameters do not appear in the conditional likelihood below
//                        const arma::mat& y, arma::mat& e, arma::mat& S) {
//
//   // more comments are in the conditional log-likelihood function below!
//   int nobs = y.n_cols;
//   int nvalid = nobs - skip;
//
//   // vec(ptr_aux_mem, number_of_elements, copy_aux_mem = true, strict = false)
//   arma::vec b0vec = arma::vec(b0.memptr(), m*m, false, true);
//   arma::vec bvec = arma::vec(b.memptr(), m*m*q, false, true);
//   arma::vec avec = arma::vec(a.memptr(), m*m*(p+1), false, true);
//   arma::vec Sigmavec = arma::vec(Sigma.memptr(), m*m, false, true);
//
//   // construct system matrices for given parameters th
//   b0vec = hb0 + Hb0 * th;
//   bvec  = hb  + Hb * th;
//   avec  = ha  + Ha * th;
//   b = solve(b0, b);
//   a = solve(b0, a);
//   Sigmavec  = hSigmaL  + HSigmaL * th;
//
//   a0 = a.head_cols(m);
//   double lndeta0;
//   double sign;
//   log_det(lndeta0, sign, a0);
//
//   // noise covariance matrix sigma_L * t(sigma_L)
//   Sigma = Sigma * Sigma.t();
//
//   double lndetSigma;
//   log_det(lndetSigma, sign, Sigma);
//
//   // compute residuals
//   solve_ARMA_cpp(b, a, y, e, 1);
//   // noise covariance matrix S
//   S = (e.tail_cols(nvalid) * e.tail_cols(nvalid).t()) / nvalid;
//
//   // THAT'S THE DIFFERENCE TO THE CONDITIONAL LIKELIHOOD!
//   // trace ( S * inv(Sigma) )
//   Sigma = solve(Sigma, S);
//   double trSS = accu( Sigma.diag() );   // http://arma.sourceforge.net/docs.html#accu
//
//   // FOR DEBUGGING
//   // Rcout << " nvalid: " << nvalid << " mln2pi: " << m*log(2*(arma::datum::pi)) << " trSS: " << trSS << " lndetSigma: "  << lndetSigma  << " lndeta0: "  << -2*lndeta0 << std::endl;
//
//   double ll = (-0.5) * (m*log(2*(arma::datum::pi)) + trSS + lndetSigma - 2*lndeta0);
//
//   return ll;
// }
//
// //' @rdname logLik_ARMA_cpp
// //' @export
// //' @keywords internal
// // [[Rcpp::export]]
// double ClogLik_RMFD_cpp(const arma::vec& params_deep,
//                         int dim_out, int dim_in, int deg_c, int deg_d, int skip,
//                         const arma::vec& hc, const arma::mat& Hc, arma::mat& c,
//                         const arma::vec& hd0, const arma::mat& Hd0, arma::mat& d0,
//                         const arma::vec& hd, const arma::mat& Hd, arma::mat& d,
//                         const arma::mat& y,
//                         arma::mat& e, arma::mat& S) {
//
//   // Number of observations: Note that each observation corresponds to one COLUMN (to optimize memory access)
//   int nobs = y.n_cols;
//   int nvalid = nobs - skip;
//
//   // Create vectors for the memory address of the system parameters (in which to map the deep parameters with an affine trafo!)
//   //   Seems to be quite efficient and elegant!
//   //
//   // http://arma.sourceforge.net/docs.html#memptr
//   // http://arma.sourceforge.net/docs.html#Col (see "Advanced constructors")
//   //   vec(ptr_aux_mem, number_of_elements, copy_aux_mem = true, strict = false)
//   arma::vec d0_vec = arma::vec(d0.memptr(), dim_out*dim_in, false, true);
//   arma::vec d_vec = arma::vec(d.memptr(), dim_out*dim_in*deg_d, false, true);
//   arma::vec c_vec = arma::vec(c.memptr(), dim_in*dim_in*(deg_c+1), false, true);
//
//   // construct system matrices for given "deep" parameters th: Apply affine map (deep_params -> lin_params = h + H * deep_params)
//   c_vec  = hc  + Hc * params_deep;
//   d0_vec = hd0 + Hd0 * params_deep;
//   d_vec  = hd  + Hd * params_deep;
//   b = solve(b0, b); // arma::solve()? http://arma.sourceforge.net/docs.html#solve
//   a = solve(b0, a);
//
//   // Calculate log(det(a_0)) which appears through the Jacobian in the density transformation
//   a0 = a.head_cols(m); // http://arma.sourceforge.net/docs.html#submat
//   double lndeta0;
//   double sign;
//   log_det(lndeta0, sign, a0); // http://arma.sourceforge.net/docs.html#log_det: determinant is equal to exp(val)*sign, where val is obtained from log_det(val, sign, A)
//
//   // compute residuals: For given data y, calculate residuals e
//   //
//   // Note that the order of coefficients in the (wide) AR matrix are reversed. The
//   // solve_rmfd_cpp(const arma::mat& poly_inv, const arma::mat& poly_fwd, arma::mat& data_in, arma::mat& data_out, int t0) {
//
//   solve_ARMA_cpp(b, a, y, e, 1); // Arguments: (AR, MA, inputs, outputs, starting time t_0). Only outputs are not "const"s and are overwritten
//   // noise covariance matrix S
//   S = (e.tail_cols(nvalid) * e.tail_cols(nvalid).t()) / nvalid;
//
//   double lndetSigma;
//   log_det(lndetSigma, sign, S);
//
//   // Rcout << "nobs: " << nobs  << " nvalid: " << nvalid << " mln2pi: " << m*log(2*(arma::datum::pi)) << " trSS: " << trSS << " lndetS: "  << lndetSigma << std::endl;
//
//   double ll = (-0.5) * (m*log(2*(arma::datum::pi)) + m + lndetSigma - 2*lndeta0); // m is the output dimension
//
//   return ll;
// }
//
//

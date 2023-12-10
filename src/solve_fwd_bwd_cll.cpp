// ARMA_methods.cpp
// new versions 9.9.20
// replaces 'residualsARMA.cpp' and 'solveARMA.cpp'

// STSP_methods.cpp
// new versions 8.9.20
// replaces 'solveSTSP.cpp' and 'solve_STP_new.cpp'

// -----------------------------------------------------------------------------

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
// [[Rcpp::depends(RcppArmadillo)]]

// To make your C++ code callable from C++ code in other packages.
// This will generate a header file, inst/include/mypackage.h that
// can be included by other packages
// [[Rcpp::interfaces(r, cpp)]]

// -----------------------------------------------------------------------------
// 1) Compute Outputs (solve_de) ----
// -----------------------------------------------------------------------------

//' Outputs of an ARMA systems
//'
//' This internal helper function computes the outputs of an ARMA system
//' \deqn{a_0 y_t + a_1 y_{t-1} + \cdots + a_p y_{t-p} = b_0 u_t + \cdots + b_q u_{t-q}}{
//'       a[0] y[t] + a[1] y[t-1] + ... + a[p] y[t-p] = b[0] u[t] + ... + b[q] u[t-q]}
//'
//' Values \eqn{y_t}{y[t]}, \eqn{u_t}{u[t]} for \eqn{t\leq 0}{t\le 0} are implicitly set to be zero.
//' However, by starting the iteration with some \eqn{t_0>1}{t0>1} we can enforce non-zero
//' initial values.
//'
//' @note
//' Use this procedure with care!
//'
//' * The procedure does \bold{not} check the input arguments. We require \eqn{m > 0},
//'   \eqn{p \geq 0}{p \ge 0}, \eqn{n(q+1) \geq 0}{n(q+1)\ge 0} and
//'   \eqn{1 \leq t_0 \leq N}{1 \le t_0  \le N}.
//' * The procedure \bold{overwrites} the input argument \code{y}.
//' * The data matrices are organized columnwise (to avoid memory shuffling)!
//' * Note also the non standard representation of the coefficient matrices.
//'
//' @param A1 \eqn{(m, mp)} matrix \eqn{-a_0^{-1}(a_p,...,a_1)}{-a[0]^{-1}(a[p],...,a[1])}.
//' @param B \eqn{(m, n(q+1))} matrix \eqn{a_0^{-1}(b_0,...,b_q}{a[0]^{-1}(b[0],...,b[q])}.
//' @param u \eqn{(n, N)} matrix with the inputs \eqn{(u_1,...,u_N}{(u[1],...,u[N])}.
//' @param y \eqn{(m, N)} matrix with the outputs \eqn{(y_1,...,y_N}{(y[1],...,y[N])}.
//' @param t0 integer, start iteration at t = t0.
//'
//' @return This RcppArmadillo routine returns \code{NULL} but \bold{overwrites}
//'         the input argument \code{y} with the computed outputs!
//'
//' @export
//'
//' @rdname outputs_ARMA_cpp
//' @name outputs_ARMA_cpp
//'
//' @seealso \code{\link{outputs_ARMA_cpp}}, \code{\link{residuals_ARMA_cpp}},  \code{\link{cll_theta_ARMA_cpp}},
//'    \code{\link{outputs_STSP_cpp}}, \code{\link{residuals_STSP_cpp}},  \code{\link{cll_theta_STSP_cpp}} and
//'    \code{\link{solve_de}}, \code{\link{solve_inverse_de}} and \code{\link{ll}}.
//'
//' @examples
//' # generate a random ARMA(2,1) model (3 outputs, 2 inputs)
//' p = 2
//' q = 1
//' m = 3
//' n = 2
//' model = test_armamod(dim = c(m, n), degrees = c(p,q), digits = 2)
//' A = unclass(model$sys$a)
//' a0 = A[,,1]
//' A1 = -A[,,(p+1):2]
//' dim(A1) = c(m, m*p)
//' A1 = solve(a0, A1)
//' B = unclass(model$sys$b)
//' dim(B) = c(m, n*(q+1))
//' B = solve(a0, B)
//'
//' # generate random noise sequence (sample size N = 10)
//' n.obs = 10
//' u = matrix(rnorm(n.obs*n), nrow = n, ncol = n.obs)
//' print(u)
//'
//' # generate matrix for the outputs
//' y = matrix(0, nrow = m, ncol = n.obs)
//'
//' # call outputs_ARMA_cpp()
//' outputs_ARMA_cpp(A1, B, t0 = 2, u, y) # start with t>=2
//' print(u)
//' print(y)  # y is overwritten with the computed outputs
// [[Rcpp::export]]
void outputs_ARMA_cpp(const arma::mat& A1, const arma::mat& B, int t0,
                      const arma::mat& u, arma::mat& y) {

  // check 1 <= t0 <= N
  int N = y.n_cols;
  if ((N == 0) || (t0 > N)) {
    // nothing to do
    return;
  }
  if (t0 < 1) {
    stop("t0 < 1!");
  }

  // check m >= 1, p >= 0
  int m = y.n_rows;
  if (m < 1) {
    stop("m < 1!");
  }
  int p = A1.n_cols / m;
  int Nm = N*m;

  // check n >= 0, (q+1) >= 0
  int n = u.n_rows;
  int nq1 = B.n_cols; // n(q+1)
  int q;
  if (n > 0) {
    q = (nq1 / n) - 1;
  } else {
    q = -1;
  }

  int k, miss, i1, i2, j1, j2;

  t0 = t0 - 1; // use the C++ indexing style !!!

  // Rcout << "outputs_ARMA_cpp m:" << m << " n:" << n << " q:" << q << " p:" << p << " Nm:" << Nm << " N:" << N << std::endl;

  // "MA" part: y[t] = b[0] u[t] + ... + q[p] u[t-q]
  if (nq1 > 0) {
    y.cols(t0, (N-1)) = B.cols(0, n - 1) * u.cols(t0, (N-1));  // y[t] = b[0]*u[t]
    k = 1;
    while ((k <= q) && (k < N)) {
      // y[t] = y[t] + b[k]*u[t-k]
      miss = std::max(0, k-t0);     // take care of missing initial values!
      if ((t0+miss) <= (N-1)) {
        y.cols(t0+miss, N-1) = y.cols(t0+miss, N-1) + B.cols(k*n, (k+1)*n-1) * u.cols(t0-k+miss, N-1-k);
      }
      ++k;
    }
  } else {
    // no inputs!
    y.tail_cols(N-t0).zeros(); // set y[t] = for t>t0
  }

  // "AR" part: solve y[t] = a[1] y[t-1] + ... + a[p] y[t-p]
  if (p > 0) {
    // http://arma.sourceforge.net/docs.html#memptr
    // vec(ptr_aux_mem, number_of_elements, copy_aux_mem = true, strict = false)
    // y_vec is the vectorised version of the matrix y (shares the memory with y)
    // Rcout << "outputs_ARMA_cpp Nm" << Nm << std::endl;
    // Rcout << y << std::endl;
    arma::vec y_vec = arma::vec(y.memptr(), Nm, false, true);
    // Rcout << y_vec << std::endl;

    i1 = t0*m;
    i2 = i1 + m - 1;
    j1 = i1 - m*p;
    j2 = i1 - 1;

    // take care of missing initial values for t <= p !!!
    while ((j1 < 0) && (i2 < Nm)) {
      miss = std::max(0, -j1);

      // Rcout << "i: " << i1 << ":" << i2 << " j: " << j1 << ":" << j2 << std::endl;
      // Rcout << "A1.tail: " << (m*p-miss) << " j: " << j1 + miss << ":" << j2 << std::endl;
      if (j2 >= 0) {
        y_vec.subvec(i1,i2) = y_vec.subvec(i1,i2) + A1.tail_cols(m*p - miss) * y_vec.subvec(j1 + miss, j2);
      }

      i1 += m;
      i2 += m;
      j1 += m;
      j2 += m;
    }

    // t > q
    while (i2 < Nm) {
      // Rcout << "i: " << i1 << ":" << i2 << " j: " << j1 << ":" << j2 << std::endl;

      y_vec.subvec(i1,i2) = y_vec.subvec(i1,i2) + A1 * y_vec.subvec(j1,j2);

      i1 += m;
      i2 += m;
      j1 += m;
      j2 += m;
    }
  }

  return;
}

//' Outputs of a statespace system
//'
//' @description
//' This internal helper function computes the outputs and states
//' for a statespace system of the form
//' \deqn{a_{t+1} = A a_t + B u_t, \; y_t   = C a_t + D u_t}{
//'       a[t+1] = A a[t] + B u[t] and y[t]   = C a[t] + D u[t]}
//'
//' @note
//' Use this procedure with care!
//'
//' * The procedure does \bold{not} check the input arguments.
//' * The procedure \bold{overwrites} the input arguments \code{a} and \code{u}.
//' * The data matrices are organized columnwise (to avoid memory shuffling)!
//'
//' @param A \eqn{(s,s)} matrix.
//' @param B \eqn{(s,n)} matrix.
//' @param C \eqn{(m,s)} matrix.
//' @param D \eqn{(m,n)} matrix.
//' @param u \eqn{(n,N)} matrix with the inputs/disturbances:
//'          \eqn{(u_1,u_2,\ldots,u_N)}{(u[1],u[2],...,u[N])}.
//' @param a \eqn{(s,N+1)} matrix. This matrix is overwritten with the (computed) states:
//'            \eqn{(a_1,a_2,\ldots,a_N,a_{N+1})}{(a[1],a[2],\ldots,a[N],a[N+1])}.
//'            On input \code{a[,1]} must hold the initial state \eqn{a_1}{a[1]}.
//' @param y \eqn{(m,N)} matrix. This matrix is overwritten with (computed) outputs:
//'           \eqn{(y_1,y_2,\ldots,y_N)}{(y[1],y[2],...,y[N])}.
//'
//' @return This RcppArmadillo routine returns \code{NULL} but \bold{overwrites} the input
//'         arguments \code{a} and \code{u}!
//' @export
//'
//' @rdname outputs_STSP_cpp
//' @name outputs_STSP_cpp
//'
//' @seealso \code{\link{outputs_ARMA_cpp}}, \code{\link{residuals_ARMA_cpp}},  \code{\link{cll_theta_ARMA_cpp}},
//'    \code{\link{outputs_STSP_cpp}}, \code{\link{residuals_STSP_cpp}},  \code{\link{cll_theta_STSP_cpp}} and
//'    \code{\link{solve_de}}, \code{\link{solve_inverse_de}} and \code{\link{ll}}.
//'
//' @examples
//' # generate a random statespace model (3 outputs, 2 inputs and 4 states)
//' m = 3
//' n = 2
//' s = 4
//' model = test_stspmod(dim = c(m, n), s = s, digits = 2)
//'
//' # generate random noise sequence (sample size N = 10)
//' n.obs = 10
//' u = matrix(rnorm(n.obs*n), nrow = n, ncol = n.obs)
//' print(u)
//'
//' # generate matrix for the state sequence
//' a = matrix(0, nrow = s, ncol = n.obs+1)
//' a[,1] = rnorm(s) # random initial state a[1]
//' print(a)
//'
//' # generate matrix for the outputs
//' y = matrix(0, nrow = m, ncol = n.obs)
//'
//' # call outputs_STSP_cpp()
//' outputs_STSP_cpp(model$sys$A, model$sys$B, model$sys$C, model$sys$D, u, a, y)
//' print(u)
//' print(a)  # a is overwritten with the computed states
//' print(y)  # y is overwritten with the computed outputs

// [[Rcpp::export]]
void outputs_STSP_cpp(const arma::mat& A, const arma::mat& B,
                       const arma::mat& C, const arma::mat& D,
                       const arma::mat& u, arma::mat& a, arma::mat& y) {

   // NOTE: if y,u have different dimensions then
   //       u cannot be overwritten with the computed outputs y.
   unsigned long int nobs = u.n_cols;
   unsigned long int s = A.n_rows;
   unsigned long int m = D.n_rows;
   unsigned long int n = D.n_cols;

   // take care of the cases s=0, m=0, n=0
   // is this really necessary?
   if (s > 0) {
     if (n > 0) {
       a.tail_cols( nobs ) = B * u;
     } else {
       a.tail_cols( nobs ).zeros();
     }
     for (unsigned long int t = 0; t < nobs; t++) {
       a.col(t+1) = A * a.col(t) + a.col(t+1);
     }
   }

   if (m > 0) {
     if (s > 0 ) {
       if (n > 0) {
         y = C * a.head_cols(nobs) + D*u;
       } else {
         y = C * a.head_cols(nobs);
       }
     } else {
       if (n > 0) {
         y = B*u;
       } else {
         y.zeros();
       }
     }
   }

   return;
 }


//' Forward-backward solution of statespace systems
//'
//' @description
//' **DEPRECATED**?
//' This internal helper function computes the outputs of an, in general **unstable**, statespace system
//' \deqn{a_{t+1} = A a_t + B u_t, \; y_t = C a_t + D u_t}{
//'       a[t+1] = A a[t] + B u[t] and y[t] = C a[t] + D u[t]}
//' by forward and backward recursion. The procedure assumes that the state transition matrix
//' \eqn{A} is block upper triangular, where the upper block \eqn{A_{11}}{A[11]} is stable (i.e.
//' all eigenvalues have moduli less than one) and the
//' lower block \eqn{A_{22}}{A[22]} is unstable (i.e. all eigenvalues have moduli larger than one).
//' This function is mainly used in the routine \code{\link{innovation_form}}.
//'
//' @note
//' Use this procedure with care!
//'
//' * The procedure does \bold{not} check the input arguments. We require \eqn{m > 0}, \eqn{n > 0}.
//'   Furthermore it is assumed that the state transition matrix \eqn{A} is block upper triangular,
//'   as explained above.
//' * The procedure \bold{overwrites} the input arguments \code{y}, \code{as} and \code{au}.
//' * The data matrices are organized columnwise (to avoid memory shuffling)!
//'
//' @param A \eqn{(s, s)} matrix.
//' @param B \eqn{(s, n)} matrix.
//' @param C \eqn{(m, s)} matrix.
//' @param D \eqn{(m, n)} matrix.
//' @param u \eqn{(n, N)} matrix with the inputs \eqn{(u_1,...,u_N}{(u[1],...,u[N])}.
//' @param au \eqn{(su,N+1)} matrix. This matrix is **overwritten** with the (computed)
//'           states of the unstable part of the system.
//'            \eqn{(a_{u1},a_{u2},\ldots,a_{uN},a_{u,N+1})}{(au[1],au[2],\ldots,au[N],au[N+1])}.
//'            On input \code{au[,N+1]} must hold the "initial" state \eqn{a_{u,N+1}}{au[N+1]}.
//' @param as \eqn{(ss,N+1)} matrix. This matrix is **overwritten** with the (computed)
//'           states of the stable part of the system.
//'            \eqn{(a_{s1},a_{s2},\ldots,a_{sN},a_{s,N+1})}{(as[1],as[2],\ldots,as[N],as[N+1])}.
//'            On input \code{as[,1]} must hold the "initial" state \eqn{a_{s1}}{as[1]}.
//' @param y \eqn{(m,N)} matrix. This matrix is **overwritten** with (computed) outputs:
//'           \eqn{(y_1,y_2,\ldots,y_N)}{(y[1],y[2],...,y[N])}.
//'
//' @return This RcppArmadillo routine returns \code{NULL} but \bold{overwrites}
//'         the input argument \code{y}, \code{au} and \code{as} with the computed outputs
//'         and states!
//'
//' @seealso \code{\link{outputs_STSP_cpp}} and \code{\link{innovation_form}}.
//'
//' @export
//' @name fbsolve_STSP_cpp
//' @rdname fbsolve_STSP_cpp
//'

// [[Rcpp::export]]
void fbsolve_STSP_cpp(const arma::mat& A, const arma::mat& B, const arma::mat& C, const arma::mat& D,
                      const arma::mat& u, arma::mat& au, arma::mat& as, arma::mat& y) {
  int s_stable = as.n_rows;
  int s_unstable = au.n_rows;
  int s = s_stable + s_unstable;
  int N = u.n_cols;

  if (N == 0) return;
  if ((u.n_rows == 0) || (y.n_rows == 0)) stop("m=0 or n=0");

  int t;
  arma::mat A11, A12, A22, C1, C2, B1, B2;

  // y[t] = C1 as[t]+ C2 au[t] + D u[t]
  y = D * u;

  if (s_stable > 0) {
    A11 = A.submat(0, 0, s_stable-1, s_stable-1);
    B1  = B.rows(0, s_stable-1);
    C1  = C.cols(0, s_stable-1);
    if (s_unstable > 0) {
      A12 = A.submat(0, s_stable, s_stable-1, s-1);
    }
  }

  if (s_unstable > 0) {
    // system is not stable

    A22 = arma::inv(   A.submat(s_stable, s_stable, s-1, s-1) );
    B2  = arma::solve(-A.submat(s_stable, s_stable, s-1, s-1), B.rows(s_stable, s-1));
    C2 = C.cols(s_stable, s-1);

    // backward recursion for the unstable part
    // au = (au[1],...,au[N],au[N+1])
    int t = N-1;
    while (t >= 0) {
      // au[t] = (A22^(-1)) au[t+1] + (-A22^(-1)*B2) u[t]
      au.col(t) = A22*au.col(t+1) + B2 * u.col(t);
      // y[t] = C1 as[t]+ C2 au[t] + D u[t]
      y.col(t) = C2 * au.col(t) + y.col(t);
      t--;
    }

    if (s_stable > 0) {
      // forward recursion for the stable part
      // as = (as[1],...,as[N],as[N+1])
      t = 0;
      while(t < N) {
        // as[t+1] = A11 as[t] + A12 au[t] + B1 u[t]
        as.col(t+1) = A11 * as.col(t) + A12 * au.col(t) + B1 * u.col(t);
        // y[t] = C1 as[t]+ C2 au[t] + D u[t]
        y.col(t) = C1 * as.col(t) +  y.col(t);
        t++;
      }
    }
  } else {
    // inverse system is stable!

    if (s_stable > 0) {
      // forward recursion for the stable part
      // as = (as[1],...,as[N],as[N+1])
      t = 0;
      while(t < N) {
        // as[t+1] = A11 as[t] + B1 u[t]
        as.col(t+1) = A11 * as.col(t) + B1 * u.col(t);
        // y[t] = C1 as[t] + D u[t]
        y.col(t) = C1 * as.col(t) +  y.col(t);
        t++;
      }
    }
  }

  return;
}

//' Simulating Output from an RMFD Model (or obtain residuals)
//'
//' This RcppArmadillo function calculates for given inputs \code{data_in}
//' of dimension \eqn{(n \times nobs)}{(n x nobs)},  where \eqn{nobs} is the sample size,
//' the outputs of dimension \eqn{m}.
//' Note that data matrices are "transposed" in the sense that every column
//' corresponds to one observation because of memory management.
//' \code{data_out} is thus of dimension \eqn{(m x nobs)}.
//' This function is intended for internal use and thus arguments are not checked.
//'
//' @param poly_inv Matrix of dimension \eqn{(n \times n p)}{(n x n p)},
//'     representing a square matrix polynomial \eqn{c(z)} with \eqn{c_0}{c[0]} equal to the identity matrix (and therefore not stored).
//'     The coefficients need to be in \strong{reverse} direction, i.e. \eqn{(c_p, ... , c_1)}{(c[p], ... , c[1])},
//'     where \eqn{p} denotes the degree of \eqn{c(z)}.
//' @param poly_fwd Matrix of dimensions \eqn{(m \times n(q+1))}{(m x n(q+1))},
//'    representing a (possibly tall) matrix polynomial \eqn{d(z)} of dimension
//'    \eqn{(m \times n)}{(m x n)}, where \eqn{m \geq n}{m \ge n}.
//'    The coefficient are stored "as usual" and including \eqn{d_0}{d[0]}, i.e. \eqn{(d_0, d_1, ... , d_{q-1}, d_{q})}{(d[0], d[1], ... , d[q-1], d[q])},
//'    where \eqn{q} denotes the degree of \eqn{d(z)}.
//' @param data_in Matrix of dimension \eqn{(n \times n_obs)}{(n x n_obs)}, i.e. \eqn{(u_1, ..., u_T)}{(u[1], ..., u[T])}.
//'     Inputs to the RMFD system.
//' @param data_out Matrix of dimension \eqn{(m \times n_obs)}{(m x n_obs)}, i.e. \eqn{(y_1, ..., y_T)}{(y[1], ..., t[T])}.
//'     Outputs of the RMFD system. Initially zero and will be \strong{overwritten}.
//' @param t0 Integer. Time index from which we should start calculating a solution.
//'     Usually equal to 1.
//'
//' @return \code{data_out} is overwritten with the outputs of the RMFD system.
//'
//' @export
//'
//' @name solve_rmfd_cpp
//' @rdname solve_rmfd_cpp

// [[Rcpp::export]]
void solve_rmfd_cpp(const arma::mat& poly_inv, const arma::mat& poly_fwd, arma::mat& data_in, arma::mat& data_out, int t0) {
  // (c,d,u,y,t0)
  // We assume that dim_out >= dim_in > 0 have been checked in the R-code which calls this function, i.e. no further checks are performed here

  // int dim_out = data_out.n_rows;
  int dim_in = data_in.n_rows;
  int n_obs = data_out.n_cols;

  int deg_inv = poly_inv.n_cols / dim_in; // c(z) is parametrized as (c[p], ..., c[2], c[1]), i.e. c[0] = I_q is not included
  int deg_fwd = poly_fwd.n_cols / dim_in - 1; // d(z) must be parametrized in normal direction and includes d[0], i.e. "d = (d[0], d[1], ... , d[q])"

  // FOR DEBUGGING:
  // Rcpp::Rcout << "n: " << dim_out << " q: " << dim_in << " deg_c: " << deg_inv << " deg_d: " << deg_fwd << " n_obs: " << n_obs << std::endl;

  // http://arma.sourceforge.net/docs.html#memptr
  // vec(ptr_aux_mem, number_of_elements, copy_aux_mem = true, strict = false)
  // needed for d(z) part
  // arma::vec data_out_vec = arma::vec(data_out.memptr(), dim_out*n_obs, false, true);
  arma::vec data_in_vec = arma::vec(data_in.memptr(), dim_in*n_obs, false, true);

  // FOR DEBUGGING:
  // Rcpp::Rcout << "data_in: \n" << data_in << "\n" << std::endl;
  // Rcpp::Rcout << "data_in_vec: \n" << data_in_vec << "\n" << std::endl;


  // --------------------------------------------------------
  // c(z) Part: "Toeplitz Inverse"
  // We solve for given u[t] and c(z) for v[t] in v_t = c(z)^{-1} u_t

  // FOR DEBUGGING:
  // Rcpp::Rcout << "poly_inv: \n" << poly_inv << std::endl;
  // Rcpp::Rcout << "poly_fwd: \n" << poly_fwd << std::endl;

  if (deg_inv > 0) {

    // Indices
    int j1 = (-deg_inv+t0-1)*dim_in; // starting index for vector of length deg_inv*dim_in, i.e. (u_{t-deg_d}',...,u_{t-1}')'
    int j2 = (t0-1)*dim_in-1; // ending index for vector of length p*n
    int i1 = (t0-1)*dim_in; // starting index for vector of length q, i.e. u_t
    int i2 = t0*dim_in-1; // ending index for vector of length q, i.e. u_t

    // Starting value
    int t = t0-1;

    // FOR DEBUGGING:
    // Rcpp::Rcout << "t0-1: " << t0-1 <<  std::endl;
    // Rcpp::Rcout << "t: " << t << ", j1: "  << j1 << ", j2: "  << j2 << ", i1: "  << i1 << ", i2: "  << i2 << std::endl;

    // Loop for remaining values (starting values are dealt with in the sense that zeros are cbind()ed in R)
    while (t < n_obs) {

      // FOR DEBUGGING:
      // Rcpp::Rcout << "t: " << t << ", j1: "  << j1 << ", j2: "  << j2 << ", i1: "  << i1 << ", i2: "  << i2 << std::endl;

      // u[t] = d0inv (y[t] - (d[deg_d],..,d[1])(u[t-deg_d]',...,u[t-1]')')
      data_in_vec.subvec(i1, i2) = data_in_vec.subvec(i1, i2) - poly_inv * data_in_vec.subvec(j1, j2);

      // Increase all indices
      t++;
      j1 = j1 + dim_in;
      j2 = j2 + dim_in;
      i1 = i1 + dim_in;
      i2 = i2 + dim_in;
    }
  }
  // FOR DEBUGGING:
  // Rcpp::Rcout << "data_in after Toeplitz inverse: \n" << data_in << "\n" << std::endl;
  // Rcpp::Rcout << "data_in_vec after Toeplitz inverse: \n" << data_in_vec << "\n" << std::endl;

  // --------------------------------------------------------
  // poly_fwd(z) Part: Easy. Toeplitz multiplication for (n x q) system

  // multiply all shifted inputs with b_i (if it makes sense). When writing this with Toeplitz matrices, we see that b_1 is applied to u_1 until u_{T-1}, and b_q to u_1 till u_{T-q}
  for (int i = 0; i <= deg_fwd; i++) {

    // y[t] = y[t] + b[i] u[t-i] for all t where it makes sense
    data_out.cols(t0-1, n_obs-1) = data_out.cols(t0-1, n_obs-1) + poly_fwd.cols(i*dim_in, (i+1)*dim_in-1) * data_in.cols(t0-1-i, n_obs-1-i);
  }

  return;
 }


// -----------------------------------------------------------------------------
// 2) Compute Residuals (solve_inverse_de) ----
// -----------------------------------------------------------------------------

//' Residuals of an ARMA system
//'
//' This internal helper function computes the residuals and the directional derivatives of the
//' residuals of an ARMA system of the form
//' \deqn{a_0 y_t + a_1 y_{t-1} + \cdots + a_p y_{t-p} = b_0 u_t + \cdots + b_q u_{t-q}}{
//'       a[0] y[t] + a[1] y[t-1] + ... + a[p] y[t-p] = b[0] u[t] + ... + b[q] u[t-q]}
//'
//' Values \eqn{y_t}{y[t]}, \eqn{u_t}{u[t]} for \eqn{t\leq 0}{t\le 0} are implicitly set to be zero.
//' However, by starting the iteration with some \eqn{t_0>1}{t0>1} we can enforce non-zero
//' initial values.
//'
//' @note
//' Use this procedure with care!
//'
//' * The procedure does \bold{not} check the input arguments. We require \eqn{m = n > 0},
//'   \eqn{p,q \geq 0}{p,q \ge 0} and \eqn{1 \leq t_0 \leq N}{1 \le t0  \le N}.
//' * The procedure \bold{overwrites} the input argument \code{u} (and \code{dU}).
//' * The data matrices are organized columnwise (to avoid memory shuffling)!
//' * Note also the non standard representation of the coefficient matrices.
//'
//' @param ib0 \eqn{(m, m)} matrix, **inverse** of the coefficient matrix \eqn{b[0]}{b[0]}.
//' @param B1 \eqn{(m, mq)} matrix, \eqn{-b_0^{-1}(b_q,...,b_1)}{-b[0]^(-1)(b[q],...,b[1])}.
//' @param A \eqn{(m, n(q+1))} matrix \eqn{b_0^{-1}(a_0,...,a_p}{b[0]^(-1)(a[0],...,a[p])}.
//' @param t0 integer, start iteration at t = t0.
//' @param y \eqn{(m, N)} matrix with the observed outputs \eqn{(y_1,...,y_N}{(y[1],...,y[N])}.
//' @param u \eqn{(m, N)} matrix. This matrix is **overwritten** with the computed
//'        residuals \eqn{(u_1,...,u_N}{(u[1],...,u[N])}.
//' @param dU \eqn{(mN, m^2(p+q+2))} matrix or an empty matrix. If non empty then this
//'        matrix is **overwritten** with the directional derivatives of the vectorized residuals.
//'        The \eqn{j}-th column of \code{dU} is the derivative of \eqn{vec(u)} with respect
//'        to the \eqn{j}-th entry of
//'        \eqn{\mathrm{vec}(a_0,a_1,\ldots,a_p,b_0,\ldots,b_q)}{vec(a[0],a[1],...,a[p],b[0],...,b[q])}
//'
//' @return This RcppArmadillo routine returns \code{NULL} but \bold{overwrites}
//'         the input arguments \code{u} (and \code{dU})!
//'
//' @export
//'
//' @rdname residuals_ARMA_cpp
//' @name residuals_ARMA_cpp
//'
//' @seealso \code{\link{outputs_ARMA_cpp}}, \code{\link{residuals_ARMA_cpp}},  \code{\link{cll_theta_ARMA_cpp}},
//'    \code{\link{outputs_STSP_cpp}}, \code{\link{residuals_STSP_cpp}},  \code{\link{cll_theta_STSP_cpp}} and
//'    \code{\link{solve_de}}, \code{\link{solve_inverse_de}} and \code{\link{ll}}.
//'
//' @examples
//' # generate a random ARMA(2,1) model (3 outputs, 2 inputs)
//' p = 2
//' q = 1
//' m = 2
//' model = test_armamod(dim = c(m, m), degrees = c(p,q), digits = 2)
//'
//' # prepare parameters for "outputs_ARMA_cpp"
//' A = unclass(model$sys$a)
//' a0 = A[,,1]
//' A1 = -A[,,(p+1):2]
//' dim(A1) = c(m, m*p)
//' A1 = solve(a0, A1)
//'
//' B = unclass(model$sys$b)
//' dim(B) = c(m, m*(q+1))
//' B = solve(a0, B)
//'
//' # generate random noise sequence (sample size N = 10)
//' n.obs = 10
//' u = matrix(rnorm(n.obs*m), nrow = m, ncol = n.obs)
//'
//' # generate matrix for the outputs
//' y = matrix(0, nrow = m, ncol = n.obs)
//'
//' # compute outputs
//' t0 = 2   # start iterations from t>=t0=2
//' outputs_ARMA_cpp(A1, B, t0, u, y)
//'
//' # recompute the disturbances/residuals from the given outputs:
//' B = unclass(model$sys$b)
//' ib0 = B[,,1]
//' B1 = -B[,,(q+1):2]
//' dim(B1) = c(m, m*q)
//' B1 = solve(ib0, B1)
//'
//' A = unclass(model$sys$a)
//' dim(A) = c(m, m*(p+1))
//' A = solve(ib0, A)
//'
//' ib0 = solve(ib0)
//'
//' uu = u + 0 # "deep copy" of the disturbances
//' uu[, t0:(n.obs)] = 0 # clear values for t >= t0
//' residuals_ARMA_cpp(ib0, B1, A, t0 = 2, y, uu, diag(0))
//' all.equal(u, uu) # check
//'
//' # compute directional derivatives of residuals
//' dU = matrix(0, nrow = n.obs*m, ncol = (m^2)*(p+q+2))
//' residuals_ARMA_cpp(ib0, B1, A, t0 = 2, y, uu, dU)

// [[Rcpp::export]]
void residuals_ARMA_cpp(const arma::mat& ib0, const arma::mat& B1, const arma::mat& A,
                        int t0, const arma::mat& y, arma::mat& u, arma::mat& dU) {

  // check 1 <= t0 <= N
  int N = y.n_cols;
  if ((N == 0) || (t0 > N)) {
    // nothing to do
    return;
  }
  if (t0 < 1) {
    stop("t0 < 1!");
  }

  // check m = n > 0
  int m = y.n_rows;
  int n = u.n_rows;
  if ((m !=n) || (m < 1)) {
    stop("only (non-empty) square systems (n = m > 0) are supported!");
  }
  int q = B1.n_cols / m;
  int p = (A.n_cols / m) - 1;
  if (p < 0) {
    stop("p < 0!");
  }
  int Nm = N*m;

  // compute residuals, by calling outputs_ARMA_cpp (and switching the roles of a(z) and b(z))
  // Rcout << "residuals_ARMA_cpp m:" << m << " n:" << n << " q:" << q << " p:" << p << " Nm:" << Nm << " N:" << N << std::endl;
  outputs_ARMA_cpp(B1, A, t0, y, u);

  // compute directional derivatives of residuals only if dU is a non-empty matrix
  if (dU.n_elem > 0) {
    // int K = m*m*(p+q+2);  // note that dU must be an (Nm, K) dimensional matrix!

    // create "working" matrices c, v, du and a vectorized version of du
    // Rcout << "residuals_ARMA_cpp Nm" << Nm << " m:" << m << " N:" << N << std::endl;

    arma::mat c = arma::mat(m, 1, arma::fill::zeros);
    arma::mat v = arma::mat(1, N, arma::fill::zeros);
    arma::mat du = arma::mat(m, N, arma::fill::zeros);
    // Rcout << du << std::endl;
    arma::vec du_vec = arma::vec(du.memptr(), Nm, false, true);
    // Rcout << du_vec << std::endl;

    // compute derivative of u[t] with respect to a[k][i,j] ******************************
    int ind = 0;  // column index of the Jacobi matrix dU
    for (int j = 0; j < m; j++ ) {
         for (int i = 0; i < m; i++ ) {
           // for k=0: b[0]*du[t] + b[1]*du[t-1] + ... + b[q]*du[-qt] = (ei * ej') * y[t],
           // where ei, ej are the canonical unit vectors in R^m.
           c = ib0.cols(i,i);
           v = y.rows(j,j);
           outputs_ARMA_cpp(B1, c, t0, v, du);

           for (int k = 0; k <= p; k++) {
             // FOR DEBUGGING
             // Rcout << "J(" << J.n_rows << "x" << J.n_cols << ") ind: " << ind+(m*m*k) << " mk: "
             // << m*k << " Nm-mk: " << Nm - m*k << std::endl;

             if (k > 0) {
               // set the first k elements equal to 0
               dU.col(ind+(m*m*k)).head(m*k).zeros();
             }
             dU.col(ind+(m*m*k)).tail(Nm-m*k) = du_vec.head(Nm-m*k);
           }
           ind++;
         }
     }

    // compute derivative of u[t] with respect to b[k][i,j] ******************************
    ind = m*m*(p+1);
    for (int j = 0; j < m; j++ ) {
      for (int i = 0; i < m; i++ ) {
        // for k=0: b[0]*du[t] + b[1]*du[t-1] + ... + b[q]*du[-qt] = -(ei * ej') * u[t],
        // where ei, ej are the canonical unit vectors in R^m.
        c = -ib0.cols(i,i);
        v = u.rows(j,j);
        outputs_ARMA_cpp(B1, c, t0, v, du);

        for (int k = 0; k <= q; k++) {
          // FOR DEBUGGING
          // Rcout << "J(" << J.n_rows << "x" << J.n_cols << ") ind: " << ind+(m*m*k) << " mk: "
          // << m*k << " Nm-mk: " << Nm - m*k << std::endl;

          if (k > 0) {
            // set the first k elements equal to 0
            dU.col(ind+(m*m*k)).head(m*k).zeros();
          }
          dU.col(ind+(m*m*k)).tail(Nm-m*k) = du_vec.head(Nm-m*k);
        }
        ind++;
      }
    }
  }

  return;
}

//' Residuals of a statespace system
//'
//' This internal helper function computes the residuals
//' (and the directional derivatives of the residuals)
//' for a statespace system of the form
//' \deqn{a_{t+1} = A a_t + B u_t, \; y_t = C a_t + D u_t}{
//'       a[t+1] = A a[t] + B u[t] and y[t] = C a[t] + D u[t]}
//' The system must be square and non-empty, i.e. \eqn{m=n>0}.
//'
//'
//' @note
//' Use this procedure with care!
//'
//' * The procedure does \bold{not} check the input arguments.
//' * The procedure \bold{overwrites} the input arguments
//'    \code{a}, \code{u} and \code{dU}.
//' * The data matrices are organized columnwise (to avoid memory shuffling)!
//'
//' @param A \eqn{(s,s)} matrix.
//' @param B \eqn{(s,m)} matrix.
//' @param C \eqn{(m,s)} matrix.
//' @param D \eqn{(m,m)} matrix, must be regular.
//' @param y \eqn{(m,N)} matrix with the outputs: \eqn{(y_1,y_2,\ldots,y_N)}{(y[1],y[2],...,y[N])}.
//' @param a \eqn{(s,N+1)} matrix. This matrix is overwritten with the (computed) states:
//'            \eqn{(a_1,a_2,\ldots,a_N,a_{N+1})}{(a[1],a[2],\ldots,a[N],a[N+1])}.
//'            On input \code{a[,1]} must hold the initial state \eqn{a_1}{a[1]}.
//' @param u \eqn{(m,N)} matrix. This matrix is overwritten with (computed) residuals:
//'           \eqn{(u_1,u_2,\ldots,u_N)}{(u[1],u[2],...,u[N])}.
//' @param dPI \eqn{((m+s)^2,K)} matrix.
//' @param dU \eqn{(mN,K)} matrix or \eqn{(0,0)} matrix. This matrix is overwritten with the
//'           directional derivatives of the residuals. However, if
//'           the matrix is empty then no derivatives are computed.
//'
//' @return This RcppArmadillo implementation returns \code{NULL} but \bold{overwrites} the input
//'         arguments \code{a}, \code{u} and \code{dU}!
//'
//' @export
//'
//' @rdname residuals_STSP_cpp
//' @name residuals_STSP_cpp
//'
//' @seealso \code{\link{outputs_ARMA_cpp}}, \code{\link{residuals_ARMA_cpp}},  \code{\link{cll_theta_ARMA_cpp}},
//'    \code{\link{outputs_STSP_cpp}}, \code{\link{residuals_STSP_cpp}},  \code{\link{cll_theta_STSP_cpp}} and
//'    \code{\link{solve_de}}, \code{\link{solve_inverse_de}} and \code{\link{ll}}.
//'
//' @examples
//' # generate a random statespace model (3 outputs, 3 inputs and 4 states)
//' m = 2
//' s = 3
//' model = test_stspmod(dim = c(m, m), s = s, digits = 2)
//'
//' # generate random noise sequence (sample size N = 10)
//' n.obs = 10
//' u = matrix(rnorm(n.obs*m), nrow = m, ncol = n.obs)
//'
//' # generate matrix for the state sequence
//' a = matrix(0, nrow = s, ncol = n.obs+1)
//' a[,1] = rnorm(s) # random initial state a[1]
//'
//' # generate matrix for the outputs
//' y = matrix(0, nrow = m, ncol = n.obs)
//'
//' # compute outputs and states
//' outputs_STSP_cpp(model$sys$A, model$sys$B, model$sys$C, model$sys$D, u, a, y)
//'
//' # recompute the states and disturbances/residuals from the given outputs:
//' uu = u + 0 # "deep copy" of the disturbances
//' aa = a + 0 # and the states
//' aa[, 2:(n.obs+1)] = 0 # clear all states a[t], t > 1
//' residuals_STSP_cpp(model$sys$A, model$sys$B, model$sys$C, model$sys$D,
//'                    y, aa, uu, diag(0), diag(0))
//' all.equal(u, uu) # check
//' all.equal(a, aa) # check
//'
//' # compute directional derivatives of residuals
//' dPI = diag((m+s)^2)
//' dU = matrix(0, nrow = n.obs*m, ncol = ncol(dPI))
//' residuals_STSP_cpp(model$sys$A, model$sys$B, model$sys$C, model$sys$D,
//'                    y, aa, uu, dPI, dU)
//'
//' # check the directional derivatives
//' eps = 1e-8
//' dU_num = matrix(0, nrow = m*n.obs, ncol = (m+s)^2)
//' dPI = matrix(0, nrow = (m+s), ncol = (m+s))
//' for (k in (1:((m+s)^2))) {
//'   dPI[] = 0
//'   dPI[k] = eps
//'   uu = u + 0 # "deep copy" of the disturbances
//'   aa = a + 0 # and the states
//'   residuals_STSP_cpp(model$sys$A + dPI[1:s,1:s],
//'                      model$sys$B + dPI[1:s,(s+1):(s+m)],
//'                      model$sys$C + dPI[(s+1):(s+m),1:s],
//'                      model$sys$D + dPI[(s+1):(s+m),(s+1):(s+m)],
//'                      y, aa, uu, diag(0), diag(0))
//'   dU_num[, k] = c(uu - u )/eps  # num. approx. of the derivative in direction "dPI"
//' }
//' # relative error of the numerical approximation
//' junk = (abs(dU)+abs(dU_num))
//' junk[junk == 0] = 1
//' 2*abs(dU_num - dU)/junk
//'
// [[Rcpp::export]]
void residuals_STSP_cpp(const arma::mat& A, const arma::mat& B,
                         const arma::mat& C, const arma::mat& D,
                         const arma::mat& y,
                         arma::mat& a, arma::mat& u,
                         const arma::mat& dPI, arma::mat& dU) {

   unsigned long int nobs = y.n_cols; // sample size
   unsigned long int s = A.n_rows;    // state dimension
   unsigned long int m = D.n_rows;    // input dimension
   unsigned long int n = D.n_cols;    // output dimension
   if ((n != m) | (m <= 0)) stop("only square, non-empty systems (n = m > 0) are supported.");

   // solve the state space system for u[t]!
   // a[t+1] = A a[t] + B u[t]                       and y[t] = C a[t] + D u[t]
   // u[t] = -D^(-1) C a[t] + D^(-1) y[t]            and a[t+1] = A a[t] + B u[t]

   arma::mat iD = inv(D);                  // iD = D^(-1)
   arma::mat iC = -solve(D, C);            // iC = -D^(-1) C

   // compute residuals/disturbances u[t]
   u = iD * y;
   if (s > 0) {
     for (unsigned long int t = 0; t < nobs; t++) {
       u.col(t)   = iC * a.col(t) + u.col(t);
       a.col(t+1) = A * a.col(t) + B * u.col(t);
     }
   }

   // compute directional derivatives of residuals only if dU is a non-empty matrix
   if (dU.n_elem > 0) {
     unsigned long int K = dU.n_cols;

     arma::mat dA(s,s);
     arma::mat dB(s,m);
     arma::mat dC(m,s);
     arma::mat dD(m,m);

     // pi corresponds to the stacked system matrices [A B]
     //                                               [C D]
     arma::mat dpi(s+m, s+m);
     // http://arma.sourceforge.net/docs.html#memptr
     // vec(ptr_aux_mem, number_of_elements, copy_aux_mem = true, strict = false)
     // dpi_vec is the vectorised version of the matrix sys (shares the memory with dpi)
     arma::vec dpi_vec = arma::vec(dpi.memptr(), (m+s)*(m+s), false, true);

     arma::mat du(m, nobs);
     arma::vec du_vec = arma::vec(du.memptr(), m*nobs, false, true);
     // du_vec is the vectorised version of the matrix du (shares the memory with du)

     // compute the direction derivative of the residuals in direction dpi (dA, dB, dC, dD)
     if (s == 0) {
       for (unsigned long int i = 0; i < K; i++) {
         dpi_vec = dPI.col(i);

         // RcppArmadillo syntax for submatrices:
         //   X.submat( first_row, first_col, last_row, last_col )
         dD = dpi.submat(s, s, s+m-1, s+n-1);
         dD = - iD * dD * iD;   // dD is overwritten with d(iD)
         du = dD * y;
         dU.col(i) = du_vec; // store into i-th column
       } // end for (i = ...)
     } else {
       arma::mat da(s, nobs + 1);

       for (unsigned long int i = 0; i < K; i++) {
         dpi_vec = dPI.col(i);

         // RcppArmadillo syntax for submatrices:
         //   X.submat( first_row, first_col, last_row, last_col )
         dA = dpi.submat(0,0,s-1,s-1);
         dB = dpi.submat(0,s,s-1,s+m-1);
         dC = dpi.submat(s,0,s+m-1,s-1);
         dD = dpi.submat(s,s,s+m-1,s+m-1);

         dC = -iD * dD *iC - iD * dC;                // dC is overwritten with d(iC)
         dD = - iD * dD * iD;                        // dD is overwritten with d(iD)

         // Rcout << i << dA << dB << dC << dD << std::endl;

         // directional derivatives of residuals in direction dpi = (dA, dB, dC, dD):
         // u[t]   = iC a[t] + iD y[t]  and
         // a[t+1] =  A a[t] +  B u[t]
         // implies
         // du[t]   = iC da[t] + d(iC) a[t] + d(iD) y[t]  and
         // da[t+1] =  A da[t] + B du[t] + dA a[t] + dB u[t]
         // where
         // d(iC) = d(-D^(-1) C)     = D^(-1) dD D^(-1) C - D^(-1) dC
         // d(iD) = d(D^(-1))        = -D^(-1) dD D^(-1)
         du = dC * a.cols(0, (nobs-1)) + dD * y;
         da.col(0).zeros(); // set da[1] = 1
         da.tail_cols(nobs) = dA * a.head_cols(nobs) + dB * u;

         for (unsigned long int t = 0; t < nobs; t++) {
           du.col(t)   = iC * da.col(t) + du.col(t);
           da.col(t+1) = A * da.col(t) + B * du.col(t) + da.col(t+1);
         }

         dU.col(i) = du_vec;  // store into i-th column

       } // end for(i = ...)
     } // end if (s > 0)

   } // end if (dU.n_elem > 0)

   return;
 }

// -----------------------------------------------------------------------------
// 3) Conditonal Log Likelihood ----
// -----------------------------------------------------------------------------

//' Compute the (concentrated) conditional log likelihood for ARMA models
//' described by a model template.
//'
//' This internal helper function computes the (concentrated) conditional log Likelihood
//' of ARMA systems of the form
//' \deqn{a_0 y_t + a_1 y_{t-1} + \cdots + a_p y_{t-p} = b_0 u_t + \cdots + b_q u_{t-q}}{
//'       a[0] y[t] + a[1] y[t-1] + ... + a[p] y[t-p] = b[0] u[t] + ... + b[q] u[t-q]}
//' The conditional likelihood is computed for **zero** initial values \eqn{u_s=y_s=0}{u[s]=y[s]=0} for
//' \eqn{s\leq 0}{s\le 0}.
//'
//' This function is mainly used by the function factory \code{\link{ll_FUN}}. For a more detailed
//' documentation of the (concentrated) conditional log Likelihood, see \code{\link{ll}}.
//'
//' The procedure first constructs the ARMA parameter matrices from the given vector
//' \code{th} of "deep" parameters.
//'
//' \itemize{
//' \item  AR parameters \code{vec((a[0],a[1],...,a[p])) = h_A + H_A * th}.
//' \item  MA parameters \code{vec(b[0]) = h_b + H_b * th} and
//'        \code{vec(-(b[q],...,b[1])) = h_B + H_B * th}
//' \item  Left square root of noise covariance matrix
//'        \eqn{\Sigma = LL'} and \code{vec(L) = h_L + H_L * th}.
//' }
//'
//' The residuals (and their directional derivatives) are computed with
//' \code{\link{residuals_ARMA_cpp}}.
//'
//' @note
//' Use this procedure with care!
//'
//' * The procedure does \bold{not} check the input arguments.
//' * The procedure \bold{overwrites} some of the input arguments
//' * The data matrices are organized columnwise (to avoid memory shuffling)!
//' * Note also the non standard representation of the coefficient matrices.
//'
//' @param th \eqn{(K)} dimensional vector of "deep" parameters.
//' @param y \eqn{(m,N)} matrix with the observed outputs:
//'        \eqn{(y_1,y_2,\ldots,y_N)}{(y[1],y[2],...,y[N])}.
//' @param skip (integer), omit the first "skip" residuals, when computing the likelihood.
//' @param concentrated (bool), if TRUE then the *concentrated*, conditional log Likelihood is computed
//' @param ib0 \eqn{(m, m)} matrix, is **overwritten** with the matrix \eqn{b_0^{-1}a_0}{b[0]^(-1)a[0]}.
//' @param H_b \eqn{(m^2, K)} matrix.
//' @param h_b \eqn{((m^2)}-dimensional vector. Note that
//'         \code{vec(b[0]) = H_b*th + h_b}.
//' @param B1 \eqn{(m, mq)} matrix, is **overwritten** with
//'        \eqn{-b_0^{-1}(b_q,...,b_1)}{-b[0]^(-1)(b[q],...,b[1])}.
//' @param H_B \eqn{((m^2)*q, K)} matrix.
//' @param h_B \eqn{((m^2)*q)}-dimensional vector. Note that
//'        \code{vec(-(b[q],...,b[1])) = H_B*th + h_B}.
//' @param a0 \eqn{(m, m)} matrix, is **overwritten** with
//'        \eqn{a_0}{a[0]}.
//' @param A \eqn{(m, m(q+1))} matrix, is **overwritten** with
//'        \eqn{b_0^{-1}(a_0,...,a_p}{b[0]^(-1)(a[0],...,a[p])}.
//' @param H_A \eqn{((m^2)*(p+1), K)} matrix.
//' @param h_A \eqn{((m^2)*(p+1))}-dimensional vector. Note that
//'        \code{vec((a[0],a[1],...,a[p])) = H_A*th + h_A}.
//' @param L \eqn{(m,m)} matrix. If (concentrated==FALSE) then \code{L} is **overwritten** with
//'        the left square \eqn{L} of the noise covariance matrix \eqn{\Sigma=LL'} corresponding
//'        to the deep parameters th. However, if (concentrated==TRUE) then
//'        L is **overwritten** with sample covariance matrix of the computed residuals!
//' @param H_L \eqn{(m^2, K)} matrix.
//' @param h_L \eqn{(m^2)}-dimensional vector. Note that
//'        \code{vec(L) = H_L*th + h_L}.
//' @param u \eqn{(m,N)} matrix. This matrix is **overwritten** with (computed) residuals:
//'           \eqn{(u_1,u_2,\ldots,u_N)}{(u[1],u[2],...,u[N])}.
//' @param dU \eqn{(mN,(m^2)(p+q+2))} matrix or \eqn{(0,0)} matrix. If non empty this
//'           matrix is **overwritten** with the
//'           directional derivatives of the residuals. However, if
//'           the matrix is empty then no derivatives are computed.
//'
//' @seealso \code{\link{outputs_ARMA_cpp}}, \code{\link{residuals_ARMA_cpp}},  \code{\link{cll_theta_ARMA_cpp}},
//'    \code{\link{outputs_STSP_cpp}}, \code{\link{residuals_STSP_cpp}},  \code{\link{cll_theta_STSP_cpp}} and
//'    \code{\link{solve_de}}, \code{\link{solve_inverse_de}} and \code{\link{ll}}.
//'
//' @return (double) log Likelihood
//'
//' @export
//'
//' @rdname cll_theta_ARMA_cpp
//' @name cll_theta_ARMA_cpp
// [[Rcpp::export]]
double cll_theta_ARMA_cpp(const arma::vec& th, const arma::mat& y, unsigned long int skip, bool concentrated,
                          arma::mat& ib0, const arma::mat& H_b, const arma::vec& h_b,
                          arma::mat& B1, const arma::mat& H_B, const arma::vec& h_B,
                          arma::mat& a0, arma::mat& A, const arma::mat& H_A, const arma::vec& h_A,
                          arma::mat& L, const arma::mat& H_L, const arma::vec& h_L,
                          arma::mat& u, arma::mat& dU) {

  // Number of observations. Note that the data matrix y is an (m, N) matrix, i.e.
  // each observation y_t corresponds to one column of the data matrix y. (to optimize memory access)
  int m = y.n_rows;
  int nobs = y.n_cols;
  int nvalid = nobs - skip;
  int p = (A.n_cols / m) - 1;
  int q = B1.n_cols / m;
  // Rcout << "m: " << m << " nobs: " << nobs << " nvalid: " << nvalid << std::endl;

  // Create vectorized versions of the (system) matrices (ib0, B1, A),
  // which share the memory with the matrices. Hence we can easily
  // access these matrices either as "matrix" or as "vectors".
  //
  // http://arma.sourceforge.net/docs.html#memptr
  // http://arma.sourceforge.net/docs.html#Col (see "Advanced constructors")
  //   vec(ptr_aux_mem, number_of_elements, copy_aux_mem = true, strict = false)
  arma::vec ib0_vec = arma::vec(ib0.memptr(), m*m, false, true);
  // Rcout << "ib0" << std::endl << ib0  << ib0_vec << std::endl;
  arma::vec B1_vec = arma::vec(B1.memptr(), m*m*q, false, true);
  // Rcout << "B1" << std::endl << B1 << B1_vec << std::endl;
  arma::vec A_vec = arma::vec(A.memptr(), m*m*(p+1), false, true);
  // Rcout << "A" << std::endl << A << A_vec << std::endl;

  // construct system matrices for given parameters th
  ib0_vec = h_b + H_b * th;      // = b[0]
  if (q > 0) {                   // note the negative sign and the reverse order!
    B1_vec  = h_B  + H_B * th;   // = -(b[q],b[q-1],...,b[1])
    B1 = solve(ib0, B1);         // = -(b[0])^(-1) * (b[q],b[q-1],...,b[1])
  }
  A_vec  = h_A  + H_A * th;       // = (a[0],a[1],...,a[p])
  a0 = A.head_cols(m);            // = a[0]
  A = solve(ib0, A);              // = (b[0])^(-1) * (a[0],a[1],...,a[p])

  ib0 = inv(ib0);                 // = (b[0])^(-1)

  // compute residuals
  // Rcout << ib0 << std::endl;
  // Rcout << B1 << std::endl;
  // Rcout << A << std::endl;
  // Rcout << y << std::endl;
  // Rcout << u << std::endl;
  // Rcout << "dU" << std::endl << dU << std::endl;

  if (dU.n_elem == 0) {
    // we  only need the residuals
    outputs_ARMA_cpp(B1, A, 1L, y, u);
  } else {
    // we also need the derivatives of the residuals
    residuals_ARMA_cpp(ib0, B1, A, 1L, y, u, dU);
    // Rcout << "dU" << std::endl << dU << std::endl;
    // Rcout << "finished residuals_ARMA_cpp" << std::endl;
  }

  double ll, sign, lndetk0, lndetSigma, trSS, mln2pi;

  if (concentrated) {
    // concentrated, conditional log Likelihood

    // sample covariance matrix of residuals
    L = ( u.head_cols(nvalid) * u.head_cols(nvalid).t() ) / nvalid;
    log_det(lndetSigma, sign, L);
    trSS = m;
  } else {
    // conditional log Likelihood

    // L_vec is the vectorised version of the matrix L (shares the memory with L)
    // Rcout << "m: " << m << "m*m: " << m*m << " L" << std::endl;
    // Rcout << L << std::endl;
    arma::vec L_vec = arma::vec(L.memptr(), m*m, false, true);
    // Rcout << L_vec << std::endl;
    L_vec = h_L + H_L * th;
    // Rcout << L << std::endl;

    log_det(lndetSigma, sign, L);
    lndetSigma = 2*lndetSigma; // log(det( sigma ));
    // compute sum_t u[t]' sigma^(-1) u[t]
    u = solve(L, u);
    trSS = accu(square(u.head_cols(nvalid))) / nvalid;
  }

  mln2pi = m*log(2 * arma::datum::pi);

  // take care of the case where the lag zero coefficient k[0] = a[0]^(-1) b[0]
  // of the impulse response is not equal to the identity!
  // lndetk0 = log(abs(det(k0))) = -log(abs(det(b[0]^-1 a[0])))
  log_det(lndetk0, sign, A.head_cols(m)); // b[0]^(-1) * a[0] = k[0]^(-1)
  lndetk0 = -lndetk0;

  ll = (-0.5) * (mln2pi + trSS + lndetSigma + 2*lndetk0);
  return(ll);
}


//' Compute the (concentrated) conditional log likelihood for a statespace system
//' described by a model template.
//'
//' This is an internal helper function, used by the function factory \code{\link{ll_FUN}}. For a more detailed
//' documentation of the conditional log Likelihood, see \code{\link{ll}}.
//' The conditional likelihood is computed for the initial state \eqn{a_1} given in the first column `a[,1]` of the
//' matrix `a`.
//'
//' @param th \eqn{(K)} dimensional vector of "deep" parameters.
//' @param y \eqn{(m,N)} matrix with the observed outputs:
//'        \eqn{(y_1,y_2,\ldots,y_N)}{(y[1],y[2],...,y[N])}.
//' @param skip (integer), skip the first residuals, when computing the sample covariance of the
//'        residuals.
//' @param concentrated (bool), if TRUE then the *concentrated*, conditional log Likelihood is computed
//' @param pi \eqn{(m+s,m+s)} matrix, is overwritten with the system matrix
//'        \eqn{[A,B | C,D]}.
//' @param H_pi \eqn{(m+s)^2, K)} matrix.
//' @param h_pi \eqn{((m+s)^2)}-dimensional vector. Note that \code{vec(pi) = H_pi*th + h_pi}.
//' @param L \eqn{(m,m)} matrix. If (concentrated==FALSE) then L is overwritten with
//'        the left square of the noise covariance matrix L corresponding
//'        to the deep parameters th. However, if (concentrated==TRUE) then
//'        L is overwritten with sample covariance matrix of the computed residuals!
//' @param H_L \eqn{(m^2, K)} matrix.
//' @param h_L \eqn{(m^2)}-dimensional vector. Note that
//'        \code{vec(L) = H_L*th + h_L}.
//' @param a \eqn{(s,N+1)} matrix. This matrix is overwritten with the (computed) states:
//'            \eqn{(a_1,a_2,\ldots,a_N,a_{N+1})}{(a[1],a[2],\ldots,a[N],a[N+1])}.
//'            On input \code{a[,1]} must hold the initial state \eqn{a_1}{a[1]}.
//' @param u \eqn{(m,N)} matrix. This matrix is overwritten with (computed) residuals:
//'           \eqn{(u_1,u_2,\ldots,u_N)}{(u[1],u[2],...,u[N])}.
//' @param dU \eqn{(mN,K)} matrix or \eqn{(0,0)} matrix. This matrix is overwritten with the
//'           directional derivatives of the residuals. However, if
//'           the matrix is empty then no derivatives are computed.
//'
//' @return (double) log Likelihood
//'
//' @export
//'
//' @seealso \code{\link{outputs_ARMA_cpp}}, \code{\link{residuals_ARMA_cpp}},  \code{\link{cll_theta_ARMA_cpp}},
//'    \code{\link{outputs_STSP_cpp}}, \code{\link{residuals_STSP_cpp}},  \code{\link{cll_theta_STSP_cpp}} and
//'    \code{\link{solve_de}}, \code{\link{solve_inverse_de}} and \code{\link{ll}}.
//'
//' @rdname cll_theta_STSP_cpp
//' @name cll_theta_STSP_cpp

// [[Rcpp::export]]
double cll_theta_STSP_cpp(const arma::vec& th, const arma::mat& y, unsigned long int skip, bool concentrated,
                          arma::mat& pi, const arma::mat& H_pi, const arma::vec& h_pi,
                          arma::mat& L, const arma::mat& H_L, const arma::vec& h_L,
                          arma::mat& a, arma::mat& u, arma::mat& dU) {

  unsigned long int s = a.n_rows;
  unsigned long int m = y.n_rows;
  unsigned long int nobs = y.n_cols;
  unsigned long int nvalid = nobs - skip;

  // http://arma.sourceforge.net/docs.html#memptr
  // vec(ptr_aux_mem, number_of_elements, copy_aux_mem = true, strict = false)
  // pi_vec is the vectorised version of the matrix pi (shares the memory with pi)
  arma::vec pi_vec = arma::vec(pi.memptr(), (m+s)*(m+s), false, true);

  pi_vec = h_pi + H_pi * th;

  // A = pi.submat(0,0,s-1,s-1);
  // B = pi.submat(0,s,s-1,s+m-1);
  // C = pi.submat(s,0,s+m-1,s-1);
  // D = pi.submat(s,s,s+m-1,s+m-1);


  // compute residuals
  // void residuals_STSP_cpp(const arma::mat& A, const arma::mat& B,
  //                         const arma::mat& C, const arma::mat& D,
  //                         const arma::mat& y,
  //                         arma::mat& a, arma::mat& u,
  //                         const arma::mat& dPI, arma::mat& dU)

  if (s > 0) {
    residuals_STSP_cpp(pi.submat(0,0,s-1,s-1),      // A
                       pi.submat(0,s,s-1,s+m-1),    // B
                       pi.submat(s,0,s+m-1,s-1),    // C
                       pi.submat(s,s,s+m-1,s+m-1),  // D
                       y, a, u, H_pi, dU);
  } else {
    residuals_STSP_cpp(arma::mat(0,0),   // A
                       arma::mat(0,m),   // B
                       arma::mat(m,0),   // C
                       pi,               // D
                       y, a, u, H_pi, dU);
  }

  double ll, sign, lndetk0, lndetSigma, trSS, mln2pi;

  if (concentrated) {
    // concentrated, conditional log Likelihood

    // sample covariance matrix of residuals
    L = ( u.head_cols(nvalid) * u.head_cols(nvalid).t() ) / nvalid;
    log_det(lndetSigma, sign, L);
    trSS = m;
  } else {
    // conditional log Likelihood

    // L_vec is the vectorised version of the matrix L (shares the memory with L)
    arma::vec L_vec = arma::vec(L.memptr(), (m)*(m), false, true);
    L_vec = h_L + H_L * th;

    log_det(lndetSigma, sign, L);
    lndetSigma = 2*lndetSigma; // log(det( sigma ));
    // compute sum_t u[t]' sigma^(-1) u[t]
    u = solve(L, u);
    trSS = accu(square(u.head_cols(nvalid))) / nvalid;
  }

  mln2pi = m*log(2 * arma::datum::pi);
  // take care of the case where the lag zero coefficient k[0] = D
  // of the impulse response is not equal to the identity!
  // lndetk0 = log(abs(det(k0)))
  log_det(lndetk0, sign, pi.submat(s,s,s+m-1,s+m-1));

  ll = (-0.5) * (mln2pi + trSS + lndetSigma + 2*lndetk0);
  return(ll);
}



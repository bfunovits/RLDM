# 1)  Solve (linear) Difference Equations ----

## 1.1) Calculate Outputs From Inputs ----

#' Solve (linear) Difference Equations
#'
#' The procedure `solve_de()` solves the difference equations associated to (V)ARMA models
#' \deqn{a_0 y_t + a_1 y_{t-1} + \cdots + a_p y_{t-p} = b_0 u_t  + b_1 u_{t-1} + ... b_1 u_{t-q}}{
#'       a[0] y[t] + a[1] y[t-1] + ... + a[p] y[t-p] = b[0] u[t]  + b[1] u[t-1] + ... b[q] u[t-q]}
#' or state space models
#' \deqn{a_{t+1} = A a_t + B u_t \mbox{ and } y_t = C a_t + D u_t.}{
#'       a[t+1] = A a[t] + B u[t] and y[t] = C a[t] + D u[t].}
#'
#' `solve_de()` computes the outputs \eqn{y_t}{y[t]} for \eqn{t=1,\ldots,N}{t=1,...,N} for
#' given disturbances \eqn{u_t}{u[t]} \eqn{t=1,\ldots,N}{t=1,...,N}.
#' The starting values  (\eqn{u_t}{u[t]} and \eqn{y_t}{y[t]} for \eqn{t\leq 0}{t \le 0} for VARMA models
#' and \eqn{a_1}{a[1]} for state space models) may be given as optional arguments.
#' The default is to use zero vectors.
#'
#' For the reverse direction, i.e. to reconstruct the disturbances if the outputs are given,
#' the function `solve_inverse_de` may be used. In this case the system must be square and
#' the matrix \eqn{D} respectively \eqn{b_0}{b[0]} must be invertible.
#'
#' These functions are mainly intended for internal use and hence only some basic checks
#' on the input parameters are performed.
#'
#' @param sys [rationalmatrices::lmfd()] or [rationalmatrices::stsp()] object
#'            which describes the difference equation.
#' @param u \eqn{(N,n)} matrix with the noise (\eqn{u_t}{u[t]}, \eqn{t=1,...,N}).
#' @param y \eqn{(N,m)} matrix with the outputs (\eqn{y_t}{y[t]}, \eqn{t=1,...,N}).
#' @param u0 \eqn{(h,n)} dimensional matrix with starting values
#'           for the disturbances \eqn{(u_{1-h}, \ldots, u_{-1}, u_0)}{(u[1-h], ..., u[-1], u[0])}.
#'           Note that the last row corresponds to \eqn{u_0}{u[0]}, the last but one row
#'           to \eqn{u_{-1}}{u[-1]} and so on. If \eqn{h>q} then only the last \eqn{q} rows of
#'           `u0` are used. In the case \eqn{h<q} the "missing" initial values are set
#'           to zero vectors. \cr
#'           The default value `u0=NULL` sets all initial values \eqn{u_t}{u[t]}, \eqn{t \leq 0}{t \le 0}
#'           equal to zero vectors.
#' @param y0 \eqn{(h,m)} dimensional matrix with starting values
#'           for the outputs \eqn{(y_{1-h}, \ldots, y_{-1}, y_0)}{(y[1-h], ..., y[-1], y[0])}.
#'           This (optional) parameter is interpreted analogously to `u0`.
#' @param a1 \eqn{m} dimensional vector with the initial state \eqn{a_1}{a[1]}.
#'           If `a1=NULL` then a zero vector is used.
#' @param ... not used.
#'
#'
#' @return List with slots
#' \item{y}{\eqn{(N,n)} matrix with the (computed) outputs.}
#' \item{u}{\eqn{(N,n)} matrix with the (computed) noise.}
#' \item{a}{\eqn{(N+1,n)} matrix with the (computed) states (\eqn{a_t}{a[t]}, \eqn{t=1,...,N+1}).
#'          Note that this matrix has (\eqn{N+1}) rows! This slot is only present for state space models.}
#' @export
#'
#' @examples
#'
#' ### generate a random ARMA(2,1) model (with two outputs) #########
#' model = test_armamod(dim = c(2,2), degrees = c(2,1),
#'                      digits = 2, bpoles = 1, bzeroes = 1)
#'
#' # generate random noise sequence (sample size N = 100)
#' u = matrix(rnorm(100*2), nrow = 100, ncol = 2)
#'
#' # generate random initial values
#' u0 = matrix(rnorm(2), nrow = 1, ncol = 2) # u[0]
#' y0 = matrix(rnorm(2), nrow = 1, ncol = 2) # y[0]
#'
#' # compute outputs "y[t]"
#' # note that y0 has only one row, thus y[-1] is set to zero!
#' data = solve_de(model$sys, u = u, y0 = y0, u0 = u0)
#'
#' # we can reconstruct the noise "u" from given outputs "y"
#' data1 = solve_inverse_de(model$sys, y = data$y, u0 = u0, y0 = y0)
#' all.equal(data$u, data1$u)
#'
#' ### generate a random state space model (3 outputs and 4 states) ##
#' model = test_stspmod(dim = c(3,3), s = 4,
#'                      digits = 2, bpoles = 1, bzeroes = 1)
#'
#' # generate random noise sequence (sample size N = 100)
#' u = matrix(rnorm(100*3), nrow = 100, ncol = 3)
#'
#' # generate random initial state a[1]
#' a1 = rnorm(4)
#'
#' # compute outputs "y[t]"
#' data = solve_de(model$sys, u = u, a1 = a1)
#'
#' # we can reconstruct the noise "u" from given outputs "y"
#' data1 = solve_inverse_de(model$sys, y = data$y, a1 = data$a[1,])
#' all.equal(data$u, data1$u)
#'
solve_de = function(sys, u, ...)  {
  UseMethod("solve_de", sys)
}


#' @rdname solve_de
#' @export
solve_de.stsp = function(sys, u, a1 = NULL, ...) {
  A = sys$A
  B = sys$B
  C = sys$C
  D = sys$D

  m = nrow(D)
  n = ncol(D)
  s = nrow(A)

  # check 'u'
  if ( (!is.matrix(u)) || (ncol(u) != n) ) stop('disturbances "u" are not compatible to the model!')
  n.obs = nrow(u)

  if (is.null(a1)) {
    a1 = double(s)
  }
  if ( (!is.vector(a1)) || ( length(a1) != s ) ) stop('initial state "a1" is not compatible!')

  if (n.obs == 0) {
    return(list(y = matrix(0, nrow = n.obs, ncol = m), u = u, a = matrix(a1, nrow = 1, ncol = s)))
  }

  # note: outputs_STPS_cpp() uses transposed data matrices!!!
  u = t(u)
  y = matrix(0, nrow = m, ncol = n.obs)
  a = matrix(0, nrow = s, ncol = n.obs + 1)
  a[,1] = a1

  if (s > 0) {
    # solve state space system
    .Call(`_RLDM_outputs_STSP_cpp`, A, B, C, D, u, a, y)
  } else {
    y = D %*% u
  }

  out = list(y = t(y), u = t(u), a = t(a))
  return(out)
}


#' @rdname solve_de
#' @export
solve_de.lmfd = function(sys, u, u0 = NULL, y0 = NULL, ...) {

  a = unclass(sys$a)
  b = unclass(sys$b)
  m = dim(a)[1]
  p = dim(a)[3] - 1
  n = dim(b)[2]
  q = dim(b)[3] - 1

  if ( (m*(p+1)) == 0 ) stop('illegal ARMA system (m=0 or p<0)')

  # check 'u'
  if ( (!is.matrix(u)) || (ncol(u) != n) ) stop('disturbances "u" are not compatible to the model!')
  n.obs = nrow(u)

  # check starting values u0
  if (is.null(u0)) {
    u0 = matrix(0, nrow = max(p,q), ncol = n)
  }
  if ( (!is.matrix(u0)) || (ncol(u0) != n) ) stop('initial disturbances "u0" are not compatible to the model!')
  # just keep the last max(p,q) rows / prepend zero rows if necessary
  qq = min(nrow(u0), max(p,q))
  u0 = rbind(matrix(0, nrow = max(p,q) - qq, ncol = n), u0[iseq(nrow(u0)-qq+1,nrow(u0)),,drop = FALSE])

  # check starting values y0
  if (is.null(y0)) {
    y0 = matrix(0, nrow = max(p,q), ncol = m)
  }
  if ( (!is.matrix(y0)) || (ncol(y0) != m) ) stop('initial outputs "y0" are not compatible to the model!')
  # just keep the last max(p,q) rows / prepend zero rows if necessary
  pp = min(nrow(y0), max(p,q))
  y0 = rbind(matrix(0, nrow = max(p,q) - pp, ncol = m), y0[iseq(nrow(y0)-pp+1, nrow(y0)),,drop = FALSE])

  if (n.obs == 0) {
    return(list(y = matrix(0, nrow = n.obs, ncol = m), u = u))
  }

  # convert ARMA system
  #    a[0] y[t] + a[1] y[t-1] + ... = b[0] u[t] + b[1] u[t-1] + ...
  # to system of the form
  #    y[t] = a[1] y[t-1] + ... + b[0] u[t] ...
  # and create parameter matrices a = (a[p], ..., a[1]) and b = (b[0], ..., b[q])
  a0 = matrix(a[,,1], nrow = m, ncol = m)

  dim(b) = c(m,n*(q+1))
  if ((n*(q+1)) > 0) {
    b = solve(a0, b)
  }

  # note for outputs_ARMA_cpp we have to reshuffle the AR parameter as:  a = (a[p],...,a[1])
  if (p>0) {
    a = a[,,(p+1):2, drop = FALSE]
    dim(a) = c(m, m*p)
    a = -solve(a0, a)
  } else {
    a = matrix(0, nrow = m, ncol = 0)
  }


  # transpose u matrix and prepend initial values
  u = t(rbind(u0,u))     # (n, max(p,q)+n.obs)

  # create y matrix (including starting values)
  y = t(rbind(y0, matrix(0, nrow = n.obs, ncol = m)))  # (m, max(p,q)+n.obs)

  # solve ARMA system
  .Call(`_RLDM_outputs_ARMA_cpp`, a, b, 1+max(p,q), u, y)

  y = t(y[, (max(p,q)+1):(n.obs+max(p,q)), drop = FALSE])
  u = t(u[, (max(p,q)+1):(n.obs+max(p,q)), drop = FALSE])

  return(list(y = y, u = u))
}


#' @rdname solve_de
#' @export
solve_de.rmfd = function(sys, u, u0 = NULL, y0 = NULL, ...) {

  # Interger-valued paramaters ####
  polm_c = sys$c
  polm_d = sys$d
  dim_out = dim(polm_d)[1] # dimension of output
  dim_in = dim(polm_c)[1] # dimension of input
  deg_c = dim(polm_c)[3] # dim.polm()! I need generic fct of polm objects later
  deg_d = dim(polm_d)[3] # dim.polm() = degree of poly

  # Check integer-values parameters
  if ( (dim_out*(deg_c+1)) == 0 ) stop('Illegal RMFD system: Output dimension "dim_out" is zero or c(z) is identically zero.')

  # Inputs et al ####
  # Check inputs 'u'
  if ( (!is.matrix(u)) || (ncol(u) != dim_in) ) stop('Disturbances "u" are not compatible with the model!')
  n_obs = nrow(u)

  # check starting values u0
  if (is.null(u0)) {
    u0 = matrix(0, nrow = max(deg_c, deg_d), ncol = dim_in)
  }
  if ( (!is.matrix(u0)) || (ncol(u0) != dim_in) ) stop('initial disturbances "u0" are not compatible with the model!')

  # just keep the last max(p,q) rows / prepend zero rows if necessary
  qq = min(nrow(u0), max(deg_c, deg_d))
  u0 = rbind(matrix(0, nrow = max(deg_c, deg_d) - qq, ncol = dim_in),
             u0[iseq(nrow(u0)-qq+1,nrow(u0)),,drop = FALSE])

  # check starting values y0
  if (is.null(y0)) {
    y0 = matrix(0, nrow = max(deg_c, deg_d), ncol = dim_out)
  }
  if ( (!is.matrix(y0)) || (ncol(y0) != dim_out) ) stop('initial outputs "y0" are not compatible to the model!')
  # just keep the last max(deg_c, deg_d) rows / prepend zero rows if necessary
  pp = min(nrow(y0), max(deg_c, deg_d))
  y0 = rbind(matrix(0, nrow = max(deg_c, deg_d) - pp, ncol = dim_out), y0[iseq(nrow(y0)-pp+1, nrow(y0)),,drop = FALSE])

  # Check sample size/number of observations
  if (n_obs == 0) {
    return(list(y = matrix(0, nrow = n_obs, ncol = dim_out), u = u))
  }

  # Finalise data_input matrix: transpose u matrix and prepend initial values
  u = t(rbind(u0,
              u))     # (dim_in, max(deg_c, deg_d)+n_obs)

  # Initialise the data_output matrix: y matrix (including starting values) will be overwritten
  y = t(rbind(y0,
              matrix(0, nrow = n_obs, ncol = dim_out)))  # (dim_out, max(deg_c, deg_d)+n_obs)


  # Transform polynomials ####
  # Convert RMFD system (Even though it should be )
  #   y[t] = d[0] v[t] + ... + d[q] v[t-q]
  #     c[0] v[t] + c[1] v[t-1] + ... + c[p] v[t-p]  = u[t]
  # such that c[0] = I_q
  # and create parameter matrices c = (c[p], ..., c[1]) and d = (d[0], ..., d[q])
  c0 = matrix(unclass(polm_c)[,,1], nrow = dim_in, ncol = dim_in)
  c0is_identity = all(c0 == diag(dim_in))

  if(dim_in>1 && !c0is_identity){
    # Post-multiply c(z) and d(z) with inverse of c[0]
    polm_c = polm_c %r% solve(c0)
    polm_d = polm_d %r% solve(c0)
  }

  polm_c_rev = unclass(polm_c)[,,(deg_c+1):1, drop = FALSE]
  dim(polm_c_rev) = c(dim_in, dim_in*(deg_c+1))
  if(deg_c>0){
    polm_c_rev = polm_c_rev[, 1:(dim_in*deg_c), drop = FALSE]
  } else {
    polm_c_rev = matrix(0, nrow = dim_in, ncol = 0)
  }

  polm_d = unclass(polm_d)
  dim(polm_d) = c(dim_out, dim_in*(deg_d+1))

  # Rcpp function call: Solve RMFD system ####
  #.Call(`_RLDM_solve_rmfd_cpp`, poly_inv, poly_fwd, data_in, data_out, t0)
  .Call(`_RLDM_solve_rmfd_cpp`, polm_c_rev, polm_d, u, y, max(deg_c,deg_d)+1)

  y = t(y[, (max(deg_c,deg_d)+1):(n_obs+max(deg_c,deg_d)), drop = FALSE])
  u = t(u[, (max(deg_c,deg_d)+1):(n_obs+max(deg_c,deg_d)), drop = FALSE])

  return(list(y = y, u = u))
}

## 1.2) Calculate Inputs from (Observed) Outputs ----

#' @rdname solve_de
#' @export
solve_inverse_de = function(sys, y, ...)  {
  UseMethod("solve_inverse_de", sys)
}


#' @rdname solve_de
#' @export
solve_inverse_de.stsp = function(sys, y, a1 = NULL, ...) {
  A = sys$A
  B = sys$B
  C = sys$C
  D = sys$D

  m = nrow(D)
  n = ncol(D)
  s = nrow(A)
  if (m != n) stop('"solve_inverse_de" is only implemented for square systems!')

  # check 'y'
  if ( (!is.matrix(y)) || (ncol(y) != m) ) stop('outputs "y" are not compatible to the model!')
  n.obs = nrow(y)

  # check 'a1'
  if (is.null(a1)) {
    a1 = double(s)
  }
  if ( (!is.vector(a1)) || ( length(a1) != s ) ) stop('initial state "a1" is not compatible!')

  if (n.obs == 0) {
    return(list(y = y, u = matrix(0, nrow = n.obs, ncol = n), a = matrix(a1, nrow = 1, ncol = s)))
  }

  # note: residuals_STSP_cpp() uses transposed data matrices!!!
  y = t(y)
  u = matrix(0, nrow = n, ncol = n.obs)
  a = matrix(0, nrow = s, ncol = n.obs + 1)
  a[,1] = a1
  dPI = matrix(0, nrow = 0, ncol = 0)
  dU = matrix(0, nrow = 0, ncol = 0)

  if (s > 0) {
    # solve (inverse) state space system
    .Call(`_RLDM_residuals_STSP_cpp`, A, B, C, D, y, a, u, dPI, dU)
  } else {
    u = solve(D, y)
  }

  out = list(y = t(y), u = t(u), a = t(a))
  return(out)
}

#' @rdname solve_de
#' @export
solve_inverse_de.lmfd = function(sys, y, u0 = NULL, y0 = NULL, ...) {

  a = unclass(sys$a)
  b = unclass(sys$b)
  m = dim(a)[1]
  p = dim(a)[3] - 1
  n = dim(b)[2]
  q = dim(b)[3] - 1

  if ( (m != n) || (min(m, p+1, q+1) <= 0) ) {
    stop('"solve_inverse_de" is only implemented for square systems (m = n > 0 and p,q >=0)!')
  }

  # check 'y'
  if ( (!is.matrix(y)) || (ncol(y) != m) ) stop('outputs "y" are not compatible to the model!')
  n.obs = nrow(y)

  # check starting values u0
  if (is.null(u0)) {
    u0 = matrix(0, nrow = max(p,q), ncol = n)
  }
  if ( (!is.matrix(u0)) || (ncol(u0) != n) ) stop('initial disturbances "u0" are not compatible to the model!')
  # just keep the last max(p,q) rows / prepend zero rows if necessary
  qq = min(nrow(u0), max(p,q))
  u0 = rbind(matrix(0, nrow = max(p,q) - qq, ncol = n), u0[iseq(nrow(u0)-qq+1,nrow(u0)),,drop = FALSE])

  # check starting values y0
  if (is.null(y0)) {
    y0 = matrix(0, nrow = max(p,q), ncol = m)
  }
  if ( (!is.matrix(y0)) || (ncol(y0) != m) ) stop('initial outputs "y0" are not compatible to the model!')
  # just keep the last max(p,q) rows / prepend zero rows if necessary
  pp = min(nrow(y0), max(p,q))
  y0 = rbind(matrix(0, nrow = max(p,q) - pp, ncol = m), y0[iseq(nrow(y0)-pp+1, nrow(y0)),,drop = FALSE])

  if (n.obs == 0) {
    return(list(y = y, u = matrix(0, nrow = n.obs, ncol = n)))
  }

  # convert ARMA system
  #    a[0] y[t] + a[1] y[t-1] + ... = b[0] u[t] + b[1] u[t-1] + ...
  # to system of the form
  #    u[t] = b[1] u[t-1] + ... + a[0] y[t] ...
  # and create parameter matrices b = (b[q], ..., b[1]) and a = (a[0], ..., a[p])
  b0 = matrix(b[,,1], nrow = m, ncol = m)

  dim(a) = c(m, m*(p+1))
  if ((m*(p+1)) > 0) {
    a = solve(b0, a)
  }

  # note for outputs_ARMA_cpp we have to reshuffle the MA parameter as:  b = (b[q],...,q[1])
  if (q > 0) {
    b = b[,,(q+1):2, drop = FALSE]
    dim(b) = c(m, m*q)
    b = -solve(b0, b)
  } else {
    b = matrix(0, nrow = m, ncol = 0)
  }

  # transpose y matrix and prepend initial values
  y = t(rbind(y0, y))     # (m, max(p,q)+n.obs)

  # create u matrix (including starting values)
  u = t(rbind(u0, matrix(0, nrow = n.obs, ncol = n)))  # (n, max(p,q)+n.obs)

  # solve ARMA system
  # we could also use "residuals_ARMA_cpp"!
  .Call(`_RLDM_outputs_ARMA_cpp`, b, a, max(p,q)+1, y, u)

  y = t(y[, (max(p,q)+1):(n.obs+max(p,q)), drop = FALSE])
  u = t(u[, (max(p,q)+1):(n.obs+max(p,q)), drop = FALSE])

  return(list(y = y, u = u))
}

#' @rdname solve_de
#' @export
solve_inverse_de.rmfd = function(sys, y, u0 = NULL, y0 = NULL, ...) {

  # Integer-valued paramaters ####
  polm_c = sys$c
  polm_d = sys$d
  dim_out = dim(polm_d)[1] # dimension of output
  dim_in = dim(polm_c)[1] # dimension of input
  deg_c = dim(polm_c)[3] # dim.polm()! I need generic fct of polm objects later
  deg_d = dim(polm_d)[3] # dim.polm() = degree of poly

  # Check integer-values parameters
  if ( (dim_out*(deg_c+1)) == 0 ) stop('Illegal RMFD system: Output dimension "dim_out" is zero or c(z) is identically zero.')

  # Inputs et al ####
  # Check inputs 'y'
  if ( (!is.matrix(y)) || (ncol(y) != dim_out) ){
    stop('Disturbances "u" are not compatible with the model!')
  }
  n_obs = nrow(y)

  # check starting values u0
  if (is.null(u0)) {
    u0 = matrix(0, nrow = max(deg_c, deg_d), ncol = dim_in)
  }
  if ( (!is.matrix(u0)) || (ncol(u0) != dim_in) ) stop('initial disturbances "u0" are not compatible with the model!')

  # just keep the last max(p,q) rows / prepend zero rows if necessary
  qq = min(nrow(u0), max(deg_c, deg_d))
  u0 = rbind(matrix(0, nrow = max(deg_c, deg_d) - qq, ncol = dim_in),
             u0[iseq(nrow(u0)-qq+1,nrow(u0)),,drop = FALSE])

  # check starting values y0
  if (is.null(y0)) {
    y0 = matrix(0, nrow = max(deg_c, deg_d), ncol = dim_out)
  }
  if ( (!is.matrix(y0)) || (ncol(y0) != dim_out) ) stop('initial outputs "y0" are not compatible to the model!')
  # just keep the last max(deg_c, deg_d) rows / prepend zero rows if necessary
  pp = min(nrow(y0), max(deg_c, deg_d))
  y0 = rbind(matrix(0, nrow = max(deg_c, deg_d) - pp, ncol = dim_out),
             y0[iseq(nrow(y0)-pp+1, nrow(y0)), , drop = FALSE])

  # Check sample size/number of observations
  if (n_obs == 0) {
    return(list(y = matrix(0, nrow = n_obs, ncol = dim_out), u = u))
  }

  # Finalise data_input matrix: transpose y matrix and prepend initial values
  y = t(rbind(y0,
              y))     # (dim_out, max(deg_c, deg_d)+n_obs)

  # Initialise the data_output matrix: y matrix (including starting values) will be overwritten
  u = t(rbind(u0,
              matrix(0, nrow = n_obs, ncol = dim_in)))  # (dim_in, max(deg_c, deg_d)+n_obs)

  # Transform polynomials ####
  # y[t] is given
  # 1) Obtain w[t] = [d0inv d(z)]^{-1} [d0inv y[t]] => we need d = d0inv*(d[q], ..., d[1])
  # 2) Obtain u[t] = c(z) w[t] => we need c = (I_q, c[1], ..., c[p])

  # __Transform c(z) to (I_q, c_1, ... , c_p) ####
  # Check whether c[0] is identity, if not => transform
  c0 = matrix(unclass(polm_c)[,,1], nrow = dim_in, ncol = dim_in)
  c0is_identity = all(c0 == diag(dim_in))

  if(dim_in>1 && !c0is_identity){
    # Post-multiply c(z) and d(z) with inverse of c[0]
    polm_c = polm_c %r% solve(c0)
    polm_d = polm_d %r% solve(c0)
  }

  # Reshuffling
  polm_c = unclass(polm_c)
  dim(polm_c) = c(dim_in, dim_in*(deg_c+1))

  # Prepare d(z): (d_q, ... , d_1) ####
  d0 = matrix(unclass(polm_d)[,,1], dim_out, dim_in)
  d0inv = MASS::ginv(d0)

  if(deg_d>0){
    polm_d = unclass(polm_d)
    polm_d = polm_d[,,(deg_d+1):2]
    dim(polm_d) = c(dim_out, dim_in*deg_d)
    polm_d = d0inv %*% polm_d
  } else {
    polm_d = matrix(0, nrow = dim_out, ncol = 0)
  }

  # Transform data (enough for the case deg_d=0) ####
  y = d0inv %*% y

  # Rcpp function call: Solve RMFD system ####
  #.Call(`_RLDM_solve_rmfd_cpp`, poly_inv, poly_fwd, data_in, data_out, t0)
  .Call(`_RLDM_solve_rmfd_cpp`, polm_d, polm_c, y, u, max(deg_c,deg_d)+1)

  y = t(y[, (max(deg_c,deg_d)+1):(n_obs+max(deg_c,deg_d)), drop = FALSE])
  u = t(u[, (max(deg_c,deg_d)+1):(n_obs+max(deg_c,deg_d)), drop = FALSE])

  return(list(y = y,
              u = u))
}

# 2) Toeplitz Calculations ----


#' Toeplitz Calculations
#'
#' Multiplication of stacked data vector with a block Toeplitz matrix (`toepl_fwd()` for MA calculations) or
#' "inversion" of a block Toeplitz matrix in order to perform calculations equivalent to
#' multiplying a given stacked data vector with the inverse of a lower-triangular banded block Toeplitz matrix (`toepl_inv()` for AR calculations).
#' Note that matrix polynomials can be mapped one-to-one to banded lower-triangular block Toeplitz matrices.
#'
#' @section MA-type Toeplitz calculations:
#' Given a polynomial matrix of degree \eqn{q} and dimension \eqn{(m \times n)}{(m x n)}, where \eqn{m \geq n}{m \ge n},
#' and given a "wide" input data matrix of dimension \eqn{(n \times nobs)}{(n x nobs)}, where \eqn{nobs} is the number of observations
#' such that each column corresponds to one observation and the number of columns is equal to the number of observations,
#' we calculate a "wide" output data matrix of dimension \eqn{(m \times nobs)}.
#'
#' The function name `toepl_fwd` stems from the multiplication of the "stacked" input data vector
#' \eqn{(u_1', \ldots , u_{nobs}')'}{(u[1]', ... , u[nobs]')'} with a banded lower-triangular block Toeplitz matrix \eqn{T}
#' of dimension \eqn{(nobs m \times nobs n)}{(nobs m x nobs n)}
#' whose block elements depend only on the difference between the row- and column-index
#' such that  \eqn{T_{i,j} = d_{i-j}}{T[i,j]=d[i-j]}.
#'
#' @section AR-type Toeplitz calcuations:
#' Given a square polynomial matrix \eqn{c(z)} for which \eqn{c_0}{c[0]} is the identity matrix and given
#' a wide data matrix `y = data_in`, obtain the solution `u`, a wide data matrix, of the Toeplitz equation
#' \deqn{T (y_1', \ldots , y_{nobs}')' =  (u_1', ... , u_{nobs}')'}{
#'       T (y[1]', \ldots , u[nobs]')' =  (u[1]', ... , u[nobs]')'}
#'
#' Note that the zero-lag coefficient is discarded and the coefficients are in reversed order
#' since this simplifies computations and implementation.
#'
#' @param polm_wide Wide matrix \eqn{( d_0, d_1, \ldots, d_q )}{(d[0], d[1], ... , d[q])}
#'     of dimension \eqn{(m \times (q+1) n)}{(m x (q+1)n)} which represents a
#'     matrix polynomial \eqn{d(z)} of degree \eqn{q}.
#' @param polm_rev Wide matrix \eqn{(c_p, ... , c_1)}{(c[p], ... , c[1])} of dimension \eqn{(n \times (q+1) n)}{(n x p n)}
#'     with coefficients ordered in reverse direction, and
#'     zero-lag coefficient matrix not included.
#'     It represents a square polynomial matrix \eqn{c(z)} with \eqn{c_0}{c[0]} equal to the identity matrix and of degree \eqn{p}.
#' @param data_in Data matrix of dimension \eqn{(dim(inputs) x nobs)}.
#'     Every column corresponds to one observation.
#' @param t0 Integer. Time index where calculations should start.
#'     Default set to 1. For AR calculations, \eqn{degree + 1} would be another smart option.
#'
#' @return data_out Data matrix of dimension \eqn{(dim_out x n_obs)}
#' @export
#'
#' @examples
#' p = test_polm(dim = c(2,2), degree = 2) %>% unclass()
#' dim(p) = c(2,2*3)
#' data = matrix(rnorm(2*100), 2, 100)
#' toepl_fwd(p, data)
toepl_fwd = function(polm_wide, data_in, t0 = 1){


  # (d_0, d_1, ... , d_q)

  # Integer-valued parameters
  dim_out = nrow(polm_wide)
  dim_in = nrow(data_in)
  n_obs = ncol(data_in)
  deg_p = ncol(polm_wide)/dim_in - 1

  # Allocate data_out
  data_out = matrix(0, nrow = dim_out, ncol = n_obs)

  # Append starting values if t0 is not big enough
  add_zeros = max(deg_p-t0+1, 0)
  if (add_zeros>0){
    data_in = cbind(matrix(0, dim_in, add_zeros),
                    data_in)
    n_obs = n_obs + add_zeros
    t0 = t0 + add_zeros
  }

  for (ix in 0:deg_p){
    data_out = data_out + polm_wide[ , (1+ix*dim_in):((ix+1)*dim_in), drop = FALSE] %*% data_in[, (t0-ix):(n_obs-ix), drop = FALSE]
  }

  return(data_out)
}


#' @aliases toeplitz_calculations
#' @rdname toepl_fwd
#'
#' @export
#' @examples
#' p = test_polm(dim = c(2,2), degree = 2) %>% unclass()
#' dim(p) = c(2,2*3)
#' data = matrix(rnorm(2*100), 2, 100)
#' toepl_inv(p, data)
toepl_inv = function(polm_rev, data_in, t0 = 1){

  # polm = (c_p, ... , c_1)

  # Integer-valued parameters
  dim_out = nrow(polm_rev)
  dim_in = nrow(data_in)
  n_obs = ncol(data_in)
  deg_p = ncol(polm_rev)/dim_in

  # Fail or return early
  stopifnot(dim_out == dim_in)
  if (deg_p == 0){
    return(data_in)
  }

  # Append starting values if t0 is not big enough
  add_zeros = deg_p-t0+1
  if (add_zeros > 0){
    data_in = cbind(matrix(0, dim_in, add_zeros),
                    data_in)
    n_obs = n_obs + add_zeros
    t0 = t0 + add_zeros
  }

  if(deg_p>0){
    for (ix_t in t0:n_obs){
      data_in[, ix_t] = data_in[, ix_t, drop = FALSE] - polm_rev %*% matrix(c(data_in[, (ix_t-deg_p):(ix_t-1)]), ncol = 1)
    }
  }

  return(data_in[, t0:n_obs, drop = FALSE])
}

# 3) RMFD specific ----


#' Solve RMFD system for given inputs
#'
#' Compute the outputs \eqn{y_t}{y[t]} of RMFD(p, q) systems of the form
#' \deqn{y_t = d_0 v_t + d_1 v_{t-1} + \cdots + d_q v_{t-q}}{y[t] = d[0] v_t + d[1] v[t-1] + ... + d[q] v[t-q]}
#' where
#' \deqn{v_t + c_1 v_{t-1} + \cdots + c_p v_{t-p} = u_t}{v[t] + c[1] v[t-1] + ... + c[p] v[t-p] = u[t]}
#' for given inputs \eqn{u_t}{u[t]} which are contained column-wis in the data matrix `data_input`, see below.
#'
#' Values \eqn{y_t}{y[n]}, \eqn{u_t}{u[t]} for \eqn{t\leq 0}{t\le 0} are implicitely set to be zero.
#' However, if we start the iteration with some \eqn{t_0>1}{t0>1} we can enforce non-zero initial values.
#'
#' The routines are used internally and hence do **not** check their arguments.
#' We require the number of outputs to be positive \eqn{m > 0}, the number of inputs to be non-negative \eqn{n \geq 0}{n \ge 0},
#' the degree of \eqn{c(z)} %>%  to be non-negative \eqn{p \geq 0}{p \ge 0}, the degree of \eqn{d(z)} is unrestricted, i.e. \eqn{(q+1) \geq 0}{(q+1) \ge 0},
#' and for the starting value and the sample size \eqn{1 \leq t_0 \leq N}{1 \le t0 \le N} holds.
#' Note also that the RcppArmadillo implementation **overwrites** the input argument `y`.
#' Use this procedure with care!
#'
#' Note the non standard arguments:
#' The polynomial matrices \eqn{c(z)} and \eqn{d(z)} are saved as "wide" matrices.
#' The order of the coefficients in \eqn{c(z)} is reversed and there is no \eqn{c_0}{c[0]} coefficient (because it is required to be the identity matrix).
#' The order of the coefficients in \eqn{d(z)} is as usual, \eqn{d_0}{d[0]} is available too.
#' The data matrices are organized columnwise (to avoid memory shuffling)!
#'
#' @param polm_c,polm_d [polm] objects. Describe jointly RMFD.
#'     \eqn{c(z)} is square, has the identity as zero-lag coefficient, and its coefficients will be reversed in the procedure: \eqn{(c_p,...,c_1)}{(c[p],...,c[1])}.
#'     \eqn{d(z)} might be tall, and its zero-lag coefficient matrix is in general free.
#' @param data_input \eqn{(n, N)} matrix with the inputs \eqn{(u_1,...,u_N}{(u[1],...,u[N])}.
#' @param t0 integer, start iteration at t=t0.
#'
#' @return The R implementation `solve_RMFD_R()` returns the matrix `y` with the computed outputs.
#'     \deqn{y_t = d(z) c(z)^{-1} u_t}{y[t] = d{+}(z) c(z)^{-1} u[t]}
#'     The internal RcppArmadillo implementation returns `NULL` but **overwrites** the input argument.
#'     Note that the RcppArmadillo implementation has a different user interface (it is intended for internal use only).
#'
#' @export
solve_RMFD_R = function(polm_c, polm_d, data_input, t0 = 1) {

  # Unclassing polm objects
  polm_c = unclass(polm_c)
  polm_d = unclass(polm_d)

  # Integer-valued parameters
  dim_input = nrow(data_input) # dimension of inputs
  dim_output = dim(polm_d)[1] # dimension of outputs
  deg_c = dim(polm_c)[3]-1 # degree of c(z); needs to be transformed to (c_p, ... , c_1)
  deg_d = dim(polm_d)[3]-1 # degree of d(z); needs to be transformed to (d_0, d_1, ... , d_q)


  # Prepare c(z): (c_p, ... , c_1)
  if(deg_c > 0){
    polm_c = polm_c[,,(deg_c+1):2, drop = FALSE]
    dim(polm_c) = c(dim_input, dim_input*deg_c)
  }

  # Calculate v_t = c(z)^{-1} u_t
  if (deg_c > 0){ # nothing to do in case deg_c == 0
    data_input = toepl_inv(polm_rev = polm_c, data_in = data_input, t0)
  }

  # Prepare d(z): (d_0, d_1, ... , d_q)
  dim(polm_d) = c(dim_output, dim_input*(deg_d+1))

  # Calculate y_t = d(z) v_t
  data_output = toepl_fwd(polm_wide = polm_d, data_in = data_input, t0)

  return(data_output)
}


#' Obtain Inputs of RMFD System for Given Data
#'
#' Compute the inputs \eqn{u_t}{u[t]} of RMFD(p, q) systems of the form
#' \deqn{y_t = d_0 v_t + d_1 v_{t-1} + \cdots + d_q v_{t-q}}{y[t] = d[0] v[t] + d[1] v[t-1] + ... + d[q] v[t-q]}
#' where
#' \deqn{v_t + c_1 v_{t-1} + \cdots + c_p v_{t-p} = u_t}{v[t] + c[1] v[t-1] + ... + c[p] v[t-p] = u[t]}
#' for given data \eqn{y_t}{y[t]} which are contained column-wise in the matrix `data_output`, see below.
#'
#' Values \eqn{y_t}{y[n]}, \eqn{u_t}{u[t]} for \eqn{t\leq 0}{t\le 0} are implicitely set to be zero.
#' However, if we start the iteration with some \eqn{t_0>1}{t0>1} we can enforce non-zero initial values.
#'
#' The routines are used internally and hence do **not** check their arguments.
#' We require the number of outputs to be positive \eqn{m > 0}, the number of inputs to be non-negative \eqn{n \geq 0}{n \ge 0},
#' the degree of \eqn{c(z)} to be non-negative \eqn{p \geq 0}{p \ge 0}, the degree of \eqn{d(z)} is unrestricted, i.e. \eqn{(q+1) \geq 0}{(q+1) \ge 0},
#' and for the starting value and the sample size \eqn{1 \leq t_0 \leq N}{1 \le t0 \le N} holds.
#' Note also that the RcppArmadillo implementation **overwrites** the input argument `y`.
#' Use this procedure with care!
#'
#' Note the non standard arguments:
#' The polynomial matrices \eqn{c(z)} and \eqn{d(z)} are saved as "wide" matrices.
#' The order of the coefficients in \eqn{c(z)} is reversed and there is no \eqn{c_0} coefficient (because it is required to be the identity matrix).
#' The order of the coefficients in \eqn{d(z)} is as usual, \eqn{d_0} is available too.
#' The data matrices are organized columnwise (to avoid memory shuffling)!
#'
#' @param polm_c,polm_d polm objects. Jointly, they reprsent an RMFD.
#'     \eqn{c(z)} is square, has the identity as zero-lag coefficient.
#'     \eqn{d(z)} might be tall, and its zero-lag coefficient matrix is in general free.
#' @param data_output \eqn{(m, N)} matrix with the outputs \eqn{(y_1,...,y_N}{(y[1],...,y[N])}.
#' @param t0 integer, start iteration at t=t0.
#'
#' @return The R implementation `solve_inverse_RMFD_R` returns the matrix `u` with the computed inputs
#'     \deqn{u_t = d^{+}(z) c(z) y_t}{u[t] = d^{+}(z) c(z) y[t]} in its columns.
#'     The internal RcppArmadillo implementation returns `NULL` but **overwrites** its input argument!
#'     Note that the RcppArmadillo implementation has a different user interface (it is intended for internal use only).
#'
#' @export
#' @seealso solve_RMFD_R
solve_inverse_RMFD_R = function(polm_c, polm_d, data_output, t0=1) {

  # Unclassing polm objects
  polm_c = unclass(polm_c)
  polm_d = unclass(polm_d)

  # Integer-valued parameters
  dim_input = ncol(polm_d) # dimension of inputs
  dim_output = nrow(polm_d) # dimension of outputs
  deg_c = dim(polm_c)[3]-1 # degree of c(z); needs to be transformed to (c_0, c_1,  ... , c_p)
  deg_d = dim(polm_d)[3]-1 # degree of d(z); needs to be transformed to (d_q, ... , d_1)

  # Prepare d(z): (d_q, ... , d_1)
  d0 = matrix(polm_d[,,1], dim_output, dim_input)
  d0inv = MASS::ginv(d0)

  if(deg_d>0){
    polm_d = polm_d[,,(deg_d+1):2]
    dim(polm_d) = c(dim_output, dim_input*deg_d)
    polm_d = d0inv %*% polm_d
  } else { # maybe not needed bc if deg_d == 0, then fct toepl_inv() is not invoked
    polm_d = matrix(0, nrow = dim_output, ncol = 0)
  }

  # Prepare data (enough for the case deg_d=0)
  data_output = d0inv %*% data_output

  # Calculate v_t in (d0inv * y_t) = (d0inv * d(z)) v_t
  if (deg_d>0){
    data_output = toepl_inv(polm_rev = polm_d, data_in = data_output, t0 = t0)
  }

  # Prepare c(z): (c_0, c_1, ... , c_p)
  dim(polm_c) = c(dim_input, dim_input*(deg_c+1))

  # Calculate u_t in v_t = c(z)^{-1} u_t
  data_output = toepl_fwd(polm_wide = polm_c, data_in = data_output, t0 = t0)

  return(data_output)
}

# 4) Checks: R functions for checking RcppArmadillo implementation ----

#' Solve ARMA system
#'
#' Compute the outputs of ARMA(p, q) systems of the form
#' \deqn{y_t = a_1 y_{t-1} + ... + a_p y_{t-p} + b_0 u_t + \cdots + b_q u_{t-q}}{
#'       y[t] = a[1] y[t-1] + ... + a[p] y[t-p] + b[0] u[t] + ... + b[q] u[t-q]}
#'
#' Values \eqn{y_t}{y[n]}, \eqn{u_t}{u[t]} for \eqn{t\leq 0}{t\le 0} are implicitly set to be zero.
#' However, if we start the iteration with some \eqn{t_0>1}{t0>1} we can enforce non-zero
#' initial values.
#'
#' The routines are used internally and hence do **not** check their arguments.
#' We require \eqn{m > 0}, \eqn{p \geq 0}{p \ge 0}, \eqn{n \geq 0}{n \ge 0}, \eqn{(q+1) \geq 0}{(q+1) \ge 0}
#' and \eqn{1 \leq t_0 \leq N}{1 \le t0 \le N}.
#' Note also that the RcppArmadillo implementation **overwrites** the input argument `y`.
#' Use this procedure with care!
#'
#' Note the non standard arguments: The order of the AR coefficients is reversed.
#' The data matrices are organized column-wise (to avoid memory shuffling)!
#'
#' @param a \eqn{(m, mp)} matrix \eqn{(a_p,...,a_1)}{(a[p],...,a[1])}.
#' @param b \eqn{(m, n(q+1))} matrix \eqn{(b_0,...,b_q}{(b[0],...,b[q])}.
#' @param u \eqn{(n, N)} matrix with the inputs \eqn{(u_1,...,u_N}{(u[1],...,u[N])}.
#' @param y \eqn{(m, N)} matrix with the outputs \eqn{(y_1,...,y_N}{(y[1],...,y[N])}.
#' @param t0 integer, start iteration at t=t0.
#'
#' @return The R implementation `solve_ARMA_R` returns the matrix `y` with the computed outputs.
#'         The RcppArmadillo implementation returns `NULL` but **overwrites** the input argument `y`!
#'
#' @export
#' @name solve_ARMA
solve_ARMA_R = function(a, b, u, y, t0) {
  m = nrow(y)
  n.obs = ncol(y)
  n = nrow(u)
  p = ncol(a)/m
  q = ncol(b)/n - 1
  if (p > 0) {
    dim(a) = c(m,m,p)
    a = a[,,p:1,drop = FALSE]
  }
  dim(b) = c(m,n,q+1)

  for (t in (t0:n.obs)) {
    y[, t] = matrix(b[,,1], nrow = m, ncol = n) %*% u[,t]
    if (q > 0) {
      for (i in (1:q)) {
        if ((t-i) > 0) y[, t] = y[, t] + matrix(b[ , , i+1], nrow = m, ncol = n) %*% u[, t-i]
      }
    }
    if (p>0) {
      for (i in (1:p)) {
        if ((t-i) > 0) y[, t] = y[, t] + matrix(a[ , , i], nrow = m, ncol = m) %*% y[, t-i]
      }
    }
  }
  return(y)
}

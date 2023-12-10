#' Constructor for LMFD (ARMA) Models
#'
#' A left-matrix fraction description (LMFD) plus parameterisation of noise covariance.
#'
#' In Hannan, Deistler (2012, page 7), RMFDs are also called dynamic adjustment forms.
#' Internally, MFDs are lists with slots `sys`, `sigma_L`, `names`, `label`.
#'
#' @param sys [rationalmatrices::lmfd()] or [rationalmatrices::rmfd()] object
#' @param sigma_L Left-factor of noise covariance,
#'   i.e. the covariance \eqn{\sigma} is obtained as `sigma_L * t(sigma_L)`.
#'   If `sigma_L` is a vector of dimension \eqn{n}, where \eqn{n} is the input dimension, only the diagonal elements are parametrized.
#'   If it is a vector of dimension \eqn{n^2}, then the elements of `sigma_L` are filled column by column.
#' @param names optional vector of character strings
#' @param label optional character string
#'
#' @return Object of class `armamod`.
#'
#' @export
#'
#' @examples
#' x = armamod(sys = lmfd(c(1, 0.5), 1), sigma_L = diag(1))
#' x
armamod = function(sys, sigma_L = NULL, names = NULL, label = NULL) {
  if (!inherits(sys, 'lmfd')) stop('"sys" must be an lmfd object')

  d = dim(sys)
  m = d[1]
  n = d[2]

  if (is.null(sigma_L)) sigma_L = diag(n)
  if (!is.numeric(sigma_L)) stop('parameter sigma_L is not numeric')
  if ( is.vector(sigma_L) ) {
    if (length(sigma_L) == n) sigma_L = diag(sigma_L, nrow = n, ncol = n)
    if (length(sigma_L) == (n^2)) sigma_L = matrix(sigma_L, nrow = n, ncol = n)
  }
  if ( (!is.matrix(sigma_L)) || any(dim(sigma_L) != n) ) {
    stop('"sigma_L" is not compatible')
  }

  x = list(sys = sys, sigma_L = sigma_L, names = names, label = label)
  x = structure(x, class = c('armamod', 'rldm'))
  return(x)
}

#' Constructor for RMFD Models
#'
#' `r lifecycle::badge("experimental")`
#'
#' A right-matrix fraction description (RMFD) plus parameterisation of noise covariance.
#' In \insertCite{Hannan.Deistler12}{RLDM}, RMFDs are also called dynamic adjustment forms.
#' Internally, MFDs are lists with slots `sys`, `sigma_L`, `names`, `label`.
#' **Many of the generic functions which construct derived objects like the autocovariance `autocov()` are not yet implemented for `rmfdmod` objects.**
#'
#' @param sys [rationalmatrices::rmfd()] object
#' @inheritParams armamod
#'
#' @return Object of class `rmfdmod`.
#'
#' @references
#' \insertRef{Hannan.Deistler12}{RLDM}
#'
#' @export
#'
#' @examples
#' y = rmfdmod(sys = test_rmfd(dim = c(3,2), degrees = c(2,2)))
#' y
rmfdmod = function(sys, sigma_L = NULL, names = NULL, label = NULL) {
  # Check sys input
  if (!inherits(sys, 'rmfd')) stop('"sys" must be an rmfd object')

  d = dim(sys)
  m = d[1]
  n = d[2]

  # Parameterisation of sigma through left-factor: Different cases
  if (is.null(sigma_L)) sigma_L = diag(n)
  if (!is.numeric(sigma_L)) stop('parameter sigma_L is not numeric')
  if ( is.vector(sigma_L) ) {
    if (length(sigma_L) == n) sigma_L = diag(sigma_L, nrow = n, ncol = n)
    if (length(sigma_L) == (n^2)) sigma_L = matrix(sigma_L, nrow = n, ncol = n)
  }
  if ( (!is.matrix(sigma_L)) || any(dim(sigma_L) != n) ) {
    stop('"sigma_L" is not compatible')
  }

  x = list(sys = sys, sigma_L = sigma_L, names = names, label = label)
  x = structure(x, class = c('rmfdmod', 'rldm'))
  return(x)
}

#' Creator for stspmod class
#'
#' @param sys [rationalmatrices::stsp()] object
#' @param sigma_L noise covariance left
#' @param names optional vector of character strings
#' @param label optional chracter string
#'
#' @return Object of class `stspmod`.
#'
#' @export
#'
#' @examples
#' x = stspmod(sys = test_stsp(dim = c(2,2), s = 2), sigma_L = diag(2))
#' x
stspmod = function(sys, sigma_L = NULL, names = NULL, label = NULL) {
  if (!inherits(sys, 'stsp')) stop('"sys" must be an stsp object')

  d = dim(sys)
  m = d[1]
  n = d[2]

  if (is.null(sigma_L)) sigma_L = diag(n)
  if (!is.numeric(sigma_L)) stop('parameter sigma_L is not numeric')
  if ( is.vector(sigma_L) ) {
    if (length(sigma_L) == n) sigma_L = diag(sigma_L, nrow = n, ncol = n)
    if (length(sigma_L) == (n^2)) sigma_L = matrix(sigma_L, nrow = n, ncol = n)
  }
  if ( (!is.matrix(sigma_L)) || any(dim(sigma_L) != n) ) {
    stop('"sigma_L" is not compatible')
  }

  x = list(sys = sys, sigma_L = sigma_L, names = names, label = label)
  x = structure(x, class = c('stspmod', 'rldm'))
  return(x)
}


# Test models ----

#' Create Test ARMA model
#'
#' This simple tool may be used to create a random ARMA model
#' \deqn{y_t + a_1 y_{t-1} + \cdots + a_p y_{t-p} = b_0 u_t + b_1 u_{t-1} + \cdots + b_q u_{t-q}}{
#'       y[t] + a[1] y[t-1] + ... + a[p] y[t-p] = b[0] u[t] + b[1] u[t-1] + ... + b[q] u[t-q]}
#' with given order \eqn{(p,q)}.
#'
#' We require \eqn{m>0} and \eqn{p\geq 0}{p\ge 0}.
#'
#' The \eqn{b_0}{b[0]} matrix defaults to a \eqn{(m,n)}-dimensional diagonal matrix
#' with ones on the diagonal (`diag(x=1, nrow = m, ncol = n)`). However, one may
#' also pass an arbitray (compatible) matrix to the procedure.
#' This matrix may contain `NA`'s, which then are replaced by random numbers.
#'
#' The \eqn{sigma_L} matrix defaults to a \eqn{(n,n)}-dimensional lower, triangular matrix
#' However, one may also pass an arbitrary (compatible) \eqn{sigma_L} matrix to the procedure.
#'
#' The user may prescribe lower bounds for the moduli of the zeroes and/or poles of the transfer function
#' \deqn{k(z) = a^{-1}(z) b(z).}
#' In this case the procedure simply generates (up to n.trials) random models until a model is found
#' which satisfies the constraint. The standard deviation of the normal distribution, which is used to
#' generate the random entries, is decreased in each step. Of course this is a very crude method and
#' it may fail or need a very large number of randomly generated matrices.
#'
#' @param dim integer vector `c(m,n)`.
#' @param degrees integer vector `c(p,q)`.
#' @param digits integer, if non NULL then the randomly generated numbers are rounded to
#'               "digits" number of decimal places.
#' @param b0     \eqn{(m,n)} dimensional matrix (or `NULL`). See the details below.
#' @param sigma_L     \eqn{(n,n)} dimensional matrix (or `NULL`). See the details below.
#' @param bpoles lower bound for the moduli of the poles of the corresponding transfer function (or NULL).
#' @param bzeroes lower bound for the moduli of the zeroes of the corresponding tranmsfer function (or NULL).
#'                This parameter is ignored for non-square matrices (m != n).
#' @param n.trials maximum number of trials.
#'
#' @return [armamod()] object, which represents the generated ARMA model.
#' @export
#'
#' @examples
#' ### generate a random ARMA(1,1) model (with two outputs)
#' ### we require that the model is stable and minimum phase
#' model = try(test_armamod(dim = c(2,2), degrees = c(1,1), digits = 2, bpoles = 1, bzeroes = 1))
#' if (!inherits(model, 'try-error')) {
#'    print(model)
#'    print(abs(poles(model$sys)))
#'    print(abs(zeroes(model$sys)))
#' }
test_armamod = function(dim = c(1,1), degrees = c(1,1), b0 = NULL, sigma_L = NULL,
                        digits = NULL, bpoles = NULL, bzeroes = NULL, n.trials = 100) {
  # check input parameter "dim"
  dim = as.integer(dim) # note: as.integer converts to vector!
  if ((length(dim) != 2) || (dim[1] <= 0) || (dim[2] < 0)) {
    stop('argument "dim" must be a vector of integers with length 2, dim[1] > 0 and dim[2] >= 0!')
  }
  # check input parameter "degrees"
  degrees = as.integer(degrees) # note: as.integer converts to vector!
  if ((length(degrees) != 2) || (degrees[1] < 0) || (degrees[2] < (-1))) {
    stop('argument "degrees" must be a vector of integers with length 2, degrees[1] >= 0 and degrees[2] >= -1!')
  }

  m = dim[1]
  n = dim[2]
  p = degrees[1]
  q = degrees[2]

  if (p == 0) bpoles = NULL
  if ( (m != n) || (q <= 0) ) bzeroes = NULL

  # check input parameter "b0"
  if (!is.null(b0)) {
    if ( (!is.numeric(b0)) || (!is.matrix(b0)) || any(dim(b0) != dim) ) {
      stop('"b0" must be a compatible, numeric matrix')
      i = is.na(b0)
      theta = stats::rnorm(sum(i))
      if (!is.null(digits)) theta = round(theta, digits)
      b0[i] = theta
    }
  }
  else {
    b0 = diag(x = 1, nrow = m, ncol = n)
  }

  # check input parameter "sigma_L"
  if (!is.null(sigma_L)) {
    if ( (!is.numeric(sigma_L)) || (!is.matrix(sigma_L)) || any(dim(sigma_L) != n) ) {
      stop('"sigma_L" must be a compatible, numeric matrix')
    }
  }
  else {
    sigma_L = matrix(stats::rnorm(n*n), nrow = n, ncol = n)
    if (n >0) {
      sigma_L = t(chol(sigma_L %*% t(sigma_L)))
    }
    if (!is.null(digits)) sigma_L = round(sigma_L, digits)
  }

  i.trial = 0
  err = TRUE
  sd = 1
  while ( (i.trial < n.trials) && (err) ) {
    a = cbind(diag(m), matrix(stats::rnorm(m*m*p, sd = sd), nrow = m, ncol = m*p))
    dim(a) = c(m,m,p+1)
    if (q >= 0) {
      b = cbind(b0, matrix(stats::rnorm(m*n*q, sd = sd), nrow = m, ncol = n*q))
      dim(b) = c(m, n, q+1)
    } else {
      b = array(0, dim = c(m, n, q+1))
    }
    if (!is.null(digits)) {
      a = round(a, digits)
      b = round(b, digits)
    }
    sys = lmfd(a, b)

    err = FALSE
    if ( !is.null(bpoles) ) {
      err = try(min(abs(poles(sys, print_message = FALSE))) <= bpoles, silent = TRUE)
      if (inherits(err, 'try-error')) err = TRUE
    }
    if ( (!err) && (!is.null(bzeroes)) ) {
      err = try((min(abs(zeroes(sys, print_message = FALSE))) <= bzeroes), silent = TRUE)
      if (inherits(err, 'try-error')) err = TRUE
    }
    i.trial = i.trial + 1
    sd = sd/1.1
  }
  if (err) {
    stop('Could not generate a suitable ARMA model with ', n.trials, ' trials!')
  }

  model = armamod(sys = sys, sigma_L = sigma_L)
  return(model)
}

#' Create Test state space Model
#'
#' This simple tool may be used to create a random, state space model
#' \deqn{a_{t+1} = A a_t + B u_t \mbox{ and } y_t = C a_t + D u_t.}{
#'       a[t+1] = A a[t] + B u[t] and y[t] = C a[t] + D u[t].}
#'
#' If the Kronecker indices (parameter `nu`)
#' are given, then a state space model in echelon canonical form is generated. This means that
#' some of the entries of the \eqn{A,B,C} matrices are fixed to be one or zero and the
#' others are considerd as "free". See also [rationalmatrices::Kronecker-Indices()].
#' The entries of the \eqn{A, B, C} matrices, which are not a priori fixed
#' are randomly generated.
#'
#' If only the state dimension \eqn{s} (parameter `s`)
#' is given, then all entries of the \eqn{A, B, C} matrices
#' are considered as "free".
#'
#' The \eqn{D} matrix defaults to a \eqn{(m,n)}-dimensional diagonal matrix
#' with ones on the diagonal (`diag(x=1, nrow = m, ncol = n)`). However, one may
#' also pass an arbitray (compatible) \eqn{D} matrix to the procedure.
#' This matrix may contain `NA`'s, which then are replaced by random numbers.
#'
#' The \eqn{sigma_L} matrix defaults to a \eqn{(n,n)}-dimensional lower, triangular matrix
#' However, one may also pass an arbitray (compatible) \eqn{sigma_L} matrix to the procedure.
#'
#' The user may prescribe lower bounds for the moduli of the zeroes and/or poles of the transfer function
#' \deqn{k(z) = C(I_m z{-1} - A)^{-1} B + D.}{k(z) = C(I z{-1} - A)^{-1} B + D.}
#' In this case the procedure simply generates (up to n.trials) random models until a model is found
#' which satisfies the constraint. The standard deviation of the normal distribution, which is used to
#' generate the random entries, is decreased in each step. Of course this is a very crude method and
#' it may fail or need a very large number of randomly generated matrices.
#'
#' Note also, that the generated model may be non-minimal.
#'
#' @inheritParams test_armamod
#' @param s     state dimension (or NULL).
#' @param nu    vector with the Kronecker indices (or `NULL`). Either the state space dimension `s`
#'              or the Kronecker indices `nu` must be non `NULL`. If both parameters are given,
#'              then the parameter `s` is ignored.
#' @param D     \eqn{(m,n)} dimensional matrix (or `NULL`). See the details below.
#'
#' @return [stspmod()] object, which represents the generated model.
#' @export
#'
#' @examples
#' ### random state space model with two outputs and state dimension s = 3
#' ### The model is required to be stable and minimum phase
#' model = try(test_stspmod(dim = c(2,2), s = 3, digits = 2, bpoles = 1, bzeroes = 1))
#' if (!inherits(model, 'try-error')) {
#'    print(model)
#'    print(min(abs(poles(model$sys))) > 1)
#'    print(min(abs(zeroes(model$sys))) > 1)
#'    print(pseries2nu(pseries(model$sys, lag.max = 10))) # Kronecker indices
#' }
#'
#' ### random state space model with three outputs and 2 inputs in echelon canonical form
#' ### D is lower triangular (with ones on the diagonal)
#' ### the model is required to stable (the transfer function has no poles within the unit circle)
#' model = try(test_stspmod(dim = c(3, 2), nu = c(2,3,0),
#'                          D = matrix(c(1,NA,NA,0,1,NA), nrow = 3, ncol = 2),
#'                          digits = 2, bpoles = 1))
#'
#' if (!inherits(model, 'try-error')) {
#'    print(model)
#'    print(min(abs(poles(model$sys))) > 1)
#'    print(pseries2nu(pseries(model$sys, lag.max = 10))) # Kronecker indices
#' }
test_stspmod = function(dim = c(1,1), s = NULL, nu = NULL, D = NULL, sigma_L = NULL,
                        digits = NULL, bpoles = NULL, bzeroes = NULL, n.trials = 100) {


  # check input parameter "dim"
  dim = as.integer(dim) # note: as.integer converts to vector!
  if ((length(dim) != 2) || (min(dim) < 0)) {
    stop('argument "dim" must be a vector of non negative integers with length 2!')
  }
  m = dim[1]
  n = dim[2]

  # check input parameter "D"
  if (!is.null(D)) {
    if ( (!is.numeric(D)) || (!is.matrix(D)) || any(dim(D) != dim) ) {
      stop('"D" must be a compatible, numeric matrix')
    }
  }
  else {
    D = diag(x = 1, nrow = m, ncol = n)
  }

  # check input parameter "sigma_L"
  if (!is.null(sigma_L)) {
    if ( (!is.numeric(sigma_L)) || (!is.matrix(sigma_L)) || any(dim(sigma_L) != n) ) {
      stop('"sigma_L" must be a compatible, numeric matrix')
    }
  }
  else {
    sigma_L = matrix(stats::rnorm(n*n), nrow = n, ncol = n)
    if (n >0) {
      sigma_L = t(chol(sigma_L %*% t(sigma_L)))
    }
    if (!is.null(digits)) sigma_L = round(sigma_L, digits)
  }

  # check input parameter "nu"
  if (!is.null(nu)) {
    nu = as.integer(nu) # as.integer converts to vector
    if ((length(nu) != m) || (min(nu) < 0)) {
      stop('"nu" must be a vector of non negative integers with length equal to dim[1]!')
    }
  } else {
    if (is.null(s)) stop('either "s" or "nu" must be specified')
    s = as.integer(s)[1]
    if (s < 0) stop('parameter "s" must be a nonnegative integer')
  }

  sys = try(test_stsp(dim = dim, s = s, nu = nu, D = D, digits = digits,
                      bpoles = bpoles, bzeroes = bzeroes, n.trials = n.trials))
  if (inherits(sys, 'try-error')) {
    stop('Could not generate a suitable state space model with ', n.trials, ' trials!')
  }

  model = stspmod(sys = sys, sigma_L = sigma_L)
  return(model)
}

#' Generate a Random Model
#'
#' This function may be used to generate random state space or VARMA models.
#' The main argument is a model template, which defines the type of model to generate,
#' see e.g. [model structures()].
#' If bounds for the poles and/or the zeroes are given,
#' then the procedure simply generates random models until a model which satisfies the constraint is found.
#' Of course this is a very crude method and it may need a large number of randomly generated model.
#'
#' @param template A model template as computed e.g. by [model2template()].
#' @param ntrials.max Maximum number of trials.
#' @param bpoles,bzeroes Lower bounds on the poles and zeroes of the model
#'        (such that stability and invertibility assumptions are satisfied).
#'        If set to `NA`, then the corresponding test is skipped.
#' @param rand.gen (optional) A function to generate the random, "free" parameters.
#' @param ... Additional parameters, passed on to `rand.gen`.
#'   In particular, if the "free" paramameters are generated by [rnorm()],
#'   then the standard deviation `sd` may be set.
#'   Choosing small values for `sd` makes it easier to find a stable and miniphase model.
#'   Of course this "trick" only works if the reference model, which is obtained
#'   with a zero parameter vector, satisfies the constraints.
#'
#' @return Model object whose class depends on the template.
#'
#' @export
#'
#' @examples
#' # Generate a random VARMA model in echelon form ############
#'
#' # Compute the appropriate model template
#' tmpl = tmpl_arma_echelon(nu = c(1,2,1))
#'
#' # Create a random model, which is stable but not necessarily miniphase
#' model = r_model(tmpl, bpoles = 1, sd = 0.5)
#' model
#'
#' # Check whether the poles satisfy the constraint
#' min(abs(poles(model)))
#' min(abs(zeroes(model)))
r_model = function(template, ntrials.max = 100, bpoles = NULL, bzeroes = NULL,
                   rand.gen = stats::rnorm, ...) {

  # Initialize variable
  n.par = template$n.par
  constraint_satisfied = FALSE
  ntrials = 0

  # take care of the non-square case
  if (template$order[1] != template$order[2]) bzeroes = NULL

  # Generate random models until the constraints are satisfied
  while ( (!constraint_satisfied) && (ntrials < ntrials.max) ) {
    ntrials = ntrials + 1
    theta = rand.gen(n.par, ...)
    model = fill_template(theta, template)

    # Check constraints on poles and zeros
    constraint_satisfied = TRUE
    if (!is.null(bpoles)) {
      poles = poles(model, print_message = FALSE)
      if (length(poles) > 0) {
        constraint_satisfied = (min(abs(poles)) > bpoles)
      }
    }
    if ( constraint_satisfied && (!is.null(bzeroes)) ) {
      zeroes = zeroes(model, print_message = FALSE)
      if (length(zeroes) > 0) {
        constraint_satisfied = (min(abs(zeroes)) > bzeroes)
      }
    }
  }

  # Throw error if constraint is not satisfied after a certain number of trials
  if (!constraint_satisfied){
    stop('Could not generate a suitable model with ', ntrials, ' trials!')
  }

  return(model)
}



# Transformation between Different Rational Models ----

#' Coerce to State Space Model
#'
#' The function [as.stsp.pseries()] calls
#' [pseries2stsp()] with default parameters.
#' Of course the [pseries()] object must contain
#' sufficiently many lags. NOT YET implemented
#'
#'
#' @param obj object
#' @param method character string
#' @param ... optional additional parameters
#'
#' @return object of class [stspmod()].
#'
#' @rdname as.stspmod
#' @export
as.stspmod = function(obj, ...){
  UseMethod("as.stspmod", obj)
}


#' @rdname as.stspmod
#' @export
as.stspmod.armamod = function(obj, ...){
  obj = unclass(obj)
  obj$sys = as.stsp(obj$sys)
  class(obj) = c('stspmod', ' rldm')
  return(obj)
}

## Transformation of State Space System to Innovation Form ----

#' Innovation Form state space Model
#'
#' `r lifecycle::badge("experimental")`
#' Convert a given state space model into innovation form, i.e. the transformed model satisfies
#' \itemize{
#' \item \eqn{D=I}
#' \item the model is stable and minimum phase.
#' }
#' The procedure "flips" bad poles and zeroes by the helper functions [reflect_zeroes()] and [reflect_poles()].
#' The transformed model is an equivalent description of the process in terms of second order moments.
#' This means that the spectral density is not changed.
#'
#' @param model A state space model, i.e. an object of type [stspmod()].
#' @param echelon_form boolean, if `TRUE` the innovation form model is in addition
#'                    transformed to echelon form.
#' @param y `NULL` or a data sample. `as.matrix(y)` should return an
#'          (N,m)-dimensional numeric matrix. If not `NULL` the noise covariance matrix is
#'          estimated from this sample.
#'
#' @return List with slots
#' \item{model}{the state space model in innovation form.}
#' \item{z}{(complex) vector with the zeroes of the innovation form model.}
#' \item{z_flipped}{(boolean) vector which indicates which zeroes have been flipped.}
#' \item{p}{(complex) vector with the poles of the innovation form model.}
#' \item{p_flipped}{(boolean) vector which indicates which poles have been flipped.}
#'
#' @export
#'
#' @examples
#' # in order to get reproducable results
#' set.seed(342)
#'
#' model = r_model(tmpl_stsp_full(3, 3, 5))
#' print(model, digits = 4)
#' # the model has two non-minimum phase zeroes and two non-stable poles.
#' z = zeroes(model, print_message = FALSE)
#' abs(z)
#' p = poles(model, print_message = FALSE)
#' abs(p)
#'
#' # convert to innnovation form, by flipping the "bad" poles and zeroes.
#' out = innovation_form(model)
#' print(out$model, digits = 4)
#' flip = function(x) {
#'   x[abs(x) < 1] = 1/x[abs(x) < 1]
#'   return(x)}
#' data.frame(poles.inno = out$p, flipped = out$p_flipped,
#'            poles.ori = p[match_vectors(out$p, p)],
#'            zeroes.inno = out$z, flipped = out$z_flipped,
#'            zeroes.ori = z[match_vectors(out$z, flip(z))])
#'
#' # check that the innovation form model describes the same process,
#' # by checking that the spectral density is not changed!
#' junk = spectrald(model, n.f = 128)
#' junk1 = spectrald(out$model, n.f = 128)
#' all.equal(junk, junk1)
#'
#' # reset seed
#' set.seed(NULL)
innovation_form = function(model, echelon_form = TRUE, y = NULL) {
  if (!inherits(model, 'stspmod')) stop('only implemented for state space models')

  sys = model$sys
  sigma_L = model$sigma_L

  d = dim(sys)
  m = d[1]
  if ((m <1) || (d[2] != m)) stop('The system must be non-empty and square (m = n > 0)!')
  s = d[3]

  if (!is.null(y)) {
    y = try(as.matrix(y))
    if ( inherits(y, 'try-error') || (!is.numeric(y)) || (!is.matrix(y)) ) {
      stop('parameter "y" must be a data object which may be coerced to a matrix with "as.matrix(y)')
    }
    if (ncol(y) != m) stop('the data "y" is not compatible to the model!')
    n.obs = nrow(y)
    if (n.obs < m) stop("sample is too small for estimating the noise covariance matrix (N < m)")
    y = t(y)
  } else {
    n.obs = 0
  }

  if (s == 0) {
    # "almost" nothing to do
    D = sys$D
    sys = sys %r% solve(D)
    if (n.obs == 0) {
      sigma_L = D %*% sigma_L
      sigma_L = t(chol(tcrossprod(sigma_L)))
    } else {
      sigma_L = t(chol(crossprod(y) / n.obs))
    }
    model = stspmod(sys = sys, sigma_L = sigma_L)
    return(list(model = model, z = complex(length.out = 0L), z_flipped = logical(0L),
                p = complex(length.out = 0L), p_flipped = logical(0L)))
  }

  # check the zeroes of the system  ############
  # are there any non-minimum phase zeroes?
  inv_sys = sys^{-1} # inverse system

  # compute Schur decomposition of inv_sys$A, arrange stable eigenvalues at the top
  out = schur(inv_sys$A, select = 'iuc')
  z = out$lambda
  z_flipped = (abs(z) > 1)  # which zeroes have to be flipped
  # check: the first k eigenvalues are stable and the remaining s-k eigenvalues are unstable
  z_stable = logical(s)
  z_stable[iseq(1,out$k)] = TRUE
  if (any(z_flipped != (!z_stable))) stop('schur decomposition of (A-BD^{-1}C) failed', out$lambda)
  s_unstable = s - out$k
  s_stable = out$k

  if (s_unstable > 0) {
    # the system is not minimum phase! flip the "bad" zeroes

    # determine the noise covariance matrix sigma
    if (n.obs > 0) {
      # transform inverse system, such that the state transition matrix is
      # block upper triangular with A11 stable and A22 is unstable:
      if ((s_unstable > 0) && (s_stable > 0)) {
        # apply state transformation
        inv_sys = stsp(A = out$S, B = t(out$U) %*% inv_sys$B,
                       C = inv_sys$C %*% out$S, D = inv_sys$D)
      }

      as = matrix(0, nrow = s_stable, ncol = n.obs + 1)
      au = matrix(0, nrow = s_unstable, ncol = n.obs + 1)
      u = matrix(0, nrow = m, ncol = n.obs)

      # fbsolve_STSP_cpp(inv_sys$A, inv_sys$B, inv_sys$C, inv_sys$D,
      #                  y, au, as, u)
      .Call(`_RLDM_fbsolve_STSP_cpp`,
            inv_sys$A, inv_sys$B, inv_sys$C, inv_sys$D,
            y, au, as, u)
      sigma_L = t(chol(tcrossprod(u) / n.obs))
    }
    # if no data is given we use the noise covariance of the model (model$sigma_L)

    # normalize noise covariance matrix sigma = I
    sys = sys %r% sigma_L
    sigma_L = diag(m)

    # flip non miniphase zeroes
    sys = reflect_zeroes(sys, 1/z[z_flipped])
    z[z_flipped] = 1/z[z_flipped]
  }

  # check the poles of the system  ############
  # are there any non-stable poles?

  # compute Schur decomposition of sys$A, arrange stable eigenvalues at the top
  out = schur(sys$A, select = 'iuc')
  p = out$lambda
  p_flipped = (abs(p) > 1)  # which zeroes have to be flipped
  # check: the first k eigenvalues are stable and the remaining s-k eigenvalues are unstable
  p_stable = logical(s)
  p_stable[iseq(1,out$k)] = TRUE
  if (any(p_flipped != (!p_stable))) stop('schur decomposition of (A) failed', out$lambda)
  s_unstable = s - out$k
  s_stable = out$k

  if (s_unstable > 0) {
    # the system is not stable! flip the "bad" poles

    # determine the noise covariance matrix sigma
    if (n.obs > 0) {
      # the system should be minimum phase, thus we use
      # residuals_STSP_cpp

      a = matrix(0, nrow = s, ncol = n.obs + 1)
      u = matrix(0, nrow = m, ncol = n.obs)

      # residuals_STSP_cpp(sys$A, sys$B, sys$C, sys$D, y, a, u, diag(0), diag(0))
      .Call(`_RLDM_residuals_STSP_cpp`, sys$A, sys$B, sys$C, sys$D, y, a, u, diag(0), diag(0))
      sigma_L = t(chol(tcrossprod(u) / n.obs))
    }
    # if no data is given we use the noise covariance of the model (model$sigma_L)

    # normalize noise covariance matrix sigma = I
    sys = sys %r% sigma_L
    sigma_L = diag(m)

    # flip non stable poles
    sys = reflect_poles(sys, 1/p[p_flipped])
    p[p_flipped] = 1/p[p_flipped]
  }

  # determine the noise covariance matrix sigma ####
  if (n.obs > 0) {
    # the system should be minimum phase, thus we use
    # residuals_STSP_cpp

    a = matrix(0, nrow = s, ncol = n.obs + 1)
    u = matrix(0, nrow = m, ncol = n.obs)

    # residuals_STSP_cpp(sys$A, sys$B, sys$C, sys$D, y, a, u, diag(0), diag(0))
    .Call(`_RLDM_residuals_STSP_cpp`, sys$A, sys$B, sys$C, sys$D, y, a, u, diag(0), diag(0))
    sigma_L = t(chol(tcrossprod(u) / n.obs))
  }
  # if no data is given we use the noise covariance of the model (model$sigma_L)

  # set D = I
  D = sys$D
  sigma_L = t(chol(tcrossprod(D %*% sigma_L)))
  sys = sys %r% solve(D)

  if (echelon_form) {
    O = obs_matrix(sys$A, sys$C)
    sys = state_trafo(sys, O[1:s, ,drop = FALSE], inverse = FALSE)
  }

  model = stspmod(sys = sys, sigma_L = sigma_L)

  return(list(model = model, z = 1/z, z_flipped = z_flipped,
              p = 1/p, p_flipped = p_flipped))
}

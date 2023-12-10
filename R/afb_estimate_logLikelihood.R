# 1) Likelihood Evaluation for a Model ----

#' Log Likelihood Methods
#'
#' Tools and methods for the computation of the (conditional or exact) Gaussian log-likelihood
#' of ARMA, RMFD, and state space models.
#' For functions which serve as input for optimizers like [optim][stats::optim],
#' see [ll_theta] and [ll_FUN] (where the latter is a *function factory* which
#' generates a closure that serves as such an input).
#'
#' The procedure three choices ...
#'
#' For an ARMA model
#' \deqn{a_0 y_t + a_1 y_{t-1} + \cdots + a_p y_{t-p} = b_0 u_t + b_1 u_{t-1} + \cdots + b_q u_{t-q}}{
#'      a[0] y[t] + a[1] y[t-1] + \dots + a[p] y[t-p] = b[0] u[t] + b[1] u[t-1] + \dots + b[q] u[t-q]}
#' with Gaussian noise \eqn{u_t \sim N(0,\Sigma)}{u[t] ~ N(0,\Sigma)} an *approximation* of the
#' **scaled** log likelihood is
#' \deqn{ll = -(1/2)(m \ln(2\pi) + \mathrm{tr}(S\Sigma^{-1}) + \ln\det \Sigma + 2 \ln\det (a_0^{-1}b_0)}{
#'      ll = -(1/2)(m ln(2\pi) + tr(S\Sigma^{-1}) + ln det \Sigma + 2 ln det (a[0]^{-1}b[0])}
#' where \eqn{S} denotes the sample covariance of the residuals of the model
#' \deqn{S=\frac{1}{N-s}\sum_{t=s+1}^N e_t e_t'}{
#'       S=(1/(N-s)) \sum_{t=s+1}^N e[t] e[t]'}
#' The residuals are computed from a sample \eqn{y_t, t=1,\ldots,N}{y[t], t=1,...,N} by solving the (inverse) ARMA system
#' \deqn{b_0 e_t = -b_1 e_{t-1} - \cdots - b_q e_{t-q} + a_0 y_t + a_1 y_{t-1} + \cdots + a_p y_{t-p}}{
#'       b[0] e[t] = -b[1] e[t-1] - \dots -b[q] e[t-q] + a[0] y[t] + a[1] y[t-1] + \dots + a[p] y[t-p]}
#' and setting all unknown initial values \eqn{y_t=0}{y[t]=0} and \eqn{e_t=0}{e[t]=0} for
#' \eqn{t\leq 0}{t\le 0} equal to zero.
#' See e.g. [solve_inverse_de()].
#' Note that the log Likelihood here is **scaled** by a factor \eqn{1/(N-s)} and
#' that the first \eqn{s} observations are **skipped** when computing the sample covariance matrix.
#'
#' The log-likelihood may be easily maximized with respect to the noise covariance matrix \eqn{\Sigma}.
#' For given \eqn{S}, the optimal value for \eqn{\Sigma} is \eqn{\Sigma=S}.
#' If we plug this maximizer into the log-likelihood function, we obtain the "concentrated" log likelihood function
#' \deqn{cll = -(1/2)(m \ln(2\pi) + m + \ln\det S + 2 \ln\det (a_0^{-1}b_0)}{
#'       cll = -(1/2)(m ln(2\pi) + m + ln det S + 2 ln det (a[0]^{-1}b[0])}
#' which only depends on the sample `y` and the ARMA parameter matrices \eqn{a_i}{a[i]} and \eqn{b_i}{b[i]}.
#'
#' For state space models the (approximate) log likelihood is computed quite analogously.
#'
#' @note
#' To be precise, the functions returns \eqn{1/(N-s)} times the (approximate) log likelihood.
#'
#' The above routines only handle the case of centered data,
#' i.e. it is assumed that the output process \eqn{(y_t)}{(y[t])} has mean zero!
#'
#' The computation of the concentrated log likelihood assumes
#' that the model structure does **not** impose restrictions on the noise covariance matrix.
#'
#' @param obj Object of type [armamod()], [rmfdmod()], or [stspmod()].
#' @param y Data sample given as \eqn{(N,m)} dimensional matrix, or
#'          a "time series" object (in the sense that `as.matrix(y)`
#'          should return an \eqn{(N,m)}-dimensional numeric matrix).
#'          Missing values (`NA`, `NaN` and `Inf`) are **not** supported.
#' @param which (character) which likelihood to compute.
#' @param skip (integer) skip initial observations. This parameter is only used, when the
#'           (concentrated) conditional likelihood is computed.
#' @param P1 \eqn{(s,s)} dimensional covariance matrix of the error of the initial state estimate.
#'           If `NULL` then the state covariance \eqn{P=APA'+B\Sigma B'} is used.
#'           This parameter is only used, when the (exact) likelihood is
#'           computed via the Kalman Filter. See [ll_kf()] for more details.
#' @param a1 \eqn{s} dimensional vector, which holds the initial estimate for the state at time t=1.
#'           If `a1=NULL`, then a zero vector is used. This parameter is only used, when the
#'           (exact) likelihood is computed via the Kalman Filter. See [ll_kf()] for more details.
#' @param tol (small) tolerance value (or zero) used by the Kalman Filter routines, see [ll_kf()].
#' @param ... Not used.
#'
#' @return (double) the (scaled) log Likelihood of the model.
#'
#' @name ll
#' @export
#'
#' @examples
#' # Generate a random model in echelon form model (m = 3)
#' tmpl = tmpl_arma_echelon(nu = c(2,1,1))
#' model = r_model(template = tmpl, bpoles = 1, bzeroes = 1, sd = 0.25)
#' diag(model$sigma_L) = 1 # scale the diagonal entries of sigma_L
#' print(model)
#' # extract the corresponding free/deep parameters
#' th = extract_theta(model, tmpl)
#'
#' # generate a sample with 50 observations
#' y = sim(model, n.obs = 50, n.burn_in = 100)$y
#'
#' # conditional log likelihood
#' # the following statements return the same ll value!
#' ll(model, y, which = 'conditional', skip = 0)
#' ll_theta(th, template= tmpl, y, which = 'conditional', skip = 0)
#' llfun = ll_FUN(tmpl, y, which = 'conditional', skip = 0)
#' llfun(th)
#'
#' # concentrated, conditional log likelihood
#' # the following statements return the same ll value!
#' ll(model, y, which = 'concentrated', skip = 0)
#' ll_theta(th, template= tmpl, y, which = 'concentrated', skip = 0)
#' llfun = ll_FUN(tmpl, y, which = 'concentrated', skip = 0)
#' llfun(th)
ll = function(obj, y, which, ...) {
  UseMethod("ll")
}




#' @name ll
#' @export
ll.armamod = function(obj, y, which = c('concentrated', 'conditional'), skip = 0L, ...) {

  # Retrieve lmfd object from armamod input ####
  sys = obj$sys

  # Retrieve integer-valued parameters ####
  order = unname(dim(sys))
  m = order[1]
  n = order[2]
  p = order[3]
  q = order[4]

  # Check dimensions of inputs ####
  if ( (m < 1) || (m != n) ) {
    # print(c(order, m,n,p,q))
    stop('only square, non empty ARMA models are supported')
  }

  # Check data input ####
  # coerce data object(s) to matrices
  y = try(as.matrix(y))
  if ( inherits(y, 'try-error') || (!is.numeric(y)) || (!is.matrix(y)) ) {
    stop('input "y" must be a data object which may be coerced to a matrix with "as.matrix(y)')
  }
  if (length(y) == 0) stop('"y" contains no data')
  n.obs = nrow(y)  # sample size
  if (ncol(y) != m) stop('"y" is not compatible with the model!')
  if (any(!is.finite(y))) stop('"y" contains NAs/NaNs/Infs!')

  # if (is.null(skip)) skip = max(p,q)
  skip = as.integer(skip)[1]
  if ( (skip < 0)) stop('"skip" must be a non negative integer.')
  n.valid = n.obs - skip
  if (n.valid < m) stop("too few observations, (n.obs-skip) < m!")

  # take care of the case where the lag zero coefficient k[0] ####
  # of the impulse response is not equal to the identity!
  k0 = matrix(unclass(rationalmatrices::pseries(sys, lag.max = 0)), nrow = m, ncol = m)
  lndetk0 = log(abs(det(k0)))

  # Obtain residuals for given system parameters ####
  u = solve_inverse_de(sys, y)$u

  # sample covariance matrix of residuals ####
  S = ( crossprod(u[(skip+1):n.obs, , drop = FALSE]) ) / n.valid

  if (which == 'concentrated') {
    lndetSigma = log(det(S))
    trSS = m
  } else {
    sigma = obj$sigma_L
    sigma = sigma %*% t(sigma)
    lndetSigma = log( abs( det(sigma ) ))  # abs() should not be necessary, but maybe helps in the case that
                                           # sigma is close to singular
    trSS = sum( diag( solve(sigma, S) ) )
  }
  mln2pi = m*log(2*pi)
  ll = (-1 / 2) * (mln2pi + trSS + lndetSigma + 2*lndetk0)

  return(ll)
}

#' @name ll
#' @export
ll.stspmod = function(obj, y, which = c('concentrated', 'conditional', 'kf', 'kf2'),
                      skip = 0L, P1 = NULL, a1 = NULL, tol = 1e-8, ...) {

  which = match.arg(which)

  sys = obj$sys
  order = unname(dim(sys))
  m = order[1]
  n = order[2]
  s = order[3]

  if ( (m < 1) || (m != n) ) {
    # print(c(order,m,n,s))
    stop('only square, non empty state space models are supported')
  }

  # Check data
  # coerce data object(s) to matrices
  y = try(as.matrix(y))
  if ( inherits(y, 'try-error') || (!is.numeric(y)) || (!is.matrix(y)) ) {
    stop('input "y" must be a data object which may be coerced to a matrix with "as.matrix(y)')
  }
  if (length(y) == 0) stop('"y" contains no data')
  n.obs = nrow(y)  # sample size
  if (ncol(y) != m) stop('"y" is not compatible with the model!')
  if (any(!is.finite(y))) stop('"y" contains NAs/NaNs/Infs!')

  if ( (which == 'conditional') || (which == 'concentrated') ) {
    # if (is.null(skip)) skip = max(p,q)
    skip = as.integer(skip)[1]
    if ( (skip < 0)) stop('"skip" must be a non negative integer.')
    n.valid = n.obs - skip
    if (n.valid < m) stop("too few observations, (n.obs-skip) < m!")

    # take care of the case where the lag zero coefficient k[0] ####
    # of the impulse response is not equal to the identity!
    k0 = sys$D
    lndetk0 = log(abs(det(k0)))

    # Obtain residuals for given system parameters ####
    u = solve_inverse_de(sys, y)$u

    # sample covariance matrix of residuals ####
    S = ( crossprod(u[(skip+1):n.obs, , drop = FALSE]) ) / n.valid

    if (which == 'concentrated') {
      lndetSigma = log(det(S))
      trSS = m
    } else {
      sigma = obj$sigma_L
      sigma = sigma %*% t(sigma)
      lndetSigma = log( abs( det(sigma ) ))  # abs() should not be necessary, but maybe helps in the case that
      # sigma is close to singular
      trSS = sum( diag( solve(sigma, S) ) )
    }
    mln2pi = m*log(2*pi)
    ll = (-1 / 2) * (mln2pi + trSS + lndetSigma + 2*lndetk0)

    return(ll)
  }
  # (concentrated) conditional log Likelihood

  # exact log Likelihood via Kalman Filter or square root Kalman filter
  if ( (which == 'kf') || (which == 'kf2') ) {

    # check 'P1'
    # note we do not check that P1 is symmetric and psd.
    if (is.null(P1)) {
      if (s > 0) {
        P1 = lyapunov(obj$sys$A, tcrossprod(obj$sys$B %*% obj$sigma_L),
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

    # check 'a1'
    if (is.null(a1)) a1 = double(s)
    a1 = as.vector(a1)
    if (length(a1) != s) stop('parameter "a1" is not compatible with the model')

    sigma = tcrossprod(obj$sigma_L)
    if (which == 'kf') {
      # call 'Kalman filter'
      R = obj$sys$D %*% sigma %*% t(obj$sys$D)
      S = obj$sys$B %*% sigma %*% t(obj$sys$D)
      Q = obj$sys$B %*% sigma %*% t(obj$sys$B)
      # return(list(obj$sys$A, obj$sys$C, Q, R, S, t(y), P1, a1))
      out = .Call(`_RLDM_ll_kf_cpp`, obj$sys$A, obj$sys$C, Q, R, S, t(y), P1, a1, tol)
      return( out )
    } else {
      # call 'square root Kalman filter'
      # H_t = sigma_R %*% cbind(t(obj$sys$D), t(obj$sys$B))
      H = rbind(obj$sys$D, obj$sys$B) %*% obj$sigma_L
      if (s >0) {
        P1_R = try(chol(P1), silent = TRUE)
        if (inherits(P1_R,'try-error')) {
          # cholesky decomposition failed, P1 is singular?
          # try eigenvalue decomposition
          evd = eigen(P1, symmetric = TRUE)
          if (min(evd$values) < -sqrt(.Machine$double.eps)) {
            stop('"P1" is not positive semidefinite!')
          }
          P1_R = diag(sqrt(pmax(evd$values,0)), nrow = s, ncol = s) %*% t(evd$vectors)
        }
      } else {
        P1_R = matrix(0, nrow = 0, ncol = 0)
      }
      ll = .Call(`_RLDM_ll_kf2_cpp`, obj$sys$A, obj$sys$C, t(H), t(y), P1_R, a1, tol)
      return( ll )
    }
  }
  # exact log Likelihood via Kalman Filter or square root Kalman filter

  stop('this should not happen.')
}


## 1.1) Same but Deprecated ----

#' Maximum Likelihood Estimation
#'
#' `r lifecycle::badge("superseded")`
#' This is a naive implementation of Maximum Likelihood Estimation.
#' Rather use [ll()] and [ll_FUN()].
#'
#' @note
#' * The optimization is computed with the general-purpose routine [stats::optim()].
#' * An initial estimate is **needed**.
#' * The procedure does **not** respect constraints like stability or minimum phase.
#' * The case of the *conditional, concentrated* likelihood is somewhat special.
#'   In this case the model template must have a particular structure: (1) The noise covariance
#'   is parametrized via the left cholesky factor. (2) The last \eqn{m(m+1)/2} components
#'   of the parameter vector \eqn{\theta} parametrize this left cholesky factor and the other components
#'   describe the system. (This implies that there is no overlap/dependency betweeen the "system parameters"
#'   and the "noise parameters".)
#'
#' @param y sample, i.e. an \eqn{(N,m)} dimensional matrix,
#'   or a "time series" object (i.e. `as.matrix(y)` should return an
#'   \eqn{(N,m)}-dimensional numeric matrix). Missing values (`NA`, `NaN` and
#'   `Inf`) are **not** supported.
#' @param tmpl a model template which describes the model class, see [model structures()].
#'   Note that only the case of (non-empty, square) state space or ARMA models is implemented.
#' @param th Initial parameter estimate.
#' @param which (character string) determines which "likelihood" to be used, see also [ll()].
#'        The option `"kf"`  is only supported for state space models.
#' @param method,hessian,control are passed on to the optimization routine [stats::optim()].
#'
#' @return A list with components
#' \item{model}{The estimated model.}
#' \item{th}{The corresponding vector of *deep* parameters.}
#' \item{ll}{The log likelihood of the estimated model.}
#' \item{which}{The type of likelihood used.}
#' \item{counts, convergence, message, hessian}{as returned by  [stats::optim()].}
#'
#' @export
#'
#' @examples
#' # Generate a random model in echelon form model (m = 3)
#' tmpl = tmpl_stsp_echelon(nu = c(2,1,1))
#' model = r_model(template = tmpl, bpoles = 1, bzeroes = 1, sd = 0.25)
#' diag(model$sigma_L) = 1 # scale the diagonal entries of sigma_L
#' print(model)
#' # extract the corresponding free/deep parameters
#' th = extract_theta(model, tmpl)
#'
#' # generate a sample with 500 observations
#' y = sim(model, n.obs = 500, n.burn_in = 100)$y
#'
#' # We are cheating here and use the true model parameters
#' # as starting values for the optimization routine:
#'
#' # estimate the model with the "exakt log likelihood"
#' out = est_ML(y, tmpl, th, which = 'kf')
#' KL_divergence(model, out$model)
#'
#' # estimate the model with "conditional log likelihood"
#' out = est_ML(y, tmpl, th, which = 'conditional')
#' KL_divergence(model, out$model)
#'
#' # estimate the model with "concentrated, conditional log likelihood"
#' out = est_ML(y, tmpl, th, which = 'concentrated')
#' KL_divergence(model, out$model)
est_ML = function(y, tmpl, th, which = c('concentrated', 'conditional', 'kf'),
                  method = c("BFGS", "Nelder-Mead", "CG", "L-BFGS-B", "SANN","Brent"),
                  hessian = FALSE, control = list()) {


  # coerce data object to matrix
  y = try(as.matrix(y))
  if ( inherits(y, 'try-error') || (!is.numeric(y)) || (!is.matrix(y)) ) {
    stop('input "y" must be a data object which may be coerced to a matrix with "as.matrix(y)')
  }
  if (length(y) == 0) stop('"y" contains no data')
  n.obs = nrow(y)  # sample size
  if (any(!is.finite(y))) stop('"y" contains NAs/NaNs/Infs!')

  # check template
  if (!((tmpl$class == 'stspmod') || (tmpl$class == 'armamod'))) {
    stop('Only state space models and ARMA models are implemented!')
  }
  if (tmpl$n.par == 0) stop('the template has no free/deep parameters!')
  order = unname(tmpl$order)
  if ( (order[1] != order[2]) || (order[1] < 1) )  {
    stop('Only square, non-empty models (m = n  > 0) are implemented!')
  }
  m = order[1]
  if (ncol(y) != m) stop('"y" is not compatible with the template "tmpl"!')

  th = as.vector(th)
  if (length(th) != tmpl$n.par) {
    stop('The initial parameter estimate "th" is not compatible with the template "tmpl"!')
  }

  # optim method
  method = match.arg(method)

  # optim control parameters
  if (is.null(control$fnscale)) control$fnscale = -1
  # make sure that fnscale is negative!
  if (control$fnscale > 0) control$fnscale = -control$fnscale

  if (is.null(control$maxit)) control$maxit = 50 + length(th)^2
  if (is.null(control$trace)) control$trace = 0

  which = match.arg(which)
  if ((tmpl$class == "armamod") && (which == 'kf')) {
    stop('The option "which = \'kf\'" is only implemented for state space models.')
  }

  # evaluate ll at the initial parameter
  ll0 = try(ll_theta(th, tmpl, y, which))
  if (inherits(ll0, 'try-error') || (!is.finite(ll0))) {
    stop('Could not evaluate the log Likelihood at the initial parameter value.!')
  }

  if (which == 'concentrated') {
    # strip sigma_L parameters from tmpl
    # we assume cholesky structure, ....
    # no checks so far!!!
    tmpl0 = tmpl
    tmpl0$n.par = tmpl$n.par - m*(m+1)/2
    if (tmpl0$n.par <= 0) stop('the (stripped) template has no free/deep parameters!')
    tmpl0$H = tmpl$H[, 1:tmpl0$n.par, drop = FALSE]
    th0 = th[1: tmpl0$n.par]

    llfn = ll_FUN(tmpl0, y, 'concentrated', err = ll0)
    llgr = ll_FUN(tmpl0, y, 'gr_concentrated')

    out = stats::optim(par = th0, fn = llfn, gr = llgr,
                       method = method, hessian = hessian, control = control)

    th0 = out$par
    ll = llfn(th0)

    model = fill_template(th0, tmpl0)
    u = solve_inverse_de(model$sys, y)$u
    sigma_L = t(chol(crossprod(u) / n.obs))
    model$sigma_L = sigma_L
    th = extract_theta(model, tmpl, on_error = 'stop')

  } else {
    llfn = ll_FUN(tmpl, y, which, err = ll0)
    llgr = NULL

    out = stats::optim(par = th, fn = llfn, gr = llgr,
                       method = method, hessian = hessian, control = control)

    th = out$par
    ll = llfn(th)
    model = fill_template(th, tmpl)
  }

  return(list(model = model, th = th, ll = ll, which = which,
              counts = out$counts, convergence = out$convergence,
              message = out$message, hessian = out$hessian))
}



# 2) Likelihood Optimisation wrt Deep Parameters ----

#' Log-likelihood Given Deep Parameters
#'
#' `r lifecycle::badge("superseded")`
#' See [ll_FUN()].
#' The template `template` is filled with the deep parameters in `th`.
#' Subsequently, the S3 method [ll()] is called for the class provided in the template
#' and the value of the **scaled** log-likelihood function is returned, see [ll()].
#'
#' @param th Vector of deep parameter
#' @param template A model template, see [model structures].
#' @inheritParams ll
#'
#' @return Value of log-likelihood for a given deep/free parameter vector `th` and
#'         a model structure defined via `template`.
#'         Note that this function simply calls `ll(fill_template(th, template), y, which, ...)`.
#' @export
ll_theta = function(th, template, y, which, ...) {
  model = try(fill_template(th, template))
  if (inherits(model, 'try-error')) {
    stop('ll_theta(): the parameter vector *th* and the template are not compatible!')
  }
  return(ll(model, y, which = which, ...))
}


#' Log Likelihood Function Factory
#'
#' Creates a function similar to [ll_theta()] but faster and more memory efficient.
#' The model structure (`template`) and the data (`y`) are encoded within the
#' generated closure (a function plus its enclosing environment).
#' The generated function calls compiled C/C++ code (see [RcppArmadillo-package][RcppArmadillo::RcppArmadillo-package]) and
#' hence is much faster than calling `ll_theta(th, template, y, ...)`.
#'
#' @return A function, `llfun(th)` say, which computes the log-likelihood for given
#' *deep* parameters `th`. This function may be used for ML estimation of the model.
#'
#' @param template A model template, see [model structures].
#' @param y sample, i.e. an \eqn{(N,m)} dimensional matrix, or a "time series"
#'   object (i.e. `as.matrix(y)` should return an \eqn{(N,m)}-dimensional
#'   numeric matrix). Missing values (`NA`, `NaN` and `Inf`) are
#'   **not** supported.
#' @param which (string) Determines the type of ll function.
#' @param skip (integer) skip initial observations. If `NULL` then
#'     `skip` is set to \eqn{0} for state space models and to
#'     \eqn{\max(p,q)}{max(p,q)} for ARMA models. This parameter
#'     is only used for the cases "concentrated", "conditional" and "gr_concentrated"
#' @param tol (double) tolerance used by [ll_kf_cpp()].
#' @param err (double) return value for the case "kf", if the computation
#'        of the initial state covariance fails.
#'
#' @return
#' Function `fn(th)`
#' @export
#'
#' @examples
#' # Generate a random model in echelon form model (m = 3)
#' tmpl = tmpl_stsp_echelon(nu = c(2,1,1))
#' model = r_model(template = tmpl, bpoles = 1, bzeroes = 1, sd = 0.25)
#' diag(model$sigma_L) = 1 # scale the diagonal entries of sigma_L
#' print(model)
#' # extract the corresponding free/deep parameters
#' th = extract_theta(model, tmpl)
#'
#' # generate a sample with 50 observations
#' y = sim(model, n.obs = 50)$y
#'
#' # conditional log likelihood
#' # the following statements return the same ll value!
#' ll(model, y, 'conditional')
#' ll_theta(th, tmpl, y, 'conditional')
#' fn = ll_FUN(tmpl, y, 'conditional')
#' fn(th)
#'
#' # concentrated conditional log likelihood
#' # the following statements return the same ll value!
#' ll(model, y, 'concentrated')
#' ll_theta(th, tmpl, y, 'concentrated')
#' fn = ll_FUN(tmpl, y, 'concentrated')
#' fn(th)
#' # for this case, we may also compute the (analytic) gradient
#' gr = ll_FUN(tmpl, y, 'gr_concentrated')
#' gr(th)
#'
#' # log likelihood (via Kalman filter)
#' # the following statements return the same ll value!
#' ll(model, y, 'kf2')
#' ll_theta(th, tmpl, y, 'kf2')
#' ll(model, y, 'kf')
#' ll_theta(th, tmpl, y, 'kf')
#' fn = ll_FUN(tmpl, y, 'kf')
#' fn(th)
ll_FUN = function(template, y,
                  which = c('concentrated', 'conditional', 'kf', 'gr_concentrated'),
                  skip = 0L, tol = 1e-8, err = NA_real_) {

  which = match.arg(which)
  force(y)
  force(skip)
  force(tol)
  force(err)

  # we should check that "template" is a valid template!
  order = template$order
  m = order[1]
  n = order[2]
  if ( (m < 1) || (m != n) ) {
    # print(c(order, m, n))
    stop('only square, non empty systems are supported')
  }

  # coerce data object to matrix
  y = try(as.matrix(y))
  if ( inherits(y, 'try-error') || (!is.numeric(y)) || (!is.matrix(y)) ) {
    stop('input "y" must be a data object which may be coerced to a matrix with "as.matrix(y)')
  }
  if (length(y) == 0) stop('"y" contains no data')
  n.obs = nrow(y)  # sample size
  if (ncol(y) != m) stop('"y" is not compatible with the template!')
  if (any(!is.finite(y))) stop('"y" contains NAs/NaNs/Infs!')
  y = t(y) # transpose data matrix!


  # Case 1: state space models ----
  if (template$class == 'stspmod') {

    s = order[3]
    pi = matrix(0, nrow = s+m, ncol = s+m)
    H_pi = template$H[1:((m+s)^2), , drop = FALSE]
    h_pi = template$h[1:((m+s)^2)]
    L = matrix(0, nrow = m, ncol = m)
    H_L = template$H[((m+s)^2+1):((m+s)^2+m^2), , drop = FALSE]
    h_L = template$h[((m+s)^2+1):((m+s)^2+m^2)]

    if (which == 'kf') {
      VAR = matrix(0, nrow = s+m, ncol = s+m)
      P1 = matrix(0, nrow = s, ncol = s)
      rm(template)

      fun = function(th) {
        ll = .Call(`_RLDM_ll_kf_theta_cpp`,
                   th, y, pi, H_pi, h_pi, L, H_L, h_L,
                   VAR, P1, tol, err)
        # return(list(ll=ll, th=th, pi=pi, VER=VAR, P1=P1, L=L))
        return(ll)
      }

      return(fun)
    }

    a = matrix(0, nrow = s, ncol = n.obs + 1)
    u = matrix(0, nrow = m, ncol = n.obs)

    if ((which == 'concentrated') || (which == 'conditional')) {
      dU = matrix(0, nrow = 0, ncol = 0)
      concentrated = (which == 'concentrated')
      rm(template)

      fun = function(th) {
        ll = .Call(`_RLDM_cll_theta_STSP_cpp`,
                   th, y, skip, concentrated,
                   pi, H_pi, h_pi, L, H_L, h_L,
                   a, u, dU)
        return(ll)
      }

      return(fun)
    }

    # which = 'gr_concentrated'
    dU = matrix(0, nrow = m*n.obs, ncol = ncol(H_pi))
    ind_D = matrix(1:((m+s)^2), nrow = s+m, ncol = s+m)[iseq(s+1,s+m), iseq(s+1,s+m)]  # indices of "D"
    H_D = H_pi[ind_D, , drop = FALSE]
    concentrated = TRUE

    rm(template, ind_D)

    fun = function(th) {
      ll = .Call(`_RLDM_cll_theta_STSP_cpp`,
                 th, y, skip, concentrated,
                 pi, H_pi, h_pi, L, H_L, h_L,
                 a, u, dU)
      # cll_theta_STSP_cpp (with concentrated = TRUE) overwrites
      #    pi with the "system matrix" rbind((A,B),(C,D))
      #    L with the sample variance of the residuals
      #    u with the residuals
      #    dU with the corresponding Jacobian
      gr =  ( c(solve(L, u[, (skip+1):n.obs, drop = FALSE])) %*%
                dU[(skip*m+1):(n.obs*m), , drop = FALSE] ) * ((-1) / (n.obs-skip))

      # take the term log(det(k(0))) = log(det(D)) into account
      gr = gr - as.vector( solve(t( pi[(s+1):(s+m), (s+1):(s+m)]) ) ) %*% H_D
      return(as.vector(gr))
    }

    return(fun)
  }

  # Case 2: ARMA models ----
  if (template$class == "armamod") {

    p = order[3]
    q = order[4]
    Pi = array(1:((m^2)*(p+q+2)), dim = c(m,m,p+q+2))
    H = template$H
    h = template$h

    # b[0]
    ib0 = matrix(0, nrow = m, ncol = m)
    i = Pi[ , , p+2]
    H_b = H[i, , drop = FALSE]
    h_b = h[i]

    # -(b[q],...,b[1])
    B1 = matrix(0, nrow = m, ncol = m*q)
    if (q > 0) {
      i = Pi[, , (p+q+2):(p+3)]
    } else {
      i = integer(0)
    }
    H_B = -H[i, , drop = FALSE]
    h_B = -h[i]

    # (a[0],...,a[p])
    a0 = matrix(0, nrow = m, ncol = m)
    A = matrix(0, nrow = m, ncol = m*(p+1))
    i = Pi[, , 1:(p+1)]
    H_A = H[i, , drop = FALSE]
    h_A = h[i]

    # sigma_L = L
    L = matrix(0, nrow = m, ncol = m)
    H_L = H[((m^2)*(p+q+2) + 1):((m^2)*(p+q+3)), , drop = FALSE]
    h_L = h[((m^2)*(p+q+2) + 1):((m^2)*(p+q+3))]

    u = matrix(0, nrow = m, ncol = n.obs)

    if ((which == 'concentrated') || (which == 'conditional')) {
      dU = matrix(0, nrow = 0, ncol = 0)
      concentrated = (which == 'concentrated')
      rm(template, h, H, Pi, i)

      fun = function(th) {
        ll = .Call(`_RLDM_cll_theta_ARMA_cpp`,
                   th, y, skip, concentrated,
                   ib0, H_b, h_b,
                   B1, H_B, h_B,
                   a0, A, H_A, h_A,
                   L, H_L, h_L,
                   u, dU)
        return(ll)
      }

      return(fun)
    }

    # which = 'gr_concentrated'
    dU = matrix(0, nrow = m*n.obs, ncol = (m^2)*(p+q+2))
    H_pi = H[1:((m^2)*(p+q+2)), , drop = FALSE]
    H_a = H[1:(m^2), , drop = FALSE]
    # ind_D = matrix(1:((m+s)^2), nrow = s+m, ncol = s+m)[iseq(s+1,s+m), iseq(s+1,s+m)]  # indices of "D"
    concentrated = TRUE

    rm(template, H, h, Pi, i)

    fun = function(th) {
      ll = .Call(`_RLDM_cll_theta_ARMA_cpp`,
                 th, y, skip, concentrated,
                 ib0, H_b, h_b,
                 B1, H_B, h_B,
                 a0, A, H_A, h_A,
                 L, H_L, h_L,
                 u, dU)
      # cll_theta_ARMA_cpp (with concentrated = TRUE) overwrites
      #    L with the sample covariance of the residuals
      #    u with the residuals
      #    dU with the corresponding Jacobian
      #    ib0 with the inverse of b[0]
      #    a0 with a[0]

      gr =  ( as.vector(solve(L, u[, (skip+1):n.obs, drop = FALSE])) %*%
                dU[(skip*m+1):(n.obs*m), , drop = FALSE] ) * ((-1) / (n.obs-skip))
      gr = as.vector( gr %*% H_pi )

      # take the term log(det(k(0))) = log(det(a[0]^{-1}b[0])) into account
      gr = gr + as.vector( as.vector(t(solve(a0))) %*% H_a -
                           as.vector(t(ib0)) %*% H_b )
      return(gr)
    }

    return(fun)
  }

  stop('only "stspmod" and "armamod" templates are supported!')
}


# 3) Kalman related ----

#' Kalman Filter
#'
#' These functions implement the "standard" Kalman filter and the "square root" Kalman filter (also called "square root covariance filter") for time invariant, linear state space systems without exogenous inputs, see e.g. \insertCite{AndersonMoore2005}{RLDM}.
#'
#' The model considered is
#' \deqn{a_{t+1} = A a_t + Bu_t}{a[t+1] = A a[t] + B u[t]}
#' \deqn{y_t = C a_t + Du_t}{y[t] = C a[t] + D u[t]}
#' with \eqn{m}-dimensional outputs \eqn{y_t}{y[t]}, \eqn{s}-dimensional states
#' \eqn{a_t}{a[t]} and \eqn{n}-dimensional disturbances \eqn{u_t}{u[t]}.
#' The disturbances are white noise with a covariance matrix
#' \eqn{\mathbf{E} u_tu_t'=\Sigma}{E u[t]u[t]'=\Sigma}.
#' Note that the disturbances and the outputs may have *different* dimensions, however,
#' only "wide" systems with (\eqn{m\leq n}{m\le  n}) are implemented.
#'
#' The Kalman filter is a recursive scheme to compute the linear, least squares predictions
#' for \eqn{a_{t+1}}{a[t+1]} and \eqn{y_{t+1}}{y[t+1]} given the observations
#' \eqn{y_t,\ldots,y_1}{y[t],\ldots,y[1]} up to time \eqn{t}. These predictions are notated with
#' \eqn{a_{t+1|t}}{a[t+1|t]} and \eqn{y_{t+1|t}}{y_[t+1|t]}, the
#' prediction error for the output \eqn{y_{t+1}}{y[t+1]} is
#' \eqn{\epsilon_{t+1|t}=(y_{t+1}-y_{t+1|t})}{\epsilon[t+1|t]=(y[t+1]-y[t+1|t])}
#' and the corresponding variances of the prediction errors are
#' \deqn{\Pi_{t+1|t}=\mathbf{E}(a_{t+1}-a_{t+1|t})
#' (a_{t+1}-a_{t+1|t})',}{P[t+1|t]=E(a[t+1]-a_[t+1|t])(a[t+1]-a_[t+1|t])',}
#' \deqn{\Sigma_{t+1|t}=\mathbf{E}(\epsilon_{t+1|t}
#' \epsilon_{t+1|t}').}{\Sigma[t+1|t]=E(\epsilon_[t+1|t]\epsilon_[t+1|t]').}
#'
#' The standard form of the Kalman filter is based on the parameter matrices \eqn{A,C}, the variance of
#' "state disturbances"
#' \eqn{Q=\mathbf{E}(Bu_t (Bu_t)')=(B\Sigma B')}{Q=E(Bu[t](Bu[t])')=(B\Sigma B')}, the variance
#' of the "measurement disturbances"
#' \eqn{R=\mathbf{E}(Du_t (Du_t)')=(D\Sigma D')}{R=(Du[t](Du[t])')=E(D\Sigma D')} and the covariance
#' \eqn{S=\mathbf{E}(Bu_t(Du_t)')=(B\Sigma D')}{S=(Bu[t](Du[t])')=E(B\Sigma D')}.
#' Furthermore we need the initial prediction
#' \eqn{a_{1|0}}{a[1|0]} and the corresponding error variance
#' \eqn{\Pi_{1|0}}{P[1|0]}.
#'
#' For the square root form of the filter we need the "square roots"
#' \eqn{\Pi_{1|0}^{1/2}}{P[1|0]^{1/2}} and \eqn{\Sigma^{1/2}}, i.e. matrices such that
#' \eqn{\Pi_{1|0} = \Pi_{1|0}^{1/2} (\Pi_{1|0}^{1/2})'}{P[1|0] = P[1|0]^{1/2} (P[1|0]^{1/2})'}
#' and \eqn{\Sigma = \Sigma^{1/2}(\Sigma^{1/2})'}. In addition, we define
#' \eqn{H=(D',B')'\Sigma^{1/2}}.
#'
#' @details
#' The routines `kf_cpp`, `kf2_cpp` are \pkg{RcppArmadillo} implementations of the standard form and
#' of the square root form of the Kalman filter. The wrapper function `kf` takes an [stspmod()] object,
#' which describes the state space model, and then calls the approriate `RcppArmadillo` function.
#'
#' Square root Kalman filter: For the square root
#' \eqn{\Pi_{1|0}^{1/2}}{P[1|0]^{1/2}} the procedure first tries the Cholesky decomposition.
#' If this fails (since \eqn{\Pi_{1|0}^{1/2}}{P[1|0]^{1/2}} is (close to) singular),
#' then `ll_kf` tries to compute a symmetric square root via the eigenvalue decomposition
#' of \eqn{\Pi_{1|0}^{1/2}}{P[1|0]^{1/2}}.
#'
#' @param model [stspmod()] object, which represents the state space model.
#' @param y sample, i.e. an \eqn{(N,m)} dimensional matrix,
#'   or a "time series" object (i.e. `as.matrix(y)` should return an
#'   \eqn{(N,m)}-dimensional numeric matrix). Missing values (`NA`, `NaN` and
#'   `Inf`) are **not** supported.
#' @param y_t     \eqn{(m,N)} transposed data matrix `y_t = t(y)`.
#' @param method Character string. If `method="kf"` then `kf` calls
#'               `kf_cpp` ("standard form" of the Kalman filter) and for
#'               `method="kf2"` the "square root" form of the Kalman filter is used,
#'               i.e. `kf2_cpp` is called. Up to numerical errors the outputs should not
#'               depend on the chosen method.
#' @param P1    \eqn{(s,s)} dimensional covariance matrix of the error of the initial state estimate,
#'              i.e. \eqn{\Pi_{1|0}}{P[1|0]}.
#'              If `NULL`, then the state covariance \eqn{P = APA'+B\Sigma B'} is used.
#'              Note that this scheme assumes that the state space model is stable,
#'              i.e. that the state transition matrix \eqn{A} is stable.
#' @param P1_R  (right) square root of `P1`, i.e. `P1 = t(P1_R) \%*\% P1_R`.
#' @param a1    \eqn{s} dimensional vector, which holds the initial estimate \eqn{a_{1|0}}{a[1|0]}
#'              for the state at time \eqn{t=1}.  If `a1=NULL`, then a zero vector is used.
#' @param A     \eqn{(s,s)} dimensional state transition matrix \eqn{A}.
#' @param C     \eqn{(m,s)} dimensional matrix \eqn{C}.
#' @param Q,R,S The variance, covariance matrices of the "state disturbances" (\eqn{Bu_t}{Bu[t]})
#'              and the "measurement disturbances" (\eqn{Du_t}{Du[t]}) as described above.
#'              These matrices must be of dimension \eqn{(s,s)},  \eqn{(m,m)} and \eqn{(s,m)} respectively.
#' @param H_t    \eqn{(n,s+m)} dimensional matrix. This parameter corresponds to the transpose \eqn{H'} of
#'              \eqn{H=(D',B')'\Sigma^{1/2}}.
#'
#' @return List with components
#' \item{e}{\eqn{(N,m)} dimensional matrix with the standardized one-step
#'          ahead prediction errors. The \eqn{t}-th row of the matrix `e`
#'          corresponds to
#'          \deqn{e_t = \Sigma_{t|t-1}^{-1/2}\epsilon_{t|t-1}.}{e[t] = \Sigma[t|t-1]^{-1/2}\epsilon[t|t-1].}
#'          If the model is correctly specified then these standardized residuals are white noise
#'          with a unit covariance matrix. So they may be used for validaton of the model.}
#' \item{a}{\eqn{(N+1,s)} dimensional matrix with the estimated states. The
#'          \eqn{t}-th row of the matrix `a` corresponds to
#'          \eqn{a_{t|t-1}}{a[t|t-1]}. Given `y` and `a`, the one step ahead predictions
#'          \eqn{y_{t|t-1}}{y[t|t-1]} may be computed with `yh = a \%*\% t(C)`.}
#' \item{ll}{(scaled) Gaussian log likelihood of the model
#'          \deqn{-\frac{1}{2N}\sum_{t=1}^{N}\left(m\log(2\pi) + \log\det\Sigma_{t|t-1} +
#'           (y_t - y_{t|t-1})' \Sigma_{t|t-1}^{-1} (y_t - y_{t|t-1}) \right).}{
#'           (-1/(2N))\sum_{t=1}^{N}[ m log(2\pi)  + log det \Sigma[t|t-1] +
#'            (y[t] - y[t|t-1])' \Sigma[t|t-1]^{-1} (y[t] - y[t|t-1]) ].}}
#' \item{P1}{\eqn{(s,s)} dimensional covariance matrix of the error of the state prediction
#'           \eqn{a_{N+1|N}}{a[N+1|N]}, i.e. this matrix corresponds to \eqn{\Pi_{N+1|N}}{P[N+1|N]}.}
#'
#' @section Notes:
#' The `RcppArmadillo` functions (`kf_cpp` and `kf2_cpp`) do not check the input parameters,
#' so these function must be used with some care.
#'
#' The procedures only accept "wide" state space systems (\eqn{m \leq n}{m \le n}), since for
#' "tall" systems (\eqn{m > n}) the variance of the prediction errors
#' (\eqn{\Sigma_{t+1|t}}{\Sigma[t+1|t]}) is singular for \eqn{t} larger than some threshold.
#'
#' Up to now, there is no support for models with exogenous inputs.
#'
#' @references
#' \insertRef{AndersonMoore2005}{RLDM}
#'
#' @seealso There exist a number of R packages which implement the Kalman filter, e.g.
#' \pkg{KFAS}. However, most of these packages do not allow for correlations between the
#' "state noise" and the "measurement noise".
#'
#' If only the likelihood is needed, then one may use [ll_kf()].
#'
#' @name kf
#' @rdname kf
#' @export
#'
#' @examples
#' s = 4  # state dimension
#' m = 2  # number of outputs
#' n = 3  # number of inputs,
#' n.obs = 100 # sample size
#'
#' # generate a (stable) state space model
#' tmpl = tmpl_stsp_full(m, n, s, sigma_L = "chol")
#' model = r_model(tmpl, bpoles = 1, sd = 0.5)
#' # generate a sample
#' data = sim(model, n.obs = n.obs, a1 = NA)
#'
#' # compute Q, R, S and P1
#' sigma_L = model$sigma_L
#' sigma = tcrossprod(sigma_L)
#' R = model$sys$D %*% sigma %*% t(model$sys$D)
#' S = model$sys$B %*% sigma %*% t(model$sys$D)
#' Q = model$sys$B %*% sigma %*% t(model$sys$B)
#' P1 = lyapunov(model$sys$A, Q)
#'
#' # call Kalman filter. Note y_t = t(y)!
#' out = kf_cpp(model$sys$A, model$sys$C, Q, R, S, t(data$y), P1, double(s))
#' # use the wrapper function
#' out_test = kf(model, data$y, method = 'kf')
#' all.equal(out, out_test)
#'
#' # compute H and square root of P1
#' H = rbind(model$sys$D, model$sys$B) %*% sigma_L
#' P1_R = chol(P1)
#'
#' # call square root Kalman filter. Note H_t = t(H) and y_t = t(y)!
#' out_test = kf2_cpp(model$sys$A, model$sys$C, t(H), t(data$y), P1_R, double(s))
#' all.equal(out, out_test)
#' # use the wrapper function
#' out_test = kf(model, data$y, method = 'kf2')
#' all.equal(out, out_test)
#'
#' # The one step ahead predictions for y[t] may be computed by
#' yh = out$a %*% t(model$sys$C)
#' # and the (non scaled) prediction errors are
#' uh = data$y - out$a[1:n.obs,] %*% t(model$sys$C)
kf = function(model, y, method = c('kf','kf2'), P1 = NULL, a1 = NULL) {

  # check input parameters

  method = match.arg(method)

  if (!inherits(model, 'stspmod')) stop('model must be an "stspmod" object!')

  d = dim(model$sys)
  m = d[1]
  n = d[2]
  s = d[3]

  if ((m == 0) || (n < m)) stop('only wide, non empty state space systems (n>=m>0) are supported.')

  # check 'y'
  # coerce data object(s) to matrices
  y = try(as.matrix(y))
  if ( inherits(y, 'try-error') || (!is.numeric(y)) || (!is.matrix(y)) ) {
    stop('input "y" must be a data object which may be coerced to a matrix with "as.matrix(y)')
  }
  if (length(y) == 0) stop('"y" contains no data')
  n.obs = nrow(y)  # sample size
  if (ncol(y) != m) stop('data "y" is not compatible with the model!')
  if (any(!is.finite(y))) stop('"y" contains NAs/NaNs/Infs!')

  # check 'P1'
  # note we do not check that P1 is symmetric and psd.
  if (is.null(P1)) {
    P1 = lyapunov(model$sys$A, tcrossprod(model$sys$B %*% model$sigma_L), non_stable = 'stop')
  }
  if (!is.numeric(P1))  stop('parameter "P1" is not numeric')
  if ( is.vector(P1) ) {
    if (length(P1) == s) P1 = diag(P1, nrow = s, ncol = s)
    if (length(P1) == (s^2)) P1 = matrix(P1, nrow = s, ncol = s)
  }
  if ( (!is.matrix(P1)) || any(dim(P1) != s) ) {
    stop('parameter "P1" is not compatible with the model')
  }

  # check 'a1'
  if (is.null(a1)) a1 = double(s)
  a1 = as.vector(a1)
  if (length(a1) != s) stop('parameter "a1" is not compatible with the model')

  sigma = tcrossprod(model$sigma_L)
  if (method == 'kf') {
    # call 'Kalman filter'
    R = model$sys$D %*% sigma %*% t(model$sys$D)
    S = model$sys$B %*% sigma %*% t(model$sys$D)
    Q = model$sys$B %*% sigma %*% t(model$sys$B)
    out = .Call(`_RLDM_kf_cpp`, model$sys$A, model$sys$C, Q, R, S, t(y), P1, a1)
    return( out )
  } else {
    # call 'square root Kalman filter'
    # H_t = sigma_R %*% cbind(t(model$sys$D), t(model$sys$B))
    H = rbind(model$sys$D, model$sys$B) %*% model$sigma_L
    P1_R = try(chol(P1), silent = TRUE)
    if (inherits(P1_R,'try-error')) {
      # cholesky decomposition failed, P1 is singular?
      # try eigenvalue decomposition
      evd = eigen(P1, symmetric = TRUE)
      if (min(evd$values) < -sqrt(.Machine$double.eps)) {
        stop('"P1" is not positive semidefinite!')
      }
      P1_R = diag(sqrt(pmax(evd$values,0)), nrow = s, ncol = s) %*% t(evd$vectors)
    }
    out = .Call(`_RLDM_kf2_cpp`, model$sys$A, model$sys$C, t(H), t(y), P1_R, a1)
    return( out )
  }
}


#' @title
#' Gaussian log Likelihood of a State Space Model
#'
#' @description
#' These routines compute the log Likelihood for time invariant, linear state space models of the form
#' \deqn{a_{t+1} = A a_t + Bu_t}{a[t+1] = A a[t] + B u[t]}
#' \deqn{y_t = C a_t + Du_t}{y[t] = C a[t] + D u[t]}
#' with \eqn{m}-dimensional outputs \eqn{y_t}{y[t]}, \eqn{s}-dimensional states
#' \eqn{a_t}{a[t]} and \eqn{n}-dimensional disturbances \eqn{u_t}{u[t]}.
#' The disturbances are white noise with a covariance matrix
#' \eqn{\mathbf{E} u_t u_t'=\Sigma}{E u[t]u[t]'=\Sigma}.
#' Note that the disturbances and the outputs may have *different* dimensions, however,
#' only "wide" systems with (\eqn{m\leq n}{m\le  n}) are implemented.
#'
#' The Gaussian log likelihood (for the case of Gaussian disturbances
#' \eqn{u_t\sim N(0,\Sigma)}{u[t] ~ N(0,\Sigma)} and
#' \eqn{a_1\sim N(a_{1|0},\Pi_{1|0})}{a[1] ~ N(a[1|0],P[1|0])}) here is computed by the
#' standard Kalman Filter or the square root Kalman filter, see [kf()].
#' The Kalman filter is a recursive
#' scheme to compute the linear, least squares predictions
#' for \eqn{a_{t+1}}{a[t+1]} and \eqn{y_{t+1}}{y[t+1]} given the observations
#' \eqn{y_t,\ldots,y_1}{y[t],\ldots,y[1]} up to time \eqn{t}. These predictions are notated with
#' \eqn{a_{t+1|t}}{a[t+1|t]} and \eqn{y_{t+1|t}}{y_[t+1|t]}, the
#' prediction error for the output \eqn{y_{t+1}}{y[t+1]} is
#' \eqn{\epsilon_{t+1|t}=(y_{t+1}-y_{t+1|t})}{\epsilon[t+1|t]=(y[t+1]-y[t+1|t])}
#' and the corresponding variances of the prediction errors are
#' \deqn{\Pi_{t+1|t}=\mathbf{E}(a_{t+1}-a_{t+1|t})
#' (a_{t+1}-a_{t+1|t})',}{P[t+1|t]=E(a[t+1]-a_[t+1|t])(a[t+1]-a_[t+1|t])',}
#' \deqn{\Sigma_{t+1|t}=\mathbf{E}(\epsilon_{t+1|t}
#' \epsilon_{t+1|t}').}{\Sigma[t+1|t]=E(\epsilon_[t+1|t]\epsilon_[t+1|t]').}
#'
#' The standard form of the Kalman filter is based on the parameter matrices \eqn{A,C}, the variance of
#' "state disturbances"
#' \eqn{Q=\mathbf{E}(Bu_t (Bu_t)')=(B\Sigma B')}{Q=E(Bu[t](Bu[t])')=(B\Sigma B')}, the variance
#' of the "measurement disturbances"
#' \eqn{R=\mathbf{E}(Du_t (Du_t)')=(D\Sigma D')}{R=(Du[t](Du[t])')=E(D\Sigma D')} and the covariance
#' \eqn{S=\mathbf{E}(Bu_t(Du_t)')=(B\Sigma D')}{S=(Bu[t](Du[t])')=E(B\Sigma D')}.
#' Furthermore we need the initial prediction
#' \eqn{a_{1|0}}{a[1|0]} and the corresponding error variance
#' \eqn{\Pi_{1|0}}{P[1|0]}.
#'
#' For the square root form of the filter we need the "square roots"
#' \eqn{\Pi_{1|0}^{1/2}}{P[1|0]^{1/2}} and \eqn{\Sigma^{1/2}}, i.e. matrices such that
#' \eqn{\Pi_{1|0} = \Pi_{1|0}^{1/2} (\Pi_{1|0}^{1/2})'}{P[1|0] = P[1|0]^{1/2} (P[1|0]^{1/2})'}
#' and \eqn{\Sigma = \Sigma^{1/2}(\Sigma^{1/2})'}. In addition, we define
#' \eqn{H=(D',B')'\Sigma^{1/2}}.
#'
#' The (scaled) Gaussian log Likelihood of this model then may be expressed as
#'          \deqn{\frac{-1}{2N}\sum_{t=1}^{N}\left(m\log(2\pi) + \log\det\Sigma_{t|t-1} +
#'           (y_t - y_{t|t-1})' \Sigma_{t|t-1}^{-1} (y_t - y_{t|t-1}) \right).}{
#'           (-1/(2N))\sum_{t=1}^{N}[ m log(2\pi)  + log det \Sigma[t|t-1] +
#'            (y[t] - y[t|t-1])' \Sigma[t|t-1]^{-1} (y[t] - y[t|t-1]) ].}
#'
#' @details
#' The core routines are  `ll_kf_cpp` and `ll_kf2_cpp` which are \pkg{RcppArmadillo} implementations
#' of the standard and the square root Kalman filter. The function `ll_kf` is a wrapper function,
#' which extracts the necessary parameters from an [stspmod()] object,
#' computes the initial covariance matrix `P1` and the initial state
#' estimate `a1` (if not provided) and then calls `ll_kf_cpp` or `ll2_kf_cpp`.
#'
#' Square root Kalman filter: For the square root
#' \eqn{\Pi_{1|0}^{1/2}}{P[1|0]^{1/2}} the procedure first tries the Cholesky decomposition.
#' If this fails (since \eqn{\Pi_{1|0}^{1/2}}{P[1|0]^{1/2}} is (close to) singular),
#' then `ll_kf` tries to compute a symmetric square root via the eigenvalue decomposition
#' of \eqn{\Pi_{1|0}^{1/2}}{P[1|0]^{1/2}}.
#'
#' @section Notes:
#'
#' The procedures only accept "wide" state space systems (\eqn{m \leq n}{m \le n}), since for
#' "tall" systems (\eqn{m > n}) the variance of the prediction errors
#' (\eqn{\Sigma_{t+1|t}}{\Sigma[t+1|t]}) is singular for \eqn{t} larger than some threshold.
#'
#'
#' @inheritParams kf
#' @param method Character string. If `method="kf"` then `ll_kf` calls
#'               `ll_kf_cpp` ("standard form" of the Kalman filter) and for
#'               `method="kf2"` the "square root" form of the Kalman filter is used,
#'               i.e. `ll_kf2_cpp` is called. Up to numerical errors the outputs should not
#'               depend on the chosen method.
#' @param tol (small) tolerance value (or zero). In order to speed up the computations,
#'              the algorithm(s) switch to a constant Kalman gain when there is no significant change in
#'              state error covariance. This behavior is controlled by the parameter `tol` and
#'              may be switched off by setting `tol=0`.
#'
#' @return (double) The Gaussian log Likelihood of the model.
#'
#'
#' @seealso [kf()] for computation of the predictions \eqn{a_{t+1|t}}{a[t+1|t]} and
#' \eqn{y_{t+1|t}}{y[t+1|t]}.
#'
#'
#' @name ll_kf
#' @rdname ll_kf
#' @export
#'
#' @examples
#' s = 4  # state dimension
#' m = 2  # number of outputs
#' n = m  # number of inputs (square case m=n)
#' n.obs = 100 # sample size
#'
#' # generate a (stable) state space model (in innovation form)
#' tmpl = tmpl_stsp_full(m, n, s, sigma_L = "chol")
#' model = r_model(tmpl, bpoles = 1, sd = 0.5)
#' # generate a sample
#' data = sim(model, n.obs = n.obs, a1 = NA)
#'
#' # compute Q, R, S and P1
#' sigma_L = model$sigma_L
#' sigma = tcrossprod(sigma_L)
#' R = model$sys$D %*% sigma %*% t(model$sys$D)
#' S = model$sys$B %*% sigma %*% t(model$sys$D)
#' Q = model$sys$B %*% sigma %*% t(model$sys$B)
#' P1 = lyapunov(model$sys$A, Q)
#'
#' # compute H and square root of P1
#' H = rbind(model$sys$D, model$sys$B) %*% sigma_L
#' P1_R = chol(P1)
#'
#' # compute logLikelihood (via Kalman Filter)
#' ll = ll_kf(model, data$y)
#'
#' # compute logLikelihood (via square root Kalman Filter)
#' ll_test = ll_kf(model, data$y, method = 'kf2')
#' all.equal(ll, ll_test)
#'
#' # directly call Rcpp routines, note H_t = t(H) and y_t = t(y)
#' ll_test = ll_kf_cpp(model$sys$A, model$sys$C, Q, R, S,
#'                     t(data$y), P1, a1 = double(s), tol=1e-8)
#' all.equal(ll, ll_test)
#' ll_test = ll_kf2_cpp(model$sys$A, model$sys$C, t(H),
#'                      t(data$y), P1_R, a1 = double(s), tol=1e-8)
#' all.equal(ll, ll_test)
#'
#' # call the "full" kf routines
#' out = kf(model, data$y)
#' all.equal(ll, out$ll)
#' out = kf(model, data$y, method = 'kf2')
#' all.equal(ll, out$ll)
ll_kf = function(model, y, method = c('kf', 'kf2'),  P1 = NULL, a1 = NULL, tol = 0) {

  method = match.arg(method)

  stopifnot("Model must be an 'stspmod' object." = inherits(model, 'stspmod'))

  d = dim(model$sys)
  m = d[1]
  n = d[2]
  s = d[3]

  if ((m == 0) || (n < m)) stop('only wide, non empty state space systems (n>=m>0) are supported.')

  # check 'y'
  # coerce data object(s) to matrices
  y = try(as.matrix(y))
  if ( inherits(y, 'try-error') || (!is.numeric(y)) || (!is.matrix(y)) ) {
    stop('input "y" must be a data object which may be coerced to a matrix with "as.matrix(y)')
  }
  if (length(y) == 0) stop('"y" contains no data')
  n.obs = nrow(y)  # sample size
  if (ncol(y) != m) stop('data "y" is not compatible with the model!')
  if (any(!is.finite(y))) stop('"y" contains NAs/NaNs/Infs!')

  # check 'P1'
  # note we do not check that P1 is symmetric and psd.
  if (is.null(P1)) {
    if (s >0) {
      P1 = lyapunov(model$sys$A, tcrossprod(model$sys$B %*% model$sigma_L), non_stable = 'stop')
    } else {
      P1 = matrix(0, nrow = 0, ncol = 0)
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

  # check 'a1'
  if (is.null(a1)) a1 = double(s)
  a1 = as.vector(a1)
  if (length(a1) != s) stop('parameter "a1" is not compatible with the model')

  sigma = tcrossprod(model$sigma_L)
  if (method == 'kf') {
    # call 'Kalman filter'
    R = model$sys$D %*% sigma %*% t(model$sys$D)
    S = model$sys$B %*% sigma %*% t(model$sys$D)
    Q = model$sys$B %*% sigma %*% t(model$sys$B)
    # return(list(model$sys$A, model$sys$C, Q, R, S, t(y), P1, a1))
    out = .Call(`_RLDM_ll_kf_cpp`, model$sys$A, model$sys$C, Q, R, S, t(y), P1, a1, tol)
    return( out )
  } else {
    # call 'square root Kalman filter'
    # H_t = sigma_R %*% cbind(t(model$sys$D), t(model$sys$B))
    H = rbind(model$sys$D, model$sys$B) %*% model$sigma_L
    if (s > 0) {
      P1_R = try(chol(P1), silent = TRUE)
      if (inherits(P1_R,'try-error')) {
        # cholesky decomposition failed, P1 is singular?
        # try eigenvalue decomposition
        evd = eigen(P1, symmetric = TRUE)
        if (min(evd$values) < -sqrt(.Machine$double.eps)) {
          stop('"P1" is not positive semidefinite!')
        }
        P1_R = diag(sqrt(pmax(evd$values,0)), nrow = s, ncol = s) %*% t(evd$vectors)
      }
    } else {
      P1_R = matrix(0, nrow = 0, ncol = 0)
    }
    out = .Call(`_RLDM_ll_kf2_cpp`, model$sys$A, model$sys$C, t(H), t(y), P1_R, a1, tol)
    return( out )
  }
}

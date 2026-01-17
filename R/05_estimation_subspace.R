# subspace estimates #############

#' Helper Functions for Order Estimation
#'
#' These helper function are used in the subspace estimation routine
#' [subspace methods]
#' for the estimation of the order of a state space model.
#' For a discussion of order estimation in the context of
#' subspace methods see e.g \insertCite{Bauer2001}{RLDM}.
#'
#' `estorder_max` simply returns the maximum order `s.max` considered.
#'
#' `estorder_rkH` estimates the order by an estimate of the rank of the Hankel matrix.
#' If the maximum singular value is smaller than or equal to `.Machine$double.eps`
#' then the estimate is set to `s=0`. Otherwise the estimated order is equal to
#' the number of singular values which are larger than or equal to `tol` times the
#' maximum singular value.
#'
#' The function `estorder_MOE`
#' searches for a "gap" in the singular values. The order is set to
#' maximum \eqn{s} which satisfies
#' \deqn{\ln(\sigma_s) > (\ln(\sigma_1)+\ln(\sigma_m))/2}{
#'       ln(\sigma[s]) > (ln(\sigma[1])+ln(\sigma[m]))/2}
#' where \eqn{\sigma_m}{\sigma[m]} is the minimum, non zero singular value.
#' This scheme is also implemented in the N4SID procedure of the
#' system identification toolbox of MATLAB (Ljung, 1991).
#'
#' The function `estorder_SVC` implements the so called *Singular Value Criterion*
#' \deqn{svc(s) = \sigma_{s+1}^2 + c(N)d(s)/N}{svc(s)=\sigma[s+1]^2 + c(N)d(s)/N}
#' (see e.g. \insertCite{Bauer2001}{RLDM}). Here \eqn{\sigma_s}{\sigma[s]} is the \eqn{s}-th
#' singular value of the weighted Hankel marix, \eqn{N} is the sample size,
#' \eqn{d(s) = 2ms} denotes the number of parameters for a state space
#' model with \eqn{s} states (and \eqn{m} outputs) and \eqn{c(N)} is a "penalty" term.
#' This term is
#' \eqn{c(N)=\ln(N)}{c(N)=ln(N)} for `penalty = "lnN"` and
#' \eqn{c(N)=fp\ln(N)}{c(N)=f p ln(N)} for `penalty = "fplnN"`. The estimate of the
#' order is the minimizer of this criterion.
#'
#' `estorder_IVC` estimates the order via an information criterion of the form
#' \deqn{ivc(s) = \ln\det\hat\Sigma_{s} + c(N)d(s)/N}{ivc(s)=ln det \Sigma[s] + c(N)d(s)/N}
#' where \eqn{\hat\Sigma_s}{\Sigma[s]} is the estimate of the noise covariace matrix
#' obtained from a model with order \eqn{s}. Here the term \eqn{c(N)} is chosen as
#' \eqn{c(N)=2} for `penalty = "AIC"` and
#' \eqn{c(N)=\ln(N)}{c(N)=f p log(N)} for `penalty = "BIC"`.
#'
#' Note also that for both routines `estorder_SVC` and `estorder_IVC`
#' one may also set `penalty` to an arbitrary numeric value! E.g.
#' setting `penalty=-1` would ensure that  `estorder_SVC`  always choses
#' the maximum possible order `s=s.max`.
#'
#' @param s.max (integer) maximum order (maximum state space dimension).
#' @param Hsv   vector with the Hankel singular values (must have at least s.max entries).
#' @param lndetSigma  (s.max+1) dimensional vector with the logarithms of the determinants of the
#'                    respective estimated noise covariance matrices.
#' @param n.par (s.max+1) dimensional vector with the respective number of parameters.
#' @param n.obs (integer) sample size \eqn{N} (or `Inf`).
#' @param Hsize two dimensional integer vector with the number of block rows and block columns
#'              of the Hankel matrix (`Hsize = c(f,p)`).
#' @param penalty determines the penalty term. See the details below.
#' @param tol (small) tolarenace used to determine the rank of the Hankel matrix.
#' @param ... optional additional parameters.
#'
#' @return Either NULL or a list with slots `$s` (the selected/estimated order)
#'         and `$criterion` (an (`s.max+1`) dimensional vector).
#'
#' @references
#' \insertRef{Bauer2001}{RLDM}
#'
#' @name subspace order estimates
#' @export
#'
estorder_SVC = function(s.max, Hsv=NULL, n.par, n.obs, Hsize, penalty = 'lnN', ...) {
  # cat('estorder_SVC start\n')
  # print(penalty)
  if (is.null(Hsv)) return(NULL)
  if (isTRUE(all.equal(penalty, 'lnN'))) penalty = log(n.obs)
  if (isTRUE(all.equal(penalty, 'fplnN'))) penalty = prod(Hsize)*log(n.obs)
  if (is.character(penalty)) stop('unknown option penalty=', penalty)
  penalty = as.numeric(penalty)[1]
  # print(Hsv)
  # print(n.par)
  # print(n.obs)
  # print(penalty)
  penalty = penalty / n.obs
  # take care of the case n.obs = Inf: log(Inf)/Inf = NaN
  if (!is.finite(penalty)) penalty = 0

  # Hsv should have at least s.max entries!
  if (length(Hsv) < (s.max+1)) Hsv = c(Hsv, 0)
  criterion = Hsv[1:(s.max+1)]^2 + n.par*penalty
  s = which.min(criterion) - 1
  # print(criterion)
  # print(s)
  # cat('estorder_SVC end\n')
  return(list(s = s, criterion = criterion))
}

#' @name subspace order estimates
#' @export
estorder_IVC = function(s.max, lndetSigma=NULL, n.par, n.obs, penalty = 'BIC', ...) {
  # cat('estorder_IVC start\n')
  # print(penalty)
  if (is.null(lndetSigma)) return(NULL)
  if (isTRUE(all.equal(penalty, 'AIC'))) penalty = 2
  if (isTRUE(all.equal(penalty, 'BIC'))) penalty = log(n.obs)
  if (is.character(penalty)) stop('unknown option penalty=', penalty)
  penalty = as.numeric(penalty)[1]
  # print(lndetSigma)
  # print(n.par)
  # print(n.obs)
  # print(penalty)
  penalty = penalty / n.obs
  # take care of the case n.obs = Inf: log(Inf)/Inf = NaN
  if (!is.finite(penalty)) penalty = 0

  criterion = lndetSigma + n.par*penalty
  s = which.min(criterion) - 1
  # print(criterion)
  # print(s)
  # cat('estorder_IVC end\n')
  return(list(s = s, criterion = criterion))
}

#' @name subspace order estimates
#' @export
estorder_max = function(s.max, ...) {
  return(list(s = s.max, criterion = c(rep(1, s.max), 0)))
}

#' @name subspace order estimates
#' @export
estorder_rkH = function(s.max, Hsv = NULL, tol = sqrt(.Machine$double.eps), ...) {
  if (is.null(Hsv)) return(NULL)
  if (Hsv[1] <= .Machine$double.eps) {
    return(list(s = 0, criterion = numeric(s.max+1)))
  }
  criterion = as.integer(c(Hsv[1:s.max],0) > tol * Hsv[1])
  s = sum(criterion)
  return(list(s = s, criterion = criterion))
}

# MOE(n), which is implemented in the N4SID procedure of the system
# identifcation toolbox of MATLAB (Ljung, 1991):
# The idea here is to formalise the search for a "gap" in
# the singular values.

#' @name subspace order estimates
#' @export
estorder_MOE = function(s.max, Hsv = NULL, ...) {
  if (is.null(Hsv)) return(NULL)
  criterion = Hsv
  if (Hsv[1] <= .Machine$double.eps) {
    return(list(s = 0, criterion = numeric(s.max+1)))
  }
  criterion = as.integer(log(Hsv) > (log(Hsv[1])+log(Hsv[length(Hsv)]))/2)
  criterion = c(criterion[1:s.max], 0)
  s = sum(criterion)
  return(list(s = s, criterion = criterion))
}


# core computations of the AOKI method #####################
aoki = function(m, s, svd.H, gamma0, Rp, Rf) {
  if (s == 0) {
    return(list(s = 0, A = matrix(0, nrow = 0, ncol = 0), B = matrix(0, nrow = 0, ncol = m),
                C = matrix(0, nrow = m, ncol = 0), D = diag(m),
                sigma = gamma0, sigma_L = t(chol(gamma0))))
  }

  sv2 = sqrt(svd.H$d[1:s])

  # (approximately) factorize H as H = U V
  U = svd.H$u[, 1:s, drop = FALSE] * matrix(sv2, nrow = nrow(svd.H$u), ncol = s, byrow = TRUE)
  V = t( svd.H$v[, 1:s, drop = FALSE] * matrix(sv2, nrow = nrow(svd.H$v), ncol = s, byrow = TRUE) )
  U = t(Rf) %*% U
  V = V %*% Rp

  C = U[1:m, , drop = FALSE]
  M = V[, 1:m, drop = FALSE]
  A = suppressWarnings(stats::lsfit(U[1:(nrow(U) - m), , drop = FALSE],
                                    U[(m+1):nrow(U), , drop = FALSE],
                                    intercept = FALSE)$coef)
  A = unname(A)
  A = matrix(A, nrow = s, ncol = s)

  # solve the Riccati equation to get the "B" matrix (Kalman gain)
  out = suppressWarnings(riccati(A, M, C, G = gamma0, only.X = FALSE))

  # construct state space model object
  sigma = out$sigma
  sigma_L = t(chol(sigma))
  return(list(s = s, A = A, B = out$B, C = C, D = diag(m), sigma = sigma, sigma_L = sigma_L))
}

#' Subspace Helper Methods
#'
#' These procedure implement two subspace algorithms for the estimation of state space models,
#' the *AOKI* method, as described in \insertCite{Aoki90}{RLDM} and the *CCA* algorithm
#' (see e.g. \insertCite{DahlenScherrer2004}{RLDM} or \insertCite{Bauer2001}{RLDM}).
#' These subspace algorithms center on the weighted Hankel matrix
#' \deqn{(R_f')^{-T} H_{fp} R_p^{-1}}{(Rf')^{-1} H Rp^{-1}}
#' where the block *Hankel* matrix \eqn{H_{fp}}{H} is the covariance between the
#' "past" \eqn{(y_{t-1}',\cdots,y_{t-p}')'}{(y[t-1]',...,y[t-p]')'} and the
#' "future" \eqn{(y_{t}',\cdots,y_{t+f-1}')'}{(y[t]',...,y[t+f-1]')'} and
#' \eqn{R_f}{Rf} and \eqn{R_p}{Rp} are the cholesky factors of the covariance matrices
#' of the "future" and the "past" respectively. The singular values of this weighted Hankel matrix
#' are the canonical correlation coefficients between the past and the future.
#' Note that the implementation here always sets \eqn{f = p+1}.
#'
#' AOKIs method is a *realization* algorithm, i.e. it reconstructs the underlying state space
#' model from the (population) autocovariance function. To this end a Riccati equation has to be solved,
#' see [riccati()].
#' If an estimated ACF is fed into this algorithm one obtains an estimate for the
#' state space model. However note that this may fail (in particular the Riccati
#' equation may have no positive definite solution) if the estimate of the ACF is *not* positive
#' definite, if the Hankel matrix is too small or if state dimension is not correct.
#'
#' The CCA method estimates the state space model by first constructing an estimate of the states.
#' Then the parameter matrices are estimated via simple LS regressions. This procedure does
#' not give the "true" model, even in the case when the population ACF is used. However,
#' the "distance" between the true model
#' and the estimated model converges to zero,
#' if the estimate of the ACF converges to the population ACF and
#' the size \eqn{p} of the Hankel matrix converges to infinity.
#'
#' There are two implementations of the CCA method:
#' 1. The routine [est_stsp_cca_sample()]
#'    operates directly on the supplied data.
#' 2. The routine [est_stsp_cca()] uses an (estimated) autocovariance function.
#'
#' These algorithms may also be used as simple "model reduction algorithms". If we want to
#' approximate a high dimensional state space model by a model of lower order, we may proceed
#' as follows. First we compute the ACF of the high dimensional model and then fed
#' this ACF into the subspace routines `est_stsp_cca` or `est_stsp_aoki`,
#' however setting the maximum order `s.max` to some value less than the true order.
#'
#' **Order Estimation**
#'
#' The order estimation is based on the Hankel singular values \eqn{\sigma_s}{\sigma[s]} and/or the
#' log det values of the estimated noise covariance matrices \eqn{\ln\det \hat{\Sigma}_s}{ln det \Sigma[s]}.
#' Using only the Hankel singular values has the advantage that only *one* model has to be estimated,
#' whereas otherwise estimates for *all* models with orders \eqn{s=0,\ldots,s_{\max}}{s=0,dots,s[max]}
#' have to be computed.
#'
#' In order to exploit this (small) advantage of singular values based criteria the
#' order estimation runs as follows:
#' First the procedures call
#' \cr
#' `  estorder(s.max, Hsv, n.par, m, n.obs, Hsize=c(f,p), ...) `
#' \cr
#' Here `Hsv` is an \eqn{pm} dimensional vector with the Hankel singular values and
#' `n.par` is an \eqn{(s_{\max}+1)}{(s[max]+1)} dimensional vector with the respective number of
#' parameters of the models with orders \eqn{s=0,\ldots,s_{\max}}{s=0,dots,s[max]}.
#' If this call returns an estimate of the order then the procedures estimate a
#' corresonding state space model.
#'
#' If this call fails (i.e returns `NULL`) then the procedures estimate
#' all models with orders up to \eqn{s_{\max}}{s[max]} and the corresponding
#' noise covariance matrices. The order then is estimated by calling
#' \cr
#' `  estorder(s.max = s.max, Hsv, lndetSigma, n.par, m, n.obs, Hsize, ...) `
#' \cr
#' where `lndetSigma` is the vector with the log det values of the estimated
#' noise covariance matrices (\eqn{\ln\det \hat{\Sigma}_s}{ln det \Sigma[s]}).
#'
#' The package offers some predefined order selection procedures
#' (see also [subspace order estimates]):
#'
#' `estorder_max(s.max, ...)` simply returns the maximum order `s.max` considered.
#'
#' `estorder_rkH(s.max, Hsv, tol, ...)` estimates the order by an estimate
#' of the rank of the Hankel matrix.
#'
#' `estorder_MOE(s.max, Hsv, ...)` estimates the order by searching
#' for a "gap" in the singular values.
#'
#' `estorder_SVC(s.max, Hsv, n.par, n.obs, Hsize, penalty, ...)`
#' implements the so called *Singular Value Criteria*, see \insertCite{Bauer2001}{RLDM}:
#' \deqn{svc(s) = \sigma_{s+1}^2 + c(N)d(s)/N}{svc(s)=\sigma[s+1]^2 + c(N)d(s)/N}
#' Here \eqn{\sigma_s}{\sigma[s]} is the \eqn{s}-th
#' singular value of the weighted Hankel marix, \eqn{N} is the sample size,
#' \eqn{d(s) = 2ms} denotes the number of parameters for a state space
#' model with \eqn{s} states (and \eqn{m} outputs) and \eqn{c(N)}) is a "penalty"
#' (depending on the sample size).
#'
#' The above order estimation procedures only use the Hankel singular values, whereas
#' the following procedure is based on the estimated noise covariances.
#'
#' `estorder_IVC(s.max, lndetSigma, n.par, n.obs, penalty, ...)` estimates the order
#' via an information criterion of the form
#' \deqn{ivc(s) = \ln\det\hat\Sigma_{s} + c(N)d(s)/N}{ivc(s)=ln det \Sigma[s] + c(N)d(s)/N}
#' where \eqn{\hat\Sigma_s}{\Sigma[s]} is the estimate of the noise covariace matrix
#' obtained from a model with order \eqn{s}, \eqn{d(s)} denotes the number of parameters
#' and \eqn{c(N)} is a "penalty" (depending on the sample size).
#'
#' For both `estorder_SVC` and `estorder_IVC` the (optional) parameter
#' `penalty` controls the penalty term \eqn{c(N)}.
#'
#' Note also that for  `keep_models==TRUE` the estimation procedures compute *all*
#' models even in the case of a Hankel singular value based selection criterion.
#'
#' @param gamma \eqn{(m,m,L+1)}-dimensional array with the (sample) autocovariance function.
#' @param y \eqn{(N,m)}-dimensional matrix or an object, which may be coerced to
#'          a matrix with \code{as.matrix{y}}.
#' @param s.max (integer) maximum possible order.
#' @param p     (integer) number of block columns of the Hankel matrix (size of the "past")
#' @param estorder function, to estimate the order of the system.
#' @param keep_models (boolean) should the function return a list with estimated system of
#'                    order `0:s.max`?
#' @param n.obs sample size \eqn{N}.
#' @param mean_estimate Character string giving the method used to estimate
#'        the mean \eqn{\mu = E y_t}{\mu = E y[t]}.
#'        Default is to use the sample mean.
#' @param ... additional parameters, passed on to the order estimation routine.
#'
#' @return List with slots:
#' \item{model}{a [stsp()] object, which represents the estimated state space model.}
#' \item{models}{either `NULL` (if `!keep_models`) or a list with the parameters
#'               of the estimated models with orders (`s=0:s.max+1`). This slot may e.g. be
#'               used to estimate the model order by some user defined model selection procedure.}
#' \item{s}{(integer) the estimate of the model order.}
#' \item{info}{list with information about the data and the design parameters of the estimation procedure.}
#' \item{stats}{((s.max+1)-by-5)-dimensional matrix with statistics of the (estimated) state space models.}
#'
#'
#' @references
#' \insertRef{Bauer2001}{RLDM}
#'
#' @export
#' @name subspace helpers
est_stsp_aoki = function(gamma, s.max, p, estorder = estorder_SVC,
                         keep_models = FALSE, n.obs = NULL, ...) {

  # check gamma
  if ( (!is.numeric(gamma)) || (!is.array(gamma)) ||
       (length(dim(gamma)) != 3) || (dim(gamma)[1] != dim(gamma)[2]) ) {
    stop('input "gamma" is not a valid 3-D array.')
  }
  d = dim(gamma)
  m = d[1]
  if (m == 0) {
    stop('the autocovariance function "gamma" is empty.')
  }
  lag.max = d[3] - 1
  if (lag.max <= 2) {
    stop('the autocovariance function must contain at least 2 lags.',
         ' lag.max=', lag.max)
  }
  if (is.null(n.obs)) {
    n.obs = Inf
  }
  n.obs = as.numeric(n.obs)[1]
  if (is.finite(n.obs)) {
    n.obs = as.integer(n.obs)[1]
    if (n.obs <= 0) stop('the sample size "n.obs" must be non negative')
  } else {
    n.obs = Inf
  }

  # p <=> past
  p = as.integer(p)[1]
  if (p < 1) stop('the number of lags "p" must be positive')
  if (lag.max < (2*p)) {
    stop('the autocovariance function must contain at least 2p lags.',
         ' lag.max=', lag.max, ' < 2*p=', 2*p)
  }
  # f <=> future
  f = p + 1

  s.max = as.integer(s.max)[1]
  if (s.max < 0) stop('maximum state dimension "s.max" must be non negative.')
  if (s.max > (p*m)) {
    stop('the number of lags "p" is too small for the desired maximum state dimension:',
         ' p*m=', p*m, ' < s.max=', s.max)
  }

  gamma0 = matrix(gamma[,,1], nrow = m, ncol = m)

  # covariance between "future" (y[t]', y[t+1]',...,y[t+f-1]')' and
  #                    "past" (y[t-1]', y[t-2]',...,y[t-p]')'
  H = bhankel(gamma[,,-1,drop=FALSE], d = c(f,p))
  # H0 = H
  # junk = svd(H, nu = 0, nv = 0)$d
  # print(junk / junk[1])

  # covariance of the "past" (y[t-1]', y[t-2]',...,y[t-p]')'
  Gp = btoeplitz(R = gamma[,,1:p,drop = FALSE])
  # covariance matrix of the "future" (y[t]', y[t+1]',...,y[t+f-1]')'
  Gf = btoeplitz(C = gamma[,,1:f,drop = FALSE])

  # cholesky factors of covariance matrices
  Rp = chol(Gp)
  Rf = chol(Gf)

  # weighted Hankel matrix: inv(t(Rf)) * H * inv(Rp)
  # note: inv(t(Rf)) * H * inv(Rp) is the covariance between the
  # "standardized" future and the "standardized" past.
  # The singular value decomposition of this matrix defines the
  # canonical correlation coefficients between future and past.
  ## junk = H
  H = backsolve(Rf, H, transpose = TRUE)
  ## testthat::expect_equivalent(junk, t(Rf) %*% H)
  ## junk = H
  H = t(backsolve(Rp, t(H), transpose = TRUE))
  ## testthat::expect_equivalent(junk, H %*% Rp)

  # compute SVD of weighted Hankel matrix
  svd.H = svd(H)
  Hsv = svd.H$d
  # print(Hsv / Hsv[1])

  # construct a matrix to collect the selection criteria for the models with
  # state dimension s = 0:s.max
  stats = matrix(NA_real_, nrow = s.max+1, 5)
  colnames(stats) = c('s', 'n.par', 'Hsv', 'lndetSigma', 'criterion')
  stats[, 's'] = 0:s.max                # state dimension
  stats[, 'n.par'] = 2*stats[, 's']*m   # number of parameters
  # Hankel singular values
  # print(stats)
  # print(c(NA_real_, Hsv[1:s.max]))
  stats[, 'Hsv'] = c(NA_real_, Hsv[iseq(1,s.max)])

  info = list(m = m, Hsize = c(f,p), n.obs = n.obs, s.max = s.max, Hsv = Hsv)

  if (s.max == 0) {
    # maximum order = 0 ==> there is nothing to compute
    sigma = gamma0
    sigma_L = t(chol(sigma))
    model= list(s = 0, A = matrix(0, nrow = 0, ncol = 0), B = matrix(0, nrow = 0, ncol = m),
                C = matrix(0, nrow = m, ncol = 0), D = diag(m), sigma = sigma, sigma_L = sigma_L)
    stats[1, 'lndetSigma'] = 2*sum(log(diag(sigma_L)))
    if (keep_models) {
      models = list(model)
    } else {
      models = NULL
    }
    model = stspmod(sys = stsp(A = model$A, B = model$B, C = model$C, D = model$D),
                    sigma_L = model$sigma_L)
    return( list(model = model, models = models,
                 s = 0, stats = stats, info = info) )
  }

  # estimate order based on the Hankel singular values
  out.estorder = do.call(estorder, list(s.max = s.max, Hsv = Hsv,
                                        n.par = stats[, 'n.par'],
                                        m = m, n.obs = n.obs, Hsize=c(f,p), ...))

  s.hat = integer(0)
  if (!is.null(out.estorder)) {
    s.hat = out.estorder$s
    stats[, 'criterion'] = out.estorder$criterion
  }
  if (is.null(out.estorder) || (keep_models)) {
    s = (0:s.max)
  } else {
    s = s.hat
  }
  # cat('first round\n')
  # print(out.estorder)
  # cat(s, 's.hat', s.hat, '\n')

  models = vector('list', length(s))
  for (i in (1:length(s))) {
    si = s[i]
    j = which(stats[, 's'] == si)
    # estimate model with order s = s[i]
    out = try(aoki(m, si, svd.H, gamma0, Rp, Rf), silent = TRUE)
    if (!inherits(out, 'try-error')) {
      models[[i]] = out
      lndetSigma = 2*sum(log(diag(out$sigma_L)))
      stats[j, 'lndetSigma'] = lndetSigma
    } else {
      # 'aoki' failed!!!
      model[[i]] = list(s = si)
      stats[j, 'lndetSigma'] = NA_real_
    }
  }

  if (is.null(out.estorder)) {
    # estimate order based on the lndetSigma values
    out.estorder = do.call(estorder, list(s.max = s.max, Hsv = Hsv,
                                          lndetSigma = stats[, 'lndetSigma'],
                                          n.par = stats[, 'n.par'],
                                          m = m, n.obs = n.obs, Hsize=c(f,p), ...))
    if (is.null(out.estorder)) stop('could not estimate the model order?!')
    # cat('second round\n')
    # print(out.estorder)
    # cat(s, 's.hat', s.hat, '#models', length(models), '\n')
    s.hat = out.estorder$s
    stats[, 'criterion'] = out.estorder$criterion
  }


  if (length(s) == 1) {
    model = models[[1]]
  } else {
    model = models[[s.hat+1]]
  }
  if (is.null(model$A)) stop('AOKI failed for the selecetd model order!')
  model = stspmod( sys = stsp(A = model$A, B = model$B, C = model$C, D = model$D),
                   sigma_L = model$sigma_L )
  if (!keep_models) models = NULL

  return( list(model = model, models = models,
               s = s.hat, stats = stats, info = info ) )

}

#' @export
#' @name subspace helpers
est_stsp_cca = function(gamma, s.max, p, estorder = estorder_SVC,
                        keep_models = FALSE, n.obs = NULL, ...) {

  # check gamma
  if ( (!is.numeric(gamma)) || (!is.array(gamma)) ||
       (length(dim(gamma)) != 3) || (dim(gamma)[1] != dim(gamma)[2]) ) {
    stop('input "gamma" is not a valid 3-D array.')
  }
  d = dim(gamma)
  m = d[1]
  if (m == 0) {
    stop('the autocovariance function "gamma" is empty.')
  }
  lag.max = d[3] - 1
  if (lag.max <= 2) {
    stop('the autocovariance function must contain at least 2 lags.',
         ' lag.max=', lag.max)
  }
  if (is.null(n.obs)) {
    n.obs = Inf
  }
  n.obs = as.numeric(n.obs)[1]
  if (is.finite(n.obs)) {
    n.obs = as.integer(n.obs)[1]
    if (n.obs <= 0) stop('the sample size "n.obs" must be non negative')
  } else {
    n.obs = Inf
  }

  # p <=> past
  p = as.integer(p)[1]
  if (p < 1) stop('the number of lags "p" must be positive')
  if (lag.max < (2*p)) {
    stop('the autocovariance function must contain at least 2p lags.',
         ' lag.max=', lag.max, ' < 2*p=', 2*p)
  }
  # f <=> future
  f = p + 1

  s.max = as.integer(s.max)[1]
  if (s.max < 0) stop('maximum state dimension "s.max" must be non negative.')
  if (s.max > (p*m)) {
    stop('the number of lags "p" is too small for the desired maximum state dimension:',
         ' p*m=', p*m, ' < s.max=', s.max)
  }

  gamma0 = matrix(gamma[,,1], nrow = m, ncol = m)

  # covariance between "future" Y[t] = (y[t]', y[t+1]',...,y[t+f-1]')'
  #                  and "past" X[t] = (y[t-1]', y[t-2]',...,y[t-p]')'
  YX = bhankel(gamma[,,-1,drop=FALSE], d = c(f,p))

  # covariance matrix of "past" and "present" (y[t], y[t-1]', y[t-2]',...,y[t-p]')'
  X1X = btoeplitz(R = gamma[,,1:(p+1),drop = FALSE])
  # cat(dim(X1X), m,p, '\n')
  # covariance matrix of "past" X[t] = (y[t-1]', y[t-2]',...,y[t-p]')'
  XX = X1X[1:(m*p), 1:(m*p), drop = FALSE]
  # covariance between y[t] and "past" X[t] = (y[t-1]', y[t-2]',...,y[t-p]')'
  yX = X1X[1:m, (m+1):(m*(p+1)), drop = FALSE]
  # covariance between y[t] and the "next past" X[t+1] = (y[t]', y[t-1]',...,y[t-p+1]')'
  yX1 = XX[1:m, 1:(m*p), drop = FALSE]
  # covariance matrix of y[t] (=gamma[0])
  yy = XX[1:m, 1:m, drop = FALSE]
  # covariance of the "next past" X[t+1] = (y[t], y[t-1]', ...,y[t+1-p]')'
  #                and the "past" X[t] = (y[t-1]', y[t-2]',...,y[t-p]')'
  X1X = X1X[1:(m*p), (m+1):(m*(p+1)), drop = FALSE]

  # covariance matrix of the "future" Y[t] = (y[t]', y[t+1]',...,y[t+f-1]')'
  YY = btoeplitz(C = gamma[,,1:f,drop = FALSE])

  # cholesky factors of covariance matrices
  Rp = chol(XX)
  Rf = chol(YY)

  # weighted Hankel matrix: inv(t(Rf)) * YX * inv(Rp)
  # note: inv(t(Rf)) * YX * inv(Rp) is the covariance between the
  # "standardized" future and the "standardized" past.
  # The singular value decomposition of this matrix defines the
  # canonical correlation coefficients between future and past.
  YX = backsolve(Rf, YX, transpose = TRUE)
  YX = t(backsolve(Rp, t(YX), transpose = TRUE))

  # compute SVD of weighted Hankel matrix
  svd.YX = svd(YX, nu = 0)
  Hsv = svd.YX$d # Hankel singular values
  # print(Hsv / Hsv[1])

  # construct a matrix to collect the selection criteria for the models with
  # state dimension s = 0:s.max
  stats = matrix(NA_real_, nrow = s.max+1, 5)
  colnames(stats) = c('s', 'n.par', 'Hsv', 'lndetSigma', 'criterion')
  stats[, 's'] = 0:s.max                # state dimension
  stats[, 'n.par'] = 2*stats[, 's']*m   # number of parameters
  # Hankel singular values
  # print(stats)
  # print(c(NA_real_, Hsv[iseq(1,s.max)]))
  stats[, 'Hsv'] = c(NA_real_, Hsv[iseq(1,s.max)])

  info = list(m = m, Hsize = c(f,p), n.obs = n.obs, s.max = s.max, Hsv = Hsv)

  if (s.max == 0) {
    # maximum order = 0 ==> there is nothing to compute
    sigma = gamma0
    sigma_L = t(chol(sigma))
    model= list(s = 0, A = matrix(0, nrow = 0, ncol = 0), B = matrix(0, nrow = 0, ncol = m),
                C = matrix(0, nrow = m, ncol = 0), D = diag(m), sigma = sigma, sigma_L = sigma_L)
    stats[1, 'lndetSigma'] = 2*sum(log(diag(sigma_L)))
    if (keep_models) {
      models = list(model)
    } else {
      models = NULL
    }
    model = stspmod(sys = stsp(A = model$A, B = model$B, C = model$C, D = model$D),
                    sigma_L = model$sigma_L)
    return( list(model = model, models = models,
                 s = 0, stats = stats, info = info) )
  }

  # estimate order based on the Hankel singular values
  out.estorder = do.call(estorder, list(s.max = s.max, Hsv = Hsv,
                                        n.par = stats[, 'n.par'],
                                        m = m, n.obs = n.obs, Hsize=c(f,p), ...))

  s.hat = integer(0)
  if (!is.null(out.estorder)) {
    s.hat = out.estorder$s
    stats[, 'criterion'] = out.estorder$criterion
  }
  if (is.null(out.estorder) || (keep_models)) {
    s = (0:s.max)
  } else {
    s = s.hat
  }
  # cat('first round\n')
  # print(out.estorder)
  # cat(s, 's.hat', s.hat, '\n')

  # transform covariances for the model with the maximum order max(s)
  if (max(s) > 0) {
    # T = V' * inv(R_p')
    T = t(backsolve(Rp, svd.YX$v[, 1:max(s), drop = FALSE]))
    # check
    # print(T %*% XX %*% t(T)) # = identity
    yX = yX %*% t(T)
    yX1 = yX1 %*% t(T)
    X1X = T %*% X1X %*% t(T)
  }

  models = vector('list', length(s))
  D = diag(m)
  for (i in (1:length(s))) {
    # estimate model with order s = s[i]
    # state x[t] = X[1:si,t]
    si = s[i]
    j = which(stats[, 's'] == si)
    if (si == 0) {
      # there is nothing to do
      A = matrix(0, nrow = 0, ncol = 0)
      B = matrix(0, nrow = 0, ncol = m)
      C = matrix(0, nrow = m, ncol = 0)
      sigma = yy
    } else {
      # C is defined from the regression of y[t] on x[t]: C = Eyx (Exx^{-1})
      C = yX[, 1:si, drop = FALSE]

      # sigma = E u[t] u[t]' = gamma(0) - C E x[t] x[t]' C'
      sigma = yy - C %*% t(C)

      # A is defined from the regression of x[t+1] onto x[t]: A = Ex1x (Exx)^{-1}
      A = X1X[1:si, 1:si, drop = FALSE]

      # B is defined from the regression of x[t+1] onto u[t] = y[t] - Cx[t]
      # B = (E x[t+1] y[t]' - E x[t+1] x{t]' C' ) sigma^{-1}
      #   = (E x[t+1] y[t]' - A C' ) sigma^{-1}
      ux1 = yX1[ , 1:si, drop = FALSE] - C %*% t(A)
      B = t(solve(sigma, ux1))
    }
    sigma_L = t(chol(sigma))
    lndetSigma = 2*sum(log(diag(sigma_L)))

    models[[i]] = list(s = si, A=A, B=B, C=C, D = D, sigma = sigma, sigma_L = sigma_L)
    stats[j, 'lndetSigma'] = lndetSigma
  }

  if (is.null(out.estorder)) {
    # estimate order based on the lndetSigma values
    out.estorder = do.call(estorder, list(s.max = s.max, Hsv = Hsv,
                                          lndetSigma = stats[, 'lndetSigma'],
                                          n.par = stats[, 'n.par'],
                                          m = m, n.obs = n.obs, Hsize=c(f,p), ...))
    if (is.null(out.estorder)) stop('could not estimate the model order?!')
    # cat('second round\n')
    # print(out.estorder)
    # cat(s, 's.hat', s.hat, '#models', length(models), '\n')
    s.hat = out.estorder$s
    stats[, 'criterion'] = out.estorder$criterion
  }

  if (length(s) == 1) {
    model = models[[1]]
  } else {
    model = models[[s.hat+1]]
  }
  model = stspmod( sys = stsp(A = model$A, B = model$B, C = model$C, D = model$D),
                   sigma_L = model$sigma_L )
  if (!keep_models) models = NULL

  return( list(model = model, models = models,
               s = s.hat, stats = stats, info = info ) )
}


#' @export
#' @name subspace helpers
est_stsp_cca_sample = function(y, s.max, p, estorder = estorder_SVC, keep_models = FALSE,
                               mean_estimate = c('sample.mean', 'zero'), ...) {

  mean_estimate = match.arg(mean_estimate)

  y = try(as.matrix(y))
  if ( inherits(y, 'try-error') || (!is.numeric(y)) || (!is.matrix(y)) ) {
    stop('could not coerce the input "y" to a numeric matrix')
  }
  m = ncol(y)      # number of "outputs"
  n.obs = nrow(y)  # sample size
  if (m*n.obs == 0) stop('"y" contains no data')
  if (mean_estimate == 'sample.mean') {
    y.mean = colMeans(y)
    y = scale(y, center = y.mean, scale = FALSE)
  } else {
    y.mean = numeric(m)
  }

  # p <=> past
  p = as.integer(p)[1]
  if (p < 1) stop('the number of lags "p" must be positive')
  n.valid = n.obs - 2*p
  if ((n.valid-1) < m*(p+1)) {
    stop('the sample size is too small for the desired number of lags:',
         ' (n.obs-2*p-1)=', n.valid-1, ' < m*(p+1)=', m*(p+1))
  }
  # f <=> future
  f = p + 1

  s.max = as.integer(s.max)[1]
  if (s.max < 0) stop('maximum state dimension "s.max" must be non negative.')
  if (s.max > (p*m)) {
    stop('the number of lags "p" is too small for the desired maximum state dimension:',
         ' p*m=', p*m, ' < s.max=', s.max)
  }

  # "past" X[t] = (y[t-1]', y[t-2]',...,y[t-p]')'
  X = matrix(0, nrow = n.valid, ncol = m*p)
  for (i in (1:p)) X[, ((i-1)*m+1):(i*m)] = y[(p+1-i):(n.valid+p-i), ]
  # "future" y[t] = (y[t]', y[t+1]',...,y[t+f-1]')'
  Y = matrix(0, nrow = n.valid, ncol = m*f)
  for (i in (0:(f-1))) Y[, (i*m+1):((i+1)*m)] = y[(p+1+i):(n.valid+p+i), ]
  y = y[(p+1):(n.valid+p), ]

  # QR decomposition of Y,X
  X = qr.Q(qr(X)) # X is semiorthogonal: X'X = I
  Y = qr.Q(qr(Y)) # Y is semiorthogonal: Y'Y = I

  # compute SVD of weighted "Hankel matrix"
  svd.YX = svd((t(Y) %*% X), nu = 0)
  Hsv = svd.YX$d # Hankel singular values

  # construct a matrix to collect the selection criteria for the models with
  # state dimension s = 0:s.max
  stats = matrix(NA_real_, nrow = s.max+1, 5)
  colnames(stats) = c('s', 'n.par', 'Hsv', 'lndetSigma', 'criterion')
  stats[, 's'] = 0:s.max                # state dimension
  stats[, 'n.par'] = 2*stats[, 's']*m   # number of parameters
  # Hankel singular values
  # print(stats)
  # print(c(NA_real_, Hsv[iseq(1,s.max)]))
  stats[, 'Hsv'] = c(NA_real_, Hsv[iseq(1,s.max)])

  info = list(m = m, Hsize = c(f,p), n.obs = n.obs, s.max = s.max, Hsv = Hsv)

  if (s.max == 0) {
    # maximum order = 0 ==> there is nothing to compute
    sigma = ( t(y) %*% y ) / n.valid
    sigma_L = t(chol(sigma))
    model= list(s = 0, A = matrix(0, nrow = 0, ncol = 0), B = matrix(0, nrow = 0, ncol = m),
                C = matrix(0, nrow = m, ncol = 0), D = diag(m), sigma = sigma, sigma_L = sigma_L)
    stats[1, 'lndetSigma'] = 2*sum(log(diag(sigma_L)))
    if (keep_models) {
      models = list(model)
    } else {
      models = NULL
    }
    model = stspmod(sys = stsp(A = model$A, B = model$B, C = model$C, D = model$D),
                    sigma_L = model$sigma_L)
    return( list(model = model, models = models, y.mean = y.mean,
                 s = 0, stats = stats, info = info) )
  }

  # estimate order based on the Hankel singular values
  out.estorder = do.call(estorder, list(s.max = s.max, Hsv = Hsv,
                                        n.par = stats[, 'n.par'],
                                        m = m, n.obs = n.obs, Hsize=c(f,p), ...))

  s.hat = integer(0)
  if (!is.null(out.estorder)) {
    s.hat = out.estorder$s
    stats[, 'criterion'] = out.estorder$criterion
  }
  if (is.null(out.estorder) || (keep_models)) {
    s = (0:s.max)
  } else {
    s = s.hat
  }
  # cat('first round\n')
  # print(out.estorder)
  # cat(s, 's.hat', s.hat, '\n')

  # transform X -> X * V (for the model with the maximum order max(s))
  if (max(s) > 0) {
    X = X %*% svd.YX$v[, 1:max(s), drop = FALSE] #  X'X = identity!'
    # print(t(X) %*% X)
  }

  models = vector('list', length(s))
  D = diag(m)
  for (i in (1:length(s))) {
    # estimate model with order s = s[i]
    # state x[t] = X[1:si,t]
    si = s[i]
    j = which(stats[, 's'] == si)
    if (si == 0) {
      # there is nothing to do
      A = matrix(0, nrow = 0, ncol = 0)
      B = matrix(0, nrow = 0, ncol = m)
      C = matrix(0, nrow = m, ncol = 0)
      sigma = t(y) %*% y / n.valid
    } else {
      # C is defined from the regression of y[t] on x[t]: C = Eyx (Exx^{-1})
      # note X'X = I
      # cat(si, dim(X), s, '\n')
      C = t(y) %*% X[, 1:si, drop = FALSE]

      # sigma = E u[t] u[t]' = gamma(0) - C E x[t] x[t]' C'
      # u[t] = y[t] - C x[t]
      u = y - X[ , 1:si, drop = FALSE] %*% t(C)
      sigma = (t(u) %*% u ) / n.valid

      # AB is defined from the regression of x[t+1] onto (x[t]', u[t}')'
      AB = stats::lsfit(cbind(X[1:(n.valid-1), 1:si, drop = FALSE],
                              u[1:(n.valid-1), ,drop = FALSE]),
                        X[2:n.valid, 1:si, drop = FALSE],
                        intercept = FALSE) $coefficients
      AB = unname(AB)
      AB = matrix(AB, nrow = si+m, ncol = si)
      A = t(AB[1:si,,drop = FALSE])
      B = t(AB[(si+1):(si+m),,drop = FALSE])
    }
    sigma_L = t(chol(sigma))
    lndetSigma = 2*sum(log(diag(sigma_L)))

    models[[i]] = list(s = si, A=A, B=B, C=C, D = D, sigma = sigma, sigma_L = sigma_L)
    stats[j, 'lndetSigma'] = lndetSigma
  }

  if (is.null(out.estorder)) {
    # estimate order based on the lndetSigma values
    out.estorder = do.call(estorder, list(s.max = s.max, Hsv = Hsv,
                                          lndetSigma = stats[, 'lndetSigma'],
                                          n.par = stats[, 'n.par'],
                                          m = m, n.obs = n.obs, Hsize=c(f,p), ...))
    if (is.null(out.estorder)) stop('could not estimate the model order?!')
    # cat('second round\n')
    # print(out.estorder)
    # cat(s, 's.hat', s.hat, '#models', length(models), '\n')
    s.hat = out.estorder$s
    stats[, 'criterion'] = out.estorder$criterion
  }

  if (length(s) == 1) {
    model = models[[1]]
  } else {
    model = models[[s.hat+1]]
  }
  model = stspmod( sys = stsp(A = model$A, B = model$B, C = model$C, D = model$D),
                   sigma_L = model$sigma_L )
  if (!keep_models) models = NULL

  return( list(model = model, models = models, y.mean = y.mean,
               s = s.hat, stats = stats, info = info ) )
}

#' Estimate State Space Models with Subspace Methods
#'
#' Estimate (respectively construct) a state space model from a given sample or a given
#' (sample) autocovariance function.
#'
#' The procedure implements three subspace algorithms for the estimation of state space models,
#' the *AOKI* method, as described in \insertCite{Aoki90}{RLDM} and the *CCA*
#' and *MEST* algorithms
#' (see e.g. \insertCite{DahlenScherrer2004}{RLDM}). All three algorithms center on the
#' weighted Hankel matrix
#' \deqn{(R_f')^{-T} H_{fp} R_p^{-1}}{(Rf')^{-1} H Rp^{-1}}
#' where the block *Hankel* matrix \eqn{H_{fp}}{H} is the covariance between the
#' "past" \eqn{(y_{t-1}',\cdots,y_{t-p}')'}{(y[t-1]',...,y[t-p]')'} and the
#' "future" \eqn{(y_{t}',\cdots,y_{t+f-1}')'}{(y[t]',...,y[t+f-1]')'} and
#' \eqn{R_f}{Rf} and \eqn{R_p}{Rp} are the cholesky factors of the covariance matrices
#' of the "future" and the "past" respectively. The singular values of this weighted Hankel matrix
#' are the canonical correlation coefficients between the past and the future.
#' Note that the implementation here always sets \eqn{f = p+1}.
#'
#' AOKIs method is a *realization* algorithm, i.e. it reconstructs the underlying state space
#' model from the (population) autocovariance function. To this end a Riccati equation has to be solved,
#' see [riccati()].
#' If an estimated ACF is fed into this algorithm one obtains an estimate for the
#' state space model. However note that this may fail (in particular the Riccati
#' equation may have no positive definite solution) if the estimate of the ACF is *not* positive
#' definite, if the Hankel matrix is too small or if state dimension is not correct.
#'
#' The CCA method estimates the state space model by first constructing an estimate of the states.
#' Then the parameter matrices are estimated via simple LS regressions. This procedure does
#' not give the "true" model, even in the case when the population ACF is used. However,
#' the "distance" between the true model
#' and the estimated model converges to zero,
#' if the estimate of the ACF converges to the population ACF and
#' the size \eqn{p} of the Hankel matrix converges to infinity.
#'
#' There are two implementations of the CCA method:
#' 1. If `obj` is a "time series" object and
#'    `sample2acf==FALSE` then the helper function
#'    [est_stsp_cca_sample()] is called. This implementation
#'    of CCA operates directly on the supplied data.
#' 2. If `obj` is an [autocov()] object or
#'    when `obj` is a "time series" object and
#'    `sample2acf==TRUE` then the helper function [est_stsp_cca()]
#'    is called. This implementation uses an (estimated) autocovariance function.
#'    For a time series object, first the sample
#'    autocovariance function is computed and then fed into the helper function.
#'
#' The key idea of the MEST algorithm is to first estimate a "long" AR model, convert this
#' AR model to a state space model and then to use a "balancing and truncation" method to obtain the
#' final estimate of the state space model. This scheme may be obtained by
#' calling `est_stsp_ss` with the option `extend_acf=TRUE`: This option
#' instructs the procedure to first estimate an AR(p) model and then to use this model
#' to "extend" the ACF, i.e.
#' to compute the values of the ACF for lags \eqn{p+1,\ldots,2p}{p+1,\dots,2p}.
#' Then this extended ACF is fed into the helper function  [est_stsp_cca()].
#'
#' Note that MEST uses the autocovariance function. So for a "time series" object one has
#' to set `sample2acf=TRUE`.
#'
#' These algorithms may be used for model reduction (i.e. to find a model with a
#' smaller state space dimension than the true model) and for estimation (by feeding a sample autocovariance function in).
#'
#' These algorithms may also be used as simple "model reduction algorithms". If we want to
#' approximate a high dimensional state space model by a model of lower order, we may proceed
#' as follows. First we compute the ACF of the high dimensional model and then fed
#' this ACF into the subspace routine `est_stsp_ss`, however setting the maximum order
#' `s.max` to some value less than the true order.
#' Note thet the AOKI procedure however, may break down, since it is not guaranteed that
#' the Riccati equation, which needs to be solved, has a positive semidefinite solution.
#'
#' **Size of the Hankel matrix**
#'
#' If the input parameter `p=NULL` then \eqn{p} is chosen as follows. The procedure
#' estimates the order of a "long" AR model with the AIC criterion. The size of
#' the "past" \eqn{p} then is set to \eqn{p = p_f\hat{p}_{AIC}}{p = p[f]phat}
#' where \eqn{p_f}{p[f]} is a factor (defaults to \eqn{2}) and
#' \eqn{\hat{p}_{AIC}}{phat} is the (AIC) estimate of the order of the "long"
#' AR model. See also [est_ar()].
#'
#' **Estimation of the Mean**
#'
#' If the input parameter `obj` is an [autocov()] object (which contains no info
#' about the mean \eqn{\mu=E y_t}{\mu=E y[t]}) the "estimate" of the mean is simply set to
#' a vector of `NA`'s.
#'
#' If the input parameter `obj` is a "time series" object, then there are two options.
#' For `mean_estimate == 'zero'` the procedure assumes that the process is centered
#' (\eqn{\mu=E y_t=0}{\mu=E y[t]=0}) and thus sets the estimate to a zero vector.
#' In the case `mean_estimate == 'sample.mean'`  the sample mean of the data is used.
#'
#' **Order Estimation**
#'
#' The input parameter `s.max` defines the maximum order considered.
#'
#' The order estimation is based on the Hankel singular values \eqn{\sigma_s}{\sigma[s]} and/or the
#' log det values of the estimated noise covariance matrices \eqn{\ln\det \hat{\Sigma}_s}{ln det \Sigma[s]}.
#' Using only the Hankel singular values has the advantage that only *one* model has to be estimated,
#' whereas otherwise estimates for *all* models with orders \eqn{s=0,\ldots,s_{\max}}{s=0,dots,s[max]}
#' have to be computed.
#'
#' In order to exploit this (small) advantage of singular values based criteria the
#' order estimation runs as follows:
#' First the procedures call
#' \cr
#' `  estorder(s.max, Hsv, n.par, m, n.obs, Hsize=c(f,p), ...) `
#' \cr
#' Here `Hsv` is an \eqn{pm} dimensional vector with the Hankel singular values and
#' `n.par` is an \eqn{(s_{\max}+1)}{(s[max]+1)} dimensional vector with the respective number of
#' parameters of the models with orders \eqn{s=0,\ldots,s_{\max}}{s=0,dots,s[max]}.
#' If this call returns an estimate of the order then the procedures estimate a
#' corresonding state space model.
#'
#' If this call fails (i.e returns `NULL`) then the procedures estimate
#' all models with orders up to \eqn{s_{\max}}{s[max]} and the corresponding
#' noise covariance matrices. The order then is estimated by calling
#' \cr
#' `  estorder(s.max = s.max, Hsv, lndetSigma, n.par, m, n.obs, Hsize, ...) `
#' \cr
#' where `lndetSigma` is the vector with the log det values of the estimated
#' noise covariance matrices (\eqn{\ln\det \hat{\Sigma}_s}{ln det \Sigma[s]}).
#'
#' The package offers some predefined order selection procedures
#' (see also [subspace order estimates]):
#'
#' * `estorder_max(s.max, ...)` simply returns
#'   the maximum order `s.max` considered.
#' * `estorder_rkH(s.max, Hsv, tol, ...)` estimates
#'   the order by an estimate of the rank of the Hankel matrix.
#' * `estorder_MOE(s.max, Hsv, ...)` estimates the order by searching
#'   for a "gap" in the singular values.
#' * `estorder_SVC(s.max, Hsv, n.par, n.obs, Hsize, penalty, ...)`
#'   implements the so called *Singular Value Criteria*, see
#'   \insertCite{Bauer2001}{RLDM}:
#'   \deqn{svc(s) = \sigma_{s+1}^2 + c(N)d(s)/N}{svc(s)=\sigma[s+1]^2 + c(N)d(s)/N}
#'   Here \eqn{\sigma_s}{\sigma[s]} is the \eqn{s}-th
#'   singular value of the weighted Hankel marix, \eqn{N} is the sample size,
#'   \eqn{d(s) = 2ms} denotes the number of parameters for a state space
#'   model with \eqn{s} states (and \eqn{m} outputs) and \eqn{c(N)}) is a "penalty"
#'   (depending on the sample size).
#'   \cr
#'   The above order estimation procedures only use the Hankel singular values, whereas
#'   the following procedure is based on the estimated noise covariances.
#' * `estorder_IVC(s.max, lndetSigma, n.par, n.obs, penalty, ...)` estimates the order
#'   via an information criterion of the form
#'   \deqn{ivc(s) = \ln\det\hat\Sigma_{s} + c(N)d(s)/N}{ivc(s)=ln det \Sigma[s] + c(N)d(s)/N}
#'   where \eqn{\hat\Sigma_s}{\Sigma[s]} is the estimate of the noise covariace matrix
#'   obtained from a model with order \eqn{s}, \eqn{d(s)} denotes the number of parameters
#'   and \eqn{c(N)} is a "penalty" (depending on the sample size).
#'
#' For both `estorder_SVC` and `estorder_IVC` the (optional) parameter
#' `penalty` controls the penalty term \eqn{c(N)}.
#'
#' @section Further Notes:
#'
#' The actual computations are done by the helper routines detailed in [subspace helpers].
#'
#' The type of the [autocov()] object is irrelevenat
#' since the function always uses the slot `obj$gamma`.
#'
#' For  `keep_models==TRUE` the estimation procedure compute *all*
#' models even in the case of a Hankel singular value based selection criterion.
#'
#' @param obj Either a "time series" object (i.e `as.matrix(obj)`
#'     returns an \eqn{(N,m)}-dimensional numeric matrix)
#'     or an [autocov()] object (with \eqn{L} lags)
#'     which represents an (estimated) autocovariance function.
#'     The type of the `autocov` object is irrelevant since `est_stsp_ss` always uses the
#'     slot `obj$gamma` which contains the autocovariance function.
#' @param method Character string giving the method used to fit the model.
#' @param s.max (integer) maximum order of the state space model. If `NULL` a default value
#'        is chosen based on the sample size \eqn{N}, respectively based on the number of lags
#'        \eqn{L} of the ACF.
#' @param p (integer) number of block columns of the Hankel matrix (size of the "past"). If `NULL`
#'        then `p` is chosen by fitting a "long" AR model.
#' @param p.ar.max (integer) maximum order of the "long" AR model. If `NULL` a default choice is made.
#'        This parameter is only needed in the case `p=NULL`.
#' @param p.factor (integer) If `p=NULL`, then the number of block columns of the
#'        Hankel matrix is set to \eqn{p = p_f\hat{p}_{AIC}}{p = p[f]phat} where
#'        \eqn{p_f}{p[f]} is this parameter `p.factor` and
#'        \eqn{\hat{p}_{AIC}}{phat} is
#'        the (AIC) estimate of the order of the "long" AR model. See also [est_ar()].
#' @param extend_acf (boolean) If TRUE then the ACF is extended via an AR(p) model (MEST).
#' @param sample2acf (boolean) If `obj` is a data object and `sample2acf` is TRUE,
#'        then first the sample autocovariance function is computed and then used for the actual computations.
#' @param estorder function, used to select the order of the state space model.
#' @param keep_models (boolean) should the function return a list with estimated system of order 0:s.max?
#' @param mean_estimate Character string giving the method used to estimate the mean
#'        \eqn{\mu = E y_t}{\mu = E y[t]}.
#'        Default is to use the sample mean. See the details below.
#' @param n.obs Optional integer which gives the sample size \eqn{N}. This parameter is only used,
#'        when `obj` is an [autocov()] object. If `n.obs=NULL` then the
#'        slot `obj$n.obs` is used. Note that `obj$n.obs=NULL` or `obj$n.obs=Inf`
#'        refers to the case of a population autocovariance function,
#'        i.e. \eqn{N=\infty}. For a "time series" object the sample size is of course set to
#'        the number of observations, i.e. `n.obs = nrow(as.matrix(obj))`.
#'        The sample size \eqn{N} controls the computation of the default (maximum) orders
#'        and the estimation of the order of the state space model.
#' @param ... additional parameters, passed on to the order estimation routine.
#'
#' @return list with slots
#' \item{model}{a [stsp()] object, which represents the estimated state space model.}
#' \item{models}{either `NULL` (if `!keep_models`) or a list with the parameters
#'               of the estimated models with orders (`s=0:s.max+1`). This slot may e.g. be
#'               used to estimate the model order by some user defined model selection procedure.}
#' \item{s}{(integer) the estimate of the model order.}
#' \item{info}{list with information about the data and the design parameters of the estimation procedure.}
#' \item{stats}{((s.max+1)-by-5)-dimensional matrix with statistics of the (estimated) state space models.}
#' \item{y.mean}{estimate of the mean \eqn{\mu}.}
#'
#' @aliases CCA MEST AOKI
#'
#' @references
#' \insertRef{Aoki90}{RLDM}
#'
#' \insertRef{Bauer2001}{RLDM}
#'
#' \insertRef{DahlenScherrer2004}{RLDM}
#'
#' @export
#' @name subspace methods
#'
#' @examples
#' set.seed(3421) # in order to get reproducible results
#'
#' # create a "random", stable and miniphase state space model
#' m = 2 # number of outputs
#' s = 3 # number of states
#' s.max = 2*s
#' lag.max = max(4*s, 25)
#' n.obs = 1000
#'
#' model = r_model(tmpl_stsp_full(m, m, s),
#'                 bpoles = 1, bzeroes = 1,  sd = 0.5)
#' # scale sigma_L
#' diag(model$sigma_L) = 1
#'
#' # compute ACF
#' gam = autocov(model, lag.max = lag.max)
#'
#' # simulate data
#' data = sim(model, n.obs)
#'
#' # sample ACF
#' gam.sample = autocov(data$y, lag.max = lag.max, demean = FALSE)
#'
#' # AOKIs method ##############################################################
#' # reconstruct the true model from the population ACF
#' # "estimate" the order by the rank of the Hankel matrix
#' out = est_stsp_ss(gam, method = 'aoki', s.max = 2*s, estorder = estorder_rkH)
#'
#' # compute the ACF of the constructed model.
#' gam.hat = autocov(out$model, lag.max = lag.max)
#'
#' # check that the constructed model is equivalent to the original model
#' all.equal(dim(model$sys), dim(out$model$sys))
#' all.equal(gam, gam.hat)
#'
#'
#' # CCA based on the sample ###################################################
#' # estimate the order by a "singular value criterion"
#' out = est_stsp_ss(data$y, method = 'cca', sample2acf = FALSE, s.max = 2*s,
#'                   estorder = estorder_SVC)
#'
#' # compute the ACF of the constructed model.
#' gam.hat = autocov(out$model, lag.max = lag.max)
#'
#' all.equal(dim(model$sys), dim(out$model$sys)) # the estimated order is correct
#' all.equal(gam$gamma, gam.hat$gamma)           # but of course the estimated model is not perfect
#'
#' # CCA based on the sample ACF ###############################################
#' # estimate the order by an "information criterion"
#' out = est_stsp_ss(gam.sample, method = 'cca', s.max = 2*s,
#'                   estorder = estorder_IVC)
#'
#' # compute the ACF of the constructed model.
#' gam.hat = autocov(out$model, lag.max = lag.max)
#'
#' cat('s.hat=', dim(out$model$sys)[3], '\n') # the estimated order is s.hat=2, the true order is s=3!
#' all.equal(gam$gamma, gam.hat$gamma)        # relative error of the TRUE and the estimated ACF
#'
#' # alternatively, we may also use
#' out2 = est_stsp_ss(data$y, sample2acf = TRUE, mean_estimate = 'zero',
#'                    method = 'cca', s.max = 2*s, estorder = estorder_IVC)
#' all.equal(out$model, out2$model)
#'
#' # MEST algorithm #############################################################
#' # estimate the order by an "information criterion"
#' out = est_stsp_ss(gam.sample, method = 'cca', extend_acf = TRUE, s.max = 2*s,
#'                   estorder = estorder_IVC)
#'
#' # compute the ACF of the constructed model.
#' gam.hat = autocov(out$model, lag.max = lag.max)
#'
#' cat('s.hat=', dim(out$model$sys)[3], '\n') # the estimated order is s.hat=2, the true order is s=3!
#' all.equal(gam$gamma, gam.hat$gamma)        # relative error of the TRUE and the estimated ACF
#'
#' # make a plot of the ACFs
#' plot(gam, list(gam.sample, gam.hat), legend = c('TRUE', 'sample', 'subspace'))
#'
#' # reset seed
#' set.seed(NULL)
est_stsp_ss = function(obj, method = c('cca', 'aoki'),
                       s.max = NULL, p = NULL, p.ar.max = NULL, p.factor = 2,
                       extend_acf = FALSE, sample2acf = TRUE,
                       estorder = estorder_SVC, keep_models = FALSE,
                       mean_estimate = c('sample.mean', 'zero'), n.obs = NULL, ...) {
  method = match.arg(method)
  mean_estimate = match.arg(mean_estimate)

  if (!inherits(obj, 'autocov')) {
    y = try(as.matrix(obj))
    if ( inherits(y, 'try-error') || (!is.numeric(y)) || (!is.matrix(y)) ) {
      stop('could not coerce "obj" to a numeric matrix')
    }
    m = ncol(y)      # number of "outputs"
    n.obs = nrow(y)  # sample size (optional parameter n.obs is ignored)
    if (m == 0) stop('the "time series object" (obj) contains no data')
    if (n.obs < 3) stop('the sample size must be larger than 2')
    if (mean_estimate == 'sample.mean') {
      y.mean = colMeans(y)
      y = scale(y, center = y.mean, scale = FALSE)
    } else {
      y.mean = numeric(m)
    }
    lag.max = Inf
    gamma = NULL
  } else {
    gamma = obj$gamma
    d = dim(gamma)
    m = d[1]
    if (m == 0) {
      stop('the "autocov" object (obj) is empty.')
    }
    lag.max = d[3] - 1
    if (lag.max <= 2) {
      stop('the "autocov" object (obj) must contain at least 2 lags.',
           ' lag.max=', lag.max)
    }

    if (is.null(n.obs)) {
      n.obs = obj$n.obs
      if (is.null(n.obs)) n.obs = Inf
    }
    n.obs = as.numeric(n.obs)[1]
    if (is.finite(n.obs)) {
      n.obs = as.integer(n.obs)[1]
      if (n.obs < 3) stop('the sample size must be larger than 2')
    } else {
      n.obs = Inf
    }

    y.mean = rep(NA_real_, m)
    y = NULL
  }

  if ((is.null(gamma)) && (!sample2acf) && (method == 'aoki')) {
    stop('the "AOKI" method is only implemented for ACFs.')
  }

  # the "orders" p, p.factor, p.ar.max, ... must satisfy the following restrictions
  #
  # if extend_acf==TRUE (extrapolate ACF via AR(p) model): set ext = 2, else ext = 1
  #
  # 2 <= lag.max <= n.obs-1   (for a sample of size N = n.obs)  (=> n.obs >= 3)
  #
  # p.factor >= 1
  #
  # regression of y[t],...,y[t+p] on y[t-1],...,y[t-p]
  #   sample: n.obs - 2*p >= m*p                       => p <= n.obs/(m+2)
  #   ACF (Hankel matrix (p+1) x p): 2*p <= lag.max    => p <= lag.max/2 = ext*lag.max/2
  #      if extend_acf (extrapolate ACF via AR(p))     => p <= lag.max   = ext*lag.max/2
  #
  # p >= 1
  # p*m >= s.max                                       => p >= max(1, s.max/m)
  #
  # if p is undefined, then p is "estimated" via a long autoregression:
  #
  # AR model regression of y[t] on intercept and y[t-1],...,y[t-p]
  #   sample: n.obs - p >= p*m + 1                     => p.ar.max <= (n.obs-1)/(m+1)
  #   ACF:    p <= lag.max                             => p.ar.max <= lag.max
  #   default value for p.ar.max is                       p.ar.max = 10*log[10](n.obs)
  #
  #
  # p = p.factor*p.ar.hat <= p.factor*p.ar.max implies
  #                                         => p.factor*p.ar.max <= n.obs/(m+2)
  #                                         => p.factor*p.ar.max <= ext*lag.max/2
  #
  #                                         => p.factor*p.ar.max >= max(1, s.max/m)
  #
  # combining these restrictions (with p.factor >=1), we require:
  #                                         => p.ar.max <= (n.obs-1)/((m+2)*p.factor)
  #                                         => p.ar.max <= ext*lag.max/(2*p.factor)
  #
  #                                         => p.ar.max >= max(1, s.max/m)/p.factor

  ext = ifelse(extend_acf, 2, 1)
  if (is.null(s.max)) s.max = NA_integer_
  s.max = as.integer(s.max)[1]

  p.upper = floor(min(n.obs/(m+2), ext*lag.max/2))
  p.lower = ceiling(max(1, s.max/m, na.rm = TRUE))
  if (p.upper < p.lower) stop('the sample size or the number of lags is too small')

  if (is.null(p)) {
    # estimate the number of future/past lags (= size of the Hankel matrix)
    # by fitting an AR model
    p.ar.upper = floor(min((n.obs-1)/(m+2), ext*lag.max/2) / p.factor)
    p.ar.lower = ceiling(max(1, s.max/m, na.rm = TRUE) / p.factor)
    if (p.ar.upper < p.ar.lower) stop('the sample size or the number of lags is too small')
    if (is.null(p.ar.max)) {
      p.ar.max = max(p.ar.lower, min(p.ar.upper, round(10*log10(n.obs))))
    }
    # check p.ar.max
    if ((p.ar.lower > p.ar.max) || (p.ar.upper < p.ar.max)) {
      stop('the maximum AR order "p.ar.max" is too small or too large: ',
           'lower bound=',p.ar.lower, ', upper bound=', p.ar.upper, 'but p.ar.max=', p.ar.max)
    }

    ar_model = NULL
    p.ar.hat = -1
    # estimate the number of future/past lags (= size of the Hankel matrix)
    # by fitting an AR model
    if ( (is.null(gamma)) && (!sample2acf) ) {
      p.ar.hat = est_ar_ols(y, p.max = p.ar.max, penalty = 2/n.obs)$p
    } else {
      if (is.null(gamma)) {
        gamma = autocov(y, lag.max = 2*p.upper/ext, demean = FALSE)$gamma
      }
      ar_model = est_ar_yw(gamma, p.max = p.ar.max, penalty = 2/n.obs)
      p.ar.hat = ar_model$p
    }
    p = max(p.lower, as.integer(p.factor * p.ar.hat))
  }

  # check p
  p = as.integer(p)[1]
  if ( (p < p.lower) || (p > p.upper)) {
    stop('the number of (past) lags "p" is to small or too large for the given data.',
         ' p=', p)
  }


  # compute ACF if needed
  if ((is.null(gamma)) && (sample2acf)) {
    gamma = autocov(y, lag.max = 2*p/ext, demean = FALSE)$gamma
  }

  # exzend ACF if needed
  if ((!is.null(gamma)) && (extend_acf)) {
    if (is.null(ar_model) || (p.ar.hat != p)) {
      ar_model = est_ar_yw(gamma, p.max = p, penalty = -1)
    }
    gamma = rbind(bmatrix(gamma[,,1:(p+1),drop = FALSE], rows = c(1,3)),
                  matrix(0, nrow = m*p, ncol = m))
    a = bmatrix(ar_model$a[,,p:1,drop = FALSE], rows = 1)
    for (lag in ((p+1):(2*p))) {
      # print(c(lag,(m*lag+1),(m*(lag+1)),(m*(lag-p)+1),(m*lag)))
      gamma[(m*lag+1):(m*(lag+1)), ] = a %*% gamma[(m*(lag-p)+1):(m*lag), ,drop = FALSE]
    }
    dim(gamma) = c(m, 2*p+1, m)
    gamma = aperm(gamma,c(1,3,2))
  }

  # check s.max
  if (is.na(s.max)) {
    s.max = p*m
  }
  s.max = as.integer(s.max)[1]
  if (s.max < 0) stop('maximum state dimension "s.max" must be non negative.')
  if (s.max > (p*m)) {
    # this should not happen, due to the restrictions on p
    stop('the number of lags "p" is too small for the desired maximum state dimension:',
         ' p*m=', p*m, ' < s.max=', s.max)
  }

  n.valid = n.obs - 2*p

  if (is.null(gamma)) {
    out = est_stsp_cca_sample(y, s.max = s.max, p = p, estorder = estorder,
                              keep_models = keep_models,
                              mean_estimate = 'zero', ...)
  } else {
    if (method == 'aoki') {
      out = est_stsp_aoki(gamma, s.max = s.max, p = p, estorder = estorder,
                          keep_models = keep_models, n.obs = n.obs, ...)
    } else {
      out = est_stsp_cca(gamma, s.max = s.max, p = p, estorder = estorder,
                         keep_models = keep_models, n.obs = n.obs, ...)
    }
  }

  out$y.mean = y.mean
  # out$gamma = gamma
  return(out)
}


# Helpers ----

#' Solve a discrete time, algebraic Riccati equation
#'
#' This function solves the discrete time, algebraic Riccati equation
#' \deqn{X = AXA' + (M -AXC')(G-CXC')^{-1}(M-AXC')'}
#' where A is a square (s-by-s) matrix, M and C' are matrices of dimension (s-by-m)
#' and G is a square (positive definite) matrix of dimension (m-by-m). Given certain
#' regularity conditions (see the discussion below) the solution X computed by
#' `riccati` is a positive definite matrix (of size (s-by-s)).
#'
#' Here (within this package) this function is mainly used to construct a state space model if we have
#' given the autocovariance function of a process \eqn{(y_t)}{(y[t])} with a rational spectral density.
#' The ACF is represented by four matrices (A,M,C,G) as
#' \eqn{\gamma(0)=G} and \eqn{\gamma(k) = CA^{k-1}M} for \eqn{k>0}. This *realization problem*
#' is related to the so called *spectral factorization problem*. If the matrix \eqn{A} is stable,
#' the pair \eqn{(A,C)} is controllable, the pair \eqn{(A,M)} is controllable
#' (i.e. (A,M,C,G) is a "minimal" realization of the ACF) and if the spectral
#' density is positive definite (i.e. has no zeros) then the Riccati equation
#' has a (unique) solution \eqn{X} which is positive definite and where the matrix
#' \eqn{A - (M-AXC')(G-CXC')^{-1}C} is stable.
#' The process  \eqn{(y_t)}{(y[t])} then has a state space representation with parameters
#' \eqn{(A,B,C,D=I)}, where \eqn{B=(M-AXC')\Sigma^{-1}} and \eqn{\Sigma=(G-CXC')} is the
#' innovation covariance.
#'
#' The solution is computed via an (2m-by-2m) dimensional generalized eigenvalue problem,
#' which in turn is solved with a QZ decomposition.
#' See [QZ::qz.dgges()] and [QZ::qz.dtgsen()].
#' The eigenvalues of modulus less than one are
#' the eigenvalues of the matrix \eqn{A-BC}.
#'
#' In addition `riccati` is also used to compute the stochastically
#' balanced realization of a state space model,
#' see [rationalmatrices::grammians()] and [rationalmatrices::balance()].
#'
#' Note that this function is mainly used as a utility function and therefore no checks
#' on the given input parameters are performed.
#'
#' @param A (s-by-s) matrix
#' @param M (s-by-m) matrix
#' @param C (m-by-s) matrix
#' @param G (m-by-m) matrix
#' @param only.X boolean
#'
#' @return if (only.X) is TRUE then `riccati` just returns the solution X, otherwise a list
#'         with slots
#' \item{X}{the solution of the Ricatti equation}
#' \item{B}{the matrix \eqn{B = (M- AXC')\Sigma^{-1}}}
#' \item{sigma}{the matrix \eqn{\Sigma=G-CXC'}}
#' \item{lambda}{the (2m) vector of eigenvalues of the associated generalized eigenvalue problem. The first \eqn{m} entries
#'               are the eigenvalues of \eqn{A-BC}.}
#'
#' @export
#'
#' @examples
#' # create a "random" state space model, which satisfies the
#' # stability and the (strict) miniphase assumption
#' m = 2 # number of outputs
#' s = 4 # number of states
#'
#' model = r_model(template = tmpl_stsp_full(m, m, s),
#'                 bpoles = 1, bzeroes = 1, sd = 0.25)
#' # scale sigma
#' model$sigma_L = model$sigma_L / sqrt(sum(diag(model$sigma_L %*% t(model$sigma_L))))
#'
#' # extract the model parameter matrices
#' A = model$sys$A
#' B = model$sys$B
#' C = model$sys$C
#' sigma = model$sigma_L %*% t(model$sigma_L)
#'
#' # compute the variance of the state P = A P A' + B sigma B'
#' P = lyapunov(A, B %*% sigma %*% t(B))
#'
#' # variance of the output y[t]: G = C P C' + sigma
#' G = C %*% P %*% t(C) + sigma
#' # covariance between s[t+1] and y[t]: M = A P C' + B sigma
#' M = A %*% P %*% t(C) + B %*% sigma
#'
#' # check that P solves the Riccati equation P = APA' + (M -APC')(G-CPC')^{-1}(M-APC')'
#' all.equal(P,
#'           A %*% P %*% t(A) +
#'             (M - A %*% P %*% t(C)) %*% solve(G - C%*% P %*% t(C), t(M - A %*% P %*% t(C))))
#'
#' # compute P from the Riccati equation: P = APA' + (M -APC')(G-CPC')^{-1}(M-APC')'
#' out = riccati(A, M, C, G, only.X=FALSE)
#'
#' # check the solution
#' all.equal(P, out$X)
#' all.equal(B, out$B)
#' all.equal(sigma, out$sigma)
#'
#' # eigenvalues of (A-BC) ( <=> reciprocals of the zeroes of the system)
#' lambda = eigen(A - B %*% C, only.values=TRUE)$values
#' all.equal(sort(lambda), sort(out$lambda[1:s]))
#'
riccati = function(A, M, C, G, only.X = TRUE) {
  m = nrow(M) # number of states
  n = ncol(M) # number of outputs
  if (is.null(C)) C = matrix(0, nrow = n, ncol = m)

  # setup generalized eigenvalue problem (of dimension (2m-by-2m))
  # AA*X = BB*X*Lambda
  AA = diag(1,2*m)
  BB = AA

  tA = A - M %*% solve(G, C)
  AA[1:m,1:m] = t(tA)
  AA[(m+1):(2*m),1:m] = - M %*% solve(G, t(M))

  BB[(m+1):(2*m),(m+1):(2*m)] = tA
  BB[1:m,(m+1):(2*m)] = -t(C) %*% solve(G, C)

  # compute QZ decomposition of (AA,BB):
  # AA = Q*S*Z^H, BB = Q*T*Z^H where S,Q are upper triangular and Q,Z are unitary
  out = QZ::qz.dgges(AA, BB, vsl = FALSE, vsr = TRUE, LWORK = NULL)

  # throw a warning if qz.dgges 'fails'
  if (out$INFO<0) {
    warning('qz.dgges was not successful!')
  }

  # eigenvalues
  lambda = out$ALPHA / out$BETA

  # select eigenvalues with modulus < 1
  select = abs(lambda)<1

  # the number of stable eigenvalues should be equal to m (number of states)
  if (sum(select)!=m) {
    warning('number of eigenvalues with modulus less than one is not equal to m!')
  }

  # reorder the QZ decomposition such that the "stable" eigenvalues are on the top!
  out = QZ::qz.dtgsen(S=out$S, T=out$T, Q = out$Z, Z=out$Z, select = select, ijob = 0L,
                      want.Q = FALSE, want.Z = TRUE, LWORK = NULL, LIWORK = NULL)

  # throw a warning if qz.dtgsen fails
  if (out$INFO<0) {
    warning('qz.dtgsen was not successful!')
  }
  # recompute (ordered) eigenvalues
  lambda = out$ALPHA / out$BETA

  # compute X from the eigenvectors Z[,1:m]
  X = solve(t(out$Z[1:m,1:m,drop=FALSE]), t(out$Z[(m+1):(2*m),1:m,drop=FALSE])) # note X should be symmetric!
  X = Re( X + t(X) )/2

  # if (only.X) {
  #   return(list(X = X, info = out$INFO, lambda = lambda, AA = AA, BB = BB, out))
  # } else {
  #   sigma = G - C %*% X %*% t(C)
  #   B = t(solve(t(sigma),t(M - A %*% X %*% t(C))))
  #   return(list(X = X, B = B, sigma = sigma, info = out$INFO, lambda = lambda, AA = AA, BB = BB, out))
  # }
  if (only.X) {
    return(X)
  } else {
    sigma = G - C %*% X %*% t(C)
    B = t(solve(t(sigma),t(M - A %*% X %*% t(C))))
    return(list(X = X, B = B, sigma = sigma, lambda = lambda))
  }

}

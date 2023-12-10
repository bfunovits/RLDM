# Model comparison ----

#' Kullback–Leibler divergence
#'
#' Compute the Kullback-Leibler divergence between a "true" state space model and
#' an estimated state space model. The function only works for square systems,
#' where the "true" model is stable and the estimate is (strictly) miniphase.
#'
#' The KL divergence is computed as follows. Suppose \eqn{y_t = k(z) u_t}{y[t] = k(z) u[t]},
#' with \eqn{\mathbf{E} u_t u_t' = \Sigma}{E u[t] u[t]' = \Sigma} is the true model, and
#' let \eqn{y_t = h(z) u_t}{y[t] = h(z) u[t]}, with \eqn{\mathbf{E} u_t u_t'=\Omega}{E u[t] u[t]'=\Omega}
#' denote the estimate. W.l.o.g. we assume that the models are in innovation form,
#' i.e. \eqn{k(0) = h(0) = I}. The procedure computes the covariance matrix,
#' \eqn{\Delta = \mathbf{E} e_t e_t'}{\Delta = E e[t] e[t]'} say,
#' of the *one-step-ahead prediction errors* \eqn{e_t = h^{-1}(z) k(z) u_t}{e[t] = h^{-1} (z) k(z) u[t]}
#' and then the KL divergence
#' \deqn{
#' \mathrm{KL} = (1/2)(\mathrm{tr}(\Omega^{-1}\Delta) - m - \ln\det(\Omega^{-1}\Delta)))
#' }{
#' KL = (1/2)(tr(\Omega^{-1}\Delta) - m - ln det(\Omega^{-1}\Delta)))
#' }
#' Note that this procedure breaks down if the transfer function \eqn{h^{-1}(z) k(z)} is not stable.
#' Therefore the true models has to be stable and the estimated model has to strictly miniphase.
#'
#' @param model [stspmod()] object, true model
#' @param model_hat  [stspmod()] object, estimated model
#'
#' @return KL divergence.
#' @export
KL_divergence = function(model, model_hat) {

  if (!inherits(model, 'stspmod')) stop('input "model" must be a "stspmod" object')
  if (!inherits(model_hat, 'stspmod')) stop('input "model_hat" must be a "stspmod" object')
  d = dim(model$sys)
  d_hat = dim(model_hat$sys)
  if ((d[1] <= 0) || (d[2] != d[1])) stop('only square non empty models are supported')
  if (any (d[1:2] != d_hat[1:2])) stop('models are not compatible')
  m = d[1]

  k0_hat = model_hat$sys$D
  sigma_hat = tcrossprod(k0_hat %*% model_hat$sigma_L)
  # print(sigma_hat)

  mm = stspmod(sys = model_hat$sys^(-1) %r% model$sys, sigma_L = model$sigma_L)
  sigma = autocov(mm, lag.max = 0)$gamma
  dim(sigma) = c(m,m)
  sigma = k0_hat %*% sigma %*% t(k0_hat)
  # print(sigma)

  s = solve(sigma_hat, sigma)
  # print(s)
  ll = (1/2)*(sum(diag(s)) - ncol(s) - log(det(s)))
  return(ll)
}


#' Portmanteau Test for Serial Correlation
#'
#' Test whether the residuals of an estimated model are serially correlated. The test statistic is
#' \deqn{Q = N^2\sum_{k=1}^{K} (N-k)^{-1}\mbox{tr} (G_k G_0^{-1} G_k' G_0^{-1})}{
#'       Q = N^2\sum_{k=1}^{K} (N-k)^{-1} tr (G[k] G[0]^{-1} G[k]' G[0]^{-1})}
#' where \eqn{G_k}{G[k]} are the sample covariances of the residuals. Under the Null of a correctly
#' specified and estimated model the test statistic is asmyptotically Chi-squared distributed with
#' \eqn{Km^2-\kappa} degrees of freedom, where \eqn{\kappa} is the number of (free) parameters of the
#' model (class).
#'
#'
#' @param u        (N-by-m) matrix of residuals (or an object which may be coerced to
#'                 a matrix with `as.matrix(u)`).
#' @param lag.max  (integer) maximum number of lags.
#' @param n.par    (integer) number of parameters of the estimated model.
#'
#' @return Matrix with four columns ("`lags`" number of lags, "`df`" degrees of freedom,
#'         "`Q`" test statistics  and "`p`" p values).
#' @export
#'
#' @examples
#' u = matrix(rnorm(100*3), nrow = 100, ncol = 3)
#' pm_test(u, 4, 0)
pm_test = function(u, lag.max, n.par) {
  # check parameters
  u = try(as.matrix(u))
  if ( inherits(u, 'try-error') || (!is.numeric(u)) || (!is.matrix(u)) ) {
    stop('input "u" must be a data object which may be coerced to a matrix with "as.matrix(u)"')
  }

  n.obs = nrow(u)
  m = ncol(u)
  if ((m*n.obs)==0) stop('input "u" contains no data.')

  lag.max = as.integer(lag.max)[1]
  if ((lag.max < 1) || (lag.max >= n.obs)) {
    stop('maximum number of lags "lag.max" must be positive and smaller than the sample size "n.obs".')
  }
  n.par = as.integer(n.par)[1]
  if ((n.par < 0) || (n.par >= ((m^2)*lag.max))) {
    stop('number of free parameters "n.par" must be positive and smaller than "m^2 * lag.max".')
  }

  gamma = autocov(u, type = 'covariance', demean = FALSE, lag.max = lag.max)$gamma
  gamma0 = matrix(gamma[,,1], nrow = m, ncol = m)
  L = t(chol(gamma0))
  gamma = gamma[ , , -1, drop = FALSE]

  dim(gamma) = c(m,m*lag.max)
  gamma = forwardsolve(L, gamma) # gamma[k] -> L^-1 * gamma[k]
  dim(gamma) = c(m,m,lag.max)

  gamma = aperm(gamma, c(2,1,3)) # gamma[k] -> gamma[k]'

  dim(gamma) = c(m, m*lag.max)
  gamma = forwardsolve(L, gamma) # gamma[k] -> L^-1 * gamma[k]
  dim(gamma) = c(m,m,lag.max)

  lags = 1:lag.max
  Q = apply(gamma^2, MARGIN = 3, FUN = sum) / (n.obs-lags)
  Q = (n.obs^2) * cumsum(Q)
  df = m^2*lags - n.par
  i = which(df > 0)
  if (length(i) <= 0) {
    # this should not happen, see checks above
    stop('lag.max is too small compared to the number of parameters.')
  }
  Q = Q[i]
  df = df[i]
  lags = lags[i]
  p = stats::pchisq(Q, df = df, lower.tail = FALSE)
  test = cbind(lags, Q, df, p)
  rownames(test) = NULL
  return(test)
}

#' Compare Estimated Models
#'
#' This utility function computes a number of statistics which may be used to compare/evaluate a set of
#' estimated models. In particular the function returns the log Likelihood (`ll`),
#' the Akaike Information Criterion (`AIC`), the Bayes Information Criterion (`BIC`),
#' the Final Prediction Error (`FPE`) and the p-values of a Portmanteau test for serial correlation
#' of the residuals.
#'
#' The concentrated, conditional (scaled) log Likelihood
#' \deqn{ll = -(1/2)(m \ln(2\pi) + m + \ln\det S + 2 \ln\det (k_0)}{
#'       ll = -(1/2)(m ln(2\pi) + m + ln det S + 2 ln det (k[0])}
#' is computed with `ll(model, y, skip = 0, concentrated = TRUE)`, see [ll()].
#' Here \eqn{S} denotes the sample covariance of the residuals of the model.
#'
#' The information criteria are
#' \deqn{AIC = -2 ll + (2/N) \kappa}
#' \deqn{BIC = -2 ll + (\ln(N)/N) \kappa}{BIC = -2 ll + (ln(N)/N) \kappa}
#' where \eqn{\kappa} denotes the respective number of free parameters.
#'
#' The Final Prediction Error is
#' \deqn{FPE = \det(S)\frac{N+\kappa}{N-\kappa}}{FPE = det(S)(N+\kappa)/(N-\kappa)}
#'
#' For the portmanteau test, see [pm_test()]. If the number of lags is not
#' specified then the procedure choses a default value based on the sample size.
#' (This values is attached to the output as attribute).
#'
#' Note that the procedure (re)evaulates these measures, even if the estimates contain
#' this information (e.g. the residuals or the log likelihood may be stored in the corresponding list).
#' The reason is to have a common data set and a common evaluation procedure for all estimates.
#'
#' Typically the data `y` is the "estimation data set",
#' i.e. the data which has been used to estimate the models. However, one may also
#' pass a new data set to the procedure.
#'
#' @param estimates (named) list of estimates. Each slot should contain a list with slots
#'                  `$model` (the estimated model) and `$n.par` the corresponding
#'                  number of (free) parameters of the model (class).
#' @param y (N-by-m)-dimensional matrix with the observed data (or an object which may
#'          be coerced to a matrix with `as.matrix(y)`).
#' @param n.lags number of lags for the Portmantaeu test for serial correlation of the
#'               residuals, see also [pm_test()].
#'
#' @return Matrix with the computed statistics for the estimated models. This matrix has attributes
#'         `m`, `n.obs` and `n.lags`.
#' @export
compare_estimates = function(estimates, y, n.lags = NULL) {
  y = try(as.matrix(y))
  if ( inherits(y, 'try-error') || (!is.numeric(y)) || (!is.matrix(y)) ) {
    stop('input "y" must be a data object which may be coerced to a matrix with "as.matrix(y)"')
  }

  n.obs = nrow(y)
  m = ncol(y)
  if ((m*n.obs)==0) stop('input "y" contains no data.')

  n.estimates = length(estimates)
  which = c('#par', 'll', 'AIC', 'BIC', 'FPE', 'PM test')
  n.stats = length(which)

  stats = matrix(NA_real_, nrow = n.estimates, ncol = n.stats)
  names.estimates = names(estimates)
  if (is.null(names.estimates)) names.estimates = paste('estimate', 1:n.estimates, sep=' ')
  rownames(stats) = names.estimates
  colnames(stats) = which

  n.par = sapply(estimates, FUN = function(x) { x$n.par } )
  stats[,'#par'] = n.par
  if (is.null(n.lags)) {
    n.lags = max(1, ceiling(10*log10(n.obs)), ceiling(max(n.par)/(m^2)))
  }
  n.lags = as.integer(n.lags)[1]
  if ( (n.lags<1) || (n.lags>= n.obs) || (((m^2)*n.lags) <= max(n.par)) ) {
    stop('illegal number of lags for the PM test')
  }

  for (i in (1:n.estimates)) {
    stats[i, 'll'] = ll(estimates[[i]]$model, y, 'concentrated', skip = 0)
    stats[i, 'AIC'] = -2*stats[i, 'll'] + (2/n.obs)*n.par[i]
    stats[i, 'BIC'] = -2*stats[i, 'll'] + (log(n.obs)/n.obs)*n.par[i]
    u = solve_inverse_de(estimates[[i]]$model$sys, y)$u # residuals
    S = crossprod(u)/n.obs
    stats[i, 'FPE'] = det(S)*(n.obs + n.par[i])/(n.obs-n.par[i])
    junk = pm_test(u, n.lags, n.par[i])
    stats[i, 'PM test'] = junk[nrow(junk),'p']
  }

  attr(stats, 'm') = n.obs
  attr(stats, 'n.obs') = n.obs
  attr(stats, 'n.lags') = n.lags
  return(stats)
}


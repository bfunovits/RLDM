#' Autocovariance, Autocorelation and Partial Autocorrelation Function
#'
#' Compute respectively estimate the autocovariance, autocorrelation or partial autocorrelation function of a stationary process.
#'
#' The class of the input parameter "`obj`" determines the S3 method called and hence what is actually computed.
#'
#' **Population ACF:**
#'
#' If "`obj`" is an [armamod()] or [stspmod()] object then `autocov(obj, ...)` computes the ACF of the corresponding stationary process.
#'
#' Note however, that the function returns nonsense, if the model does not satisfy the *stability* condition.
#'
#' **Change the type of an ACF:**
#'
#' Calling `autocov(obj, type)`, where "`obj`" is an `autocov` object returns an ACF of the desired type.
#' E.g. if "`obj`" holds a partial autocorrelation function then `autocov(obj, type = 'covariance')` may be used to retrieve the corresponding autocovariance function.
#'
#' This is possible since the `autocov` object stores the "original" autocovariances in a slot named `gamma`.
#'
#' **Sample ACF:**
#'
#' The default S3 method estimates the ACF from given data.
#' It assumes that "`obj`" is a univariate or multivariate numeric time series object,
#' a (numeric) data frame or a (numeric) vector, respectively matrix
#' and then simply calls the function [stats::acf()] in the \pkg{stats} package to compute the sample autocovariance function.
#' If needed, then the corresponding sample autocorrelation, respectively sample partial autocorrelation function is computed (and returned).
#'
#' The syntax is quite analogous to [stats::acf()], so please consider the documentation of [stats::acf()] for more details.
#'
#' Note that \pkg{stats} stores autocovariance/autocorrelation functions as `(lag.max+1,m,m)` dimensional arrays, whereas \pkg{RLDM} uses `(m,m,lag.max+1)` dimensional arrays.
#'
#' The definition of partial autocorrelations used by [stats::acf()] differs from the definition used here, see e.g. \insertCite{Reinsel1997}{RLDM}.
#' Furthermore [stats::acf()] skips the lag zero partial autocorrelation coefficient
#' and thus the pacf computed by  [stats::acf()] is (lag.max,n,n) dimensional.
#'
#' The default choice for the number of lags is \eqn{10*log10(N/m)} where \eqn{N}
#' is the number of observations and \eqn{m} the number of series.
#' This number will be automatically limited to one less than the number
#' of observations in the series.
#'
#' @param obj either a [armamod()], [stspmod()], [autocov()] object
#'            or a "data" object.
#' @param type   character string giving the type of acf to be computed. Allowed values are
#'               "covariance" (the default), "correlation", or "partial". Will be partially matched.
#'               Note that the default value here is "covariance" whereas [stats::acf()]
#'               uses "correlation" as default.
#' @param ...    not used.
#' @param lag.max (integer) maximum lag.
#' @param na.action function to be called to handle missing values. [stats::na.pass()] can be used.
#' @param demean logical. Should the covariances be about the sample means?
#'
#' @return `autocov` object, i.e. a list with slots
#' \item{acf}{[rationalmatrices::pseries()] object, which stores the covariances (correlations).}
#' \item{type}{character string which indicates the type of the ACF.}
#' \item{gamma}{(m,m,lag.max+1) dimensional array which stores the autocovariance function.}
#' \item{names}{(m)-dimensional character vector or NULL. This optional slot stores the names
#'              for the components of the time series/process.}
#' \item{label}{character string or NULL.}
#' \item{n.obs}{integer or NULL. This slot stores the sample size.}
#'
#'
#' @seealso The autocovariance function of (V)ARMA processes may also be
#' computed by [stats::ARMAacf()] in the scalar case and by
#' [MTS::VARMAcov()] in the multivariate case (m > 1).
#'
#' As noted above the sample ACF is computed via the [stats::acf()] routine in the \pkg{stats} package.
#'
#' @references
#' \insertRef{Reinsel1997}{RLDM}
#'
#' @export
#'
#' @examples
#' model = stspmod(sys = stsp(A = c(0,0.2,1,-0.5), B = c(1,1,1,-1),
#'                            C = c(1,0,0,1)), sigma_L = diag(c(4,1)),
#'                 names = c('y1','y2'), label = 'test model')
#' g = autocov(model, lag.max=10)       # ACF
#' r = autocov(model, lag.max=10, type = 'correlation')  # autocorrelation function
#' r = autocov(g, type = 'correlation')                  # this gives the same result!
#' c = autocov(r, type = 'partial')     # partial autocorrelation function
#'
#' \dontrun{
#' # consider an equivalent VARMA model
#' model2 = impresp2varma(irf(model, lag.max = 20))$model
#' g2 = autocov(model2, lag.max = 10)
#' all.equal(g,g2) # of course both return the same ACF
#'
#' autocov(matrix(rnorm(100*2), nrow = 100))
#' autocov(stspmod(test_stsp(dim = c(2,2), s = 2), sigma_L = diag(2)))
#'
#' # generate a random sample with 100 observations and 3 outputs/series.
#' x = matrix(rnorm(100*3), nrow = 100, ncol = 3)
#'
#' # the covariance estimates are of course identical
#' stats_acfobj = stats::acf(x, type = 'covariance', demean = TRUE, plot = FALSE)
#' rldm_acfobj  = acf_estimate(x, type = 'covariance', demean = TRUE)
#' testthat::expect_equivalent(rldm_acfobj$gamma, aperm(stats_acfobj$acf,c(2,3,1)))
#'
#' # also the correlation estimates are identical
#' stats_acfobj = stats::acf(x, type = 'correlation', demean = TRUE, plot = FALSE)
#' rldm_acfobj  = acf_estimate(x, type = 'correlation', demean = TRUE)
#' testthat::expect_equivalent(rldm_acfobj$gamma, aperm(stats_acfobj$acf,c(2,3,1)))
#'
#' # However, the partial correlations dont match!
#' stats_acfobj = stats::acf(x, type = 'partial', demean = TRUE, plot = FALSE)
#' rldm_acfobj  = acf_estimate(x, type = 'partial', demean = TRUE)
#' testthat::expect_equivalent(rldm_acfobj$gamma[,,-1,drop=FALSE], aperm(stats_acfobj$acf,c(2,3,1)))
#' }
#'
autocov = function(obj, type, ...) {
  UseMethod("autocov", obj)
}

#' @rdname autocov
#' @export
autocov.default = function(obj, type=c('covariance','correlation','partial'), lag.max = NULL,
                           na.action = stats::na.fail, demean = TRUE, ...) {

  if (!is.null(lag.max)) {
    lag.max = as.integer(lag.max)[1]
    if (lag.max < 0) stop('negative lag.max')
  }
  type = match.arg(type)

  out_acf = try(stats::acf(obj, lag.max = lag.max, type = 'covariance', plot = FALSE,
                           na.action = na.action, demean = demean))

  if (inherits(out_acf, 'try-error')) stop('stats:acf failed, the input "obj" may not be supported.')

  gamma = aperm(out_acf$acf, c(2,3,1))
  n.obs = out_acf$n.used

  if (type == 'covariance') {
    acf = gamma
  }
  if (type == 'correlation') {
    # compute correlations
    m = dim(gamma)[1]
    # standard deviations
    s = sqrt(diag(matrix(gamma[,,1], nrow = m, ncol = m)))
    s = matrix(s, nrow = m, ncol = m)
    s2 = s * t(s)
    acf = gamma
    for (i in (1: dim(gamma)[3])) acf[,,i] = acf[,,i] / s2
  }
  if (type == 'partial') {
    acf = est_ar_dlw(gamma, penalty = -1)$partial
  }
  acf = structure(acf, class = c('pseries', 'ratm') )

  out = structure( list(acf = acf, type = type, gamma = gamma, names = NULL, label = NULL,  n.obs = n.obs),
                   class = c('autocov', 'rldm') )

  return(out)
}


#' @rdname autocov
#' @export
autocov.autocov = function(obj, type=c('covariance','correlation','partial'), ...) {
  type = match.arg(type)
  out = unclass(obj)
  gamma = out$gamma

  out$type = type

  if (type == 'covariance') {
    acf = gamma
  }
  if (type == 'correlation') {
    # compute correlations
    m = dim(gamma)[1]
    # standard deviations
    s = sqrt(diag(matrix(gamma[,,1], nrow = m, ncol = m)))
    s = matrix(s, nrow = m, ncol = m)
    s2 = s * t(s)
    acf = gamma
    for (i in (1: dim(gamma)[3])) acf[,,i] = acf[,,i] / s2
  }
  if (type == 'partial') {
    acf = est_ar_dlw(gamma, penalty = -1)$partial
  }
  acf = structure(acf, class = c('pseries', 'ratm') )
  out$acf = acf

  out = structure( out, class = c('autocov', 'rldm') )

  return(out)
}


#' @rdname autocov
#' @export
autocov.armamod = function(obj, type=c('covariance','correlation','partial'), lag.max = 12, ...) {
    lag.max = as.integer(lag.max[1])
    if (lag.max<0) stop('negative lag.max')
    type = match.arg(type)

    d = unname(dim(obj$sys))
    m = d[1]
    n = d[2]
    p = d[3]
    q = d[4]

    a = unclass(obj$sys$a)
    b = unclass(obj$sys$b)
    sigma = obj$sigma_L
    sigma = sigma %*% t(sigma)
    k = unclass(pseries(obj$sys, lag.max = q)) # impulse response

    # the ACF is computed via the Yule Walker equations (for j = 0, ... , p)
    # a[0] gamma[j] + a[1] gamma[j-1] + ... + a[p] gamma[j-p] =
    #         = E (b[0] eps[t] + b[1] eps[t-q] + ... + b[q] eps[t-q]) y[t-j]'

    # compute "right hand side" of the YW equations
    # E (b[0] eps[t] + b[1] eps[t-q] + ... + b[q] eps[t-q]) y[t-j]' =
    #    b[j] sigma k[0]' + ... + b[q] sigma k[q-j]'
    Ewy = array(0, dim = c(m, m, (max(p , q, lag.max)+1)) )
    for (j in (0:q)) {
      for (i in (j:q)) Ewy[,,j+1] = Ewy[,,j+1] + b[,,i+1] %*% sigma %*% t(k[,,i-j+1])
    }

    # construct the Yule-Walker equations in vectorized form:
    # note: vec(a[i] * gamma[j]) = (I x a[i]) vec(gamma[j]), where x stands for the Kronecker product
    # and:  vec(a[i] * gamma[-j]) = (I x a[i]) vec(gamma[-j]) = (I x a[i]) P vec(gamma[j])
    #       where P is a permutation matrix!

    # construct permutation matrix P
    junk = as.vector(t(matrix(1:(m^2), nrow = m, ncol = m)))
    P = numeric(m^2)
    P[junk] = 1:(m^2)

    # construct an array with the Kronecker products (I x a[i])
    A = array(0, dim = c(m^2, m^2, p+1))
    for (i in iseq(0, p)) A[,,i+1] = diag(1, nrow = m) %x% a[ , , i + 1]

    # left hand side of equation system
    L = array(0, dim = c(m^2, m^2, p+1, p+1))
    # right hand side of equation system
    R = array(0, dim = c(m, m, p+1))
    for (j in (0:p)) {
      # E y[t] t(y[t-j]) = ...
      R[ , , j+1] = Ewy[ , , j+1]

      L[ , , j+1, j+1] = A[ , , 1]
      for (i in (0:p)) {
        lag = j - i
        #      cat(j, i, lag,'\n')
        if (lag >= 0) {
          L[ , , j+1, lag+1] = A[ , , i+1]
        } else {
          L[ , , j+1, 1-lag] = L[ , , j+1, 1-lag] + A[, P, i+1]
        }
      }
    }
    # print(L)
    # print(R)
    L = bmatrix(L, rows = c(1,3), cols = c(2,4))
    R = as.vector(R)
    # print(cbind(L,R))

    gamma0 = solve(L, R)
    #  print(gamma0)
    gamma0 = array(gamma0, dim=c(m, m, p+1))

    if (p >= lag.max) {
      gamma = gamma0[ , , 1:(lag.max+1), drop=FALSE]
    } else {
      # extend acf for lags p+1,...lag.max
      # simply by two nested for loops!
      # however we have to transform the AR coefficients
      # a[0] y_t + a[1] y[t-1] + ...  ==> y[t] = (-a[0]^{-1} a[1])y[t-1] + ...
      for (i in iseq(1,p)) {
        a[,,i+1] = solve(a[,,1], -a[,,i+1])
      }

      gamma = dbind(d = 3, gamma0, Ewy[ , , (p+2):(lag.max+1), drop = FALSE])
      for (lag in ((p+1):lag.max)) {
        for (i in iseq(1,p)) {
          gamma[ , , lag+1] = gamma[ , , lag+1] + a[ , ,i+1] %*% gamma[ , , lag+1-i]
        }
      }
    }
    #  print(gamma)

    if (type == 'covariance') {
      acf = gamma
    }
    if (type == 'correlation') {
      # compute correlations
      m = dim(gamma)[1]
      # standard deviations
      s = sqrt(diag(matrix(gamma[,,1], nrow = m, ncol = m)))
      s = matrix(s, nrow = m, ncol = m)
      s2 = s * t(s)
      acf = gamma
      for (i in (1: dim(gamma)[3])) acf[,,i] = acf[,,i] / s2
    }
    if (type == 'partial') {
      acf = est_ar_dlw(gamma, penalty = -1)$partial
    }
    acf = structure(acf, class = c('pseries', 'ratm') )

    out = structure( list(acf = acf, type = type, gamma = gamma, names = NULL, label = NULL,  n.obs = NULL),
                     class = c('autocov', 'rldm') )
    return(out)
}

#' @rdname autocov
#' @export
autocov.stspmod = function(obj, type=c('covariance','correlation','partial'), lag.max = 12, ...) {

  lag.max = as.integer(lag.max)[1]
  if (lag.max < 0) stop('negative lag.max')
  type = match.arg(type)

  A = obj$sys$A
  B = obj$sys$B
  C = obj$sys$C
  D = obj$sys$D
  m = nrow(C)
  n = ncol(B)
  s = ncol(A)

  sigma = obj$sigma_L %*% t(obj$sigma_L)

  gamma = array(0, dim = c(m, m, lag.max+1))
  if (s == 0) {
    gamma[,,1] = D %*% sigma %*% t(D)
  } else {
    P = lyapunov(A, B %*% sigma %*% t(B), non_stable = 'stop')
    M = A %*% P %*% t(C) + B %*% sigma %*% t(D)

    gamma[,,1] = C %*% P %*% t(C) + D %*% sigma %*% t(D)
    for (i in iseq(1, lag.max)) {
      gamma[,,i+1] = C %*% M
      M = A %*% M
    }
  }

  if (type == 'covariance') {
    acf = gamma
  }
  if (type == 'correlation') {
    # compute correlations
    m = dim(gamma)[1]
    # standard deviations
    s = sqrt(diag(matrix(gamma[,,1], nrow = m, ncol = m)))
    s = matrix(s, nrow = m, ncol = m)
    s2 = s * t(s)
    acf = gamma
    for (i in (1: dim(gamma)[3])) acf[,,i] = acf[,,i] / s2
  }
  if (type == 'partial') {
    acf = est_ar_dlw(gamma, penalty = -1)$partial
  }
  acf = structure(acf, class = c('pseries', 'ratm') )

  out = structure( list(acf = acf, type = type, gamma = gamma, names = NULL, label = NULL,  n.obs = NULL),
                   class = c('autocov', 'rldm') )
  return(out)
}




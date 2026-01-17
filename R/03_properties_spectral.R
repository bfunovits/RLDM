# spectrald_methods.R ######################################
# defines the main methods for the 'spectrald' class


#' Spectral Density
#'
#' Compute the spectral density of an ARMA process or a process defined by a state space model.
#'
#' The spectral density of a stationary process with an absolutely summable
#' autocovariance function \eqn{(\gamma_j)}{(\gamma[j])} is given by
#' \deqn{
#' \Gamma(\lambda) = \frac{1}{2\pi}\sum_{j=-\infty}^{\infty} \gamma_j e^{-i\lambda j}.
#' }{
#' \Gamma(\lambda) = (1/2\pi) sum_{j=-\infty}^{\infty} \gamma[j] exp(-i\lambda j).
#' }
#'
#' For an ARMA process, or process defined by a state space model the
#' spectral density is equal to
#' \deqn{
#' \Gamma(\lambda) = \frac{1}{2\pi} K(\lambda) \Sigma K^*(\lambda)
#' }{
#' \Gamma(\lambda) = (1/2\pi) K(\lambda) \Sigma K^*(\lambda)
#' }
#' where \eqn{\Sigma} is the noise covariance, \eqn{K()} is the
#' frequency response of the model and \eqn{K^*(\lambda)} is the
#' Hermitean transpose of \eqn{K(\lambda)}.
#' See also [autocov()] and [freqresp()].
#'
#' Note that \eqn{\Gamma()} is (up to a factor \eqn{2\pi}) the discrete-time Fourier
#' transform (DTFT) of the autocovariance function and therefore the ACF \eqn{\gamma_j}{\gamma[j]}
#' may be reconstructed from the spectral density via the inverse DTFT
#' \deqn{
#' \gamma_j = \int_{-\pi}^{\pi} \Gamma(\lambda) e^{i\lambda j} d\lambda
#' }{
#' \gamma[j] = int_{-\pi}^{\pi} \Gamma(\lambda) e^{i\lambda j} d\lambda
#' }
#'
#' The S3 methods `spectrald.*` evaluate the spectral density function
#' on a grid of angular frequencies \eqn{\lambda_j = 2\pi j/N}{\lambda[j] = 2\pi j/N},
#' \eqn{j=0,\ldots,N-1}{j=0,\dots,N-1} and store the result
#' in a **spectrald** object.
#'
#'
#' There are several possible ways to specify the process. One may provide the ARMA (`armamod`)
#' respectively the state space model (`stspmod`), the autocovariance function (`autocov`) or
#' the impulse response function (`impresp`) which maps the noise to the outputs.
#' Note however, that if we have only given an `autocov` or `impresp` object then
#' the computed spectral density is only an approximation of the true spectral density since
#' only a finite number of covariances respectively impulse response coefficients are given.
#' The type of the autocovariance function ("covariances", "correlations" or "partial correlations")
#' is irrelevenat since the procedure alwayas uses the slot "`gamma`" which contains the
#' covariances.
#'
#' The default method `spectrald.default` assumes that `obj` is a "time series" object and
#' tries to coerce this object to a data matrix via `y = as.matrix(obj)`. In this case the
#' procedure computes the *periodogram* which is a simple estimate of the spectral density.
#' The periodgram may also be computed by the call `spectrald(autocov(obj, max.lag = n.obs-1))`,
#' i.e. by first computing the sample auto covariance function and then computing the corresponding
#' spectral density.
#'
#' Note that we use a different scaling than the `stats::[spectrum][stats::spectrum]` routine.
#'
#' @param obj [armamod()], [stspmod()], [autocov()],
#'            [impresp()] object or a "time series" object, i.e. an
#'            object which may be coerced to a data matrix by `y = as.matrix(obj)`.
#' @param n.f number of frequencies.
#' @param demean (logical) should the data be demeaned before computing the periodogram?
#' @param ... not used.
#'
#' @return `freqresp` object, i.e. a list with slots
#' \item{spd}{[rationalmatrices::zvalues()] object.}
#' \item{names}{(m)-dimensional character vector or NULL. This optional slot stores the names
#'              for the components of the time series/process.}
#' \item{label}{character string or NULL.}
#' \item{n.obs}{(optional) integer or NULL.}
#'
#' @export
#'
#' @examples
#' #' ### generate random 3-dimensional ARMA(1,1) model
#' # "bpoles = 1.1" implies that the poles have moduli larger than 1.1
#' # and therefore the impulse response coefficients decay with a rate (1.1)^k
#' arma_model = test_armamod(dim = c(3,3), degrees = c(1,1), bpoles = 1.1)
#'
#' # spectral density
#' spd = spectrald(arma_model)
#'
#' # compute the spectral density via the impulse response
#' spd1 = spectrald(impresp(arma_model, lag.max = 100))
#'
#' # since the impulse response quickly decays
#' # the "truncated" spectral density should be close to the true one
#' all.equal(spd, spd1)
#'
#' # compute the spectral density via the autocovariance function
#' spd1 = spectrald(autocov(arma_model, lag.max = 100))
#'
#' # since the ACF quickly decays
#' # the "truncated" spectral density should be close to the true one
#' all.equal(spd, spd1)
#'
#' # create an equivalent state space model
#' stsp_model = as.stspmod(arma_model)
#'
#' # of course the state space model gives the same spectrum
#' # as the original ARMA model
#' spd1 = spectrald(stsp_model)
#' all.equal(spd, spd1)
#'
spectrald = function(obj, n.f, ...) {
  UseMethod("spectrald", obj)
}

# Internal helper for spectrald methods using frequency response
spectrald_from_freqresp = function(obj, n.f) {
  n.f = as.integer(n.f)[1]
  if (n.f < 0) stop('the number of frequencies "n.f" must be a non negative integer')

  out = unclass(freqresp(obj, n.f))
  frr = out$frr
  frr = frr %r% (out$sigma_L / sqrt(2*pi))
  spd = frr %r% Ht(frr)

  structure(list(spd = spd, names = out$names, label = out$label),
            class = c('spectrald','rldm'))
}

#' @rdname spectrald
#' @export
spectrald.armamod = function(obj, n.f = 128, ...) {
  spectrald_from_freqresp(obj, n.f)
}

#' @rdname spectrald
#' @export
spectrald.stspmod = function(obj, n.f = 128, ...) {
  spectrald_from_freqresp(obj, n.f)
}


#' @rdname spectrald
#' @export
spectrald.autocov = function(obj, n.f = 128, ...) {
  n.f = as.integer(n.f)[1]
  if (n.f < 0) stop('the number of frequencies "n.f" must be a non negative integer')

  gamma = obj$gamma
  gamma[,,1] = gamma[,,1] / 2
  gamma = gamma/(2*pi)

  spd = dft_3D(gamma, n.f)
  spd = spd + Ht(spd)
  out = list(spd = spd, names = obj$names, label = obj$label)
  if (!is.null(obj$n.obs)) out$n.obs = obj$n.obs
  class(out) = c('spectrald', 'rldm')
  return(out)
}


#' @rdname spectrald
#' @export
spectrald.impresp = function(obj, n.f = 128, ...) {
  spectrald_from_freqresp(obj, n.f)
}


#' @rdname spectrald
#' @export
spectrald.default = function(obj, n.f = NULL, demean = TRUE, ...) {
  y = try(as.matrix(obj))
  if (inherits(y, 'try-error')) stop('could not coerce input parameter "obj" to matrix')

  n.obs = nrow(y)
  m = ncol(y)

  if (is.null(n.f)) n.f = n.obs
  n.f = as.integer(n.f)[1]
  if (n.f < 0) stop('the number of frequencies "n.f" must be a non negative integer')

  if (demean) {
    y.mean = apply(y, MARGIN = 2, FUN = mean)
    y = y - matrix(y.mean, nrow = n.obs, ncol = m, byrow = TRUE)
  }
  dim(y) = c(n.obs, 1, m)
  y = aperm(y, c(3,2,1))

  spd = dft_3D(y / sqrt(2*pi*n.obs), n.f)
  spd = spd %r% Ht(spd)

  out = list(spd = spd, names = colnames(y), label = NULL, n.obs = n.obs)
  class(out) = c('spectrald', 'rldm')

  return(out)
}


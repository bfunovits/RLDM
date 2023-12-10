# str.____ methods ##############################################################

#' Display the Structure of Objects
#'
#' @param object an object
#' @param ... not used
#'
#' @return invisible(NULL)
#'
#' @rdname str
#' @name str methods
#' @export
str.armamod = function(object, ...) {
  d = attr(object$sys, 'order')
  cat('ARMA model [',d[1],',',d[2],'] with orders p = ', d[3], ' and q = ', d[4], '\n', sep = '')
  return(invisible(NULL))
}

#' @rdname str
#' @name str methods
#' @export
str.stspmod = function(object, ...) {
  d = attr(object$sys, 'order')
  cat('state space model [',d[1],',',d[2],'] with s = ', d[3], ' states\n', sep = '')
  return(invisible(NULL))
}

#' @rdname str
#' @name str methods
#' @export
str.impresp = function(object, ...) {
  d = dim(object$irf)
  orth = FALSE
  if (d[2] > 0) {
    sigma = object$sigma_L
    sigma = sigma %*% t(sigma)
    orth = isTRUE(all.equal(sigma, diag(d[2])))
  }
  if (orth) {
    cat('Orthogonalized impulse response [',d[1],',',d[2],'] with ', d[3], ' lags\n', sep = '')
  } else {
    cat('Impulse response [',d[1],',',d[2],'] with ', d[3], ' lags\n', sep = '')
  }
  return(invisible(NULL))
}

#' @rdname str
#' @name str methods
#' @export
str.autocov = function(object, ...) {
  d = dim(object$acf)
  type = object$type
  if (type == 'covariance') {
    cat('Autocovariance function [',d[1],',',d[2],'] with ', d[3], ' lags\n', sep = '')
  }
  if (type == 'correlation') {
    cat('Autocorrelation function [',d[1],',',d[2],'] with ', d[3], ' lags\n', sep = '')
  }
  if (type == 'partial') {
    cat('Partial autocorrelation function [',d[1],',',d[2],'] with ', d[3], ' lags\n', sep = '')
  }
  return(invisible(NULL))
}

#' @rdname str
#' @name str methods
#' @export
str.fevardec = function(object, ...) {
  d = dim(object$vd)
  cat('Forecast error variance decompositon [',d[1],',',d[2],'] maximum horizon = ', d[3], '\n', sep = '')
  return(invisible(NULL))
}


#' @rdname str
#' @name str methods
#' @export
str.freqresp = function(object, ...) {
  d = dim(object$frr)
  cat('Frequency response [',d[1],',',d[2],'] with ', d[3], ' frequencies\n', sep = '')
  return(invisible(NULL))
}


#' @rdname str
#' @name str methods
#' @export
str.spectrald = function(object, ...) {
  d = dim(object$spd)
  cat('Spectral density [',d[1],',',d[2],'] with ', d[3], ' frequencies\n', sep = '')
  return(invisible(NULL))
}



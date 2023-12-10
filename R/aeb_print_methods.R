# RLDM - print.____() methods ##############################################################

#' Print Methods
#'
#' @param x RLDM object, i.e. a [armamod()], [rmfdmod()], [stspmod()],
#'          [impresp()], [autocov()], [freqresp()],
#'          [spectrum()] or [fevardec()] object.
#' @param digits (integer) if non `NULL` then correspondingly rounded numbers are printed,
#'        see [round()].
#' @param format (character string) selects specific output formats. Note that
#'        [stsp()] and [fevardec()] objects have no format option.
#'        The option `'character'` is only implemented for (V)ARMA models.
#' @param ... Further parameters are ignored.
#'
#' @return `invisible(x)`
#'
#' @rdname print
#' @name print methods
#'
#' @examples
#' # for VARMA models six different print formats are implemented ###################
#' m = armamod(test_lmfd(dim = c(2,2), degrees = c(1,1)), sigma_L = diag(2))
#' print(m, digits = 2, format = "i|jz")
NULL

#' @rdname print
#' @export
print.armamod = function(x, digits = NULL,
                         format = c('i|jz', 'i|zj', 'iz|j', 'zi|j', 'i|j|z', 'character'), ...) {
  if (!is.null(digits)) {
    digits = as.vector(as.integer(digits))[1]
  }

  format = match.arg(format)
  if ((format == 'character') && (is.complex(unclass(x$sys)))) {
    stop('the format option "character" is only implemented for ARMA models with real coefficients.')
  }

  label = ''
  if (!is.null(x$label)) label = paste(x$label, ': ', sep = '')

  d = attr(x$sys, 'order')
  m = d[1]
  n = d[2]
  p = d[3]
  q = d[4]

  cat(label, 'ARMA model [',d[1],',',d[2],'] with orders p = ', d[3], ' and q = ', d[4], '\n', sep = '')

  if ((m*m*(p+1)) > 0) {
    cat('AR polynomial a(z):\n')

    a = unclass(x$sys$a)

    # use the function print_3D (contained in rationalmatrices)
    dimnames(a) = list(paste('[', 1:m, ',]', sep = ''),
                       paste('[,', 1:m, ']', sep = ''),
                       paste('z^',0:p, sep = ''))
    print_3D(a, digits, format)
  }

  if ((m*n*(q+1)) > 0) {
    cat('MA polynomial b(z):\n')

    a = unclass(x$sys$b)

    # use the above defined internal function print_3D
    dimnames(a) = list(paste('[', 1:m, ',]', sep = ''),
                       paste('[,', 1:n, ']', sep = ''),
                       paste('z^',0:q, sep = ''))
    print_3D(a, digits, format)
  }

  if (n > 0) {
    cat('Left square root of noise covariance Sigma:\n')

    a = x$sigma_L
    dimnames(a) = list(paste('u[',1:n,']',sep = ''),paste('u[',1:n,']',sep = ''))
    print(a)
  }

  return(invisible(x))
}

#' @rdname print
#' @export
print.rmfdmod = function(x, digits = NULL,
                         format = c('i|jz', 'i|zj', 'iz|j', 'zi|j', 'i|j|z', 'character'), ...) {
  if (!is.null(digits)) {
    digits = as.vector(as.integer(digits))[1]
  }

  format = match.arg(format)
  if ((format == 'character') && (is.complex(unclass(x$sys)))) {
    stop('the format option "character" is only implemented for ARMA models with real coefficients.')
  }

  label = ''
  if (!is.null(x$label)) label = paste(x$label, ': ', sep = '')

  d = attr(x$sys, 'order')
  m = d[1]
  n = d[2]
  p = d[3]
  q = d[4]

  cat(label, 'RMFD model [',d[1],',',d[2],'] with orders p = ', d[3], ' and q = ', d[4], '\n', sep = '')

  if ((n*n*(p+1)) > 0) {
    cat('right factor polynomial c(z):\n')

    c = unclass(x$sys$c)

    # use the above defined internal function print_3D
    dimnames(c) = list(paste('[', 1:n, ',]', sep = ''),
                       paste('[,', 1:n, ']', sep = ''),
                       paste('z^',0:p, sep = ''))
    print_3D(c, digits, format)
  }

  if ((m*n*(q+1)) > 0) {
    cat('left factor polynomial d(z):\n')

    a = unclass(x$sys$d)

    # use the above defined internal function print_3D
    dimnames(a) = list(paste('[', 1:m, ',]', sep = ''),
                       paste('[,', 1:n, ']', sep = ''),
                       paste('z^',0:q, sep = ''))
    print_3D(a, digits, format)
  }

  if (n > 0) {
    cat('Left square root of noise covariance Sigma:\n')

    a = x$sigma_L
    dimnames(a) = list(paste('u[',1:n,']',sep = ''),paste('u[',1:n,']',sep = ''))
    print(a)
  }

  return(invisible(x))
}

#' @rdname print
#' @export
print.stspmod = function(x, digits = NULL, ...) {
  if (!is.null(digits)) {
    digits = as.vector(as.integer(digits))[1]
  }

  d = attr(x$sys, 'order')
  m = d[1]
  n = d[2]
  s = d[3]

  label = ''
  if (!is.null(x$label)) label = paste(x$label, ': ', sep = '')

  cat(label, 'state space model [', m, ',', n, '] with s = ', s, ' states\n', sep = '')

  a = unclass(x$sys)
  attr(a, 'order') = NULL
  if (length(a) == 0) {
    return(invisible(x))
  }

  # rounding digits
  if (!is.null(digits)) {
    a = round(a, digits)
  }

  snames = character(s)
  if (s > 0) snames = paste('s[',1:s,']',sep = '')
  xnames = character(m)
  if (m > 0) {
    if ( !is.null(x$names) && is.character(x$names) && is.vector(x$names) && (length(x$names) == m) ) {
      xnames = x$names
    } else {
      xnames = paste('x[',1:m,']',sep = '')
    }
  }
  unames = character(n)
  if (n > 0) unames = paste('u[',1:n,']',sep = '')

  rownames(a) = c(snames, xnames)
  colnames(a) = c(snames, unames)
  print(a)

  if (n > 0) {
    cat('Left square root of noise covariance Sigma:\n')

    a = x$sigma_L
    dimnames(a) = list(paste('u[',1:n,']',sep = ''),paste('u[',1:n,']',sep = ''))
    print(a)
  }

  invisible(x)
}

#' @rdname print
#' @export
print.impresp = function(x, digits = NULL,
                         format = c('i|jz', 'i|zj', 'iz|j', 'zi|j', 'i|j|z'), ...) {
  if (!is.null(digits)) {
    digits = as.vector(as.integer(digits))[1]
  }
  format = match.arg(format)

  d = dim(x$irf)
  m = d[1]
  n = d[2]
  lag.max = d[3]

  label = ''
  if (!is.null(x$label)) label = paste(x$label, ': ', sep = '')

  orth = FALSE
  if (n > 0) {
    sigma = x$sigma_L
    sigma = sigma %*% t(sigma)
    orth = isTRUE(all.equal(x$sigma_L, diag(n)))
  }

  if (orth) {
    cat(label, 'Orthogonalized impulse response [', m, ',', n, '] with ', lag.max, ' lags\n', sep = '')
  } else {
    cat(label, 'Impulse response [', m, ',', n, '] with ', lag.max, ' lags\n', sep = '')
  }

  if ((m*n*(lag.max+1)) > 0) {
    a = unclass(x$irf)

    # use the above defined internal function print_3D
    dimnames(a) = list(paste('[', 1:m, ',]', sep = ''),
                       paste('[,', 1:n, ']', sep = ''),
                       paste('lag=', 0:lag.max, sep = ''))
    print_3D(a, digits, format)
  }

  invisible(x)
}

#' @rdname print
#' @export
print.autocov = function(x, digits = NULL,
                         format = c('i|jz', 'i|zj', 'iz|j', 'zi|j', 'i|j|z'), ...) {
  if (!is.null(digits)) {
    digits = as.vector(as.integer(digits))[1]
  }
  format = match.arg(format)

  d = dim(x$acf)
  m = d[1]
  lag.max = d[3]

  label = ''
  if (!is.null(x$label)) label = paste(x$label, ': ', sep = '')

  if (x$type == 'covariance') {
    label = paste(label, 'Autocovariance function [',d[1],',',d[2],'] with ', d[3], ' lags', sep = '')
  }
  if (x$type == 'correlation') {
    label = paste(label, 'Autocorrelation function [',d[1],',',d[2],'] with ', d[3], ' lags', sep = '')
  }
  if (x$type == 'partial') {
    label = paste(label, 'Partial autocorrelation function [',d[1],',',d[2],'] with ', d[3], ' lags', sep = '')
  }

  n.obs = x$n.obs
  if (!is.null(n.obs)) {
    n.obs = as.integer(n.obs)[1]
  } else {
    n.obs = Inf
  }
  if ( is.finite(n.obs) ) {
    label = paste(label, ', sample size is ', n.obs, sep = '')
  }

  cat(label, '\n')

  if ((m*(lag.max+1)) > 0) {
    a = unclass(x$acf)

    # use the above defined internal function print_3D
    dimnames(a) = list(paste('[', 1:m, ',]', sep = ''),
                       paste('[,', 1:m, ']', sep = ''),
                       paste('lag=', 0:lag.max, sep = ''))
    print_3D(a, digits, format)
  }
  invisible(x)
}

#' @rdname print
#' @export
print.fevardec = function(x, digits = NULL,
                          format = c('i|jz', 'i|zj', 'iz|j', 'zi|j', 'i|j|z'), ...) {
  if (!is.null(digits)) {
    digits = as.vector(as.integer(digits))[1]
  }
  format = match.arg(format)

  vd = x$vd
  d = dim(vd)
  n = d[1]
  h.max = d[3]

  label = ''
  if (!is.null(x$label)) label = paste(x$label, ': ', sep = '')
  cat(label, 'Forecast error variance decompositon [', n, ',', n,'] maximum horizon = ', h.max, '\n', sep = '')

  if ((n*h.max)> 0) {
    names = x$names
    if (is.null(names)) {
      names = paste('y[', 1:n, ']', sep = '')
    }
    unames = paste('u[', 1:n, ']', sep = '')
    hnames = paste('h=', 1:h.max, sep='')

    dimnames(vd) = list(names, unames, hnames)
    print_3D(vd, digits, format)
  }

  invisible(x)
}

#' @rdname print
#' @export
print.freqresp = function(x, digits = NULL,
                          format = c('i|jz', 'i|zj', 'iz|j', 'zi|j', 'i|j|z'), ...) {
  if (!is.null(digits)) {
    digits = as.vector(as.integer(digits))[1]
  }
  format = match.arg(format)

  d = dim(x$frr)
  m = d[1]
  n = d[2]
  n.f = d[3]

  label = ''
  if (!is.null(x$label)) label = paste(x$label, ': ', sep = '')

  cat(label, 'Frequency response [',d[1],',',d[2],'] with ', d[3], ' frequencies\n', sep = '')

  if ((m*n*n.f) > 0) {
    a = unclass(x$frr)
    attr(a, 'z') = NULL
    f = (0:(n.f-1))/n.f

    # use the above defined internal function print_3D
    if ((format == 'i|jz') || (format == 'i|zj')) {
      dimnames(a) = list(paste('[', 1:m, ',]', sep = ''),
                         paste('[,', 1:n, ']', sep = ''),
                         paste('f[',1:n.f, ']', sep = ''))
    } else {
      dimnames(a) = list(paste('[', 1:m, ',]', sep = ''),
                         paste('[,', 1:n, ']', sep = ''),
                         paste('f=', round(f,3), sep = ''))
    }
    print_3D(a, digits, format)
  }

  invisible(x)
}

#' @rdname print
#' @export
print.spectrald = function(x, digits = NULL,
                           format = c('i|jz', 'i|zj', 'iz|j', 'zi|j', 'i|j|z'), ...) {
  if (!is.null(digits)) {
    digits = as.vector(as.integer(digits))[1]
  }
  format = match.arg(format)

  d = dim(x$spd)
  m = d[1]
  n = d[2]
  n.f = d[3]

  label = ''
  if (!is.null(x$label)) label = paste(x$label, ': ', sep = '')

  label = paste(label, 'Frequency response [',d[1],',',d[2],'] with ', d[3], ' frequencies', sep = '')

  n.obs = x$n.obs
  if (!is.null(n.obs)) {
    n.obs = as.integer(n.obs)[1]
  } else {
    n.obs = Inf
  }
  if ( is.finite(n.obs) ) {
    label = paste(label, ', sample size is ', n.obs, sep = '')
  }

  cat(label, '\n')

  if ((m*n.f) > 0) {
    a = unclass(x$spd)
    attr(a, 'z') = NULL
    f = (0:(n.f-1))/n.f

    # use the above defined internal function print_3D
    if ((format == 'i|jz') || (format == 'i|zj')) {
      dimnames(a) = list(paste('[', 1:m, ',]', sep = ''),
                         paste('[,', 1:n, ']', sep = ''),
                         paste('f[',1:n.f, ']', sep = ''))
    } else {
      dimnames(a) = list(paste('[', 1:m, ',]', sep = ''),
                         paste('[,', 1:n, ']', sep = ''),
                         paste('f=', round(f,3), sep = ''))
    }
    print_3D(a, digits, format)
  }

  invisible(x)
}


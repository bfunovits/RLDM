# impresp_methods.R ######################################
# defines the main methods for the 'impresp' class

#' Impulse Response Function
#'
#' Compute the (orthogonalized) impulse response function of a VARMA model or
#' a state space model. The impulse response coefficients are also called
#' *Power series parameters* of the system.
#'
#' The impulse response coefficients \eqn{(k_j \,|\, j \geq 0)}{(k[j], j \ge 0)}
#' define the map between the noise and the output process. If the model is stable then
#' the stationary solution of the ARMA system, respectively state space system, is given by
#' \deqn{
#' y_t = \sum_{j \geq 0} k_j u_{t-j}.
#' }{
#' y[t] = sum_{j \ge 0} k[j] u[t-j].
#' }
#' For a state space system the impulse response coefficients
#' are
#' \deqn{k_0 = D \mbox{ and }}{k[0] = D and}
#' \deqn{k_j = CA^{j-1}B \mbox{ for }j >0.}{k[j] = CA^{j-1}B for j>0.}
#' For an ARMA model the coefficients are (recursively) computed
#' by a comparison of coefficients in the equation
#' \deqn{
#' (a_0 + a_1 z + \cdots + a_p z^p)(k_0 + k_1 z + k_2 z^2 + \cdots ) = b_0 + b_1 z + \cdots + b_q z^q
#' }{
#' (a[0] + a[1] z + \dots + a[p] z^p)(k[0] + k[1] z + k[2] z^2 + \dots ) = b[0] + b[1] z + \dots + b[q] z^q
#' }
#'
#' The S3 methods `impresp.*` compute the coefficients \eqn{k_j}{k[j]} for
#' \eqn{j = 0,\cdots,N}{j = 0,\dots,N} and store the result, together with the left square root
#' (`sigma_L`) of the noise covariance \eqn{\Sigma}, in a **impresp** object.
#' `impresp` objects contain the complete information
#' about the underlying model, provided that the maximum lag \eqn{N} is large enough.
#' This means that one may reconstruct the underlying model from an impulse response object.
#'
#' @param obj [armamod()], [stspmod()] object or [impresp()] object.
#'            The last case may be used to transform the impulse response function to
#'            a different orthogonalization scheme.
#' @param lag.max Maximum lag of the impulse response coefficients. This parameter is ignored in
#'                the case that `obj` is an impresp object.
#' @param H An (n x n) (non singular) matrix which specifies a transformation of the noise.
#'   The noise \eqn{u_t}{u[n]} is transformed to \eqn{H^{-1}u_t}{H^{-1}u[t]} and
#'   the impulse response coefficients (\eqn{k_j \rightarrow k_j H}{k[j] -> k[j] H}) and the (left)
#'   square root of the noise covariance matrix (\eqn{L \rightarrow H^{-1}L}{L -> H^{-1}L}) are
#'   transformed correspondingly.
#'   \cr
#'   The default case `H=NULL` corresponds to the identity matrix (i.e. no transformation).
#'   \cr
#'   For `H='chol'`, the transformation matrix `H = t(chol(Sigma))`
#'   is determined from the Choleski decomposition of the noise covariance \eqn{\Sigma}.
#'   For `H='eigen'` the symmetric square root of
#'   \eqn{\Sigma} (obtained from the the eigenvalue decomposition of \eqn{\Sigma})
#'   is used. For `H='sigma_L'` the left square root of the noise covariance,
#'   which is stored in the object `obj`, is used. In these cases
#'   one obtains an *orthogonalized* impulse response function.
#'   Other orthogonalization schemes may be obtained by setting \eqn{H} to a
#'   suitable square root of \eqn{\Sigma}.
#'
#' @return
#' `impresp` object, i.e. a list with components
#' \item{irf}{`pseries` object.}
#' \item{sigma_L}{(n,n)-dimensional matrix which contains left square of noise covariance matrix.}
#' \item{names}{(n)-dimensional character vector or NULL}
#' \item{label}{character string or NULL}
#' @export
#'
#' @references
#' \insertRef{luet05}{RLDM}
#'
#' @examples
#' # IRF from state space model ################################################
#' model = stspmod(stsp(A = c(0,0.2,1,-0.5), B = c(1,1,1,-1),
#'                      C = c(1,0,0,1)),
#'                 sigma_L = matrix(c(4, 1, 1, 3), 2, 2),
#'                 names = c('y1','y2'), label = 'test model')
#'
#' # IRF
#' irf = impresp(model, lag.max=10)
#' irf
#'
#' # Orthogonalized IRF: Cholesky
#' irf_chol = impresp(model, lag.max = 10, H = 'chol')
#' irf_chol
#' print(irf_chol$sigma_L) # Sigma is (approximately equal to) the identity matrix
#'
#'
#' # IRF from VARMA model ################################################
#' model = armamod(sys = test_lmfd(dim = c(2,2), degrees = c(2,1)))
#'
#' irf = impresp(model)
#' print(irf, digits = 2, format = 'iz|j')
impresp = function(obj, lag.max, H) {
  UseMethod("impresp", obj)
}

# internal helper function which computes the transformation matrix "H"
make_H = function(type = c('chol','eigen','sigma_L'), sigma_L) {
  n = nrow(sigma_L)
  type = match.arg(type)

  # H = sigma_L
  if ( type == 'sigma_L' ) {
    H = sigma_L
  }

  # cholesky decomposition
  if ( type == 'chol' ) {
    sigma = sigma_L %*% t(sigma_L)
    H = t(chol(sigma))
  }

  # eigenvalue decomposition
  if ( type == 'eigen' ) {
    sigma = sigma_L %*% t(sigma_L)
    ed = eigen(sigma, symmetric = TRUE)
    H = ed$vectors %*% diag(x = sqrt(ed$values), nrow = n) %*% t(ed$vectors)
  }

  return(H)
}


#' @rdname impresp
#' @export
impresp.armamod = function(obj, lag.max = 12, H = NULL) {

  lag.max = as.integer(lag.max)[1]
  if (lag.max < 0) {
    stop('lag.max must be a non-negative integer.')
  }

  sys = obj$sys
  sigma_L = obj$sigma_L

  irf = pseries(sys, lag.max = lag.max)

  # Compute orthogonalized IRF #
  if ((!is.null(H)) && (nrow(sigma_L) > 0)) {
    if (is.character(H)) {
      H = make_H(H, sigma_L)
    }
    if ( (!is.numeric(H)) || (!is.matrix(H)) || (any(dim(H) != nrow(sigma_L))) ) {
      stop('H must be a numeric matrix of the same dimension as the noise.')
    }
    irf = irf %r% H
    sigma_L = solve(H, sigma_L)
  }

  out = structure(list(irf = irf, sigma_L = sigma_L, names = obj$names, label = obj$label),
                  class = c('impresp','rldm'))

  return(out)
}

#' @rdname impresp
#' @export
impresp.rmfdmod = function(obj, lag.max = 12, H = NULL) {

  lag.max = as.integer(lag.max)[1]
  if (lag.max < 0) {
    stop('lag.max must be a non-negative integer.')
  }

  sys = obj$sys
  sigma_L = obj$sigma_L

  irf = pseries(sys, lag.max = lag.max)

  # Compute orthogonalized IRF #
  if ((!is.null(H)) && (nrow(sigma_L) > 0)) {
    if (is.character(H)) {
      H = make_H(H, sigma_L)
    }
    if ( (!is.numeric(H)) || (!is.matrix(H)) || (any(dim(H) != nrow(sigma_L))) ) {
      stop('H must be a numeric matrix of the same dimension as the noise.')
    }
    irf = irf %r% H
    sigma_L = solve(H, sigma_L)
  }

  out = structure(list(irf = irf, sigma_L = sigma_L, names = obj$names, label = obj$label),
                  class = c('impresp','rldm'))

  return(out)
}

#' @rdname impresp
#' @export
impresp.stspmod = function(obj, lag.max = 12, H = NULL) {

  # code is completey identical to the code of 'impresp.armamod #

  lag.max = as.integer(lag.max)[1]
  if (lag.max < 0) {
    stop('lag.max must be a non-negative integer.')
  }

  sys = obj$sys
  sigma_L = obj$sigma_L

  irf = pseries(sys, lag.max = lag.max)

  # Compute orthogonalized IRF #
  if ((!is.null(H)) && (nrow(sigma_L) > 0)) {
    if (is.character(H)) {
      H = make_H(H, sigma_L)
    }
    if ( (!is.numeric(H)) || (!is.matrix(H)) || (any(dim(H) != nrow(sigma_L))) ) {
      stop('H must be a numeric matrix of the same dimension as the noise.')
    }
    irf = irf %r% H
    sigma_L = solve(H, sigma_L)
  }

  out = structure(list(irf = irf, sigma_L = sigma_L, names = obj$names, label = obj$label),
                  class = c('impresp','rldm'))

  return(out)
}


#' @rdname impresp
#' @export
impresp.impresp = function(obj, lag.max = NULL, H = NULL) {
  # code is almost identical to the code of 'impresp.armamod #

  sigma_L = obj$sigma_L
  irf = obj$irf

  # Compute orthogonalized IRF #
  if ((!is.null(H)) && (nrow(sigma_L) > 0)) {
    if (is.character(H)) {
      H = make_H(H, sigma_L)
    }
    if ( (!is.numeric(H)) || (!is.matrix(H)) || (any(dim(H) != nrow(sigma_L))) ) {
      stop('H must be a numeric matrix of the same dimension as the noise.')
    }
    irf = irf %r% H
    sigma_L = solve(H, sigma_L)
  }

  out = structure(list(irf = irf, sigma_L = sigma_L, names = obj$names, label = obj$label),
                  class = c('impresp','rldm'))

  return(out)
}


# Forecast Error Variance Decomposition given IRF ----

#' Forecast Error Variance Decomposition
#'
#' Computes the Forecast Errors Variance Decomposition from a given (orthogonalized)
#' impulse response function.
#'
#' @param obj [impresp()] object which represents the (orthogonalized) impulse response function.
#' @param h.max maximum forecast horizon. The default is one plus the number of lags of the `impresp`
#'              object.
#' @param H An (n x n) (non singular) matrix which renders the impulse response to a
#'   *orthogonalized* impulse response. The noise \eqn{u_t}{u[n]} is transformed to
#'   \eqn{H^{-1}u_t}{H^{-1}u[t]} and the impulse response coefficients
#'   (\eqn{k_j \rightarrow k_j H}{k[j] -> k[j] H}) and the (left)
#'   square root of the noise covariance matrix (\eqn{L \rightarrow H^{-1}L}{L -> H^{-1}L}) are
#'   transformed correspondingly.
#'   \cr
#'   The default case `H=NULL` corresponds to the identity matrix (i.e. no transformation).
#'   \cr
#'   For `H='chol'`, the transformation matrix `H = t(chol(Sigma))`
#'   is obtained from the Choleski decomposition of the noise covariance \eqn{\Sigma}.
#'   For `H='eigen'` the symmetric square root of
#'   \eqn{\Sigma} (obtained from the the eigenvalue decomposition of \eqn{\Sigma})
#'   is used. For `H='sigma_L'` the left square root of the noise covariance,
#'   which is stored in the object `obj`, is used.
#'   Other orthogonalization schemes may be obtained by setting \eqn{H} to a
#'   suitable square root of \eqn{\Sigma}.
#'   \cr
#'   The procedure checks whether the transformation yields an orthogonalized
#'   impulse response. If not, an error is thrown.
#'
#' @return
#' `fevardec` object, i.e. a list with components
#' \item{vd}{n-by-n-by-h.max array which contains the forecast error variance
#'           decomposition: `vd[i,j,h]` is the percentage of the variance of the
#'           h-step ahead forecast error of the i-th component due to the j-th
#'           orthogonalized shock.}
#' \item{v}{n-by-h.max matrix which contains the forecast error variances:
#'          `v[i,h]` is the variance of the h-step ahead forecast error for the
#'          i-th component.}
#' \item{names}{(m)-dimensional character vector}
#' \item{label}{character string or NULL}
#'
#' @export
#'
#' @examples
#' model = stspmod(sys = stsp(A = c(0,0.2,1,-0.5), B = c(1,1,1,-1),
#'                            C = c(1,0,0,1)),
#'                 sigma_L = t(chol(matrix(c(4,2,2,3),nrow=2))),
#'                 names = c('y1','y2'), label = 'test model')
#' fevardec(impresp(model, lag.max=10), H = 'chol', h.max = 5) %>% print(digits = 2, format = 'iz|j')
#' fevardec(impresp(model, lag.max=4, H = 'eigen'))            %>% print(digits = 2, format = 'iz|j')
fevardec = function(obj, h.max = NULL, H = NULL) {
  if (!inherits(obj,'impresp')) stop('"obj" must be an "impresp" object!')

  irf = obj$irf
  m = dim(irf)[1]
  n = dim(irf)[2]
  lag.max = dim(irf)[3]
  if ( (m*n*(lag.max+1)) == 0 ) stop('empty impulse response (m*n*(lag.max+1) = 0)')
  if (m != n) stop('non square impulse response')

  sigma_L = obj$sigma_L

  if (is.null(h.max)) h.max = (lag.max+1)
  h.max = as.integer(h.max)[1]
  if (h.max < 1) stop('h.max is zero!')
  if (h.max > (lag.max+1)) stop('"h.max" is too large!')

  # Compute orthogonalized IRF #
  if ((!is.null(H)) && (nrow(sigma_L) > 0)) {
    if (is.character(H)) {
      H = make_H(H, sigma_L)
    }
    if ( (!is.numeric(H)) || (!is.matrix(H)) || (any(dim(H) != nrow(sigma_L))) ) {
      stop('H must be a numeric matrix of the same dimension as the noise.')
    }
    irf = irf %r% H
    sigma_L = solve(H, sigma_L)
  }

  sigma = sigma_L %*% t(sigma_L)
  if (!isTRUE(all.equal(sigma, diag(n)))) stop('orthogonalization failed!')

  irf = unclass(irf)[ , , 1:h.max, drop = FALSE]

  vd = apply(irf^2, MARGIN = c(1,2), FUN = cumsum)  # vd[h,i,j] = sum[l=1:h] (irf[i,j,l]^2)
  v  = apply(vd, MARGIN = c(1,2), FUN = sum)        # v[h,i] = sum[j=1:n] (vd[h,i,j]) = variance of the i-th component
  #                                   of the h-step ahead error
  vd = vd / array(v, dim=c(h.max,n,n))
  vd = aperm(vd, c(2,3,1))
  v = t(v)

  fevardec = list(vd = vd, v = v,
                  names = obj$names, label = obj$label)

  class(fevardec) = 'fevardec'
  return(fevardec)
}





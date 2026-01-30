


#' sigma_L Structure
#'
#' Create templates for the left square root \eqn{L} of the noise covariance
#' matrix \eqn{\Sigma = LL'}. This means that \eqn{L} is parametrized as
#' \deqn{\mbox{vec}(L) = h + H \theta}{vec(L) = h + H \theta} with a
#' (\eqn{k}-dimensional) parameter vector \eqn{\theta}.
#'
#' The parameter `structure` has the following meaning
#' \describe{
#' \item{as_given}{Use the given parameter `sigma_L` to construct the template:
#'                 `NA` entries are considered as free and all other entries as fixed.}
#' \item{chol}{Set all entries of `sigma_L` above the diagonal to zero and then proceed as above.}
#' \item{symm}{First make `sigma_L` symmetric (`sigma_L = (sigma_L + t(sigma_L))/2`) and then use this `sigma_L` as template.
#'             However, `h`, `H` are constructed such that \eqn{h + H\theta} gives a symmetric matrix!
#'             Note that NAs overwrite fixed values, see examples.}
#' \item{identity}{Use the identity matrix as template. In this case there are no free parameters, i.e. \eqn{\theta} is an empty vector (vector with zero length).}
#' \item{full_normalized}{Ones on the diagonal, otherwise all parameters are free.} }
#'
#' @param sigma_L numeric (n x n) matrix, where the free entries are coded with NAs
#' @param structure character string, determines the "structure" of sigma_L, see the examples.
#'
#' @return List with slots
#'   \itemize{
#'   \item `h` (\eqn{n^2}-dimensional vector),
#'   \item `H` (\eqn{(n^2, k)}-dimensional matrix, where \eqn{k} denotes the number of  free/deep parameters) and
#'   \item `n.par` (integer) is the number of free/deep parameters (\eqn{=k}).
#'   }
#'
#' @export
#'
#' @examples
#' sigma_L = matrix(c(0, NA, 1, 0, 2, 3, NA, 1, 1), nrow = 3, ncol = 3)
#' sigma_L
#'
#' tmpl = tmpl_sigma_L(sigma_L, structure = 'as_given')
#' th = -(1:tmpl$n.par)
#' matrix(tmpl$h + tmpl$H %*% th, nrow = 3, ncol = 3)
#'
#' tmpl = tmpl_sigma_L(sigma_L, structure = 'chol')
#' th = -(1:tmpl$n.par)
#' matrix(tmpl$h + tmpl$H %*% th, nrow = 3, ncol = 3)
#'
#' tmpl = tmpl_sigma_L(sigma_L, structure = 'symm')
#' th = -(1:tmpl$n.par)
#' matrix(tmpl$h + tmpl$H %*% th, nrow = 3, ncol = 3)
#'
#' tmpl = tmpl_sigma_L(sigma_L, structure = 'identity')
#' tmpl$n.par # = 0
#' matrix(tmpl$h, nrow = 3, ncol = 3)
#'
#' tmpl = tmpl_sigma_L(sigma_L, structure = 'full_normalized')
#' th = -(1:tmpl$n.par)
#' matrix(tmpl$h + tmpl$H %*% th, nrow = 3, ncol = 3)
#'
tmpl_sigma_L = function(sigma_L, structure = c('as_given', 'chol', 'symm', 'identity', 'full_normalized')) {
  if ( (!is.numeric(sigma_L)) || (!is.matrix(sigma_L)) || (ncol(sigma_L) != nrow(sigma_L)) ) {
    stop('"sigma_L" is not a square, numeric matrix')
  }
  structure = match.arg(structure)

  n = ncol(sigma_L)
  if (n == 0) {
    h = numeric(0)
    H = matrix(0, nrow =0, ncol = 0)
    return(list(h = h, H = H, n.par = 0))
  }
  u = upper.tri(sigma_L, diag = FALSE)

  if (structure == 'identity') {
    sigma_L = diag(n)
  }
  if (structure == 'full_normalized') {
    sigma_L = diag(n)
    sigma_L[sigma_L == 0] = NA
  }
  if (structure == 'chol') {
    sigma_L[u] = 0
  }
  if (structure == 'symm') {
    sigma_L = ( sigma_L + t(sigma_L) ) /2
    sigma_L[u] = 0
  }

  # print(sigma_L)
  h = as.vector(sigma_L)
  ix = which(is.na(h))
  h[ix] = 0
  n.par = length(ix)
  # cat(n.par, ix)
  H = matrix(0, nrow = n^2, ncol = n.par)
  # print(H)
  if (n.par > 0) H[ix, ] = diag(n.par)

  if ((structure == 'symm') && (n > 1)) {
    i = matrix(1:(n*n), ncol = n, nrow = n)
    iU = i[u]
    iL = t(i)[u]
    h[iU] = h[iL]
    H[iU, ] = H[iL, ]
  }

  return(list(h = h, H = H, n.par = n.par))
}

#' Model Structures
#'
#' These tools define and implement **model structures** where the **model parameters** are represented
#' by an affine function of some **free** (**deep**) parameters. As an example consider
#' multivariate ARMA models. The AR coefficients \eqn{a_k}{a[k]}, the MA coefficients \eqn{b_k}{b[k]}
#' and the (left) square root of the noise covariance matrix, \eqn{L} say, are vectorized and stacked
#' into a (long) parameter vector as
#' \deqn{\pi = (\mbox{vec}(a_1)',\ldots,\mbox{vec}(a_p)',
#'              \mbox{vec}(b_1)',\ldots,\mbox{vec}(b_q)',\mbox{vec}(L)')'}{\pi =
#'              (vec(a[1])',\ldots,vec(a[p])',
#'               vec(b[1])',\ldots,vec(b[q])',vec(L)')'}
#' This parameter vector then is written as
#' \deqn{\pi = h + H\theta}
#' where \eqn{\theta} represents the **free** parameters. Of course the matrix \eqn{H} is assumed to have
#' full column rank. This parameterization scheme is quite flexible. In particular, ARMA and
#' state space models in **echelon form** may be represented by this scheme.
#' \cr
#' Templates and the related tools are mainly used for estimation and for the generation of (random) models
#' for simulations and testing.
#'
#' The functions `model2template` and `tmpl_***` generate model "templates" which
#' represent certain model structures, where the model parameters are affine functions of some
#' **free**, respectively **deep**, parameters.
#'
#' The template contains information about the model structure in the following slots
#' \enumerate{
#' \item `h, H` represent the vector \eqn{h} and the marix \eqn{H} as described above.
#'                   (The vector \eqn{\pi} of stacked model parameters is represented as
#'                   \eqn{\pi = h +H \theta} where \eqn{\theta} is a vector of *deep* parameters.)
#' \item `class=["armamod"|"rmfdmod"|"stspmod"]`: determines whether the template parametrizes
#'                   ARMA, RMFD or state space models.
#' \item `order`: an integer vector which contains the dimensions and orders of the model.
#'                   For ARMA and RMFD models `order = c(m,n,p,q)` and
#'                   for state space models `order = c(m,n,s)`.
#' \item `n.par`: the number of *free* parameters, i.e. the dimension of the
#'                   vector \eqn{\theta}.
#' \item `nu`: This optional slot contains the Kronecker indices \eqn{\nu}.
#' }
#'
#' **model2template:**
#'
#' The function [model2template()] takes an [armamod()], [rmfdmod()], or [stspmod()] object
#' where the free parameters are coded as `NA`'s, `NaN`'s or `Inf`'s and
#' constructs a corresponding model template.
#'
#' For the parametrization of the (left) square root, \eqn{L} say, of the noise covariance \eqn{\Sigma = LL'}
#' the following choices are possible: In the case `sigma_L = "as_given"` the
#' slot `model$sigma_L` of the given `model` is used to construct the template: `NA` entries
#' are considered as free and all other entries as fixed. For the choice `sigma_L = "chol"`
#' first all entries of `model$sigma_L` above the diagonal are set to zero and then
#' the template is constructed as above. In the case `sigma_L = "symm"`
#' the matrix  `model$sigma_L` is first replaced by a symmetric one and then the template
#' is constructed (according to the `NA`'s) such that the square root `L=sigma_L` is always
#' symmetric. The choice `sigma_L = "identity"` sets the matrix `L = sigma_L` to the
#' identity matrix.
#' Finally, the choice `sigma_L = "full_normalized"` sets the diagonal elements equal to ones and all other elements to NAs in `L = sigma_L`.
#'
#' **tmpl_***:**
#'
#' The functions `tmpl_***` implement the following model structures:
#' \describe{
#' \item{tmpl_arma_pq}{ARMA models ([armamod()]) with prescribed orders \eqn{(p,q}).}
#' \item{tmpl_arma_echelon}{ARMA models ([armamod()]) in echelon form, with given
#'                              Kronecker indices \eqn{\nu}.}
#' \item{tmpl_rmfd_pq}{RMFD models ([rmfdmod()]) with prescribed orders \eqn{(p,q}).}
#' \item{tmpl_rmfd_echelon}{RMFD models ([rmfdmod()]) in echelon form, with
#'                              Kronecker indices \eqn{\nu}. Note that for RMFD models the
#'                              Kronecker indices refer to the basis of the *column space*
#'                              of the Hankel matrix of the impulse response coefficients. }
#' \item{tmpl_stsp_full}{Fully parametrized state space models ([stspmod()])
#'                       with given state space dimension \eqn{s},
#'                       i.e. a models where each entry in the matrices \eqn{A,B,C}
#'                       is considered non-fixed.}
#' \item{tmpl_stsp_echelon}{State space models ([stspmod()]) in echelon form,
#'                          with given Kronecker indices \eqn{\nu}.}
#' \item{tmpl_state space_ar}{State space model representations ([stspmod()])
#'                           of AR models with given order \eqn{p}.
#'                           Here only the "square" case \eqn{m=n} is implemented.}
#' }
#'
#' For these model structures the impulse response (transfer function) is scaled such that the
#' \eqn{(m,n)}-dimensional lag zero coefficient, \eqn{k_0}{k[0]} say, is of the form
#' \describe{
#' \item{\eqn{m=n}}{\eqn{k_0}{k[0]} is the \eqn{m}-dimensional identity matrix.}
#' \item{\eqn{m<n}}{The first \eqn{m} columns of \eqn{k_0}{k[0]} form the
#'                  \eqn{m}-dimensional identity matrix and the remaining columns are zero.}
#' \item{\eqn{m>n}}{The first \eqn{n} rows of \eqn{k_0}{k[0]} form the
#'                  \eqn{n}-dimensional identity matrix and the remaining rows are *free*.}
#' }
#' For the parametrization of the (left) square root \eqn{L} of the noise covariance \eqn{\Sigma = LL'}
#' the following choices are possible: For `sigma_L="chol"` the matrix \eqn{L} is lower triangular
#' (all entries on and below the main diagonal are considered as free and the entries above the diagonal are zero).
#' For `sigma_L="symm"` the matrix \eqn{L} is symmetric (all entries on and below the main diagonal are considered
#' as free and the entries above the diagonal are such that \eqn{L=L'} holds).
#' For `sigma_L="identity"` the matrix \eqn{L} is *fixed* to the \eqn{n}-dimensional identity matrix.
#' For `sigma_L="full_normalized"` the diagonal elements of the matrix \eqn{L} are *fixed* to ones and all other elements are free.
#'
#' @param model [armamod()], [rmfdmod()] or [stspmod()] object.
#' @param m output dimension
#' @param n input dimension (= number of shocks = dimension of the noise process)
#' @param p order of the AR polynomial for ARMA models, respectively the order of the
#'          right factor polynomial \eqn{c(z)} in an RMFD model.
#' @param q order of the MA polynomial for ARMA models, respectively the order of the
#'          left factor polynomial \eqn{d(z)} in an RMFD model.
#' @param s state dimension for state space models.
#' @param nu vector of Kronecker indices. For ARMA models the Kronecker indices
#'           describe the basis rows and for RMFD models the basis columns
#'           of the Hankel matrix of the impulse response coefficients.
#' @param sigma_L (character string) determines the form of the (left) square root
#'           of the noise covariance \eqn{\Sigma}. The choice `"chol"` gives a
#'           lower triangular matrix, `"symm"` gives a symmetric matrix and
#'           `"identity"` corresponds to a fixed (identity) matrix.
#'           The  procedure `model2template` has an additional option `"as_given"`
#'           which means that the structure of the square root `sigma_L` is
#'           completely determined by the `sigma_L` slot of the given model.
#'
#' @return The functions `model2template` and `tmpl_***` return a model template.
#'
#' @seealso [r_model()], [fill_template()], [ll()], [ll_theta()] and other estimation procedures
#'
#' @examples
#'
#' # ######################################################
#' # construct a template from a model
#' # ######################################################
#'
#' # Let us consider scalar ARMA(5,1) models
#' # for quarterly data with a strong seasonal component.
#' # In order to have parsimonious models we want a[2]=a[3]=0:
#' model = armamod(lmfd(a = c(1,NA,0,0,NA,NA), b = c(1,NA)))
#' tmpl = model2template(model)
#'
#' # Let's see how the "free" parameters are mapped to the model parameters
#' print(cbind(tmpl$h, tmpl$H))
#' th = -(1:tmpl$n.par)
#' fill_template(th, tmpl)
#'
#' # Generate a random model with this structure
#' th0 = rnorm(tmpl$n.par, sd = 0.1)
#' model = fill_template(th0, tmpl)
#'
#' # Extract the "free" parameters from the model
#' th = extract_theta(model, tmpl)
#' all.equal(th, th0)
#'
#' # This model structure fixes sigma_L = 1.
#' # If we change sigma_L = 2 then the model does not fit to the template.
#' model$sigma_L = 2
#' # the default choice on_error = 'ignore', tells
#' # extract_theta to ignore this misfit:
#' th = extract_theta(model, tmpl, on_error = 'ignore')
#' # with on_error = 'warn' we get a warning and
#' # with on_error = 'stop' would throw an error.
#' th = extract_theta(model, tmpl, on_error = 'warn')
#' # We may also "ignore" sigma_L
#' th = extract_theta(model, tmpl, on_error = 'stop', ignore_sigma_L=TRUE)
#'
#' # If the orders/class of template and model does not fit
#' \dontrun{
#' model = armamod(lmfd(a = c(1,1), b = c(1,1)))
#' extract_theta(model, tmpl)
#' model = stspmod(stsp(D = 1))
#' extract_theta(model, tmpl)
#' }
#'
#' # ######################################################
#' # the parameter "sigma_L"
#' # ######################################################
#'
#' # consider a state space model (with 1 state) for a 3-dimensional process
#' model = stspmod(stsp(A = 1, B = c(1,0,0), C = c(1,1,1), D = diag(3)))
#'
#' # We may specify an arbitrary structure for the left square root (L = sigma_L)
#' # of the noise covariance Sigma. Any NA entry is considered as a "free" parameter.
#' L = matrix(c(0, NA, 1, 0, 2, 3, NA, 1, 1), nrow = 3, ncol = 3)
#' L
#' # L has 2 NA entries and thus we get a model structure with 2 free parameters.
#' model$sigma_L = L
#'
#' tmpl = model2template(model, sigma_L = 'as_given')
#' th = -(1:tmpl$n.par)
#' fill_template(th, tmpl)
#'
#' # The choice sigma_L = 'chol' forces L to be lower triangular.
#' # In the case considered here, we get a model structure with 1 free parameter.
#' tmpl = model2template(model, sigma_L = 'chol')
#' th = -(1:tmpl$n.par)
#' fill_template(th, tmpl)
#'
#' # The choice sigma_L = 'symm' forces L = sigma_L to be symmetric.
#' # In the case considered here we thus get a model structure with 2 free parameters.
#' tmpl = model2template(model, sigma_L = 'symm')
#' th = -(1:tmpl$n.par)
#' fill_template(th, tmpl)
#'
#' # The choice sigma_L = 'identity' set L equal to the identity matrix,
#' # i.e. sigma_L is fixed.
#' tmpl = model2template(model, sigma_L = 'identity')
#' th = numeric(0)
#' fill_template(th, tmpl)
#' tmpl$n.par # there are no free parameters: tmpl$n.par = 0
#'
#' # The choice sigma_L = 'full_normalized' sets the diagonal elements of L equal to ones,
#' # and leaves all other elements free.
#' tmpl = model2template(model, sigma_L = 'full_normalized')
#' th = -(1:tmpl$n.par)
#' fill_template(th, tmpl)
#'
#' @name model structures
#' @rdname model_structures
#' @export
model2template = function(model, sigma_L = c("as_given", "chol", "symm", "identity", "full_normalized")) {

  # Check inputs and obtain integer-valued parameters
  if ( !( inherits(model, 'armamod') || inherits(model, 'stspmod') || inherits(model, 'rmfdmod') ) ) {
    stop('model is not an "armamod", "rmfdmod", or "stspmod" object.')
  }
  order = unname(dim(model$sys))
  n = order[2]

  # Affine parametrisation
  h = as.vector(unclass(model$sys))
  ix = which(!is.finite(h))
  h[ix] = 0
  n.par = length(ix)
  H = matrix(0, nrow = length(h), ncol = n.par)
  H[ix,] = diag(n.par)

  # Set noise parameters
  sigma_L_structure = match.arg(sigma_L)
  sigma_L = model$sigma_L

  junk = tmpl_sigma_L(sigma_L, structure = sigma_L_structure)
  h = c(h, junk$h)
  H = bdiag(H, junk$H)

  return(list(h = h, H = H, class = class(model)[1],
              order = order, n.par = ncol(H)))
}

## LMFD Templates ----

#' @examples
#'
#' # ######################################################
#' # ARMA(p,q) models
#' # ######################################################
#'
#' m = 2 # output dimension
#' p = 1 # AR order
#' q = 1 # MA order
#'
#' # model structure with lower triangular sigma_L
#' tmpl = tmpl_arma_pq(m, n = m, p, q, sigma_L = "chol")
#' th = rnorm(tmpl$n.par)
#' th = -(1:tmpl$n.par)
#' fill_template(th, tmpl)
#'
#' # model structure with symmetric sigma_L
#' tmpl = tmpl_arma_pq(m, n = m, p, q, sigma_L = "symm")
#' fill_template(th, tmpl)
#'
#' # model structure with sigma_L = I
#' tmpl = tmpl_arma_pq(m, n = m, p, q, sigma_L = "identity")
#' # here the number of free paramaters is of course (by 3) smaller
#' # than for the above model structures!
#' fill_template(th[1:(length(th)-3)], tmpl)
#'
#' @rdname model_structures
#' @export
tmpl_arma_pq = function(m, n, p, q, sigma_L = c("chol", "symm", "identity", "full_normalized")) {

  m = as.integer(m)[1]
  n = as.integer(n)[1]
  p = as.integer(p)[1]
  q = as.integer(q)[1]
  if (min(m-1, n-1, p, q) < 0) stop('illegal dimensions/orders: m,n >= 1 and p,q >= 0 must hold')

  sigma_L = match.arg(sigma_L)

  order = c(m, n, p, q)

  sys = matrix(NA_real_, nrow = m, ncol = m*(p+1) + n*(q+1))
  # set a[0] to the identity matrix
  sys[ , 1:m] = diag(m)
  # the first min(m,n) rows of b[0] and a[0] coincide
  sys[1:min(m,n), (m*(p+1) + 1):(m*(p+1) + n)] = diag(x = 1, nrow = min(m,n), ncol = n)

  sL = matrix(NA_real_, nrow = n, ncol = n)

  # create a helper model
  sys = structure(sys, order = order, class = c('lmfd','ratm'))
  model = list(sys = sys, sigma_L = sL, names = NULL, label = NULL)
  model = structure(model, class = c('armamod', 'rldm'))

  # create template
  tmpl = model2template(model, sigma_L = sigma_L)

  return(tmpl)
}

#' @rdname model_structures
#' @export
tmpl_arma_echelon = function(nu, n = length(nu), sigma_L = c("chol", "symm", "identity", "full_normalized")) {
#' @examples
#' # Create a template
#' tmpl <- tmpl_arma_echelon()
#' tmpl
#' 
#' # Use the template with fill_template()
#' # filled <- fill_template(tmpl, theta = rnorm(tmpl$n.par))

  sigma_L = match.arg(sigma_L)
  nu = as.integer(nu)
  m = length(nu)
  if ( (n < 1) || (m < 1) ) stop('illegal dimension, (m,n) must be positive')
  if (min(nu) < 0) stop('Kronecker indices must be non negative')

  p = max(nu)
  order = c(m, n, p, p)

  # code the position of the basis rows of the Hankel matrix
  basis = rationalmatrices::nu2basis(nu)

  # coefficients of a(z) in reverse order!!!
  # a = [a[p],...,a[0]]
  a = matrix(0, nrow = m, ncol = (p + 1)*m)
  # b = [b[0],...,b[p]]
  b = matrix(0, nrow = m, ncol = (p + 1)*n)

  # code free entries with NA's
  for (i in (1:m)) {
    shift = (p-nu[i])*m
    k = nu[i]*m + i
    basis_i = basis[basis < k]
    # cat(i, shift, k, basis_i, '\n')
    a[i, k + shift] = 1
    a[i, basis_i + shift] = NA_real_
    b[i, iseq(n + 1, (nu[i]+1)*n)] = NA_real_   # i-th row has degree nu[i]
    if (i > n) b[i, 1:n] = NA_real_             # the last (m-n) rows of b[0] are free
  }
  # print(a)
  # reshuffle a -> a = [a[0],...,a[p]]
  dim(a) = c(m, m, p+1)
  a = a[, , (p+1):1, drop = FALSE]
  dim(a) = c(m, m*(p+1))

  # print(cbind(a, b))

  sL = matrix(NA_real_, nrow = n, ncol = n)

  # create a helper model
  sys = structure(cbind(a, b), order = order, class = c('lmfd','ratm'))
  model = list(sys = sys, sigma_L = sL, names = NULL, label = NULL)
  model = structure(model, class = c('armamod', 'rldm'))
  # print(model$sys)

  # create template
  tmpl = model2template(model, sigma_L = sigma_L)
  # add Kronecker indices
  tmpl$nu = nu

  # the left, upper (min(n,m), min(n,m)) blocks of a[0] and b[0] are equal!
  # matrix of linear indices
  i = matrix(1:(m*m), nrow = m, ncol = m)
  ia = as.vector(i[1:min(n,m), 1:min(n,m)])
  ib = ia + (m*m*(p+1))
  # print(rbind(ia,ib))
  tmpl$h[ib] = tmpl$h[ia]
  tmpl$H[ib, ] = tmpl$H[ia, ]

  return(tmpl)
}

## RMFD Templates ----

#' @examples
#'
#' # ######################################################
#' # RMFD(p,q) models y[t] = d(z) c(z)^-1 e[t]
#' # ######################################################
#'
#' m = 2 # output dimension
#' p = 1 # order of c(z)
#' q = 1 # order of d(z)
#'
#' # model structure with lower triangular sigma_L
#' tmpl = tmpl_rmfd_pq(m, n = m, p, q, sigma_L = "chol")
#' th = rnorm(tmpl$n.par)
#' th = -(1:tmpl$n.par)
#' fill_template(th, tmpl)
#'
#' # model structure with symmetric sigma_L
#' tmpl = tmpl_rmfd_pq(m, n = m, p, q, sigma_L = "symm")
#' fill_template(th, tmpl)
#'
#' # model structure with sigma_L = I
#' tmpl = tmpl_rmfd_pq(m, n = m, p, q, sigma_L = "identity")
#' # here the number of free paramaters is of course (by 3) smaller
#' # than for the above model structures!
#' fill_template(th[1:(length(th)-3)], tmpl)
#'
#' @rdname model_structures
#' @export
tmpl_rmfd_pq = function(m, n, p, q, sigma_L = c("chol", "symm", "identity", "full_normalized")) {

  m = as.integer(m)[1]
  n = as.integer(n)[1]
  p = as.integer(p)[1]
  q = as.integer(q)[1]
  if (min(m-1, n-1, p, q) < 0) stop('illegal dimensions/orders: m,n >= 1 and p,q >= 0 must hold')

  sigma_L = match.arg(sigma_L)

  order = c(m, n, p, q)

  sys = matrix(NA_real_, nrow = n*(p+1) + m*(q+1), ncol = n)
  # set c[0] to the identity matrix
  sys[1:n, ] = diag(n)
  # set d[0] to a "quasi identity" matrix
  sys[(n*(p+1) + 1):(n*(p+1) + min(m,n)), 1:n] = diag(x = 1, nrow = min(m,n), ncol = n)

  sL = matrix(NA_real_, nrow = n, ncol = n)

  # create a helper model
  sys = structure(sys, order = order, class = c('rmfd','ratm'))
  model = list(sys = sys, sigma_L = sL, names = NULL, label = NULL)
  model = structure(model, class = c('rmfdmod', 'rldm'))

  # create template
  tmpl = model2template(model, sigma_L = sigma_L)

  return(tmpl)
}





#' @rdname model_structures
#' @export
tmpl_rmfd_echelon = function(nu, m = length(nu), sigma_L = c("chol", "symm", "identity", "full_normalized")) {
#' @examples
#' # Create a template
#' tmpl <- tmpl_rmfd_echelon()
#' tmpl
#' 
#' # Use the template with fill_template()
#' # filled <- fill_template(tmpl, theta = rnorm(tmpl$n.par))

  # (m,n) transfer function k = d(z) c^(-1)(z) with degrees deg(c) = p, deg(d) = q
  sigma_L = match.arg(sigma_L)
  nu = as.integer(nu)
  n = length(nu)
  if ( (n < 1) || (m < 1) ) stop('illegal dimension, (m,n) must be positive')
  if (min(nu) < 0) stop('Kronecker indices must be non negative')

  p = max(nu)
  order = c(m, n, p, p)

  # code the position of the basis columns of the Hankel matrix
  basis = rationalmatrices::nu2basis(nu)

  # coefficients of c(z) in reverse order!!!
  # c = [c[p]',...,c[0]']'
  c = matrix(0, nrow = n*(p+1), ncol = n)
  # d = [d[0]',...,d[p]']'
  d = matrix(0, nrow = m*(p+1), ncol = n)

  # code free entries with NA's
  for (i in (1:n)) {
    shift = (p-nu[i])*n
    k = nu[i]*n + i
    basis_i = basis[basis < k]
    # cat(i, shift, k, basis_i, '\n')
    c[k + shift, i] = 1
    c[basis_i + shift, i] = NA_real_
    d[iseq(m + 1, (nu[i] + 1)*m), i] = NA_real_   # i-th column has degree nu[i]
    d[iseq(n + 1, m), i] = NA_real_               # the last (m-n) rows of b[0] are free
  }
  # print(c)
  # reshuffle c -> c = [c[0]',...,c[p]']'
  dim(c) = c(n, p+1, n)
  c = c[, (p+1):1, , drop = FALSE]
  dim(c) = c(n*(p+1), n)

  # print(rbind(c, d))

  sL = matrix(NA_real_, nrow = n, ncol = n)

  # create a helper model
  sys = structure(rbind(c, d), order = order, class = c('rmfd','ratm'))
  model = list(sys = sys, sigma_L = sL, names = NULL, label = NULL)
  model = structure(model, class = c('rmfdmod', 'rldm'))
  # print(model$sys)

  # create template
  tmpl = model2template(model, sigma_L = sigma_L)
  # add Kronecker indices
  tmpl$nu = nu

  # the first min(n,m) rows of d[0] and c[0] are equal!
  # matrix of linear indices
  i = matrix(1:((n+m)*(p+1)*n), nrow = (n+m)*(p+1), ncol = n)
  ic = as.vector(i[1:min(n,m), 1:n])
  id = as.vector(i[(n*(p+1)+1):(n*(p+1)+min(n,m)), 1:n])
  # print(rbind(ia,ib))
  tmpl$h[id] = tmpl$h[ic]
  tmpl$H[id, ] = tmpl$H[ic, ]

  return(tmpl)
}

## State Space Templates ----

#' @rdname model_structures
#' @export
tmpl_stsp_full = function(m, n, s, sigma_L = c("chol", "symm", "identity", "full_normalized")) {
#' @examples
#' # Create a template
#' tmpl <- tmpl_stsp_full()
#' tmpl
#' 
#' # Use the template with fill_template()
#' # filled <- fill_template(tmpl, theta = rnorm(tmpl$n.par))

  m = as.integer(m)[1]
  n = as.integer(n)[1]
  s = as.integer(s)[1]
  if (min(m-1, n-1, s) < 0) stop('illegal dimensions/orders: m,n >= 1 and s >= 0 must hold')

  sigma_L = match.arg(sigma_L)

  order = c(m, n, s)

  sys = matrix(NA_real_, nrow = (m+s), ncol = (n+s))
  sys[(s+1):(s+min(m,n)),(s+1):(s+n)] = diag(x = 1, nrow = min(m,n), ncol = n)

  sL = matrix(NA_real_, nrow = n, ncol = n)

  # create a helper model
  sys = structure(sys, order = order,  class = c('stsp','ratm'))
  model = list(sys = sys, sigma_L = sL, names = NULL, label = NULL)
  model = structure(model, class = c('stspmod', 'rldm'))

  # create template
  tmpl = model2template(model, sigma_L = sigma_L)

  return(tmpl)
}

#' @rdname model_structures
#' @export
tmpl_stsp_ar = function(m, p, sigma_L = c("chol", "symm", "identity", "full_normalized")) {
#' @examples
#' # Create a template
#' tmpl <- tmpl_stsp_ar()
#' tmpl
#' 
#' # Use the template with fill_template()
#' # filled <- fill_template(tmpl, theta = rnorm(tmpl$n.par))

  m = as.integer(m)[1]
  p = as.integer(p)[1]
  if (min(m-1, p) < 0) stop('illegal dimensions/orders: m >= 1 and p >= 0 must hold')

  sigma_L = match.arg(sigma_L)

  n = m
  s = m*p
  order = c(m, n, s)

  if (s == 0) {
    sys = diag(m)
  } else {
    C = matrix(NA_real_, nrow = m, ncol = s)
    A = rbind(matrix(0, nrow = m, ncol = s),
              diag(x = 1, nrow = m*(p-1), ncol = s))
    B = diag(x = 1, nrow = s, ncol = m)
    D = diag(m)
    sys = rbind( cbind(A,B), cbind(C, D))
  }

  sL = matrix(NA_real_, nrow = n, ncol = n)

  # create a helper model
  sys = structure(sys, order = order,  class = c('stsp', 'ratm'))
  model = list(sys = sys, sigma_L = sL, names = NULL, label = NULL)
  model = structure(model, class = c('stspmod', 'rldm'))

  # create template
  tmpl = model2template(model, sigma_L = sigma_L)

  # take care of the restriction A[1:m,] = C
  if (p > 0) {
    i = matrix(1:((m+s)*(m+s)), nrow = m+s, ncol = m+s)
    iA = i[1:m, 1:s]
    iC = i[(s+1):(s+m), 1:s]
    tmpl$H[iA, ] = tmpl$H[iC, ]
  }

  return(tmpl)
}

#' @examples
#'
#' # ######################################################
#' # state space models in echelon form
#' # ######################################################
#' nu = c(3,2,4)   # Kronecker indices
#' m = length(nu)  # number of outputs/inputs
#' tmpl = tmpl_stsp_echelon(nu = nu)
#'
#' # generate a random vector of parameters.
#' # Note that "tmpl$n.par" contains the number free parameters.
#' th = rnorm(tmpl$n.par)
#'
#' # generate a model according to this structure with the parameters th
#' model = fill_template(th, tmpl)
#' print(model)
#'
#' # we can extract the free parameters from this given model
#' all.equal(th, extract_theta(model, tmpl, on_error = 'stop'))
#'
#' # check the impulse response
#' k = impresp(model, lag.max = 2*sum(nu) + 1)
#'
#' # the lag zero coeffcient k[0] is equal to the identity
#' all.equal(unclass(k$irf)[,,1], diag(m))
#'
#' # check the Kronecker indices
#' all.equal(rationalmatrices::pseries2nu(k$irf), nu)
#'
#' @rdname model_structures
#' @export
tmpl_stsp_echelon = function(nu, n = length(nu), sigma_L = c("chol", "symm", "identity", "full_normalized")) {

  # (m,n) transfer function C(I z^-1 -A)^-1 + D
  sigma_L = match.arg(sigma_L)
  nu = as.integer(nu)
  m = length(nu)
  if ( (n < 1) || (m < 1) ) stop('illegal dimension, (m,n) must be positive')
  if (min(nu) < 0) stop('Kronecker indices must be non negative')
  s = sum(nu) # state space dimension

  D = diag(x = 1, nrow = m, ncol = n)
  if (m > n) D[(n+1):m, ] = NA_real_

  if (s == 0) {
    sys = structure(D, order = c(m, n, s),  class = c('stsp','ratm'))
  } else {
    basis = nu2basis(nu)
    AC = matrix(0, nrow = s + m, ncol = s)
    dependent = c(basis + m, 1:m)
    for (i in (1:length(dependent))) {
      d = abs(basis-dependent[i])
      if (min(d) == 0) {
        # dependent[i]-th row is in basis
        j = which(d == 0)
        AC[i, j] = 1
      } else {
        j = which(basis < dependent[i])
        AC[i, j] = NA_real_
      }
    }
    B = matrix(NA_real_, nrow = s, ncol = n)
    sys = structure(cbind( AC, rbind(B,D)), order = c(m, n, s),  class = c('stsp','ratm'))
  }

  sL = matrix(NA_real_, nrow = n, ncol = n)

  # create a helper model
  model = list(sys = sys, sigma_L = sL, names = NULL, label = NULL)
  model = structure(model, class = c('stspmod', 'rldm'))
  # print(model$sys)

  # create template
  tmpl = model2template(model, sigma_L = sigma_L)
  # add Kronecker indices
  tmpl$nu = nu

  return(tmpl)
}

# Local Model Structures ----

# DDLC parametrization
#


#' Local Model Structures
#'
#' Parametrization for "local" model classes, in particular, "Data Driven Local Coordinates"
#' as detailed in \insertCite{McKelveyHelmerssonRibarits2004}{RLDM} and
#' \insertCite{RibaritsDeistlerHanzon2005}{RLDM}.
#'
#' The function `tmpl_DDLC` and `tmpl_GRAM` construct model templates which describe
#' models in a neighborhood of a given reference model.
#' In a first step the reference state space model is transformed to \eqn{D=I} and eventually
#' (depending on the parameter `"balance"`) balanced.
#'
#' state space models are described by a quadruple \eqn{(A,B,C,D=I)} of matrices which may be
#' embedded into an \eqn{(s^2+2ms)}-dimensional euclidean space. Note that the parameter matrices
#' are not uniqely determined from the ACF or the spectral density of the process, i.e. there
#' is an inherent non identifiablity problem. For minimal models the "equivalence class"
#' of models, which represent the same ACF is given by the set of all models which
#' may be obtained by a state transformation
#' \eqn{(A,B,C,D) \rightarrow (TAT^{-1}, TB, CT^{-1}, D)}{(A,B,C,D) -> (TAT^{-1}, TB, CT^{-1}, D)}.
#'
#' The *DDLC* parametrization now considers models, \eqn{(A,B,C,D=I)}, which
#' are contained in the \eqn{2ms}-dimensional subspace, which is *orthogonal* to
#' the \eqn{s^2}-dimensional tangent space of the set of equivalent models.
#'
#' The routine `tmpl_GRAM` considers the \eqn{2ms}-dimensional subspace, where
#' models close to the reference models are "approximately" balanced.
#'
#' Both schemes may fail for "non-generic" models. `tmpl_DDLC` issues a
#' warning message and `tmpl_GRAM` throws an error, in cases where the
#' \eqn{2ms}-dimensional subspace is not well defined.
#'
#' Note that also the parametrization of the left square root `L=sigma_L` of the
#' noise covariance is "local", i.e. `th = 0` corresponds to the (balanced) reference model.
#'
#' @param model [stspmod()] object, which represents a state space model. Only the case
#'           \eqn{m = n > 0} is implemented, i.e. the output process and the noise process must
#'           be of the same dimension.
#' @param balance (character string) For `balance = "lyapunov"` or
#'           `balance = "minimum phase"` the reference model is first balanced by
#'           the respective scheme.
#' @param sigma_L (character string) determines the form of the (left) square root
#'           of the noise covariance \eqn{\Sigma}. The choice `"chol"` gives a
#'           lower triangular matrix, `"symm"` gives a symmetric matrix and
#'           `"identity"` corresponds to the (fixed) identity matrix.
#'
#'
#'
#' @seealso For the computation of Grammians and for balancing of state space models,
#'          see [rationalmatrices::balance()].
#'
#' @references
#' \insertRef{McKelveyHelmerssonRibarits2004}{RLDM}
#'
#' \insertRef{RibaritsDeistlerHanzon2005}{RLDM}.
#'
#'
#'
#' @return Model template, i.e. a list with slots
#' \item{`h`}{\eqn{((m+s)^2 + m^2)}-dimensional vector,}
#' \item{`H`}{\eqn{((m+s)^2 + m^2, k)}-dimensional matrix,}
#' \item{`class = "stspmod"`}{}
#' \item{`order = c(m,m,s)`}{ and }
#' \item{`n.par`}{number of free parameters \eqn{=k}.}
#'
#' See also [model structures()] for more details.
#'
#' @name local model structures
#' @rdname local_model_structures
#' @export
#'
#' @examples
#' # create a random state space model with m outputs and s states
#' m = 3
#' s = 6
#' tmpl = tmpl_stsp_full(m, n = m, s, sigma_L = 'symm')
#' model = r_model(tmpl, bpoles = 1.1, bzeroes = 1.1, sd = 1/s)
#' model                              # note that sigma_L is symmetric
#' model$sigma_L %*% t(model$sigma_L) # noise covariance Sigma
#'
#' # tmpl_DDLC #############################################
#'
#' # create a DDLC parametrization of a neighborhood of this model
#' tmpl = tmpl_DDLC(model, balance = 'lyapunov', sigma_L = 'chol')
#' # for th = 0, we get the original model (in balanced form)
#' model = fill_template(numeric(tmpl$n.par), tmpl)
#' model                                # note that sigma_L is lower triangular
#' model$sigma_L %*% t(model$sigma_L)  # however Sigma is the same as above
#'
#' #' apply a "small" state transformation T = (diag(s)+eps*X)
#' eps = sqrt(.Machine$double.eps)
#' sys = model$sys
#' d_sys = state_trafo(sys, diag(s) + matrix(rnorm(s^2, sd = eps), nrow = s, ncol = s))
#' d_pi = (as.vector(unclass(d_sys) - unclass(sys)))/eps
#' # The vector d_pi is (close to) an element of the tangent space
#' # of the set of models, which are generated by a state transformation
#' # of the reference model
#'
#' # by construction d_pi is (close to) orthogonal to tmpl$H
#' max(abs(d_pi %*% tmpl$H[1:((m+s)^2), , drop = FALSE]))
#'
#' # the tmpl_DDLC routine may fail in some exceptional cases
#' m = 1
#' s = 3
#' model = stspmod(sys = stsp(A = matrix(0, nrow = s, ncol = s),
#'                            B = matrix(rnorm(m*s), nrow = s, ncol = m),
#'                            C = matrix(rnorm(m*s), nrow = m, ncol = s),
#'                            D = diag(m)),
#'                 sigma_L = diag(m))
#'
#' # For this model "tmpl_DLLC" issues a warning.
#' junk = tmpl_DDLC(model, sigma_L = 'chol', balance = 'none')
#'
#' # tmpl_GRAM #############################################
#' model = fill_template(numeric(tmpl$n.par), tmpl)
#'
#' tmpl = tmpl_GRAM(model, sigma_L = 'chol')
#' model = fill_template(numeric(tmpl$n.par), tmpl)
#'
#' # check grammians
#' gr = grammians(model$sys, 'lyapunov')
#' P = gr$P
#' Q = gr$Q
#' # P=Q=diag() should hold!
#' print(round(cbind(P, P-Q), 6))
#'
#' # now consider a model close to the reference model
#' d_th = rnorm(tmpl$n.par, sd = eps)
#' d_model = fill_template(d_th, tmpl)
#' d_sys = d_model$sys
#' gr = grammians(d_sys, 'lyapunov')
#' d_P = gr$P - P
#' d_Q = gr$Q - Q
#'
#' # the "disturbed" system should still be approximately balanced!
#' print(round(cbind(d_P, d_P - d_Q)/eps, 6) )
tmpl_DDLC = function(model,
                     balance = c('none', 'lyapunov', 'minimum phase'),
                     sigma_L = c("chol", "symm", "identity")) {
  if (!inherits(model, 'stspmod')) stop('model is not an "stspmod" object.')
  sigma_L_structure = match.arg(sigma_L)
  balance = match.arg(balance)

  sys = model$sys
  sigma_L = model$sigma_L
  order = dim(sys)
  m = order[1]
  n = order[2]
  s = order[3]

  if  ((m < 1) || (m != n)) stop('only the case m = n > 0 is implemented!')

  # normalize D -> D = I
  D = sys$D
  sigma_L = D %*% sigma_L
  D = try(solve(D))
  if (inherits(D, 'try-error')) stop('"D" is singular!')
  sys = sys %r% D

  if ((balance != 'none') && (s > 0)) {
    # balance state space system
    gr = grammians(sys, which = balance)
    sys = balance(sys, gr, truncate = FALSE)$obj
  }

  sys = unclass(sys)
  A = sys[iseq(1,s), iseq(1,s), drop = FALSE]
  B = sys[iseq(1,s), iseq(s+1,s+m), drop = FALSE]
  C = sys[iseq(s+1,s+m), iseq(1,s), drop = FALSE]
  sys[iseq(s+1,s+m), iseq(s+1,s+m)] = diag(m)

  h = as.vector(sys)

  if (s > 0) {
    # action of a state space trafo
    # A -> (I + E) A (I + E)^-1 ~ A + EA - AE
    # B -> (I + E) B            ~ B + EB
    # C -> C (I + E)^-1         ~ C - CE
    # construct a ((s^2+ 2*s*m) x (s^2)) matrix T, which spans the tangent space
    # of the equivalence classes
    T = rbind( (t(A) %x% diag(s)) - (diag(s) %x% A),
               t(B) %x% diag(s),
               -diag(s) %x% C )
    junk = svd(T, nu = nrow(T), nv = 0)
    # check rank of the matrix T
    d = junk$d
    if ( d[s^2] < (1e-8)*d[1] ) {
      warning('The tangent space of the equivalence class does not have dimension s^2=',
              s^2,' (sv[1]=', d[1], ', sv[', s^2, ']=', d[s^2], ')')
    }
    H = junk$u[, (s^2 + 1):nrow(T), drop = FALSE]
    H = rbind(H, matrix(0, nrow = m*m, ncol = ncol(H)))  # D = I is fixed!
    # cat(c(m, n, s), ' ABCD:', dim(ABCD), 'T: ', dim(T), 'U: ', dim(junk$u), 'H:', dim(H), '\n')
  } else {
    H = matrix(0, nrow = length(h), ncol = 0)
  }

  # we still have to reshuffle the rows of H
  i = matrix(1:((m+s)*(m+s)), nrow = m+s, ncol = m+s)
  #     vec(A), vec(B), vec(C), vec(D)
  i = c(i[iseq(1,s), iseq(1,s)], i[iseq(1,s), iseq(s+1,s+m)],
        i[iseq(s+1,s+m), iseq(1,s)], i[iseq(s+1,s+m), iseq(s+1,s+m)])
  # cat(i, '\n')
  H[i, ] = H

  # print(cbind(h, H))

  # add sigma_L parameters

  if (sigma_L_structure == 'identity') {
    sigma_L = diag(m)
  }
  if (sigma_L_structure == 'chol') {
    sigma_L = t(chol( sigma_L %*% t(sigma_L) ))
  }
  if (sigma_L_structure == 'symm') {
    junk = eigen( sigma_L %*% t(sigma_L), symmetric = TRUE )
    sigma_L = junk$vectors %*% diag(sqrt(junk$values), nrow = m, ncol = m) %*% t(junk$vectors)
  }
  junk = tmpl_sigma_L(matrix(NA_integer_, nrow = m, ncol = m), structure = sigma_L_structure)
  junk$h = as.vector(sigma_L)
  # print(cbind(junk$h, junk$H))
  h = c(h, junk$h)
  H = bdiag(H, junk$H)
  # print(cbind(h, H))
  # print(sigma_L)

  return(list(h = h, H = H, class = 'stspmod',
              order = order, n.par = ncol(H)))
}

#' @rdname local_model_structures
#' @export
tmpl_GRAM = function(model,
#' @examples
#' # Create a template
#' tmpl <- tmpl_GRAM()
#' tmpl
#' 
#' # Use the template with fill_template()
#' # filled <- fill_template(tmpl, theta = rnorm(tmpl$n.par))
                     balance = c('lyapunov', 'minimum phase'),
                     sigma_L = c("chol", "symm", "identity")) {
  if (!inherits(model, 'stspmod')) stop('model is not an "stspmod" object.')
  sigma_L_structure = match.arg(sigma_L)
  balance = match.arg(balance)
  if (balance == 'minimum phase') stop('only the case "balance=lyapunov" has been implemented so far.')

  sys = model$sys
  sigma_L = model$sigma_L
  order = dim(sys)
  m = order[1]
  n = order[2]
  s = order[3]

  if  ((m < 1) || (m != n)) stop('only the case m = n > 0 is implemented!')

  # normalize D -> D = I
  D = sys$D
  sigma_L = D %*% sigma_L
  D = try(solve(D))
  if (inherits(D, 'try-error')) stop('"D" is singular!')
  sys = sys %r% D

  if (s > 0) {
    # balance state space system
    gr = grammians(sys, which = 'lyapunov')
    sys = balance(sys, gr, truncate = FALSE)$obj
  }

  sys = unclass(sys)
  A = sys[iseq(1,s), iseq(1,s), drop = FALSE]
  B = sys[iseq(1,s), iseq(s+1,s+m), drop = FALSE]
  C = sys[iseq(s+1,s+m), iseq(1,s), drop = FALSE]
  sys[iseq(s+1,s+m), iseq(s+1,s+m)] = diag(m)

  h = as.vector(sys)

  if (s > 0) {
    # linearize grammians!
    i = matrix(1:((m+s)^2), nrow = (m+s), ncol = (m+s))

    dX = matrix(0, nrow = s^2, ncol = (s+m)^2)
    dYY = matrix(0, nrow = s^2, ncol = (s+m)^2)

    dX_ = diag(s^2) # dA
    dX[, i[iseq(1,s), iseq(1,s)]] = dX_
    # testthat::expect_equivalent(as.vector(A), as.vector(dX %*% h))

    # d(B B') = dB B' + B dB'
    dYY_ = diag(m*s) # dB
    dYY_[t(matrix(1:(m*s), nrow = m, ncol = s)), ] = dYY_ # dB'
    # print(t(B))
    # print(as.vector(dYY_ %*% h[i[iseq(1,s), iseq(s+1,s+m)]]))
    # testthat::expect_equivalent(as.vector(t(B)), as.vector(dYY_ %*% h[i[iseq(1,s), iseq(s+1,s+m)]]))
    dim(dYY_) = c(m, s*m*s)
    dYY_ = B %*% dYY_
    dim(dYY_) = c(s*s, m*s) # B dB'
    dYY_ = dYY_ + dYY_[t(matrix(1:(s*s), nrow = s, ncol = s)), ] # B dB' + dB B'
    dYY[, i[iseq(1,s), iseq(s+1,s+m)]] = dYY_

    P_J = lyapunov_Jacobian(A, tcrossprod(B), dX, dYY)$J

    dX[] = 0
    dYY[] = 0

    dX_ = diag(s^2) # dA
    dX_[t(matrix(1:(s^2), nrow = s, ncol = s)), ] = dX_ # dA'
    dX[, i[iseq(1,s), iseq(1,s)]] = dX_
    # testthat::expect_equivalent(as.vector(t(A)), as.vector(dX %*% h))

    # d(C' C) = dC' C + C' dC
    dYY_ = diag(m*s) # dC
    # testthat::expect_equivalent(as.vector(C), as.vector(dYY_ %*% h[i[iseq(s+1,s+m), iseq(1,s)]]))
    dim(dYY_) = c(m, s*m*s)
    dYY_ = t(C) %*% dYY_
    dim(dYY_) = c(s*s, m*s) # C' dC
    dYY_ = dYY_ + dYY_[t(matrix(1:(s*s), nrow = s, ncol = s)), ] # C' dC + dC' C
    dYY[, i[iseq(s+1,s+m), iseq(1,s)]] = dYY_

    Q_J = lyapunov_Jacobian(t(A), crossprod(C), dX, dYY)$J
    # browser()
    u1 = upper.tri(diag(s), diag = TRUE)
    u2 = upper.tri(diag(s), diag = FALSE)

    # restrictions: upper diagonal elements of P are zero
    #               upper and diagonal elements of P and Q are equal
    R = rbind(P_J[u2, ], P_J[u1, ] - Q_J[u1, ])

    # compute right kernel of R
    junk = svd(R, nu = 0, nv = ncol(R))
    # print(junk$d)

    # check rank of the matrix R
    d = junk$d
    if ( d[s^2] < (1e-8)*d[1] ) {
      stop('The matrix of "restrictions" does not have dimension s^2=',
           s^2,' (sv[1]=', d[1], ', sv[', s^2, ']=', d[s^2], ')')
    }

    H = junk$v[, (s^2 + 1):((s+m)^2)]
  } else {
    H = matrix(0, nrow = length(h), ncol = 0)
  }

  # add sigma_L parameters
  if (sigma_L_structure == 'identity') {
    sigma_L = diag(m)
  }
  if (sigma_L_structure == 'chol') {
    sigma_L = t(chol( sigma_L %*% t(sigma_L) ))
  }
  if (sigma_L_structure == 'symm') {
    junk = eigen( sigma_L %*% t(sigma_L), symmetric = TRUE )
    sigma_L = junk$vectors %*% diag(sqrt(junk$values), nrow = m, ncol = m) %*% t(junk$vectors)
  }
  junk = tmpl_sigma_L(matrix(NA_integer_, nrow = m, ncol = m), sigma_L_structure)
  junk$h = as.vector(sigma_L)
  # print(cbind(junk$h, junk$H))
  h = c(h, junk$h)
  H = bdiag(H, junk$H)
  # print(cbind(h, H))
  # print(sL)

  return(list(h = h, H = H, class = 'stspmod',
              order = order, n.par = ncol(H)))
}


# Structural Time Series Models ----

#' Structural Time Series Models
#'
#' Tools for Structural Time Series Models, as described e.g. in
#' \insertCite{Harvey94}{RLDM}.
#'
#' **Local Level Model (LLM):** `tmpl_llm()`
#' \deqn{a_{t+1} = a_t + u_t,\quad y_t = a_t}{a[t+1] = a[t] + u[t], y[t] = a[t]}
#' where \eqn{(u_t)}{(u[t])} is white noise with variance \eqn{\sigma_u^2}{\sigma[u]^2}.
#' The model has one free parameter \eqn{\theta = \sigma_u}{\theta = \sigma[u]}.
#' The output process \eqn{(y_t)}{(y[t])} is a random walk.
#'
#' **Local Linear Trend Model (LLTM):** `tmpl_lltm()`
#' \deqn{a_{t+1} = a_t + b_t + u_t,\quad b_{t+1} = b_t + v_t,\quad y_t = a_t}{
#'       a[t+1] = a[t] + b[t] + u[t], b[t+1] = b[t] + v[t], y[t] = a[t]}
#' where \eqn{(u_t)}{(u[t])}, \eqn{(v_t)}{(v[t])} are two independent white noise processes
#' with variance  \eqn{\sigma_u^2}{\sigma[u]^2} and  \eqn{\sigma_v^2}{\sigma[v]^2}.
#' The model has two free parameter \eqn{\theta_1 = \sigma_u}{\theta[1] = \sigma[u]}
#' and  \eqn{\theta_2 = \sigma_v}{\theta[2] = \sigma[v]}. In general the output process is
#' integrated of order two (\eqn{I(2)}). For \eqn{sigma[v]^2=0}{sigma_v^2=0}
#' the model generates a random walk with drift and for \eqn{sigma[u]^2=0}{sigma_u^2=0}
#' one gets an integrated random walk.
#'
#' **Cyclical Models:** `tmpl_cycle(fr, rho)`
#'
#' `tmpl_cycl(fr. rho)` generates a template for scalar AR(2) models, where
#' the AR polynomial has two roots at
#' \deqn{z = \rho^{-1}\exp((\pm i 2\pi f)}{z = \rho^{-1} exp(± i 2\pi f)}
#' If the "damping factor" \eqn{\rho} is close to one then the model generates
#' processes with a strong "cyclical component" with frequency \eqn{f}.
#' For \eqn{\rho <1} the AR(2) model satisfies the stability condition, i.e.
#' the forward solution converges to a stationary process. For \eqn{\rho > 1}
#' the trajectories of the forward solution diverge exponentially.
#' The template has one free parameter, the standard deviation of the driving white noise:
#' \eqn{\theta = \sigma_u}{\theta = \sigma[u]}.
#'
#' **Seasonal Models:** `tmpl_season(s)`
#'
#' `tmpl_season(s)` generates a template for scalar seasonal models, i.e.
#' for models which generate trajectories which are "almost" periodic
#' with a given period, \eqn{s} say. The template has one free parameter,
#' the standard deviation of the driving white noise:
#' \eqn{\theta = \sigma_u}{\theta = \sigma[u]}.
#'
#' **Combine Models** `cbind_templates(...)`
#'
#' The utility `cbind_templates(...)` may be used to construct
#' models from simple "bulding blocks". Suppose e.g. that the observed process is
#' described as the sum of two (unobserved) components
#' \deqn{y_t = k_1(z) u_t + k_2(z) v_t}{y[t] = k[1](z) u[t] + k[2](z) v[t]}
#' where \eqn{(u_t)}{(u[t])}, \eqn{(v_t)}{(v[t])} are two independent white noise processes.
#' If both components are described by the templates `tmpl1` and `tmpl2`
#' then we may construct a template for the combined model simply
#' by `cbind_templates(tmpl1, tmpl2)`.
#'
#' The function `cbind_templates` only deals with state space models and
#' of course all templates must describe outputs with the same dimension.
#'
#' The functions `tmpl_llm()`, ..., `tmpl_season()` generate templates
#' for scalar time series. However, the utility `cbind_templates(...)`
#' also handles the multivariate case.
#'
#' @seealso See [model structures] and
#' [local model structures] for more details on model templates.
#'
#' @param fr,rho frequency and damping factor for cyclical components
#' @param s (integer > 1) period of seasonal component.
#' @param ... compatible (state space model) templates. The output dimensions of the state space models must be
#'            the same for all templates.
#'
#' @return Model template, i.e. a list with slots
#' \item{`h`}{\eqn{((m+s)(n+s) + m^2)}-dimensional vector,}
#' \item{`H`}{\eqn{((m+s)(n+s) + m^2, k)}-dimensional matrix,}
#' \item{`class`}{`class = "stspmod"`, only state space models are implemented}
#' \item{`order`}{`order = c(m,n,s)` (output, noise and state dimensions),}
#' \item{`n.par`}{number of free parameters \eqn{=k} and}
#' \item{`idx`}{a list with slots `state`, `noise` and `par`. These
#'                   indices code which states, noise components and parameters are
#'                   associated to the respective components. See the example(s) below.}
#'
#' @references
#' \insertRef{Harvey94}{RLDM}
#'
#' @name STSmodels
#' @rdname STSmodel
#' @export
#'
#' @examples
#' # build a structural times series model (see Harve 94) with
#' #   a "local linear trend component",
#' #   a cyclical component with period 50 (frequency 1/50),
#' #   a seasonal component with period 6 and
#' #   an AR(1) component.
#' tmpl = cbind_templates(tmpl_lltm(), tmpl_cycle(1/50,1), tmpl_season(6),
#'                        tmpl_stsp_ar(1, 1, sigma_L = 'identity'))
#' # set some "reasonable" values for the standard deviations
#' # of the respective noise and for the AR(1) coefficient.
#' model = fill_template(c(0.0, 0.1,  # parameters for trend (lltm) component
#'                             0.1,       # parameter for cyclical component
#'                             0.1,       # parameter for seasonal component
#'                            -0.5        # AR(1) coefficient
#'                            ), tmpl)
#' print(model)
#'
#' # simulate the time series (with initial states)
#' out = sim(model, n.obs = 100,
#'           a1 = c(100, 1,     # initial states for the trend component
#'                  3, 0,       # initial states for the cyclical component
#'                  5, 10, 10, -10, -10,   # ... for the seasonal component
#'                  0           # initial state for the AR(1) component
#'           ))
#'
#' # extract the contribution of the respective components
#' X = cbind(out$y,
#'  out$a[1:100,tmpl$idx$state == 1, drop = FALSE] %*% model$sys$C[1, tmpl$idx$state == 1] +
#'   out$u[,tmpl$idx$noise == 1, drop = FALSE] %*% model$sys$D[1, tmpl$idx$noise == 1],
#'  out$a[1:100,tmpl$idx$state == 2, drop = FALSE] %*% model$sys$C[1, tmpl$idx$state == 2] +
#'   out$u[,tmpl$idx$noise == 2, drop = FALSE] %*% model$sys$D[1, tmpl$idx$noise == 2],
#'  out$a[1:100,tmpl$idx$state == 3, drop = FALSE] %*% model$sys$C[1, tmpl$idx$state == 3] +
#'   out$u[,tmpl$idx$noise == 3, drop = FALSE] %*% model$sys$D[1, tmpl$idx$noise == 3],
#'  out$a[1:100,tmpl$idx$state == 4, drop = FALSE] %*% model$sys$C[1, tmpl$idx$state == 4] +
#'   out$u[,tmpl$idx$noise == 4, drop = FALSE] %*% model$sys$D[1, tmpl$idx$noise == 4])
#'
#' matplot(X, ylab = 'y', xlab = 't',
#'         type = 'l', lty = 1, col = 1:5)
#' grid()
#' legend('topleft', legend = c('y','trend','cycle','season','AR(1) noise'),
#'        lwd = 2, col = 1:5, bty = 'n')
#'
#' \dontrun{
#' # the following examples throw errors
#' # 1 is not a template
#' cbind_templates(1, tmpl_season(4))
#' # the respective output dimensions are not equal
#' cbind_templates(tmpl_season(4), tmpl_stsp_ar(2, 2))
#' # the third argument is a "VARMA template"
#' cbind_templates(tmpl_lltm(), tmpl_cycle(1/20,1), tmpl_arma_pq(1, 1, 1, 1))
#' }
tmpl_llm = function() {
  model = stspmod(sys = stsp(A = 1, B = 1, C= 1, D = 0), sigma_L = NA_real_)
  tmpl = model2template(model, sigma_L = 'as_given')
  return(tmpl)
}

#' @name STSmodels
#' @rdname STSmodel
#' @export
tmpl_lltm = function() {
#' @examples
#' # Create a template
#' tmpl <- tmpl_lltm()
#' tmpl
#' 
#' # Use the template with fill_template()
#' # filled <- fill_template(tmpl, theta = rnorm(tmpl$n.par))
  model = stspmod(sys = stsp(A = matrix(c(1,0,1,1), nrow = 2, ncol = 2), B = diag(2), C = matrix(c(1,0), nrow = 1),
                             D = matrix(0, nrow = 1, ncol = 2)), sigma_L = diag(x=NA_real_, nrow = 2, ncol = 2))
  # print(model)
  tmpl = model2template(model, sigma_L = 'as_given')
  return(tmpl)
}

#' @name STSmodels
#' @rdname STSmodel
#' @export
tmpl_cycle = function(fr, rho) {
#' @examples
#' # Create a template
#' tmpl <- tmpl_cycle()
#' tmpl
#' 
#' # Use the template with fill_template()
#' # filled <- fill_template(tmpl, theta = rnorm(tmpl$n.par))
  lambda = fr*2*pi
  z0 = exp(complex(imaginary = lambda))*rho
  # a(z) = (1 - z0*z)(1 - Conj(z0)*z)
  a = c(1, -2*Re(z0), Mod(z0)^2)
  model = stspmod(sys = stsp(A = rbind(-a[2:3], c(1,0)),
                             B = c(1, 0),
                             C = -a[2:3],
                             D = 1),
                  sigma_L = NA_real_)
  tmpl = model2template(model, sigma_L = 'as_given')
  return(tmpl)
}

#' @name STSmodels
#' @rdname STSmodel
#' @export
tmpl_season = function(s) {
#' @examples
#' # Create a template
#' tmpl <- tmpl_season()
#' tmpl
#' 
#' # Use the template with fill_template()
#' # filled <- fill_template(tmpl, theta = rnorm(tmpl$n.par))
  s = as.integer(s)[1]
  if (s<=1) stop('period "s" must be an integer > 1')
  model = stspmod(sys = stsp(A = rbind(-rep(1, s-1), diag(x = 1, nrow = s-2, ncol = s-1)),
                             B = c(1, numeric(s-2) ),
                             C = -rep(1, s-1),
                             D = 1),
                  sigma_L = NA_real_)
  tmpl = model2template(model, sigma_L = 'as_given')
  return(tmpl)
}


#' @name STSmodels
#' @rdname STSmodel
#' @export
cbind_templates = function(...) {
#' @examples
#' # Basic example
#' result <- cbind_templates()
#' result
  templates = list(...)
  n_templates = length(templates)

  if (n_templates == 0) return(NULL)

  ok = sapply(templates, FUN = is.template, 'stspmod')
  # print(ok)
  if (!all(ok)) stop('illegal arguments')

  if (n_templates == 1) return(templates[[1]])

  # collect order and n.par of all templates into a matrix
  order_tmpl = function(tmpl) {
    c(tmpl$order, tmpl$n.par)
  }
  order = sapply(templates, FUN = order_tmpl)
  # print(order)
  if ( any(dim(order) != c(4, n_templates)) || (min(order[1, ]) != max(order[1, ])) ) {
    stop('illegal/non-compatible arguments')
  }
  order = t(order)
  order = cbind(order, 1:n_templates)
  colnames(order) = c('m','n','s','n.par','k')
  # print(order)

  m = order[1, 'm']
  N = sum(order[, 'n'])
  S = sum(order[, 's'])
  N.par = sum(order[, 'n.par'])

  h = numeric((m+S)*(N+S) + N^2)
  H = matrix(0, nrow = (m+S)*(N+S) + N^2, ncol = N.par)

  # convert a matrix into a list
  rows2list = function(a) {
    n = nrow(a)
    x = vector(mode = 'list', n)
    for (i in (1:n)) x[[i]] = a[i,]
    return(x)
  }

  A = sapply(rows2list(order[, c('s','k'), drop = FALSE]),
             FUN = function(x) {matrix(x[2], nrow = x[1], ncol = x[1])})
  A = do.call(bdiag, A)
  # print(A)

  B = sapply(rows2list(order[, c('s','n','k'), drop = FALSE]),
             FUN = function(x) {matrix(x[3], nrow = x[1], ncol = x[2])})
  B = do.call(bdiag, B)
  # print(B)

  C = sapply(rows2list(order[, c('m','s','k'), drop = FALSE]),
             FUN = function(x) {matrix(x[3], nrow = x[1], ncol = x[2])})
  C = do.call(cbind, C)
  # print(C)

  D = sapply(rows2list(order[, c('m','n','k'), drop = FALSE]),
             FUN = function(x) {matrix(x[3], nrow = x[1], ncol = x[2])})
  D = do.call(cbind, D)
  # print(D)

  sigma_L = sapply(rows2list(order[, c('n','k'), drop = FALSE]),
                   FUN = function(x) {matrix(x[2], nrow = x[1], ncol = x[1])})
  sigma_L = do.call(bdiag, sigma_L)
  # print(sigma_L)

  par = sapply(rows2list(order[, c('n.par','k'), drop = FALSE]),
               FUN = function(x) {rep(x[2], x[1])})
  par = do.call(c, par)
  # print(par)

  ABCD = rbind( cbind(A, B), cbind(C, D))
  # print(ABCD)
  sys = c(ABCD, sigma_L)

  idx = list(state = C[1,], noise = D[1,], par =par)

  for (k in (1:n_templates)) {
    i = which(sys == k)
    j = which(par == k)
    h[i] = templates[[k]]$h
    H[i, j] = templates[[k]]$H
  }

  return(list(h = h, H = H, class = 'stspmod',
              order = c('m'=nrow(D), 'n'=ncol(D), 's'=nrow(A)),
              n.par = N.par, idx = idx))
}





# Filling a template and extracting the free parameters ----

#' Connect Deep Parameters with a Model
#'
#'
#' `fill_template` fills a given `template`,
#' i.e. in essence an affine mapping from the free parameters to the linear parameters of a model class specified in the template,
#' with the given free (deep) parameters `th` and returns a model (e.g. an [armamod()] or a [stspmod()]).
#' The procedure can be used to generate random models, see e.g. [r_model()], or for *M-estimates*,
#' i.e. for estimation procedures where the estimate is obtained by minimizing (maximizing) a criterion (e.g. ML or GMM estimation).
#' \cr\cr
#' The "inverse function" `extract_theta` extracts the free/deep parameters from the linear parameters of a model,
#' by using the information provided in the template.
#' To this end the procedure first constructs the vector \eqn{\pi} of the stacked (linear) model parameters and then determines the deep parameters \eqn{\theta} as the least squares solution of the equation system
#' \deqn{(\pi - h) =  H\theta}
#' If the residuals are non zero, then the model does not (exactly) fit to the model structure.
#' The threshold `tol` is used to decide whether the model fits to the template or not.
#' The parameter `on_error` decides what to do in the case of a "significant" misfit.
#' For `"ignore"` the procedure ignores the misfit, for `"warn"` the procedure issues a warning,
#'and for `"stop"` the procedure stops with an appropriate error message.
#' \cr
#' In many cases the noise covariance is not explicitly parametrized, since \eqn{\Sigma} is estimated by another route.
#' This may be accomplished by fixing sigma_L to the identity matrix, with the option `sigma_L = "identity"` for the `tmpl_***` functions.
#' In order to **extract the system parameters** (e.g. the AR and MA parameters for an ARMA model) from a model where `sigma_L` is not equal to the identity, one may use the option `ignore_sigma_L = TRUE`.
#' This ignores a possible mismatch for `sigma_L` but still checks whether the system parameters are in accordance to the model structure.
#'
#' @section Connection to Likelihood Estimation:
#'
#' These functions are important for likelihood estimation where the following instances of functionality are necessary.
#' \itemize{
#' \item When an initial estimate is given through a model (together with a template), one may use `extract_theta` to extract the deep parameters. This vector of initial free/deep parameter values needs to be supplied to an optimizer.
#' \item The optimized deep parameter values need to be filled into the model by using the structure provided by a template. This is done with `fill_template`
#' }
#'
#' @param th Vector containing free (deep) parameters.
#' @param template A template like listed in [model structures()], or a template explicitly specified by the user with [model2template()]. Essentially, a template is an affine mapping parametrised as a vector `h` and matrix `H` which connect the free (deep) parameters to the linear parameters of the model.
#' @param model A model object, i.e. an [armamod()], [stspmod()], or [rmfdmod()] object, from which deep parameters should be extracted.
#' @param tol In `extract_theta`, small double specifying the tolerance for an acceptable distance between the linear parameters and `H` times the deep parameters.
#'            Default is st to `sqrt(.Machine$double.eps)`.
#' @param on_error In `extract_theta`, character string with possible choices `ignore`, `warn`, `stop`.
#'                 Specifies what should happen when the distance between the linear parmameters and `H` times the deep parameters, is larger than the specified `tol`.
#'                 Default is `ignore`
#' @param ignore_sigma_L Boolean, default set to `FALSE`.
#'                       If TRUE, the linear and free parameters pertaining
#'                       the left square root `sigma_L` of the error covariance matrix are ignored.
#'                       See also [tmpl_sigma_L()] and [model structures()] for more detail.
#'
#' @return `fill_template` returns a model object, i.e. an [armamod()], [stspmod()], or
#'         [rmfdmod()] object, according to the class of the template (and the given parameters `th`).
#'         The function `extract_theta` returns the vector of *free* parameters for a given model and template.
#'
#' @seealso [model structures()], [ll()], [ll_theta()] and [ll_FUN()].
#'
#' @export
#' @examples
#' # Extract deep parameter from ARMA object with ARMA(p,q) template ##########
#' (armamod_obj = test_armamod(dim = c(2,2), degree = c(3,1)))
#' (tmpl_obj = tmpl_arma_pq(2, 2, 3, 1))
#'
#' extract_theta(armamod_obj, tmpl_obj)
#'
#'
#' # Fill template with deep parameters #################
#' (tmpl_obj = tmpl_arma_echelon(nu = c(3,2)))
#' # Number of columns of matrix H in affine mapping = number of free parameters
#' (n_par_deep = dim(tmpl_obj$H)[2])
#'
#' fill_template(rnorm(n_par_deep), tmpl_obj)
fill_template = function(th, template) {
  if (template$class == 'armamod') {
    order = template$order
    m = order[1]
    n = order[2]
    p = order[3]
    q = order[4]
    P = template$h + template$H %*% th
    sys = matrix(P[1:(length(P) - n*n)], nrow = m, ncol = m*(p+1)+n*(q+1))
    sys = structure(sys, order = order, class = c('lmfd','ratm'))

    sigma_L = matrix(P[((length(P) - n*n + 1)):length(P)],
                     nrow = n, ncol = n)
    model = list(sys = sys, sigma_L = sigma_L, names = NULL, label = NULL)
    model = structure(model, class = c('armamod', 'rldm'))
    return(model)
  }
  if (template$class == 'rmfdmod') {
    order = template$order
    m = order[1]
    n = order[2]
    p = order[3]
    q = order[4]
    P = template$h + template$H %*% th
    sys = matrix(P[1:(length(P) - n*n)], nrow = n*(p+1) + m*(q+1), ncol = n)
    sys = structure(sys, order = order, class = c('rmfd','ratm'))

    sigma_L = matrix(P[((length(P) - n*n + 1)):length(P)],
                     nrow = n, ncol = n)
    model = list(sys = sys, sigma_L = sigma_L, names = NULL, label = NULL)
    model = structure(model, class = c('rmfdmod', 'rldm'))
    return(model)
  }
  if (template$class == 'stspmod') {
    order = template$order
    m = order[1]
    n = order[2]
    s = order[3]
    P = template$h + template$H %*% th
    sys = matrix(P[1:(length(P) - n*n)], nrow = (m+s), ncol = (n+s))
    sys = structure(sys, order = order, class = c('stsp','ratm'))

    sigma_L = matrix(P[((length(P) - n*n + 1)):length(P)],
                     nrow = n, ncol = n)
    model = list(sys = sys, sigma_L = sigma_L, names = NULL, label = NULL)
    model = structure(model, class = c('stspmod', 'rldm'))
    return(model)
  }

  stop('illegal template')
}


#' @rdname fill_template
#' @export
extract_theta = function(model, template,
#' @examples
#' # Basic example
#' result <- extract_theta()
#' result
                         tol = sqrt(.Machine$double.eps), on_error = c('ignore','warn','stop'), ignore_sigma_L = FALSE) {

  on_error = match.arg(on_error)

  if ( !( (class(model)[1] == template$class) && all(dim(model$sys) == template$order) ) ) {
    stop('model and template are not compatible.')
  }

  P = c(as.vector(unclass(model$sys)), as.vector(model$sigma_L))
  if (ncol(template$H) == 0) return(double(0))

  out = stats::lsfit(template$H, P - template$h, intercept = FALSE)

  if (on_error != 'ignore') {
    res = out$residuals
    if (ignore_sigma_L) {
      n = template$order[2]
      res = res[iseq(1,length(res) - n^2)]   # ignore sigma_L
    }
    # is length(res)=0 possible?
    if (length(res) > 0) {
      if (max(abs(res)) > tol) {
        if (on_error == 'warn') warning(paste0('model does not match template. max(abs(res)) = ', max(abs(res))))
        if (on_error == 'stop') stop(paste0('model does not match template. max(abs(res)) = ', max(abs(res))))
      }
    }
  }

  th = unname(out$coef)
  return(th)
}

# Check Template ----

#' Check templates
#
#' @param tmpl object to be tested
#' @param class test for a specific class
#
#' @return TRUE/FALSE
#' @export
#
#' @examples
#' is.template(1)
#' is.template(tmpl_llm(), 'armamod')
#' is.template(tmpl_llm(), 'stspmod')
is.template = function(tmpl, class = c('any','stspmod','armamod','rmfdmod')) {
  class = match.arg(class)
  if (class == 'any') class = c('stspmod','armamod','rmfdmod')

  ok = try({
    cl = tmpl$class
    if (is.na(match(cl, class))) stop('unknown/wrong class')

    n.par = tmpl$n.par
    if ( (!is.vector(n.par)) || (!is.integer(n.par)) || (length(n.par)!=1) || (n.par < 0)) {
      stop('illegal n.par')
    }

    n_order = ifelse(cl == 'stspmod', 3, 4)

    order = tmpl$order
    if ( (!is.vector(order)) || (!is.integer(order)) || (length(order)!=n_order) || (min(order) < 0) ) {
      stop('illegal order')
    }

    h = tmpl$h
    if ( (!is.numeric(h)) || (!is.vector(h)) || (any(is.na(h))) ) stop('illegal h')

    H = tmpl$H
    if ( (!is.numeric(H)) || (!is.matrix(H)) || (any(is.na(H))) || (any(dim(H) != c(length(h), n.par))) ) {
      stop('illegal H')
    }

    if (cl == 'stspmod') {
      m = order[1]
      n = order[2]
      s = order[3]
      if (length(h) != ((m+s)*(n+s) + n^2)) stop('order and h are not compatible')
    } else {
      m = order[1]
      n = order[2]
      p = order[3]
      q = order[4]
      if (length(h) != ((n^2)*(p+1) + m*n*(q+1) + n^2)) stop('order and h are not compatible')
    }
    TRUE
  }, silent = TRUE)

  if (inherits(ok, 'try-error')) ok = FALSE
  return(ok)

}




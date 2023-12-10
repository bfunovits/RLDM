# poles.___ and zeroes.___ method ############################################################
# 
# 
#' Poles and Zeroes
#'
#' Compute the poles and zeroes of VARMA and Statespace models. Note that these models describe 
#' the corresponding processes as 
#' \deqn{x_t = k(B) u_t}{x[t] = k(B) + u[t]}
#' where \eqn{(u_t)}{(u[t])} is a white noise process and \eqn{k(B)} is a rational filter 
#' (\eqn{B} denotes the lag- or backward shift operator). 
#' The poles and zeroes are the poles and zeroes of the rational transfer function 
#' \eqn{k(z)} of this filter. 
#' 
#' @seealso For more details we refer to the discussion about the computation of poles and zeros 
#'          of rational matrices in the companion package, 
#'          see [rationalmatrices::poles and zeroes][rationalmatrices::poles].
#' 
#' @param x an object which represents a VARMA, RMFD or statespace model
#'          (i.e. an [armamod()], [rmfdmod()] or [stspmod()] object).
#' @param tol Double. Default set to `sqrt(.Machine$double.eps)`.
#'   Required to decide on when a root is to be considered "at infinity".
#' @param print_message Boolean. Default set to TRUE.
#'   Prints a message if roots "at infinity " are discarded.
#' @param ... not used.
#'
#' @return Vector of poles, respectively zeroes.
#'
#' @name poles and zeroes
#' @rdname poles_and_zeroes
#' @export
zeroes.armamod = function(x, tol = sqrt(.Machine$double.eps), print_message = TRUE, ...) {
  z = zeroes(x$sys, tol = tol, print_message = print_message)
  return(z)
}


#' @rdname poles_and_zeroes
#' @export
poles.armamod = function(x, tol = sqrt(.Machine$double.eps), print_message = TRUE, ...) {
  z = poles(x$sys, tol = tol, print_message = print_message)
  return(z)
}

#' @rdname poles_and_zeroes
#' @export
zeroes.rmfdmod = function(x, tol = sqrt(.Machine$double.eps), print_message = TRUE, ...) {
  z = zeroes(x$sys, tol = tol, print_message = print_message)
  return(z)
}


#' @rdname poles_and_zeroes
#' @export
poles.rmfdmod = function(x, tol = sqrt(.Machine$double.eps), print_message = TRUE, ...) {
  z = poles(x$sys, tol = tol, print_message = print_message)
  return(z)
}

#' @rdname poles_and_zeroes
#' @export
zeroes.stspmod = function(x, tol = sqrt(.Machine$double.eps), print_message = TRUE, ...) {
  z = zeroes(x$sys, tol = tol, print_message = print_message)
  return(z)
}


#' @rdname poles_and_zeroes
#' @export
poles.stspmod = function(x, tol = sqrt(.Machine$double.eps), print_message = TRUE, ...) {
  z = poles(x$sys, tol = tol, print_message = print_message)
  return(z)
}

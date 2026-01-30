# Importing the pipe ####
# According to https://github.com/sckott/analogsea/issues/32

#' Pipe operator
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL

# Common functions to interaction-version and non-interaction-version of RLS algorithms ####
# __Rcpp roxygen stuff ####
#' @useDynLib RLDM
#' @importFrom Rcpp sourceCpp
NULL


#' RLS function
#'
#' This function implements the Recursive Least Squares (RLS) algorithm with exponentially down-weighted past.
#'
#' @param y Vector of doubles. Response variable in regression. There must not be NAs in this vector.
#' @param X Matrix of doubles. Dimension = (length(y) x maximal number of regressors). NAs are not handled separately (i.e. they must be checked before this function is called).
#' @param r Matrix of doubles (a column vector of dimension \eqn{( \text{forgetting\_factors} \times 1)} containing the forgetting factors.
#' @param n_init Integer. Number of observations used for initial estimate of beta.
#'   If no value provided, at least 21 observations (3 weeks) or 3 times the number of regressors is used.
#' @param start_of_eval Integer. Starting value of evaluation period for honest prediction error
#'   (if we start too early, there are bad initial estimations involved)
#' @param end_of_train Integer. Index specifying the end of the training set. (Afterwards the forecasting period starts.)
#'   Important for calculating the honest prediction error (otherwise it would not be out-of-sample...).
#'   comment_bf: If one wants to use arx_rls_core() outside the automated forecasting framework,
#'   set \code{end_of_train} to \code{length(y)}.
#' @param enhance_conv Boolean. Indicates whether the convergence enhancing factor as described in Young (2011) page 55.
#'
#' @return List containing
#'   \itemize{
#'     \item\code{y_pred}: vector of predictions (vector of same length as input y)
#'     \item \code{fev_honest}: double. honest prediction error (calculated from \code{end_of_train - start_of_eval + 1} observations)
#'     \item \code{forgetting}: double. forgetting factor which produced the minimal honest prediction error
#'   }
#'
#' @export
#' @examples
#' # Generate simple time series data
#' set.seed(123)
#' n <- 100
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' X <- cbind(x1, x2)
#' y <- 2*x1 + 3*x2 + rnorm(n, sd = 0.5)
#' y <- matrix(y, ncol = 1)  # Convert to matrix
#'
#' # Run RLS with default forgetting factors
#' result <- arx_rls_core(y, X, n_init = 20, start_of_eval = 25, end_of_train = 80)
#' str(result)
arx_rls_core <- function(y,
                         X,
                         r = matrix(c(0.9, 0.9500000, 0.9750000, 0.9875000, 0.9937500, 0.9968750, 0.9984375, 1), ncol = 1),
                         n_init = NULL,
                         start_of_eval = NULL,
                         end_of_train = NULL,
                         enhance_conv = TRUE){

  # Number of observations for initial estimator should be at least 3 weeks or 3 times the number of regressors
  if (is.null(n_init)) {
    n_init <- max(3*dim(X)[2], 21)
  }

  # Check: End of training period
  if (is.null(end_of_train) || n_init >= end_of_train ) {
    stop("arx_rls_core(): The start and the end of the training data set must be specified and
         larger than the integer specifying the initial estimation period.")
  }

  if (is.null(start_of_eval)) {
    start_of_eval <- n_init + 1
    message("Start of evaluation period for honest prediction error set to one day after the initial period.
            Note that this might be sub-optimal because the forecast error is usually larger and not representative
            if it has not converged to a steady-state.")
  }

  # initial estimator for beta
  y_init <- y[1:n_init, , drop = FALSE]
  X_init <- X[1:n_init, , drop = FALSE]

  P <- MASS::ginv(t(X_init) %*% X_init)
  beta <- P %*% t(X_init) %*% y_init


  # RLS
  n_obs <- length(y)
  y_pred <- vector("numeric", n_obs)
  y_pred[1:n_init] <- NA
  fe_honest <- vector("numeric", n_obs)

  fev_honest <- vector("numeric", length(r))

  for (i in 1:nrow(r)) {

    for (j in (n_init + 1):n_obs) {
      x <- X[j, , drop = FALSE]
      g <- P %*% t(x) %*% (r[i] + x %*% P %*% t(x))^(-1)

      y_pred[j] <- x %*% beta
      fe_honest[j] <- y[j] - y_pred[j]

      beta <- beta + g %*% fe_honest[j]
      P <- 1/r[i] * (P - g %*% x %*% P )
    }

    fev_honest[i] <- mean(fe_honest[start_of_eval:end_of_train] ^ 2)

    if (i == 1) {
      fev_honest_final <- fev_honest[i]
      y_pred_final <- y_pred
      r_final <- r[i]
    } else if (fev_honest[i] < fev_honest[i - 1]) {
      fev_honest_final <- fev_honest[i]
      y_pred_final <- y_pred
      r_final <- r[i]
    }
  }


  return(list( y_pred = y_pred_final,
               fev_honest = fev_honest_final,
               forgetting = r_final))
  }

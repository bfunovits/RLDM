#' Model Predictions
#'
#' Compute the forecasts based on a VARMA or state space model. The procedure implements a simplified approach.
#' It uses the formulas for the prediction from an infinite past and sets the unknown initial values
#' (prior to \eqn{t < 1}) simply to zero. This simple approach assumes that the model is *stable* and
#' *strictly miniphase* and thus the disturbances \eqn{u_t}{u[t]} are the innovations of the process.
#' Note also that the forecasts for *known* exogenous inputs
#' are calculated, i.e. the "conditional forecasts". For an *honest* prediction,
#' forecasts of the exogenous inputs should be used.
#' \cr
#' The forecast error covariance matrix computed assumes that the true model is used.
#' The error which stems from the estimation of the model is not taken into account.
#' \cr
#' \cr
#' The utility function `evaluate_prediction` may be used to assess the quality of predictions.
#'
#' The utility function `evaluate_prediction` may be used to asses the quality of some given predictions.
#' (E.g. as computed by `predict`). The evaluation criteria to be used are passed as parameter `criteria`
#' to the function. This parameter is a list with components which are either character strings (for selection of
#' one of the implemented quality measures) or a user defined function. Such a function takes three arguments
#' `fun(uhat, y, utilde)` where `uhat` is a matrix of prediction errors,
#' `y` is a matrix of the (corresponding) true values and `ytilde` is a matrix with
#' the predictions of a benchmark prediction procedures. This allows to compute relative error measures.
#'
#' The benchmark predictions are passed on via the parameter `benchmark` to the procedure.
#' If this input parameter is missing then the naive \eqn{h}-step ahead predictions are used a benchmark.
#' (Therefore the user also has to specify the respective forecast horizons via the paramater `h`.)
#'
#' The following evaluation criteria are implemented (\eqn{\hat{y}_{it}}{uhat[it]} denotes
#'   the prediction error for \eqn{y_{it}}{y[it]} and \eqn{\tilde{u}_{it}}{utilde[it]} is the corresponding
#'   error of the benchmark procedure.)
#' \describe{
#' \item{MSE}{Mean Square Error}
#' \item{RMSE}{Root Mean Square Error}
#' \item{MAE}{Mean Absolute Error}
#' \item{MdAE}{Median Absolute Error}
#' \item{MAPE}{Mean Absolute Percentage Error
#'             \eqn{100 mean(|\hat{u}_{it}/y_{it}|)}{100 mean(|uhat[it]/y[it]|)}}
#' \item{MdAPE}{Median Absolute Percentage Error
#'             \eqn{100 median(|\hat{u}_{it}/y_{it}|)}{100 median(|uhat[it]/y[it]|)}}
#' \item{RMdSPE}{Root Median Square Percentage Error
#'             \eqn{100 \sqrt{median(\hat{u}^2_{it}/y^2_{it})}}{100 \sqrt{median(uhat[it]^2 / y[it]^2)}}}
#' \item{RelRMSE}{Relative Root Mean Square Error
#'    \eqn{\sqrt{mean(\hat{u}^2_{it})}/\sqrt{mean(\tilde{u}^2_{it})}}{\sqrt{mean(uhat[it]^2)}/\sqrt{mean(utilde[it]^2)}}}
#' \item{RelMAE}{Relative Mean Absolute Error
#'    \eqn{mean(|\hat{u}_{it}|)/mean(|\tilde{u}^2_{it}|)}{mean(|uhat[it]|)/mean(|utilde[it]|)}}
#' \item{PB}{Percentage better
#'   \eqn{mean(100 I(|\hat{u}_{it}|< |\tilde{u}_{it}|))}{mean(100 I(|uhat[it]| < |utilde[it]|))}}
#' \item{HR}{Hit Rate
#'    \eqn{100 mean( I((\tilde{u}_{it} - \hat{u}_{it})\tilde{u}_{it} \geq 0))}{100 mean( I((utilde[it] - uhat[it])utilde[it] \ge 0))}.
#'    To be precise this measure computes the hit rate only if the naive prediction is the benchmark.}
#' }
#'
#' The procedure also supports the evaluation on different (sub) samples.
#' The parameter `samples` is simply list of integer vectors, where each vector defines a sub sample.
#' E.g. for daily data, one could evaluate the predictions for different weekdays.
#'
#'
#'
#'
#' @param object `rldm_varma` or `rldm_ss` object which represents the model.
#' @param y \eqn{(T,n)} matrix with the observed outputs (\eqn{y_t}{y[t]}, \eqn{t=1,...,T}).
#' @param x \eqn{(T+T_0,r)}{(T+T0,r)} matrix with the exogenous inputs
#'          (\eqn{x_t}{x[t]}, \eqn{t=1,...,T+T_0}{t=1,...,T+T0}). This input parameter is ignored,
#'          if the model has no exogenous inputs. Note that the condition forecasts
#'          are computed and hence (for a model with exogenous inputs) we need the
#'          values of the inputs up to time \eqn{t=T+T_0}{t=T+T0}.
#' @param h (integer) vector of forecast horizons.
#' @param n.ahead (integer) number of time steps to look ahead (out of sample). This number is
#'           also denoted with \eqn{T_0}{T0}.
#' @param yhat \eqn{(T,n,l)} dimensional array of forecasts. The entries `yhat[t,,i]`
#'             should contain the prediction for \eqn{y_t}{y[t]}.
#' @param criteria is a list with "evaluation criteria". See below for more details.
#' @param benchmark \eqn{(T,n,l)} dimensional array with "benchmark" forecasts. If `NULL`
#'            then the naive (\eqn{h}-step ahead) forecasts are used as benchmark.
#' @param samples is a list with "(sub) samples". See below for more details.
#' @param ... not used.
#'
#' @return The function `predict` returns a list with components
#' \item{yhat}{\eqn{(T,n,l)} dimensional array for the h-step ahead forecast. \eqn{T}
#'             is the sample size, \eqn{n} is the dimension of the outputs \eqn{y_t}{y[t]}
#'             and \eqn{l} is the number of forecasts made, i.e. the length of the
#'             vector `h`. The entries `yhat[t,,i]` are the `h[i]`-step
#'             ahead forecast for \eqn{y_t}{y[t]}.}
#' \item{sigmahat}{\eqn{(n,n,l)} dimensional array, where `sigmahat[,,i]` contains the theoretical
#'             covariance matrix of the \eqn{h}-step ahead prediction error for `h=h[i]`.}
#' \item{h}{the (integer) vector of forecasts horizons considered.}
#' \item{yhat.ahead}{\eqn{(T_0,n)}{(T0,n)} dimensional matrix, which contains the "out-of-sample" forecasts
#'             for \eqn{t=T+1, t=T+2,...,t=T+T_0}{t=T+1, t=T+2, ..., t=T+T0}.}
#' \item{sigmahat.ahead}{\eqn{(n,n,T_0)}{(n,n,T0)} dimensional array, where `sigmahat.ahead[,,h]`
#'            contains the theoretical covariance matrix of the \eqn{h}-step ahead prediction
#'            error.}
#' \item{y,x}{the original data.}
#' The function `evaluate_prediction`  returns a 4-dimensional array where the dimensions refer to
#' the evaluation criteria, the (sub) samples, the predictors and the components of the output \eqn{y_t}{y[t]}.
#' Note that the evaluation criteria are applied to the (\eqn{n}) individual components as well as to
#' the joint vector and hence the 4-th dimension of the array has size \eqn{n+1}.
#' E.g. if we consider the "RMSE" and the "MAE" of the forecast errors, two samples (a training sample and a
#' test sample), 5 forecast horizons \eqn{h=1,\ldots,5}{h=1,...,5} and a process
#' \eqn{(y_t)}{(y[t])} with 2 components, then the result will be a \eqn{(2,2,5,3)}-dimensional array.
#'
#' @export
#'
#' @rdname predict
#' @name predict
#'
#' @examples
#'
#' # create a "random" ARMA(1,1) model (stable and miniphase)
#' model = test_armamod(dim = c(2,2), degrees = c(1,1), bpoles = 1, bzeroes = 1)
#'
#' # generate data (sample size n.obs = 200, "burn_in" phase has length 100.)
#' data = sim(model, n.obs = 200, n.burn_in = 100)
#'
#' # predict with true model
#' pred_true = predict(model, data$y, h = c(1,5))
#'
#' # estimate AR model, order selection by AIC
#' n.train = 110   # use the first 110 observation for estimation
#' n.test = 90     # the last 90 observations are used for (a fair) comparison
#' model_ar = est_ar(data$y[1:n.train,],  mean_estimate = "zero",
#'                   ic = 'AIC', method = 'ols')$model
#'
#' # predict with AR model
#' pred_ar = predict(model_ar, data$y, h = c(1,5))
#'
#' # estimate AR1 model (Yule-Walker)
#' model_ar1 = est_ar(data$y[1:n.train,],  mean_estimate = "zero",
#'                   penalty = -1, p.max = 1, method = 'yule-walker')$model
#'
#' # predict with AR1 model
#' pred_ar1 = predict(model_ar1, data$y, h = c(1,5))
#'
#' # evaluate prediction of the AR model (with the AR1 prediction as benchmark)
#' stats = evaluate_prediction(data$y, pred_ar$yhat, h = pred_ar$h,
#'                             criteria = list('RMSE', 'MAE', 'PB'),
#'                             samples = list(train = 21:n.train, test = (n.train+1):(n.train+n.test)),
#'                             benchmark = pred_ar1$yhat)
#'
#' # use array2data.frame for "tabular" display of the results
#' print(array2data.frame(stats, rows = 1:3, cols = 4))
#'
#' # evaluate all predictions
#' # join predictions
#' yhat  = dbind(3, pred_true$yhat, pred_ar1$yhat, pred_ar$yhat)
#'
#' # define a function to compute the "Median Relative Absolute Error"
#' MdRAE_ = function(u.hat, y, u.bench){
#'    stats::median(abs(u.hat/u.bench), na.rm = TRUE)
#' }
#' stats = evaluate_prediction(data$y, yhat,
#'                             h = c(pred_true$h, pred_ar1$h, pred_ar$h),
#'                             criteria = list('RMSE', 'MAE', MdRAE = MdRAE_),
#'                             samples = list(train = 21:n.train, test = (n.train+1):(n.train+n.test)))
#'
#' # split prediction method and forecast horizon
#' dimnames.stats = dimnames(stats)
#' stats = stats[,,c(1,3,5,2,4,6),]
#' dim(stats) = c(3,2,3,2,3)
#' dimnames(stats) = list(criterion = dimnames.stats[[1]], sample = dimnames.stats[[2]],
#'                       model = c('true','AR1','AR'), h = paste('h=',c(1,5),sep=''),
#'                       data = dimnames.stats[[4]])
#'
#' # use array2data.frame for "tabular" display of the results
#' print(array2data.frame(stats, cols = 5, rows = c(3,4,1,2)))

predict.armamod = function(object, y, h = 1, n.ahead = 0, ...) {
  # check input parameters
  d = unname(dim(object$sys))
  m = d[1]
  n = d[2]
  p = d[3]
  q = d[4]

  if ((min(c(m,n,p+1,q+1)) <= 0) || (m != n)) stop('illegal arma model (we should have m=n>0 and p,q>=0)')

  # check 'y'
  y = try(as.matrix(y))
  if (inherits(y, 'try-error')) stop('could not coerce parameter "y" to a matrix')

  if (ncol(y) != m) stop('"y" is not compatible with the ARMA model')
  n.obs = nrow(y)

  # check 'h'
  h = as.integer(h)
  if ((length(h) > 0) && (min(h) <= 0)) stop('the forecast horizon(s) "h" must be positive integers!')

  # check 'n.ahead'
  n.ahead = as.integer(n.ahead)[1]
  if (n.ahead < 0) stop('"n.ahead" must be a non negative integer!')

  if( max(h,n.ahead)==0) {
    stop('nothing to do, since max(h,n.ahead)=0!')
  }


  names = object$names
  if (is.null(names)) names = paste('y[', 1:m, ']', sep ='')


  # compute "innovations/noise"
  u = solve_inverse_de(object$sys, y = y)$u

  # h-step ahead prediction error = k[0] u[t] + k[1] u[t-1] + ... + k[h-1] u[t-h+1]
  # where k[i] is the impulse response

  # compute forecast error variances
  sigma = object$sigma_L %*% t(object$sigma_L)
  sigmahat0 = array(0, dim = c(n, n, max(h,n.ahead)))
  irf = unclass(impresp(object, lag.max = max(h, n.ahead)-1)$irf)
  sigmahat0[,,1] = irf[,,1] %*% sigma %*% t(irf[,,1])
  for (i in iseq(2, max(h, n.ahead))) {
    sigmahat0[,,i] = matrix(sigmahat0[,,i-1], nrow = m, ncol = m) +
      irf[,,i] %*% sigma %*% t(irf[,,i])
  }
  dimnames(sigmahat0) = list(u1 = names, u2 = names, h = 1:(max(h, n.ahead)))

  if (length(h) > 0) {
    yhat = array(0, dim = c(n.obs, n, max(h)))
    yhat[,,1] = y - u %*% t(matrix(irf[,,1], nrow = m, ncol = m))
    for (i in iseq(2, max(h))) {
      if (i <= n.obs) {
        yhat[i:n.obs,,i] = matrix(yhat[i:n.obs,,i-1], nrow = (n.obs-i+1), ncol = m) -
          u[1:(n.obs-i+1),,drop = FALSE] %*%  t(matrix(irf[,,i], nrow = m, ncol = m))
      }
    }
    # extract the horizon's needed
    yhat = yhat[ , , h, drop = FALSE]
    dimnames(yhat) = list(t = (1:n.obs), y = names, h = paste('h=', h, sep=''))
    sigmahat = sigmahat0[ , , h, drop = FALSE]
  } else {
    yhat = array(0, dim= c(n.obs, n, 0))
    dimnames(yhat) = list(t = (1:n.obs), y = names, h = character(0))
    sigmahat = sigmahat0[ , , h, drop = FALSE]
  }

  if (n.ahead > 0) {
    yhat.ahead = solve_de(object$sys, u = matrix(0, nrow = n.ahead, ncol = m),
                          u0 = u, y0 = y)$y
  } else {
    yhat.ahead = matrix(0, nrow = 0, ncol = m)
  }
  colnames(yhat.ahead) = names
  rownames(yhat.ahead) = iseq(n.obs+1, n.obs+n.ahead)
  sigmahat.ahead = sigmahat0[ , , iseq(1,n.ahead), drop = FALSE]

  return(list(yhat = yhat, sigmahat = sigmahat, h = h,
              yhat.ahead = yhat.ahead, sigmahat.ahead = sigmahat.ahead,
              y = y))
}


#' @rdname predict
#' @export
predict.stspmod = function(object, y, x, h = 1, n.ahead = 0, ...) {
  # this is essentially the same as predict.armamod

  # check input parameters
  d = unname(dim(object$sys))
  m = d[1]
  n = d[2]
  s = d[3]

  if ((min(c(m,n)) <= 0) || (m != n)) stop('illegal statespace model (we should have m=n>0)')

  # check 'y'
  y = try(as.matrix(y))
  if (inherits(y, 'try-error')) stop('could not coerce parameter "y" to a matrix')

  if (ncol(y) != m) stop('"y" is not compatible with the ARMA model')
  n.obs = nrow(y)

  # check 'h'
  h = as.integer(h)
  if ((length(h) > 0) && (min(h) <= 0)) stop('the forecast horizon(s) "h" must be positive integers!')

  # check 'n.ahead'
  n.ahead = as.integer(n.ahead)[1]
  if (n.ahead < 0) stop('"n.ahead" must be a non negative integer!')

  if( max(h,n.ahead)==0) {
    stop('nothing to do, since max(h,n.ahead)=0!')
  }

  names = object$names
  if (is.null(names)) names = paste('y[', 1:m, ']', sep ='')

  # compute "innovations/noise" and! the states a[t]
  data = solve_inverse_de(object$sys, y = y)
  u = data$u

  # h-step ahead prediction error = k[0] u[t] + k[1] u[t-1] + ... + k[h-1] u[t-h+1]
  # where k[i] is the impulse response

  # compute forecast error variances
  sigma = object$sigma_L %*% t(object$sigma_L)
  sigmahat0 = array(0, dim = c(n, n, max(h,n.ahead)))
  irf = unclass(impresp(object, lag.max = max(h, n.ahead)-1)$irf)
  sigmahat0[,,1] = irf[,,1] %*% sigma %*% t(irf[,,1])
  for (i in iseq(2, max(h, n.ahead))) {
    sigmahat0[,,i] = matrix(sigmahat0[,,i-1], nrow = m, ncol = m) +
      irf[,,i] %*% sigma %*% t(irf[,,i])
  }
  dimnames(sigmahat0) = list(u1 = names, u2 = names, h = 1:(max(h, n.ahead)))

  if (length(h) > 0) {
    yhat = array(0, dim = c(n.obs, n, max(h)))
    yhat[,,1] = y - u %*% t(matrix(irf[,,1], nrow = m, ncol = m))
    for (i in iseq(2, max(h))) {
      if (i <= n.obs) {
        yhat[i:n.obs,,i] = matrix(yhat[i:n.obs,,i-1], nrow = (n.obs-i+1), ncol = m) -
          u[1:(n.obs-i+1),,drop = FALSE] %*%  t(matrix(irf[,,i], nrow = m, ncol = m))
      }
    }
    # extract the horizon's needed
    yhat = yhat[ , , h, drop = FALSE]
    dimnames(yhat) = list(t = (1:n.obs), y = names, h = paste('h=', h, sep=''))
    sigmahat = sigmahat0[ , , h, drop = FALSE]
  } else {
    yhat = array(0, dim= c(n.obs, n, 0))
    dimnames(yhat) = list(t = (1:n.obs), y = names, h = character(0))
    sigmahat = sigmahat0[ , , h, drop = FALSE]
  }

  if (n.ahead > 0) {
    # this line is the only difference to predict.rldm_varma!
    yhat.ahead = solve_de(object$sys, u = matrix(0, nrow = n.ahead, ncol = n),
                          a1 = data$a[(n.obs+1),])$y
  } else {
    yhat.ahead = matrix(0, nrow = 0, ncol = n)
  }
  colnames(yhat.ahead) = names
  rownames(yhat.ahead) = iseq(n.obs+1, n.obs+n.ahead)
  sigmahat.ahead = sigmahat0[ , , iseq(1,n.ahead), drop = FALSE]

  return(list(yhat = yhat, sigmahat = sigmahat, h = h,
              yhat.ahead = yhat.ahead, sigmahat.ahead = sigmahat.ahead,
              y = y))
}


# expand_letters is an internal function
# #' Creating letters
# #'
# #' Create a 'character' vector of length n out of the letters 'l' if 'l' has less than n entries,
# #' then combinations are created.
# #'
# #' @note
# #' Should eventually be internal.
# #'
# #' @param n Integer. Number of letters to be created.
# #' @param l Vector of characters. Must be at least of length 2.
# #'
# #' @return Vector of characters of length n.
# #' @export
# #'
# #' @examples
# #' expand_letters(4, l = c('a','b'))
# #' expand_letters(6, l = c('a','b'))
expand_letters = function(n, l = letters) {
  if (n==0) return(character(0))
  if ((n>1) && (length(l) <= 1)) stop('"l" must have at least 2 entries!')
  l0 = l
  while (length(l) < n) {
    l = outer(l,l0,FUN = function(a,b) {paste(b,a,sep='')})
  }
  l[1:n]
}


MSE_ = function(u.hat, y, u.bench) {mean((u.hat)^2, na.rm = TRUE)}
RMSE_ = function(u.hat, y, u.bench) {sqrt(mean((u.hat)^2, na.rm = TRUE))}
MAE_ = function(u.hat, y, u.bench) {mean(abs(u.hat), na.rm = TRUE)}
MdAE_ = function(u.hat, y, u.bench) {stats::median(abs(u.hat), na.rm = TRUE)}
MAPE_ = function(u.hat, y, u.bench) {100*mean(abs(u.hat/y), na.rm = TRUE)}
MdAPE_ = function(u.hat, y, u.bench) {100*stats::median(abs(u.hat/y), na.rm = TRUE)}
RMdSPE_ = function(u.hat, y, u.bench) {100*sqrt(stats::median((u.hat/y)^2, na.rm = TRUE))}
RelRMSE_ = function(u.hat, y, u.bench) {sqrt(mean((u.hat)^2, na.rm = TRUE)) / sqrt(mean((u.bench)^2, na.rm = TRUE))}
RelMAE_ = function(u.hat, y, u.bench) {mean(abs(u.hat), na.rm = TRUE) / mean(abs(u.bench), na.rm = TRUE)}
PB_ = function(u.hat, y, u.bench) {100*mean(abs(u.hat) < abs(u.bench), na.rm = TRUE)}
HR_ = function(u.hat, y, u.bench) {100*mean(((u.bench-u.hat)*u.bench) >= 0, na.rm = TRUE)}

#' @rdname predict
#' @export
evaluate_prediction = function(y, yhat, h, criteria = list('RMSE'),
                               benchmark = NULL, samples = list(1:nrow(y))) {


  # check evluation criteria
  criteria = as.list(criteria)
  n.criteria = length(criteria)
  if (n.criteria == 0) stop('no evaluation criteria provided!')
  criteria.names = names(criteria)
  if (is.null(criteria.names)) {
    criteria.names = paste('crit', expand_letters(n.criteria,l = LETTERS),sep='')
  }
  for (i in (1:n.criteria)) {
    cr = criteria[[i]]
    if ( (is.character(cr)) && (length(cr) == 1) ) {
      criteria.names[[i]] = cr
      fun = switch(cr,
                   MSE = MSE_,
                   RMSE = RMSE_,
                   MAE = MAE_,
                   MdAE = MdAE_,
                   MAPE = MAPE_,
                   MdAPE = MdAPE_,
                   RMdSPE = RMdSPE_,
                   RelRMSE = RelRMSE_,
                   RelMAE = RelMAE_,
                   PB = PB_,
                   HR = HR_)
      if (is.null(fun)) stop('criterion "',cr,'" is not implemented!')
    } else {
      fun = cr
      if (!is.function(fun)) {
        stop('"criteria" has to be list of functions "fun(u.hat,y,u.bench)" or character strings!')
      }
    }
    criteria[[i]] = fun
  }

  # check 'y'
  if (!is.numeric(y)) stop('"y" is not numeric!')
  y = as.matrix(y)
  n = ncol(y)
  n.obs = nrow(y)
  if ((n.obs*n)==0) stop('"y" is empty!')
  y.names = colnames(y)
  if (is.null(y.names)) {
    y.names = paste('y[',1:n,']',sep='')
  }
  y.names = c(y.names, 'total')

  # check forecasts 'yhat'
  if (!is.numeric(yhat)) stop('"yhat" is not numeric!')
  if (is.vector(yhat)) {
    yhat = matrix(yhat, nrow = length(yhat), ncol = 1)
  }
  if (is.matrix(yhat)) {
    dim(yhat) = c(nrow(yhat),ncol(yhat),1)
  }
  if ( (!is.array(yhat)) || (length(dim(yhat))!=3) ||
       (dim(yhat)[1] != n.obs) || (dim(yhat)[2] != n) ||
       (dim(yhat)[3]<=0) ) stop('"yhat" is not a compatible array!')
  n.predictors = dim(yhat)[3]

  # check forecast horizons 'h'
  h = as.integer(as.vector(h))
  if ((length(h) != n.predictors) || (min(h) <= 0)) stop('"h" is not compatible!')

  predictor.names = dimnames(yhat)[[3]]
  if (is.null(predictor.names)) {
    predictor.names = paste('pred',expand_letters(n.predictors,l = LETTERS),'_h', h, sep = '')
  }

  # check 'samples'
  if (!is.list(samples)) samples = list(samples)
  n.samples = length(samples)
  if (n.samples == 0) stop('no samples provided!')
  for (i in (1:n.samples)) {
    samples[[i]] = as.integer(as.vector(samples[[i]]))
    if ((min(samples[[i]])<1) || (max(samples[[i]])>n.obs) ||
        (length(samples[[i]]) == 0) ) stop('illegal sample(s)!')
  }
  sample.names = names(samples)
  if (is.null(sample.names)) {
    sample.names = paste('smpl',expand_letters(n.samples,l = LETTERS),sep='')
  }

  # create benchmark prediction = naive prediction (if necessary)
  if (is.null(benchmark)) {
    benchmark = array(NA+0, dim = c(n.obs, n, n.predictors))
    for (j in (1:n.predictors)) {
      benchmark[(1+h[j]):n.obs,,j] = y[1:(n.obs-h[j]),]
    }
  }

  # check benchmark
  if (!is.numeric(benchmark)) stop('"benchmark" is not numeric!')
  if (is.vector(benchmark)) {
    benchmark = matrix(benchmark, nrow = length(benchmark), ncol = 1)
  }
  if (is.matrix(benchmark)) {
    dim(benchmark) = c(nrow(benchmark),ncol(benchmark),1)
  }
  if ( (!is.array(benchmark)) || (length(dim(benchmark))!=3) ||
       (dim(benchmark)[1] != n.obs) || (dim(benchmark)[2] != n) ||
       (dim(benchmark)[3] != n.predictors) ) stop('"benchmark" is not a compatible array!')

  # compute forecast errors
  for (j in (1:n.predictors)) {
    yhat[,,j] = y - yhat[,,j]
    benchmark[,,j] = y - benchmark[,,j]
  }

  # create stats array
  stats = array(0, dim = c(n.criteria, n.samples, n.predictors, n+1))

  for (i in (1:n.criteria)) {
    fun = criteria[[i]]
    for (j in (1:n.samples)) {
      sample = samples[[j]]
      for (k in (1:n.predictors)) {
        for (l in (1:n)) {
          stats[i,j,k,l] = fun(yhat[sample,l,k], y[sample,l], benchmark[sample,l,k])
        }
        stats[i,j,k,n+1] = fun(matrix(yhat[sample,,k], nrow = length(sample), ncol = n),
                               y[sample,,drop = FALSE],
                               matrix(benchmark[sample,,k], nrow = length(sample), ncol = n))
      }
    }
  }
  # print(str(stats))
  # print(list(y = y.names, predictor = predictor.names, criterion = criteria.names, sample = sample.names))

  dimnames.stats = list(criterion = criteria.names, sample = sample.names,
                        predictor = predictor.names, data = y.names)
  dimnames(stats) = dimnames.stats

  return(stats)
}


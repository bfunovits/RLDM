# Model Predictions

Compute the forecasts based on a VARMA or state space model. The
procedure implements a simplified approach. It uses the formulas for the
prediction from an infinite past and sets the unknown initial values
(prior to \\t \< 1\\) simply to zero. This simple approach assumes that
the model is *stable* and *strictly miniphase* and thus the disturbances
\\u_t\\ are the innovations of the process. Note also that the forecasts
for *known* exogenous inputs are calculated, i.e. the "conditional
forecasts". For an *honest* prediction, forecasts of the exogenous
inputs should be used.  
The forecast error covariance matrix computed assumes that the true
model is used. The error which stems from the estimation of the model is
not taken into account.  
  
The utility function `evaluate_prediction` may be used to assess the
quality of predictions.

## Usage

``` r
# S3 method for class 'armamod'
predict(object, y, h = 1, n.ahead = 0, ...)

# S3 method for class 'stspmod'
predict(object, y, x, h = 1, n.ahead = 0, ...)

evaluate_prediction(
  y,
  yhat,
  h,
  criteria = list("RMSE"),
  benchmark = NULL,
  samples = list(1:nrow(y))
)
```

## Arguments

- object:

  `rldm_varma` or `rldm_ss` object which represents the model.

- y:

  \\(T,n)\\ matrix with the observed outputs (\\y_t\\, \\t=1,...,T\\).

- h:

  (integer) vector of forecast horizons.

- n.ahead:

  (integer) number of time steps to look ahead (out of sample). This
  number is also denoted with \\T_0\\.

- ...:

  not used.

- x:

  \\(T+T_0,r)\\ matrix with the exogenous inputs (\\x_t\\,
  \\t=1,...,T+T_0\\). This input parameter is ignored, if the model has
  no exogenous inputs. Note that the condition forecasts are computed
  and hence (for a model with exogenous inputs) we need the values of
  the inputs up to time \\t=T+T_0\\.

- yhat:

  \\(T,n,l)\\ dimensional array of forecasts. The entries `yhat[t,,i]`
  should contain the prediction for \\y_t\\.

- criteria:

  is a list with "evaluation criteria". See below for more details.

- benchmark:

  \\(T,n,l)\\ dimensional array with "benchmark" forecasts. If `NULL`
  then the naive (\\h\\-step ahead) forecasts are used as benchmark.

- samples:

  is a list with "(sub) samples". See below for more details.

## Value

The function `predict` returns a list with components

- yhat:

  \\(T,n,l)\\ dimensional array for the h-step ahead forecast. \\T\\ is
  the sample size, \\n\\ is the dimension of the outputs \\y_t\\ and
  \\l\\ is the number of forecasts made, i.e. the length of the vector
  `h`. The entries `yhat[t,,i]` are the `h[i]`-step ahead forecast for
  \\y_t\\.

- sigmahat:

  \\(n,n,l)\\ dimensional array, where `sigmahat[,,i]` contains the
  theoretical covariance matrix of the \\h\\-step ahead prediction error
  for `h=h[i]`.

- h:

  the (integer) vector of forecasts horizons considered.

- yhat.ahead:

  \\(T_0,n)\\ dimensional matrix, which contains the "out-of-sample"
  forecasts for \\t=T+1, t=T+2,...,t=T+T_0\\.

- sigmahat.ahead:

  \\(n,n,T_0)\\ dimensional array, where `sigmahat.ahead[,,h]` contains
  the theoretical covariance matrix of the \\h\\-step ahead prediction
  error.

- y,x:

  the original data.

The function `evaluate_prediction` returns a 4-dimensional array where
the dimensions refer to the evaluation criteria, the (sub) samples, the
predictors and the components of the output \\y_t\\. Note that the
evaluation criteria are applied to the (\\n\\) individual components as
well as to the joint vector and hence the 4-th dimension of the array
has size \\n+1\\. E.g. if we consider the "RMSE" and the "MAE" of the
forecast errors, two samples (a training sample and a test sample), 5
forecast horizons \\h=1,\ldots,5\\ and a process \\(y_t)\\ with 2
components, then the result will be a \\(2,2,5,3)\\-dimensional array.

## Details

The utility function `evaluate_prediction` may be used to asses the
quality of some given predictions. (E.g. as computed by `predict`). The
evaluation criteria to be used are passed as parameter `criteria` to the
function. This parameter is a list with components which are either
character strings (for selection of one of the implemented quality
measures) or a user defined function. Such a function takes three
arguments `fun(uhat, y, utilde)` where `uhat` is a matrix of prediction
errors, `y` is a matrix of the (corresponding) true values and `ytilde`
is a matrix with the predictions of a benchmark prediction procedures.
This allows to compute relative error measures.

The benchmark predictions are passed on via the parameter `benchmark` to
the procedure. If this input parameter is missing then the naive
\\h\\-step ahead predictions are used a benchmark. (Therefore the user
also has to specify the respective forecast horizons via the paramater
`h`.)

The following evaluation criteria are implemented (\\\hat{y}\_{it}\\
denotes the prediction error for \\y\_{it}\\ and \\\tilde{u}\_{it}\\ is
the corresponding error of the benchmark procedure.)

- MSE:

  Mean Square Error

- RMSE:

  Root Mean Square Error

- MAE:

  Mean Absolute Error

- MdAE:

  Median Absolute Error

- MAPE:

  Mean Absolute Percentage Error \\100 mean(\|\hat{u}\_{it}/y\_{it}\|)\\

- MdAPE:

  Median Absolute Percentage Error \\100
  median(\|\hat{u}\_{it}/y\_{it}\|)\\

- RMdSPE:

  Root Median Square Percentage Error \\100
  \sqrt{median(\hat{u}^2\_{it}/y^2\_{it})}\\

- RelRMSE:

  Relative Root Mean Square Error
  \\\sqrt{mean(\hat{u}^2\_{it})}/\sqrt{mean(\tilde{u}^2\_{it})}\\

- RelMAE:

  Relative Mean Absolute Error
  \\mean(\|\hat{u}\_{it}\|)/mean(\|\tilde{u}^2\_{it}\|)\\

- PB:

  Percentage better \\mean(100 I(\|\hat{u}\_{it}\|\<
  \|\tilde{u}\_{it}\|))\\

- HR:

  Hit Rate \\100 mean( I((\tilde{u}\_{it} -
  \hat{u}\_{it})\tilde{u}\_{it} \geq 0))\\. To be precise this measure
  computes the hit rate only if the naive prediction is the benchmark.

The procedure also supports the evaluation on different (sub) samples.
The parameter `samples` is simply list of integer vectors, where each
vector defines a sub sample. E.g. for daily data, one could evaluate the
predictions for different weekdays.

## Examples

``` r
# create a "random" ARMA(1,1) model (stable and miniphase)
model = test_armamod(dim = c(2,2), degrees = c(1,1), bpoles = 1, bzeroes = 1)

# generate data (sample size n.obs = 200, "burn_in" phase has length 100.)
data = sim(model, n.obs = 200, n.burn_in = 100)

# predict with true model
pred_true = predict(model, data$y, h = c(1,5))

# estimate AR model, order selection by AIC
n.train = 110   # use the first 110 observation for estimation
n.test = 90     # the last 90 observations are used for (a fair) comparison
model_ar = est_ar(data$y[1:n.train,],  mean_estimate = "zero",
                  ic = 'AIC', method = 'ols')$model

# predict with AR model
pred_ar = predict(model_ar, data$y, h = c(1,5))

# estimate AR1 model (Yule-Walker)
model_ar1 = est_ar(data$y[1:n.train,],  mean_estimate = "zero",
                  penalty = -1, p.max = 1, method = 'yule-walker')$model

# predict with AR1 model
pred_ar1 = predict(model_ar1, data$y, h = c(1,5))

# evaluate prediction of the AR model (with the AR1 prediction as benchmark)
stats = evaluate_prediction(data$y, pred_ar$yhat, h = pred_ar$h,
                            criteria = list('RMSE', 'MAE', 'PB'),
                            samples = list(train = 21:n.train, test = (n.train+1):(n.train+n.test)),
                            benchmark = pred_ar1$yhat)

# use array2data.frame for "tabular" display of the results
print(array2data.frame(stats, rows = 1:3, cols = 4))
#>    criterion sample predictor      y[1]        y[2]      total
#> 1       RMSE  train       h=1  1.581590  0.02788022  1.1185268
#> 2        MAE  train       h=1  1.292886  0.02163501  0.6572603
#> 3         PB  train       h=1 51.111111 75.55555556 63.3333333
#> 4       RMSE   test       h=1  1.611489  0.02871524  1.1396757
#> 5        MAE   test       h=1  1.265299  0.02190208  0.6436004
#> 6         PB   test       h=1 45.555556 87.77777778 66.6666667
#> 7       RMSE  train       h=5  1.749987  0.50162780  1.2872618
#> 8        MAE  train       h=5  1.392897  0.40512846  0.8990125
#> 9         PB  train       h=5 52.222222 46.66666667 49.4444444
#> 10      RMSE   test       h=5  1.612239  0.50362675  1.1943520
#> 11       MAE   test       h=5  1.213223  0.40452278  0.8088731
#> 12        PB   test       h=5 44.444444 50.00000000 47.2222222

# evaluate all predictions
# join predictions
yhat  = dbind(3, pred_true$yhat, pred_ar1$yhat, pred_ar$yhat)

# define a function to compute the "Median Relative Absolute Error"
MdRAE_ = function(u.hat, y, u.bench){
   stats::median(abs(u.hat/u.bench), na.rm = TRUE)
}
stats = evaluate_prediction(data$y, yhat,
                            h = c(pred_true$h, pred_ar1$h, pred_ar$h),
                            criteria = list('RMSE', 'MAE', MdRAE = MdRAE_),
                            samples = list(train = 21:n.train, test = (n.train+1):(n.train+n.test)))

# split prediction method and forecast horizon
dimnames.stats = dimnames(stats)
stats = stats[,,c(1,3,5,2,4,6),]
dim(stats) = c(3,2,3,2,3)
dimnames(stats) = list(criterion = dimnames.stats[[1]], sample = dimnames.stats[[2]],
                      model = c('true','AR1','AR'), h = paste('h=',c(1,5),sep=''),
                      data = dimnames.stats[[4]])

# use array2data.frame for "tabular" display of the results
print(array2data.frame(stats, cols = 5, rows = c(3,4,1,2)))
#>    model   h criterion sample      y[1]       y[2]     total
#> 1   true h=1      RMSE  train 1.6934124 0.02952683 1.1976054
#> 2    AR1 h=1      RMSE  train 1.6932010 0.08461760 1.1987681
#> 3     AR h=1      RMSE  train 1.5815900 0.02788022 1.1185268
#> 4   true h=5      RMSE  train 1.8237975 0.52615332 1.3422136
#> 5    AR1 h=5      RMSE  train 1.8054364 0.51155715 1.3268933
#> 6     AR h=5      RMSE  train 1.7499873 0.50162780 1.2872618
#> 7   true h=1       MAE  train 1.3542487 0.02379324 0.6890210
#> 8    AR1 h=1       MAE  train 1.3396762 0.06517521 0.7024257
#> 9     AR h=1       MAE  train 1.2928857 0.02163501 0.6572603
#> 10  true h=5       MAE  train 1.4402109 0.43114789 0.9356794
#> 11   AR1 h=5       MAE  train 1.4176281 0.41286399 0.9152460
#> 12    AR h=5       MAE  train 1.3928966 0.40512846 0.8990125
#> 13  true h=1     MdRAE  train 0.6470609 0.05429154 0.2378054
#> 14   AR1 h=1     MdRAE  train 0.6079586 0.14506506 0.3293400
#> 15    AR h=1     MdRAE  train 0.6392057 0.05059007 0.2203919
#> 16  true h=5     MdRAE  train 0.8171496 0.98500188 0.8872385
#> 17   AR1 h=5     MdRAE  train 0.8191178 1.01087683 0.8861237
#> 18    AR h=5     MdRAE  train 0.8127587 0.94625998 0.8565500
#> 19  true h=1      RMSE   test 1.5292605 0.02802385 1.0815320
#> 20   AR1 h=1      RMSE   test 1.5527279 0.08241224 1.0994898
#> 21    AR h=1      RMSE   test 1.6114890 0.02871524 1.1396757
#> 22  true h=5      RMSE   test 1.5968055 0.48483869 1.1800119
#> 23   AR1 h=5      RMSE   test 1.5987971 0.49913896 1.1843335
#> 24    AR h=5      RMSE   test 1.6122387 0.50362675 1.1943520
#> 25  true h=1       MAE   test 1.1691349 0.02105252 0.5950937
#> 26   AR1 h=1       MAE   test 1.2020720 0.06628620 0.6341791
#> 27    AR h=1       MAE   test 1.2652987 0.02190208 0.6436004
#> 28  true h=5       MAE   test 1.1771857 0.39086052 0.7840231
#> 29   AR1 h=5       MAE   test 1.1854835 0.40537031 0.7954269
#> 30    AR h=5       MAE   test 1.2132234 0.40452278 0.8088731
#> 31  true h=1     MdRAE   test 0.6272284 0.04769589 0.1762662
#> 32   AR1 h=1     MdRAE   test 0.6488267 0.16345762 0.3293592
#> 33    AR h=1     MdRAE   test 0.7003001 0.05752692 0.2239963
#> 34  true h=5     MdRAE   test 0.7191348 0.62282870 0.6443024
#> 35   AR1 h=5     MdRAE   test 0.7371333 0.62450028 0.6636491
#> 36    AR h=5     MdRAE   test 0.6912901 0.65825976 0.6656079
```

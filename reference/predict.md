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
#>    criterion sample predictor       y[1]       y[2]      total
#> 1       RMSE  train       h=1  1.1880098  0.8340406  1.0263993
#> 2        MAE  train       h=1  0.8994266  0.6381087  0.7687677
#> 3         PB  train       h=1 52.2222222 57.7777778 55.0000000
#> 4       RMSE   test       h=1  1.0029791  0.8260473  0.9187822
#> 5        MAE   test       h=1  0.8048857  0.6557681  0.7303269
#> 6         PB   test       h=1 53.3333333 62.2222222 57.7777778
#> 7       RMSE  train       h=5  1.5345725  1.4710960  1.5031694
#> 8        MAE  train       h=5  1.1968908  1.0760500  1.1364704
#> 9         PB  train       h=5 51.1111111 47.7777778 49.4444444
#> 10      RMSE   test       h=5  1.4687151  1.6055353  1.5386468
#> 11       MAE   test       h=5  1.1476862  1.2845105  1.2160983
#> 12        PB   test       h=5 50.0000000 44.4444444 47.2222222

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
#>    model   h criterion sample      y[1]      y[2]     total
#> 1   true h=1      RMSE  train 1.1410587 0.7948097 0.9832948
#> 2    AR1 h=1      RMSE  train 1.3865012 0.9150052 1.1746532
#> 3     AR h=1      RMSE  train 1.1880098 0.8340406 1.0263993
#> 4   true h=5      RMSE  train 1.5523979 1.4718143 1.5126428
#> 5    AR1 h=5      RMSE  train 1.6091143 1.5098116 1.5602532
#> 6     AR h=5      RMSE  train 1.5345725 1.4710960 1.5031694
#> 7   true h=1       MAE  train 0.8761662 0.6219701 0.7490681
#> 8    AR1 h=1       MAE  train 1.0311656 0.7141773 0.8726714
#> 9     AR h=1       MAE  train 0.8994266 0.6381087 0.7687677
#> 10  true h=5       MAE  train 1.2334374 1.0979894 1.1657134
#> 11   AR1 h=5       MAE  train 1.2805796 1.1249307 1.2027551
#> 12    AR h=5       MAE  train 1.1968908 1.0760500 1.1364704
#> 13  true h=1     MdRAE  train 0.5067438 0.4245879 0.4617085
#> 14   AR1 h=1     MdRAE  train 0.4981103 0.4867328 0.4867328
#> 15    AR h=1     MdRAE  train 0.4515703 0.3713669 0.4242064
#> 16  true h=5     MdRAE  train 0.7023582 0.6816619 0.6838157
#> 17   AR1 h=5     MdRAE  train 0.6905275 0.6953111 0.6953111
#> 18    AR h=5     MdRAE  train 0.6555941 0.7287681 0.6911765
#> 19  true h=1      RMSE   test 0.8761931 0.7631273 0.8216075
#> 20   AR1 h=1      RMSE   test 1.1073192 0.8712328 0.9962937
#> 21    AR h=1      RMSE   test 1.0029791 0.8260473 0.9187822
#> 22  true h=5      RMSE   test 1.3909232 1.5019212 1.4474866
#> 23   AR1 h=5      RMSE   test 1.4247032 1.5239341 1.4751533
#> 24    AR h=5      RMSE   test 1.4687151 1.6055353 1.5386468
#> 25  true h=1       MAE   test 0.7094357 0.5876944 0.6485651
#> 26   AR1 h=1       MAE   test 0.9068906 0.7030084 0.8049495
#> 27    AR h=1       MAE   test 0.8048857 0.6557681 0.7303269
#> 28  true h=5       MAE   test 1.0913610 1.1799927 1.1356769
#> 29   AR1 h=5       MAE   test 1.0930442 1.1557806 1.1244124
#> 30    AR h=5       MAE   test 1.1476862 1.2845105 1.2160983
#> 31  true h=1     MdRAE   test 0.4229676 0.4209370 0.4227645
#> 32   AR1 h=1     MdRAE   test 0.5621913 0.4001367 0.4842738
#> 33    AR h=1     MdRAE   test 0.4739936 0.3579215 0.3908032
#> 34  true h=5     MdRAE   test 0.7248938 0.7631031 0.7401132
#> 35   AR1 h=5     MdRAE   test 0.7090822 0.6644633 0.7049224
#> 36    AR h=5     MdRAE   test 0.6970866 0.7656649 0.7380062
# Basic example
result <- evaluate_prediction()
#> Error in evaluate_prediction(): argument "y" is missing, with no default
result
#> Error: object 'result' not found
```

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
#>    criterion sample predictor       y[1]      y[2]     total
#> 1       RMSE  train       h=1  0.6539811  2.287527  1.682330
#> 2        MAE  train       h=1  0.5234412  1.744002  1.133722
#> 3         PB  train       h=1 56.6666667 52.222222 54.444444
#> 4       RMSE   test       h=1  0.7260734  2.516711  1.852163
#> 5        MAE   test       h=1  0.5982744  2.037086  1.317680
#> 6         PB   test       h=1 58.8888889 50.000000 54.444444
#> 7       RMSE  train       h=5  2.9084169  4.081186  3.543654
#> 8        MAE  train       h=5  2.3003365  3.214216  2.757276
#> 9         PB  train       h=5 52.2222222 54.444444 53.333333
#> 10      RMSE   test       h=5  4.0055572  5.269170  4.680205
#> 11       MAE   test       h=5  3.2446598  4.165137  3.704898
#> 12        PB   test       h=5 44.4444444 41.111111 42.777778

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
#> 1   true h=1      RMSE  train 0.6719019 2.3451931 1.7250193
#> 2    AR1 h=1      RMSE  train 0.8927488 2.4775867 1.8621810
#> 3     AR h=1      RMSE  train 0.6539811 2.2875270 1.6823304
#> 4   true h=5      RMSE  train 2.9639855 4.1420597 3.6015183
#> 5    AR1 h=5      RMSE  train 2.9458916 4.0902008 3.5642685
#> 6     AR h=5      RMSE  train 2.9084169 4.0811860 3.5436541
#> 7   true h=1       MAE  train 0.5419697 1.7596616 1.1508157
#> 8    AR1 h=1       MAE  train 0.6761560 1.8638281 1.2699921
#> 9     AR h=1       MAE  train 0.5234412 1.7440020 1.1337216
#> 10  true h=5       MAE  train 2.2678501 3.2323421 2.7500961
#> 11   AR1 h=5       MAE  train 2.2625911 3.2567509 2.7596710
#> 12    AR h=5       MAE  train 2.3003365 3.2142164 2.7572764
#> 13  true h=1     MdRAE  train 0.3714473 0.7019594 0.5446233
#> 14   AR1 h=1     MdRAE  train 0.4768955 0.7407391 0.5369586
#> 15    AR h=1     MdRAE  train 0.3718159 0.6523235 0.5372316
#> 16  true h=5     MdRAE  train 0.6353331 0.6557029 0.6557029
#> 17   AR1 h=5     MdRAE  train 0.6863424 0.6090035 0.6386169
#> 18    AR h=5     MdRAE  train 0.5802512 0.5521782 0.5767738
#> 19  true h=1      RMSE   test 0.6970850 2.4257264 1.7846675
#> 20   AR1 h=1      RMSE   test 0.9516563 2.6068528 1.9623113
#> 21    AR h=1      RMSE   test 0.7260734 2.5167108 1.8521631
#> 22  true h=5      RMSE   test 3.8720195 5.0745797 4.5135293
#> 23   AR1 h=5      RMSE   test 3.6869466 4.9815370 4.3823102
#> 24    AR h=5      RMSE   test 4.0055572 5.2691697 4.6802050
#> 25  true h=1       MAE   test 0.5570119 1.9958262 1.2764190
#> 26   AR1 h=1       MAE   test 0.7689192 2.1064866 1.4377029
#> 27    AR h=1       MAE   test 0.5982744 2.0370860 1.3176802
#> 28  true h=5       MAE   test 3.1962482 4.0687861 3.6325171
#> 29   AR1 h=5       MAE   test 2.9970296 3.9463451 3.4716873
#> 30    AR h=5       MAE   test 3.2446598 4.1651370 3.7048984
#> 31  true h=1     MdRAE   test 0.3119790 0.5488354 0.3990913
#> 32   AR1 h=1     MdRAE   test 0.3860908 0.6555889 0.4772513
#> 33    AR h=1     MdRAE   test 0.3212286 0.6320373 0.4341362
#> 34  true h=5     MdRAE   test 0.5076402 0.4327344 0.4643516
#> 35   AR1 h=5     MdRAE   test 0.4910610 0.3828607 0.4365129
#> 36    AR h=5     MdRAE   test 0.5129465 0.3895001 0.4723429
```

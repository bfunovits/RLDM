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
#>    criterion sample predictor      y[1]      y[2]     total
#> 1       RMSE  train       h=1  2.086576  2.284865  2.187968
#> 2        MAE  train       h=1  1.605899  1.821225  1.713562
#> 3         PB  train       h=1 60.000000 60.000000 60.000000
#> 4       RMSE   test       h=1  2.145560  2.332489  2.240974
#> 5        MAE   test       h=1  1.655692  1.820500  1.738096
#> 6         PB   test       h=1 51.111111 68.888889 60.000000
#> 7       RMSE  train       h=5  3.926034  6.024364  5.084619
#> 8        MAE  train       h=5  3.060212  4.707886  3.884049
#> 9         PB  train       h=5 55.555556 54.444444 55.000000
#> 10      RMSE   test       h=5  5.143429  7.940752  6.689933
#> 11       MAE   test       h=5  4.098242  6.317676  5.207959
#> 12        PB   test       h=5 53.333333 51.111111 52.222222

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
#> 1   true h=1      RMSE  train 2.1004481 2.2410891 2.1719073
#> 2    AR1 h=1      RMSE  train 2.1511088 3.5571022 2.9394085
#> 3     AR h=1      RMSE  train 2.0865760 2.2848650 2.1879680
#> 4   true h=5      RMSE  train 3.8872815 6.1180098 5.1254756
#> 5    AR1 h=5      RMSE  train 4.2226266 7.3215919 5.9764657
#> 6     AR h=5      RMSE  train 3.9260343 6.0243641 5.0846194
#> 7   true h=1       MAE  train 1.5911145 1.7118734 1.6514939
#> 8    AR1 h=1       MAE  train 1.7707579 2.6322062 2.2014820
#> 9     AR h=1       MAE  train 1.6058989 1.8212248 1.7135618
#> 10  true h=5       MAE  train 2.9814276 4.7912724 3.8863500
#> 11   AR1 h=5       MAE  train 3.2264289 5.4652213 4.3458251
#> 12    AR h=5       MAE  train 3.0602124 4.7078857 3.8840490
#> 13  true h=1     MdRAE  train 0.2102793 0.1273445 0.1805562
#> 14   AR1 h=1     MdRAE  train 0.2894324 0.1793023 0.2442877
#> 15    AR h=1     MdRAE  train 0.2351318 0.1570311 0.1972723
#> 16  true h=5     MdRAE  train 0.6447836 0.7140907 0.6620257
#> 17   AR1 h=5     MdRAE  train 0.6401162 0.8265275 0.7209322
#> 18    AR h=5     MdRAE  train 0.6116586 0.7971210 0.7208321
#> 19  true h=1      RMSE   test 2.1487276 2.3258331 2.2390321
#> 20   AR1 h=1      RMSE   test 2.1820723 3.4744623 2.9011487
#> 21    AR h=1      RMSE   test 2.1455597 2.3324893 2.2409744
#> 22  true h=5      RMSE   test 5.1315914 7.5050204 6.4287853
#> 23   AR1 h=5      RMSE   test 5.4420293 8.2481764 6.9874207
#> 24    AR h=5      RMSE   test 5.1434292 7.9407518 6.6899329
#> 25  true h=1       MAE   test 1.7040837 1.8336330 1.7688583
#> 26   AR1 h=1       MAE   test 1.7482740 2.8552381 2.3017561
#> 27    AR h=1       MAE   test 1.6556918 1.8205004 1.7380961
#> 28  true h=5       MAE   test 4.0923265 6.0219790 5.0571528
#> 29   AR1 h=5       MAE   test 4.4109470 6.5999849 5.5054660
#> 30    AR h=5       MAE   test 4.0982421 6.3176756 5.2079588
#> 31  true h=1     MdRAE   test 0.1675068 0.1563427 0.1598174
#> 32   AR1 h=1     MdRAE   test 0.1653810 0.2558399 0.2210482
#> 33    AR h=1     MdRAE   test 0.1732917 0.1466292 0.1488637
#> 34  true h=5     MdRAE   test 0.6740428 0.8374683 0.7425257
#> 35   AR1 h=5     MdRAE   test 0.6921266 0.9550540 0.8202070
#> 36    AR h=5     MdRAE   test 0.6882451 0.8602644 0.7571071
# Basic example
result <- evaluate_prediction()
#> Error in evaluate_prediction(): argument "y" is missing, with no default
result
#> Error: object 'result' not found
```

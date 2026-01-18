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
#> 1       RMSE  train       h=1  0.7464332  0.3151248  0.5729163
#> 2        MAE  train       h=1  0.6062208  0.2463317  0.4262762
#> 3         PB  train       h=1 51.1111111 53.3333333 52.2222222
#> 4       RMSE   test       h=1  0.8366101  0.3777231  0.6490729
#> 5        MAE   test       h=1  0.6926596  0.3086385  0.5006490
#> 6         PB   test       h=1 51.1111111 58.8888889 55.0000000
#> 7       RMSE  train       h=5  1.0158984  1.2603229  1.1446535
#> 8        MAE  train       h=5  0.8077797  0.9870565  0.8974181
#> 9         PB  train       h=5 58.8888889 61.1111111 60.0000000
#> 10      RMSE   test       h=5  1.3325549  1.6976224  1.5260446
#> 11       MAE   test       h=5  1.1201251  1.4008479  1.2604865
#> 12        PB   test       h=5 58.8888889 57.7777778 58.3333333

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
#> 1   true h=1      RMSE  train 0.7477685 0.3209361 0.5753945
#> 2    AR1 h=1      RMSE  train 0.7834196 0.3220445 0.5989403
#> 3     AR h=1      RMSE  train 0.7464332 0.3151248 0.5729163
#> 4   true h=5      RMSE  train 1.0201199 1.2700194 1.1518667
#> 5    AR1 h=5      RMSE  train 1.0610705 1.3159654 1.1953317
#> 6     AR h=5      RMSE  train 1.0158984 1.2603229 1.1446535
#> 7   true h=1       MAE  train 0.6067455 0.2514219 0.4290837
#> 8    AR1 h=1       MAE  train 0.6289440 0.2544288 0.4416864
#> 9     AR h=1       MAE  train 0.6062208 0.2463317 0.4262762
#> 10  true h=5       MAE  train 0.8103448 0.9980458 0.9041953
#> 11   AR1 h=5       MAE  train 0.8415633 1.0457366 0.9436500
#> 12    AR h=5       MAE  train 0.8077797 0.9870565 0.8974181
#> 13  true h=1     MdRAE  train 0.3716322 0.1601011 0.2134180
#> 14   AR1 h=1     MdRAE  train 0.4362900 0.1380195 0.2473105
#> 15    AR h=1     MdRAE  train 0.3802291 0.1474111 0.2204019
#> 16  true h=5     MdRAE  train 0.7182880 0.7885423 0.7817148
#> 17   AR1 h=5     MdRAE  train 0.8268266 0.9668042 0.8919780
#> 18    AR h=5     MdRAE  train 0.7301805 0.8001168 0.7869697
#> 19  true h=1      RMSE   test 0.8363547 0.3743463 0.6479292
#> 20   AR1 h=1      RMSE   test 0.8721171 0.3993802 0.6782672
#> 21    AR h=1      RMSE   test 0.8366101 0.3777231 0.6490729
#> 22  true h=5      RMSE   test 1.3292630 1.6959861 1.5236976
#> 23   AR1 h=5      RMSE   test 1.4364223 1.8304275 1.6452620
#> 24    AR h=5      RMSE   test 1.3325549 1.6976224 1.5260446
#> 25  true h=1       MAE   test 0.6920242 0.3055231 0.4987736
#> 26   AR1 h=1       MAE   test 0.7095644 0.3273742 0.5184693
#> 27    AR h=1       MAE   test 0.6926596 0.3086385 0.5006490
#> 28  true h=5       MAE   test 1.1179249 1.3852854 1.2516052
#> 29   AR1 h=5       MAE   test 1.2084904 1.5380311 1.3732607
#> 30    AR h=5       MAE   test 1.1201251 1.4008479 1.2604865
#> 31  true h=1     MdRAE   test 0.3303844 0.1110197 0.2029135
#> 32   AR1 h=1     MdRAE   test 0.3618924 0.1324323 0.2163938
#> 33    AR h=1     MdRAE   test 0.3462281 0.1021306 0.2055096
#> 34  true h=5     MdRAE   test 0.8782822 0.8480689 0.8562370
#> 35   AR1 h=5     MdRAE   test 0.8941141 0.9805714 0.9360010
#> 36    AR h=5     MdRAE   test 0.8840215 0.8551495 0.8651375
```

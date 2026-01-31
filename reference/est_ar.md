# Estimate Autoregressive Models

The function `est_ar` estimates (V)AR models \$\$(y_t - \mu) = a_1
(y\_{t-1} - \mu) + \cdots + a_p (y\_{t-p} - \mu) + u_t\$\$ from a given
sample or a given (sample) autocovariance function. The model order
\\p\\ is chosen by an information criterion, like *AIC* or *BIC*. The
"helper" functions `est_ar_ols`, `est_ar_yw` and `est_ar_dlw` implement
the three available estimation methods: estimation by ordinary least
squares, the Yule-Walker estimates and the Durbin-Levinson-Whittle
method.

## Usage

``` r
est_ar(
  obj,
  p.max = NULL,
  penalty = NULL,
  ic = c("AIC", "BIC", "max"),
  method = c("yule-walker", "ols", "durbin-levinson-whittle"),
  mean_estimate = c("sample.mean", "intercept", "zero"),
  n.obs = NULL
)

est_ar_yw(gamma, p.max = (dim(gamma)[3] - 1), penalty = -1)

est_ar_dlw(gamma, p.max = (dim(gamma)[3] - 1), penalty = -1)

est_ar_ols(
  y,
  p.max = NULL,
  penalty = -1,
  mean_estimate = c("sample.mean", "intercept", "zero"),
  p.min = 0L
)
```

## Arguments

- obj:

  either a "time series" object (i.e `as.matrix(obj)` returns an
  \\(N,m)\\-dimensional numeric matrix) or an
  [`autocov()`](https://bfunovits.github.io/RLDM/reference/autocov.md)
  object which represents an (estimated) autocovariance function. The
  type of the `autocov` object is irrelevant since `est_ar` always uses
  the slot `obj$gamma` which contains the autocovariance function.

- p.max:

  (integer or `NULL`) Maximum order of the candidate AR models. For the
  default choice see below.

- penalty:

  is a scalar (or NULL) which determines the "penalty" per parameter of
  the model. Note that this parameter (if not `NULL`) overrides the
  paramater `ic`.

- ic:

  (character string) Which information criterion shall be used to find
  the optimal order. Note that `ic="max"` means that an AR(p) model with
  `p=p.max` is estimated. Default is `ic="AIC"`.

- method:

  Character string giving the method used to fit the model. Note that
  'yule-walker' and 'durbin-levinson-whittle' are (up to numerical
  errors) equivalent and that the choice 'ols' is only available for a
  "time-series" object `obj`.

- mean_estimate:

  Character string giving the method used to estimate the mean \\\mu\\.
  Default is `mean_estimate = "sample.mean"`. See the details below.

- n.obs:

  Optional integer which gives the sample size \\N\\. This parameter is
  only used, when `obj` is an `autocov` object. If `n.obs=NULL` then the
  slot `obj$n.obs` is used. Note that `obj$n.obs=NULL` or
  `obj$n.obs=Inf` refers to the case of a population autocovariance
  function, i.e. \\N=\infty\\.  
  For a "time series" object the sample size is of course set to the
  number of observations, i.e. `n.obs = nrow(as.matrix(obj))`.  
  The sample size \\N\\ controls the computation of the default maximum
  order `p.max` and the computation of the information criterion.

- gamma:

  \\(m,m,lag.max+1)\\-dimensional array, which contains the (sample)
  autocovariance function.

- y:

  \\(N,m)\\-dimensional matrix, which contains the sample.

- p.min:

  (non negative integer) Minimum order of the candidate AR models. Only
  used by `est_ar_ols`.

## Value

The function `est_ar` returns a list with components

- model:

  [`armamod()`](https://bfunovits.github.io/RLDM/reference/armamod.md)
  object which represents the estimated AR model.

- p:

  optimal model order.

- stats:

  (p.max+1,4) dimensional matrix which stores the \\\ln\det(\Sigma_p)\\
  values, the number of parameters and the IC values. See the details
  below.

- y.mean:

  estimate of the mean \\\mu\\.

- ll:

  The log likelihood of the estimated model.

The "helper" functions `est_ar_yw`, `est_ar_dlw` and `est_ar_ols` return
a list with components

- a:

  `(m,m,p)`-dimensional array with the estimated AR coefficients
  \\a_i\\.

- sigma:

  `(m,m)`-dimensional matrix with the estimated noise covariance
  \\\Sigma\\.

- p:

  estimate of the AR order.

- stats:

  (p.max+1,4) dimensional matrix which stores the \\\ln\det(\Sigma_p)\\
  values, the number of parameters and the IC values. See the details
  below.

- y.mean:

  (`est_ar_ols` only) estimate of the mean \\\mu\\.

- residuals:

  (`est_ar_ols` only) `(n.obs,m)` dimensional matrix with the OLS
  residuals.

- partial:

  (`est_ar_dlw` only) `(m,m,p.max+1)` dimensional array with the
  (estimated) partial autocorrelation coefficients.

## OLS method

The helper function `est_ar_ols` implements three schemes to estimate
the mean \\\mu\\ and the AR parameters. The choice
`mean_estimate = "zero"` assumes \\\mu=0\\ and thus the AR parameters
are determined from the regression: \$\$y_t = a_1 y\_{t-1} + \cdots +
a_p y\_{t-p} + u_t \mbox{ for } t=p+1,\ldots,N\$\$ In the case
`mean_estimate = "sample.mean"` the mean \\\mu\\ is estimated by the
sample mean and the AR parameters are determined by the LS estimate of
the regression \$\$(y_t - \mu) = a_1 (y\_{t-1} - \mu) + \cdots + a_p
(y\_{t-p} - \mu) + u_t \mbox{ for } t=p+1,\ldots,N\$\$ In the last case
`mean_estimate = "intercept"`, a regression with intercept \$\$y_t = d +
a_1 y\_{t-1} + \cdots + a_p y\_{t-p} + u_t \mbox{ for }
t=p+1,\ldots,N\$\$ is considered. The estimate for \\\mu\\ then is
obtained as \$\$\mu = (I_m - a_1 - \cdots - a_p)^{-1} d\$\$ This
estimate of the mean \\\mu\\ fails if the estimated AR model has a *unit
root*, i.e. if \\(I_m - a_1 - \cdots - a_p)\\ is singular.

The sample covariance of the corresponding residuals (scaled with
\\1/(N-p)\\) serves as an estimate for the noise covariance \\\Sigma\\.

For the actual computations the routine
[`stats::lsfit()`](https://rdrr.io/r/stats/lsfit.html) in the stats
package is used.

## Yule-Walker estimates

Both `est_ar_yw` and `est_ar_dlw` use the *Yule-Walker* equations to
estimate the AR coefficients \\(a_i)\\ and the noise covariance matrix
\\\Sigma\\. However, they use a different numerical scheme to solve
these equations. The function `est_ar_dlw` uses the
*Durbin-Levinson-Whittle* recursions and in addition returns (estimates
of) the partial autocorrelation coefficients.

If `obj` is a "time series" object, then first the ACF is estimated with
a call to
[`autocov()`](https://bfunovits.github.io/RLDM/reference/autocov.md).
The option `mean_estimate = "zero"` implies that the mean is assumed to
be zero (\\\mu = 0\\) and therefore `autocov` is called with the option
`demean = FALSE`.

For `mean_estimate = "sample.mean"` or `mean_estimate = "intercept"` the
mean \\\mu\\ is estimated by the sample mean and the ACF is computed
with `demean = TRUE`.

## Estimation of the AR order

The order \\p\\ of the AR model is chosen by minimizing an information
criterion of the form \$\$IC(p) = \ln\det\Sigma_p + c(p)r(N) \mbox{ for
} p = 0,\ldots,p\_{\max}\$\$ where \\\Sigma_p\\ is the estimate of the
noise (innovation) covariance, \\c(p)\\ counts the number of parameters
of the model, and \\r(N)\\ is the "penalty" per parameter of the model.
Note that \\\log\det\Sigma\\ is up to a constant and a scaling factor
\\-(N-p)/2\\ equal to the (scaled, approximate) Gaussian log likelihood
of the model \$\$ll = -(1/2)(m \ln(2\pi) + m + \ln\det \Sigma_p)\$\$ See
also [`ll()`](https://bfunovits.github.io/RLDM/reference/ll.md). Note
that the value \\ll\\, which is returned by this routine, is the
(approximate) log Likelihood **scaled** by a factor \\1/(N-p)\\.

For an AR(p) model with intercept, the number of parameters is \\c(p) =
p m^2 + m\\ and for the AR model without intercept \\c(p) = p m^2\\.

The Akaike information criterion (AIC) corresponds to \\r(N)=2/N\\ and
the Bayes information criterion (BIC) uses the penatlty
\\r(N)=\log(N)/N\\.

For the helper routines, the user has to set the penalty term \\r(N)\\
explicitly via the input parameter "`penalty`". The default choice
`penalty = -1` means that the maximum possible order `p=p.max` is
chosen.

The function `est_ar` offers the parameter "`ic`" which tells the
routine to set the penalty accordingly. Note that the choice `ic="max"`
sets \\r(N) = -1\\ and thus again the model with maximum possible order
is fitted.

The default maximum order `p.max` is chosen as follows. The helper
functions `est_ar_yw` and `est_ar_dlw` simply chose the maximum accordng
to maximum lag of the given autocovariances,
`p.max = dim(gamma)[3] - 1`. The routine `est_ar_ols` uses the minimum
of \\12\\, \\(N-1)/(m+1)\\ and \\10\*log10(N)\\ as default. The function
`est_ar` uses the same value. However, if "`obj`" is an `autocov` object
then `p.max` is in addition bounded by the number of lags contained in
this object.

## Notes

The Yule-Walker estimates offer an easy way to reconstruct the "true"
model if the population autocovariance function is given. The noise
covariance (and thus the likelihood values) should not improve when a
model with an order larger than the true model order is "estimated".
However due to numerical errors this may not be true. As a simple trick
one may call `est_ar` (`est_ar_yw` or `est_ar_dlw`) with a very small
positive `penalty`. See the example below.

The functions are essentially equivalent to the stats routines. They are
(re) implemented for convenience, such that the input and output
parameters (models) fit to the RLDM conventions.

The AIC values of RLDM routines are equivalent to the AIC values
computed by the stats routines up to a constant and up to scaling by
\\N\\.

It seems that the Yule-Walker estimate `stats::[ar.yw][stats::ar.yw]`
uses a scaling factor \\(N - m(p+1))/N\\ for the noise covariance
\\\Sigma\\.

Finally note that `est_ar_ols`, `est_ar_yw` and `est_ar_dlw` are mainly
intended as "internal helper" functions. Therefore, these functions do
not check the validity of the input parameters.

## Examples

``` r
# set seed, to get reproducable results
set.seed(5436)

###############################################################
# generate a (bivariate) random, stable AR(3) model

m = 2
p = 3
n.obs = 100
p.max = 10
tmpl = tmpl_arma_pq(m = m, n = m, p = p, q = 0)
model = r_model(tmpl, bpoles = 1, sd = 0.25)
# make sure that the diagonal entries of sigma_L are non negative
model$sigma_L = model$sigma_L %*% diag(sign(diag(model$sigma_L)))

###############################################################
# reconstruct the true AR model from the population ACF

true_acf = autocov(model, lag.max = 12, type = 'covariance')
ARest = est_ar(true_acf, p.max = p.max, method = 'yule-walker', penalty = 1e-6)
all.equal(model, ARest$model)
#> [1] TRUE

###############################################################
# simulate a sample

y = sim(model, n.obs = n.obs, start = list(s1 = NA))$y

###############################################################
# estimate the AR(p) model with the true order p

# OLS
ARest = est_ar(y, ic = 'max', p.max = p, method = 'ols', mean_estimate = "zero")
# check the log Likelihood
p.opt = ARest$p
all.equal(ll(ARest$model, y, 'conditional', skip = p.opt), ARest$ll)
#> [1] TRUE

# Yule-Walker and Durbin-Levinson-Whittle are equivalent (up to numerical errors)
ARest = est_ar(y, ic = 'max', p.max = p, method = 'yule-walker', mean_estimate = "zero")
junk = est_ar(y, ic = 'max', p.max = p, method = 'durbin-levinson-whittle', mean_estimate = "zero")
all.equal(ARest$model, junk$model)
#> [1] TRUE

# alternatively we may first estimate the sample autocovariance function
# note that the 'type' of the ACF is irrelevant
sample_acf = autocov(y, type = 'correlation', demean = FALSE)
junk = est_ar(sample_acf, ic = 'max', p.max = p, method = 'yule-walker')
all.equal(ARest$model, junk$model)
#> [1] TRUE

###############################################################
# estimate the order of the AR model with 'AIC', estimate a model with intercept
ARest = est_ar(y, ic = 'AIC', p.max = p.max, method = 'ols', mean_estimate = "intercept")
print(ARest$model)
#> ARMA model [2,2] with orders p = 3 and q = 0
#> AR polynomial a(z):
#>      z^0 [,1]  [,2]  z^1 [,1]       [,2]  z^2 [,1]       [,2]  z^3 [,1]
#> [1,]        1     0 0.1758752 -0.2469383 0.4064368 -0.2595963 0.7758289
#> [2,]        0     1 0.3113271 -0.3092643 0.3315128 -0.1661434 0.3077410
#>            [,2]
#> [1,] -0.2624946
#> [2,] -0.2492206
#> MA polynomial b(z):
#>      z^0 [,1]  [,2]
#> [1,]        1     0
#> [2,]        0     1
#> Left square root of noise covariance Sigma:
#>          u[1]      u[2]
#> u[1] 0.144723 0.0000000
#> u[2] 0.188681 0.1501011

# compare with the stats::ar function
ARest2 = stats::ar(y, aic = TRUE, order.max = p.max, method = 'ols', intercept = TRUE)
# the estimated coefficients are equal
all.equal(unclass(ARest$model$sys$a)[,,-1], -aperm(ARest2$ar, c(2,3,1)), check.attributes = FALSE)
#> [1] TRUE
# also the AIC values are up to scaling equivalent
all.equal( ARest$stats[,'ic'] - min(ARest$stats[,'ic']), ARest2$aic/n.obs, check.attributes = FALSE)
#> [1] TRUE

# reset seed
set.seed(NULL)
set.seed(123)
# Generate example data
model <- test_stspmod(dim = c(2,2), s = 2, bpoles = 1, sigma_L = diag(2))
y <- sim(model, n.obs = 100)$y

# Compute autocovariance
acf_obj <- autocov(model, lag.max = 10)
gamma <- acf_obj$gamma

# Run estimation
result <- est_ar_yw(gamma)
result
#> $a
#> , , 1
#> 
#>           [,1]       [,2]
#> [1,]  0.593593  0.4029055
#> [2,] -1.634472 -0.4848091
#> 
#> , , 2
#> 
#>            [,1]       [,2]
#> [1,] 0.03294047 -0.1223931
#> [2,] 0.08994749  0.5193903
#> 
#> , , 3
#> 
#>             [,1]        [,2]
#> [1,] -0.04003777 -0.04905823
#> [2,]  0.13755402  0.11528428
#> 
#> , , 4
#> 
#>             [,1]        [,2]
#> [1,] -0.00855073  0.00186596
#> [2,]  0.01397270 -0.02196229
#> 
#> , , 5
#> 
#>              [,1]         [,2]
#> [1,]  0.001743794  0.004002599
#> [2,] -0.008319666 -0.012063094
#> 
#> , , 6
#> 
#>               [,1]          [,2]
#> [1,]  0.0009120883  0.0004635627
#> [2,] -0.0023173558 -0.0001580165
#> 
#> , , 7
#> 
#>              [,1]          [,2]
#> [1,] 5.626719e-06 -0.0002334706
#> [2,] 2.806182e-04  0.0008908700
#> 
#> , , 8
#> 
#>               [,1]         [,2]
#> [1,] -6.832243e-05 -0.000070400
#> [2,]  2.179598e-04  0.000146714
#> 
#> , , 9
#> 
#>               [,1]          [,2]
#> [1,] -1.065933e-05  7.103105e-06
#> [2,]  1.149316e-05 -4.403497e-05
#> 
#> , , 10
#> 
#>               [,1]          [,2]
#> [1,]  3.777965e-06  5.830430e-06
#> [2,] -1.240687e-05 -1.583706e-05
#> 
#> 
#> $sigma
#>              [,1]         [,2]
#> [1,]  1.00000e+00 -1.09468e-12
#> [2,] -1.09468e-12  1.00000e+00
#> 
#> $p
#> [1] 10
#> 
#> $stats
#>        p n.par   lndetSigma         ic
#>  [1,]  0     0 1.611048e+00   1.611048
#>  [2,]  1     4 5.164767e-01  -3.483523
#>  [3,]  2     8 1.482537e-02  -7.985175
#>  [4,]  3    12 1.188096e-03 -11.998812
#>  [5,]  4    16 1.702456e-04 -15.999830
#>  [6,]  5    20 2.166470e-06 -19.999998
#>  [7,]  6    24 1.087133e-06 -23.999999
#>  [8,]  7    28 2.535586e-08 -28.000000
#>  [9,]  8    32 3.970023e-09 -32.000000
#> [10,]  9    36 3.657306e-10 -36.000000
#> [11,] 10    40 7.634338e-12 -40.000000
#> 
set.seed(123)
# Generate example data
model <- test_stspmod(dim = c(2,2), s = 2, bpoles = 1, sigma_L = diag(2))
y <- sim(model, n.obs = 100)$y

# Compute autocovariance
acf_obj <- autocov(model, lag.max = 10)
gamma <- acf_obj$gamma

# Run estimation
result <- est_ar_dlw(gamma)
result
#> $a
#> , , 1
#> 
#>           [,1]       [,2]
#> [1,]  0.593593  0.4029055
#> [2,] -1.634472 -0.4848091
#> 
#> , , 2
#> 
#>            [,1]       [,2]
#> [1,] 0.03294047 -0.1223931
#> [2,] 0.08994749  0.5193903
#> 
#> , , 3
#> 
#>             [,1]        [,2]
#> [1,] -0.04003777 -0.04905823
#> [2,]  0.13755402  0.11528428
#> 
#> , , 4
#> 
#>             [,1]        [,2]
#> [1,] -0.00855073  0.00186596
#> [2,]  0.01397270 -0.02196229
#> 
#> , , 5
#> 
#>              [,1]         [,2]
#> [1,]  0.001743794  0.004002599
#> [2,] -0.008319666 -0.012063094
#> 
#> , , 6
#> 
#>               [,1]          [,2]
#> [1,]  0.0009120883  0.0004635627
#> [2,] -0.0023173558 -0.0001580165
#> 
#> , , 7
#> 
#>              [,1]          [,2]
#> [1,] 5.626719e-06 -0.0002334706
#> [2,] 2.806182e-04  0.0008908700
#> 
#> , , 8
#> 
#>               [,1]         [,2]
#> [1,] -6.832243e-05 -0.000070400
#> [2,]  2.179598e-04  0.000146714
#> 
#> , , 9
#> 
#>               [,1]          [,2]
#> [1,] -1.065933e-05  7.103105e-06
#> [2,]  1.149316e-05 -4.403497e-05
#> 
#> , , 10
#> 
#>               [,1]          [,2]
#> [1,]  3.777965e-06  5.830430e-06
#> [2,] -1.240687e-05 -1.583706e-05
#> 
#> 
#> $sigma
#>               [,1]          [,2]
#> [1,]  1.000000e+00 -1.095188e-12
#> [2,] -1.095188e-12  1.000000e+00
#> 
#> $p
#> [1] 10
#> 
#> $stats
#>        p n.par   lndetSigma         ic
#>  [1,]  0     0 1.611048e+00   1.611048
#>  [2,]  1     4 5.164767e-01  -3.483523
#>  [3,]  2     8 1.482537e-02  -7.985175
#>  [4,]  3    12 1.188096e-03 -11.998812
#>  [5,]  4    16 1.702456e-04 -15.999830
#>  [6,]  5    20 2.166470e-06 -19.999998
#>  [7,]  6    24 1.087133e-06 -23.999999
#>  [8,]  7    28 2.535586e-08 -28.000000
#>  [9,]  8    32 3.970023e-09 -32.000000
#> [10,]  9    36 3.657310e-10 -36.000000
#> [11,] 10    40 7.635448e-12 -40.000000
#> 
#> $partial
#> , , 1
#> 
#>            [,1]       [,2]
#> [1,]  1.0000000 -0.4553269
#> [2,] -0.4553269  1.0000000
#> 
#> , , 2
#> 
#>            [,1]        [,2]
#> [1,]  0.2531218  0.35658618
#> [2,] -0.6359387 -0.09203067
#> 
#> , , 3
#> 
#>            [,1]       [,2]
#> [1,]  0.1239848 -0.2030783
#> [2,] -0.3266716  0.6017541
#> 
#> , , 4
#> 
#>              [,1]        [,2]
#> [1,]  0.010626244 -0.04270662
#> [2,] -0.005414358  0.08320880
#> 
#> , , 5
#> 
#>             [,1]         [,2]
#> [1,] -0.00493339  0.005507714
#> [2,]  0.01907464 -0.031297476
#> 
#> , , 6
#> 
#>              [,1]         [,2]
#> [1,] -0.001546243  0.004045939
#> [2,]  0.003282471 -0.011009282
#> 
#> , , 7
#> 
#>               [,1]         [,2]
#> [1,]  0.0001401522 0.0001942242
#> [2,] -0.0009497122 0.0007089075
#> 
#> , , 8
#> 
#>               [,1]          [,2]
#> [1,]  0.0001384386 -0.0002774523
#> [2,] -0.0003924843  0.0009420754
#> 
#> , , 9
#> 
#>              [,1]          [,2]
#> [1,] 1.029896e-05 -5.667255e-05
#> [2,] 1.269272e-05  8.865205e-05
#> 
#> , , 10
#> 
#>               [,1]          [,2]
#> [1,] -8.942113e-06  1.247439e-05
#> [2,]  3.168432e-05 -5.804427e-05
#> 
#> , , 11
#> 
#>               [,1]          [,2]
#> [1,] -2.132853e-06  6.183752e-06
#> [2,]  3.827159e-06 -1.549790e-05
#> 
#> 
set.seed(123)
# Generate example data
model <- test_stspmod(dim = c(2,2), s = 2, bpoles = 1, sigma_L = diag(2))
y <- sim(model, n.obs = 100)$y

# Run estimation
result <- est_ar_ols(y)
result
#> $a
#> , , 1
#> 
#>            [,1]       [,2]
#> [1,]  0.5216462  0.4730037
#> [2,] -1.4383693 -0.7654621
#> 
#> , , 2
#> 
#>             [,1]      [,2]
#> [1,] -0.03058286 0.0201262
#> [2,] -0.57701683 0.1919959
#> 
#> , , 3
#> 
#>            [,1]      [,2]
#> [1,] 0.03399671 0.2621479
#> [2,] 0.16838782 0.2851212
#> 
#> , , 4
#> 
#>           [,1]         [,2]
#> [1,] 0.5028358  0.191399797
#> [2,] 0.1098366 -0.009784647
#> 
#> , , 5
#> 
#>             [,1]        [,2]
#> [1,]  0.06114687 -0.14561419
#> [2,] -0.03570935 -0.04731226
#> 
#> , , 6
#> 
#>             [,1]       [,2]
#> [1,] -0.16753069 -0.1340659
#> [2,]  0.03861591 -0.1099899
#> 
#> , , 7
#> 
#>            [,1]       [,2]
#> [1,] -0.2222810 0.06017286
#> [2,]  0.1619097 0.03008335
#> 
#> , , 8
#> 
#>             [,1]        [,2]
#> [1,] -0.09662243  0.09499039
#> [2,]  0.04101294 -0.08523160
#> 
#> , , 9
#> 
#>              [,1]        [,2]
#> [1,] -0.005804168  0.03833865
#> [2,] -0.131761868 -0.09399973
#> 
#> , , 10
#> 
#>             [,1]        [,2]
#> [1,] -0.12243313 -0.01736716
#> [2,] -0.05457868 -0.07349516
#> 
#> 
#> $sigma
#>             Y1          Y2
#> Y1  0.72666511 -0.05559743
#> Y2 -0.05559743  0.66870024
#> 
#> $p
#> [1] 10
#> 
#> $stats
#>        p n.par  lndetSigma         ic
#>  [1,]  0     0  1.20445544   1.204455
#>  [2,]  1     4 -0.08870917  -4.088709
#>  [3,]  2     8 -0.30410303  -8.304103
#>  [4,]  3    12 -0.37670466 -12.376705
#>  [5,]  4    16 -0.43565670 -16.435657
#>  [6,]  5    20 -0.52199874 -20.521999
#>  [7,]  6    24 -0.61678161 -24.616782
#>  [8,]  7    28 -0.67075827 -28.670758
#>  [9,]  8    32 -0.70783818 -32.707838
#> [10,]  9    36 -0.72776027 -36.727760
#> [11,] 10    40 -0.72809053 -40.728091
#> 
#> $y.mean
#> [1]  0.02165726 -0.09609867
#> 
#> $residuals
#>               [,1]         [,2]
#>   [1,]          NA           NA
#>   [2,]          NA           NA
#>   [3,]          NA           NA
#>   [4,]          NA           NA
#>   [5,]          NA           NA
#>   [6,]          NA           NA
#>   [7,]          NA           NA
#>   [8,]          NA           NA
#>   [9,]          NA           NA
#>  [10,]          NA           NA
#>  [11,]  0.75257328  2.198162516
#>  [12,] -0.91760186 -1.032808908
#>  [13,]  0.76632455 -0.606786757
#>  [14,] -1.06579517  0.251759752
#>  [15,] -0.17979718 -0.245658361
#>  [16,]  0.27724991  0.080232787
#>  [17,] -0.37302269  0.341120830
#>  [18,] -0.50015831  0.072468130
#>  [19,] -0.01569907  1.142532618
#>  [20,]  0.96697433  0.751541300
#>  [21,] -0.55902771  1.216187146
#>  [22,]  0.80665816  1.202537289
#>  [23,] -0.02411356 -0.555569966
#>  [24,]  1.03350698 -0.651661352
#>  [25,]  2.17157329  1.089114268
#>  [26,]  0.16547614 -0.764793440
#>  [27,]  0.21409777  0.008302741
#>  [28,] -0.04536647 -0.641723277
#>  [29,] -1.01896803 -0.537498856
#>  [30,] -0.35098949 -1.625324644
#>  [31,]  0.25812793  0.262770086
#>  [32,] -0.15655903 -0.246093154
#>  [33,] -0.80472319  0.549936862
#>  [34,]  0.27195315  0.971885677
#>  [35,] -0.59972480 -0.389902569
#>  [36,] -1.55592123 -0.948105060
#>  [37,] -0.39385286 -0.280085604
#>  [38,] -1.12087954 -0.391763102
#>  [39,]  1.25698575 -0.188953877
#>  [40,] -0.01136385  0.384438885
#>  [41,] -1.17009648  0.369025848
#>  [42,]  1.62920200  0.656971580
#>  [43,] -0.51412770 -0.320057101
#>  [44,] -2.22229902  1.673256300
#>  [45,] -1.48699674  1.268315893
#>  [46,]  1.12619553 -1.653433004
#>  [47,]  0.17349536 -1.006982326
#>  [48,] -1.25774385 -0.890962103
#>  [49,] -0.64212986 -1.011236836
#>  [50,] -1.35736696  0.398847416
#>  [51,]  1.13742110 -0.078472453
#>  [52,]  0.05139381  0.886701723
#>  [53,]  0.33062103 -1.002996119
#>  [54,]  0.44116995 -0.784823762
#>  [55,]  0.35244538  0.324631055
#>  [56,]  0.70706873  0.299612232
#>  [57,]  0.68436188 -1.210781876
#>  [58,] -1.15860514  2.101469622
#>  [59,] -0.25347648  0.815921861
#>  [60,]  0.68745570 -0.005964677
#>  [61,]  0.02590503  0.246083243
#>  [62,] -0.04637627 -0.045042374
#>  [63,]  0.59568567  1.562904923
#>  [64,] -0.45633772 -0.593718728
#>  [65,] -0.22115052  0.444769819
#>  [66,] -0.26259563 -0.375306481
#>  [67,] -1.01410842  1.149819586
#>  [68,] -0.35747275 -0.329157192
#>  [69,] -0.20354839  0.141739167
#>  [70,]  0.41248466 -0.064779494
#>  [71,]  0.84761359  0.065686675
#>  [72,]  0.35799828 -0.374802188
#>  [73,]  0.40488563 -0.750182143
#>  [74,] -1.33754181  1.684072263
#>  [75,]  0.15009271 -0.068410240
#>  [76,] -0.81437248 -1.435180746
#>  [77,]  1.69648596  0.776044634
#>  [78,]  0.26896221  0.175707511
#>  [79,]  0.35362148  0.127658253
#>  [80,] -0.37286734 -0.463701928
#>  [81,]  0.97945371 -0.262326806
#>  [82,] -0.21886162 -0.327893141
#>  [83,]  1.78950474  0.120786512
#>  [84,] -0.59661544  0.524661302
#>  [85,] -0.06882649 -0.161662442
#>  [86,] -1.32609204 -1.274967367
#>  [87,] -1.03635755  0.947246201
#>  [88,]  0.92884943  0.159311594
#>  [89,] -0.15720891 -0.020625734
#>  [90,] -0.41962457 -0.329405759
#>  [91,]  0.92280620 -0.870835579
#>  [92,]  1.51821116 -0.423700676
#>  [93,]  0.36918410 -0.326766709
#>  [94,] -0.18481198 -1.162690239
#>  [95,] -0.06802760 -0.553582958
#>  [96,]  0.09457118 -0.903869599
#>  [97,] -0.74957548 -0.176414690
#>  [98,]  1.89028643 -1.196545275
#>  [99,]  0.42116926  0.235887216
#> [100,]  0.49989413 -0.953022199
#> 
```

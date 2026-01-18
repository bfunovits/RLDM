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
```

# Autocovariance, Autocorelation and Partial Autocorrelation Function

Compute respectively estimate the autocovariance, autocorrelation or
partial autocorrelation function of a stationary process.

## Usage

``` r
autocov(obj, type, ...)

# Default S3 method
autocov(
  obj,
  type = c("covariance", "correlation", "partial"),
  lag.max = NULL,
  na.action = stats::na.fail,
  demean = TRUE,
  ...
)

# S3 method for class 'autocov'
autocov(obj, type = c("covariance", "correlation", "partial"), ...)

# S3 method for class 'armamod'
autocov(
  obj,
  type = c("covariance", "correlation", "partial"),
  lag.max = 12,
  ...
)

# S3 method for class 'stspmod'
autocov(
  obj,
  type = c("covariance", "correlation", "partial"),
  lag.max = 12,
  ...
)
```

## Arguments

- obj:

  either a
  [`armamod()`](https://bfunovits.github.io/RLDM/reference/armamod.md),
  [`stspmod()`](https://bfunovits.github.io/RLDM/reference/stspmod.md),
  `autocov()` object or a "data" object.

- type:

  character string giving the type of acf to be computed. Allowed values
  are "covariance" (the default), "correlation", or "partial". Will be
  partially matched. Note that the default value here is "covariance"
  whereas [`stats::acf()`](https://rdrr.io/r/stats/acf.html) uses
  "correlation" as default.

- ...:

  not used.

- lag.max:

  (integer) maximum lag.

- na.action:

  function to be called to handle missing values.
  [`stats::na.pass()`](https://rdrr.io/r/stats/na.fail.html) can be
  used.

- demean:

  logical. Should the covariances be about the sample means?

## Value

`autocov` object, i.e. a list with slots

- acf:

  `rationalmatrices::pseries()` object, which stores the covariances
  (correlations).

- type:

  character string which indicates the type of the ACF.

- gamma:

  (m,m,lag.max+1) dimensional array which stores the autocovariance
  function.

- names:

  (m)-dimensional character vector or NULL. This optional slot stores
  the names for the components of the time series/process.

- label:

  character string or NULL.

- n.obs:

  integer or NULL. This slot stores the sample size.

## Details

The class of the input parameter "`obj`" determines the S3 method called
and hence what is actually computed.

**Population ACF:**

If "`obj`" is an
[`armamod()`](https://bfunovits.github.io/RLDM/reference/armamod.md) or
[`stspmod()`](https://bfunovits.github.io/RLDM/reference/stspmod.md)
object then `autocov(obj, ...)` computes the ACF of the corresponding
stationary process.

Note however, that the function returns nonsense, if the model does not
satisfy the *stability* condition.

**Change the type of an ACF:**

Calling `autocov(obj, type)`, where "`obj`" is an `autocov` object
returns an ACF of the desired type. E.g. if "`obj`" holds a partial
autocorrelation function then `autocov(obj, type = 'covariance')` may be
used to retrieve the corresponding autocovariance function.

This is possible since the `autocov` object stores the "original"
autocovariances in a slot named `gamma`.

**Sample ACF:**

The default S3 method estimates the ACF from given data. It assumes that
"`obj`" is a univariate or multivariate numeric time series object, a
(numeric) data frame or a (numeric) vector, respectively matrix and then
simply calls the function
[`stats::acf()`](https://rdrr.io/r/stats/acf.html) in the stats package
to compute the sample autocovariance function. If needed, then the
corresponding sample autocorrelation, respectively sample partial
autocorrelation function is computed (and returned).

The syntax is quite analogous to
[`stats::acf()`](https://rdrr.io/r/stats/acf.html), so please consider
the documentation of [`stats::acf()`](https://rdrr.io/r/stats/acf.html)
for more details.

Note that stats stores autocovariance/autocorrelation functions as
`(lag.max+1,m,m)` dimensional arrays, whereas RLDM uses
`(m,m,lag.max+1)` dimensional arrays.

The definition of partial autocorrelations used by
[`stats::acf()`](https://rdrr.io/r/stats/acf.html) differs from the
definition used here, see e.g. (Reinsel 1997) . Furthermore
[`stats::acf()`](https://rdrr.io/r/stats/acf.html) skips the lag zero
partial autocorrelation coefficient and thus the pacf computed by
[`stats::acf()`](https://rdrr.io/r/stats/acf.html) is (lag.max,n,n)
dimensional.

The default choice for the number of lags is \\10\*log10(N/m)\\ where
\\N\\ is the number of observations and \\m\\ the number of series. This
number will be automatically limited to one less than the number of
observations in the series.

## References

Reinsel GC (1997). *Elements of Multivariate Time Series Analysis*.
Springer.

## See also

The autocovariance function of (V)ARMA processes may also be computed by
[`stats::ARMAacf()`](https://rdrr.io/r/stats/ARMAacf.html) in the scalar
case and by `MTS::VARMAcov()` in the multivariate case (m \> 1).

As noted above the sample ACF is computed via the
[`stats::acf()`](https://rdrr.io/r/stats/acf.html) routine in the stats
package.

## Examples

``` r
model = stspmod(sys = stsp(A = c(0,0.2,1,-0.5), B = c(1,1,1,-1),
                           C = c(1,0,0,1)), sigma_L = diag(c(4,1)),
                names = c('y1','y2'), label = 'test model')
g = autocov(model, lag.max=10)       # ACF
r = autocov(model, lag.max=10, type = 'correlation')  # autocorrelation function
r = autocov(g, type = 'correlation')                  # this gives the same result!
c = autocov(r, type = 'partial')     # partial autocorrelation function

if (FALSE) { # \dontrun{
# consider an equivalent VARMA model
model2 = impresp2varma(irf(model, lag.max = 20))$model
g2 = autocov(model2, lag.max = 10)
all.equal(g,g2) # of course both return the same ACF

autocov(matrix(rnorm(100*2), nrow = 100))
autocov(stspmod(test_stsp(dim = c(2,2), s = 2), sigma_L = diag(2)))

# generate a random sample with 100 observations and 3 outputs/series.
x = matrix(rnorm(100*3), nrow = 100, ncol = 3)

# the covariance estimates are of course identical
stats_acfobj = stats::acf(x, type = 'covariance', demean = TRUE, plot = FALSE)
rldm_acfobj  = acf_estimate(x, type = 'covariance', demean = TRUE)
testthat::expect_equivalent(rldm_acfobj$gamma, aperm(stats_acfobj$acf,c(2,3,1)))

# also the correlation estimates are identical
stats_acfobj = stats::acf(x, type = 'correlation', demean = TRUE, plot = FALSE)
rldm_acfobj  = acf_estimate(x, type = 'correlation', demean = TRUE)
testthat::expect_equivalent(rldm_acfobj$gamma, aperm(stats_acfobj$acf,c(2,3,1)))

# However, the partial correlations dont match!
stats_acfobj = stats::acf(x, type = 'partial', demean = TRUE, plot = FALSE)
rldm_acfobj  = acf_estimate(x, type = 'partial', demean = TRUE)
testthat::expect_equivalent(rldm_acfobj$gamma[,,-1,drop=FALSE], aperm(stats_acfobj$acf,c(2,3,1)))
} # }
```

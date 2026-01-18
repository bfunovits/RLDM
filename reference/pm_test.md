# Portmanteau Test for Serial Correlation

Test whether the residuals of an estimated model are serially
correlated. The test statistic is \$\$Q = N^2\sum\_{k=1}^{K}
(N-k)^{-1}\mbox{tr} (G_k G_0^{-1} G_k' G_0^{-1})\$\$ where \\G_k\\ are
the sample covariances of the residuals. Under the Null of a correctly
specified and estimated model the test statistic is asmyptotically
Chi-squared distributed with \\Km^2-\kappa\\ degrees of freedom, where
\\\kappa\\ is the number of (free) parameters of the model (class).

## Usage

``` r
pm_test(u, lag.max, n.par)
```

## Arguments

- u:

  (N-by-m) matrix of residuals (or an object which may be coerced to a
  matrix with `as.matrix(u)`).

- lag.max:

  (integer) maximum number of lags.

- n.par:

  (integer) number of parameters of the estimated model.

## Value

Matrix with four columns ("`lags`" number of lags, "`df`" degrees of
freedom, "`Q`" test statistics and "`p`" p values).

## Examples

``` r
u = matrix(rnorm(100*3), nrow = 100, ncol = 3)
pm_test(u, 4, 0)
#>      lags         Q df         p
#> [1,]    1  8.023592  9 0.5317705
#> [2,]    2 15.156435 18 0.6511991
#> [3,]    3 19.834452 27 0.8376311
#> [4,]    4 27.842198 36 0.8327771
```

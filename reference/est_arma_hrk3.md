# Different version of HRK Procedure

See
[est_arma_hrk](https://bfunovits.github.io/RLDM/reference/est_arma_hrk.md).
Stage III of the Hannan-Rissanen-Kavalieris procedure is implemented
here as well. This function returns the best model (since the iterations
might not always improve the log-likelihood value) and allows for
returning results at different stages of the HRK procedure.  
One notable differences is that the data needs to be demeaned because we
use Yule-Walker estimation in the first stage to ensure stability and
because the implementation of stage III would otherwise be even more
cumbersome.

## Usage

``` r
est_arma_hrk3(
  y,
  tmpl,
  maxit_stage2 = 5,
  tol_stage2 = 0.001,
  maxit_stage3 = 5,
  tol_stage3 = 0.001,
  info = TRUE,
  trace = FALSE,
  p.max = NULL,
  ic = c("AIC", "BIC", "max"),
  mean_estimate = c("zero", "sample.mean", "intercept"),
  tol = sqrt(.Machine$double.eps)
)
```

## Arguments

- y:

  sample, i.e. an \\(N,m)\\ dimensional matrix, or a "time series"
  object (i.e. `as.matrix(y)` should return an \\(N,m)\\-dimensional
  numeric matrix). Missing values (`NA`, `NaN` and `Inf`) are **not**
  supported.

- tmpl:

  a model template, see
  [`model structures()`](https://bfunovits.github.io/RLDM/reference/model_structures.md).
  Note that only the case is implemented, where \\a_0=b_0\\ holds, the
  diagonal entries of \\a_0=b_0\\ are equal to one and all other fixed
  elements are equal to zero. Furthermore the square root `sigma_L` of
  the noise covariance matrix is asssumed to be a lower triangular
  matrix without any further restrictions.  
  The given template is coerced to a template of this kind. If the given
  template does not comply to these restrictions, then a warning message
  is issued.

- maxit_stage2, maxit_stage3:

  Integers. Default for both stages is 5.

- tol_stage2, tol_stage3:

  Default set to `1e-3`. Maximal absolute distance of the deep
  parameters of adjacent iterations

- info:

  Boolean. Indicates whether the slot `extra_info` of the returned list
  should contain a tibble with additional info about the model at
  different iterations, e.g., min and max absolute value of zeros,
  poles, value of objective function etc.

- trace:

  (boolean) if `trace=TRUE`, then some tracing information on the
  iterations is printed.

- p.max:

  (integer or `NULL`) Maximum order of the candidate AR models. For the
  default choice see below.

- ic:

  (character string) Which information criterion shall be used to find
  the optimal AR order. Note that `ic="max"` means that an AR(p) model
  with `p=p.max` is fitted. Default is `ic="AIC"`.

- mean_estimate:

  Character string giving the method used to estimate the mean \\\mu\\.
  Default is `mean_estimate = "sample.mean"`. See the details below.

- tol:

  Small tolerance, used to check whether the mean of the supplied data
  matrix is indeed zero.

## Value

See
[`est_arma_hrk()`](https://bfunovits.github.io/RLDM/reference/est_arma_hrk.md).
The list contains the additional slots

- stage_opt:

  Since the best stable and miniphase model is returned, we also
  indicate whether this happened at stage II or stage III.

- info_tibble:

  Tibble containing all relevant info about the outcome at different
  stages and iterations.

and does not contain the slots `y.mean` (since mean is required to be
zero) and `converged`.

## Examples

``` r
data = BQdata_xts
for (pp in 0:2){
  for (qq in 0:2){
    if (pp + qq == 0){next}
    cat(paste0("p = ", pp, " q = ", qq, "\n"))
    tmpl = tmpl_arma_pq(m = 2, n = 2,
                        p = pp, q = qq)
    ff = est_arma_hrk3(data, tmpl = tmpl)
  }
}
#> p = 0 q = 1
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> p = 0 q = 2
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> p = 1 q = 0
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> p = 1 q = 1
#> p = 1 q = 2
#> p = 2 q = 0
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> p = 2 q = 1
#> p = 2 q = 2
ff$info_tibble %>% print(n=100)
#> # A tibble: 21 × 10
#>    stage iteration flipped min_pole max_pole min_zero max_zero trace log_lik
#>    <dbl>     <dbl> <lgl>      <dbl>    <dbl>    <dbl>    <dbl> <dbl>   <dbl>
#>  1     1         1 FALSE       1.23     3.13   Inf      Inf    0.945  Inf   
#>  2     2         1 FALSE       1.36     2.56     1.49     5.10 0.948   -1.33
#>  3     2         1 TRUE        1.36     2.56     1.49     5.10 0.948   -1.37
#>  4     2         2 FALSE       1.28     1.98     1.51     2.38 0.957   -1.34
#>  5     2         2 TRUE        1.28     1.98     1.51     2.38 0.957   -1.36
#>  6     2         3 FALSE       1.43    25.0      1.37     5.67 0.955   -1.34
#>  7     2         3 TRUE        1.43    25.0      1.37     5.67 0.955   -1.34
#>  8     2         4 FALSE       1.39     2.93     1.40     8.26 0.948   -1.33
#>  9     2         4 TRUE        1.39     2.93     1.40     8.26 0.948   -1.35
#> 10     2         5 FALSE       1.35     2.92     1.17     3.60 0.959   -1.34
#> 11     2         5 TRUE        1.35     2.92     1.17     3.60 0.959   -1.34
#> 12     3         1 FALSE       1.29     2.09     1.47     6.75 0.948   -1.34
#> 13     3         1 TRUE        1.29     2.09     1.47     6.75 0.963   -1.36
#> 14     3         2 FALSE       1.29     4.23     1.65     3.21 0.948   -1.34
#> 15     3         2 TRUE        1.29     4.23     1.65     3.21 0.958   -1.34
#> 16     3         3 FALSE       1.39     2.45     1.64     5.11 0.951   -1.34
#> 17     3         3 TRUE        1.39     2.45     1.64     5.11 0.958   -1.35
#> 18     3         4 FALSE       1.34     5.34     2.02     3.61 0.950   -1.34
#> 19     3         4 TRUE        1.34     5.34     2.02     3.61 0.952   -1.34
#> 20     3         5 FALSE       1.52     3.33     2.91     5.90 0.950   -1.34
#> 21     3         5 TRUE        1.52     3.33     2.91     5.90 0.952   -1.34
#> # ℹ 1 more variable: lndetSigma <dbl>
```

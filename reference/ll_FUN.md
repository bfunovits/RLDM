# Log Likelihood Function Factory

Creates a function similar to
[`ll_theta()`](https://bfunovits.github.io/RLDM/reference/ll_theta.md)
but faster and more memory efficient. The model structure (`template`)
and the data (`y`) are encoded within the generated closure (a function
plus its enclosing environment). The generated function calls compiled
C/C++ code (see
[RcppArmadillo-package](https://rdrr.io/pkg/RcppArmadillo/man/RcppArmadillo-package.html))
and hence is much faster than calling `ll_theta(th, template, y, ...)`.

## Usage

``` r
ll_FUN(
  template,
  y,
  which = c("concentrated", "conditional", "kf", "gr_concentrated"),
  skip = 0L,
  tol = 1e-08,
  err = NA_real_
)
```

## Arguments

- template:

  A model template, see [model
  structures](https://bfunovits.github.io/RLDM/reference/model_structures.md).

- y:

  sample, i.e. an \\(N,m)\\ dimensional matrix, or a "time series"
  object (i.e. `as.matrix(y)` should return an \\(N,m)\\-dimensional
  numeric matrix). Missing values (`NA`, `NaN` and `Inf`) are **not**
  supported.

- which:

  (string) Determines the type of ll function.

- skip:

  (integer) skip initial observations. If `NULL` then `skip` is set to
  \\0\\ for state space models and to \\\max(p,q)\\ for ARMA models.
  This parameter is only used for the cases "concentrated",
  "conditional" and "gr_concentrated"

- tol:

  (double) tolerance used by
  [`ll_kf()`](https://bfunovits.github.io/RLDM/reference/ll_kf.md).

- err:

  (double) return value for the case "kf", if the computation of the
  initial state covariance fails.

## Value

A function, `llfun(th)` say, which computes the log-likelihood for given
*deep* parameters `th`. This function may be used for ML estimation of
the model.

Function `fn(th)`

## Examples

``` r
# Generate a random model in echelon form model (m = 3)
tmpl = tmpl_stsp_echelon(nu = c(2,1,1))
model = r_model(template = tmpl, bpoles = 1, bzeroes = 1, sd = 0.25)
diag(model$sigma_L) = 1 # scale the diagonal entries of sigma_L
print(model)
#> state space model [3,3] with s = 4 states
#>             s[1]       s[2]        s[3]        s[4]        u[1]        u[2]
#> s[1]  0.00000000  0.0000000  0.00000000  1.00000000  0.07023323 -0.13607754
#> s[2]  0.09457491 -0.4687901 -0.03092246 -0.02045825 -0.25487373  0.05440225
#> s[3]  0.02047591  0.1742272 -0.06872531  0.42712106 -0.11331474 -0.14467849
#> s[4] -0.72399729 -0.1521027 -0.32503412 -0.05007725  0.05133291  0.42469459
#> x[1]  1.00000000  0.0000000  0.00000000  0.00000000  1.00000000  0.00000000
#> x[2]  0.00000000  1.0000000  0.00000000  0.00000000  0.00000000  1.00000000
#> x[3]  0.00000000  0.0000000  1.00000000  0.00000000  0.00000000  0.00000000
#>             u[3]
#> s[1]  0.10928421
#> s[2] -0.23300159
#> s[3]  0.46926641
#> s[4]  0.08568751
#> x[1]  0.00000000
#> x[2]  0.00000000
#> x[3]  1.00000000
#> Left square root of noise covariance Sigma:
#>            u[1]       u[2] u[3]
#> u[1]  1.0000000 0.00000000    0
#> u[2] -0.3129280 1.00000000    0
#> u[3] -0.2477061 0.05976526    1
# extract the corresponding free/deep parameters
th = extract_theta(model, tmpl)

# generate a sample with 50 observations
y = sim(model, n.obs = 50)$y

# conditional log likelihood
# the following statements return the same ll value!
ll(model, y, 'conditional')
#> [1] -4.256211
ll_theta(th, tmpl, y, 'conditional')
#> [1] -4.256211
fn = ll_FUN(tmpl, y, 'conditional')
fn(th)
#> [1] -4.256211

# concentrated conditional log likelihood
# the following statements return the same ll value!
ll(model, y, 'concentrated')
#> [1] -4.217044
ll_theta(th, tmpl, y, 'concentrated')
#> [1] -4.217044
fn = ll_FUN(tmpl, y, 'concentrated')
fn(th)
#> [1] -4.217044
# for this case, we may also compute the (analytic) gradient
gr = ll_FUN(tmpl, y, 'gr_concentrated')
gr(th)
#>  [1] -0.582604306 -0.230527610 -0.076719644 -0.163494913 -0.052758775
#>  [6] -0.048151745 -0.430755236  0.172244470 -0.182756693  0.180308000
#> [11]  0.151876817 -0.880783663 -0.005223902 -0.614676471  0.145094694
#> [16] -0.060682828  0.403265512 -0.053112623  0.137722860  0.271914889
#> [21]  0.115266821  0.683769868 -0.713761188  0.683946574  0.000000000
#> [26]  0.000000000  0.000000000  0.000000000  0.000000000  0.000000000

# log likelihood (via Kalman filter)
# the following statements return the same ll value!
ll(model, y, 'kf2')
#> [1] -4.322701
ll_theta(th, tmpl, y, 'kf2')
#> [1] -4.322701
ll(model, y, 'kf')
#> [1] -4.322701
ll_theta(th, tmpl, y, 'kf')
#> [1] -4.322701
fn = ll_FUN(tmpl, y, 'kf')
fn(th)
#> [1] -4.322701
```

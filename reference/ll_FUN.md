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
  [`ll_kf_cpp()`](https://bfunovits.github.io/RLDM/reference/ll_kf.md).

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
#>             s[1]         s[2]        s[3]        s[4]        u[1]       u[2]
#> s[1]  0.00000000  0.000000000  0.00000000  1.00000000 -0.47579669  0.1941772
#> s[2]  0.17834782 -0.001085635 -0.28582940 -0.01480361  0.09693782 -0.2800847
#> s[3]  0.09380898  0.013628061 -0.07131278  0.27512558  0.23146018 -0.1110135
#> s[4] -0.10879216 -0.278774721 -0.17177456  0.33897889  0.02478557  0.5165047
#> x[1]  1.00000000  0.000000000  0.00000000  0.00000000  1.00000000  0.0000000
#> x[2]  0.00000000  1.000000000  0.00000000  0.00000000  0.00000000  1.0000000
#> x[3]  0.00000000  0.000000000  1.00000000  0.00000000  0.00000000  0.0000000
#>             u[3]
#> s[1]  0.09142169
#> s[2] -0.10341404
#> s[3]  0.04197707
#> s[4]  0.10486197
#> x[1]  0.00000000
#> x[2]  0.00000000
#> x[3]  1.00000000
#> Left square root of noise covariance Sigma:
#>             u[1]     u[2] u[3]
#> u[1]  1.00000000 0.000000    0
#> u[2] -0.02632044 1.000000    0
#> u[3]  0.06769402 0.110399    1
# extract the corresponding free/deep parameters
th = extract_theta(model, tmpl)

# generate a sample with 50 observations
y = sim(model, n.obs = 50)$y

# conditional log likelihood
# the following statements return the same ll value!
ll(model, y, 'conditional')
#> [1] -4.283411
ll_theta(th, tmpl, y, 'conditional')
#> [1] -4.283411
fn = ll_FUN(tmpl, y, 'conditional')
fn(th)
#> [1] -4.283411

# concentrated conditional log likelihood
# the following statements return the same ll value!
ll(model, y, 'concentrated')
#> [1] -4.248823
ll_theta(th, tmpl, y, 'concentrated')
#> [1] -4.248823
fn = ll_FUN(tmpl, y, 'concentrated')
fn(th)
#> [1] -4.248823
# for this case, we may also compute the (analytic) gradient
gr = ll_FUN(tmpl, y, 'gr_concentrated')
gr(th)
#>  [1] -6.884443e-01 -3.628422e-02  7.316056e-01  2.632317e-01 -2.791062e-03
#>  [6] -6.040969e-02 -8.197604e-02  3.233422e-05  9.889647e-02 -6.908969e-01
#> [11]  1.347149e-02  3.425038e-01  6.087337e-02  3.500822e-01  1.623533e-01
#> [16] -3.689552e-02 -2.536515e-01 -4.107069e-01  1.207367e-01  7.516551e-02
#> [21] -3.309819e-02 -2.035432e-01 -2.322694e-01  9.862495e-02  0.000000e+00
#> [26]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00

# log likelihood (via Kalman filter)
# the following statements return the same ll value!
ll(model, y, 'kf2')
#> [1] -4.283794
ll_theta(th, tmpl, y, 'kf2')
#> [1] -4.283794
ll(model, y, 'kf')
#> [1] -4.283794
ll_theta(th, tmpl, y, 'kf')
#> [1] -4.283794
fn = ll_FUN(tmpl, y, 'kf')
fn(th)
#> [1] -4.283794
```

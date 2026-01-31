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
#>            s[1]        s[2]       s[3]        s[4]        u[1]        u[2]
#> s[1]  0.0000000  0.00000000  0.0000000  1.00000000 -0.17320930  0.09682150
#> s[2] -0.3306545  0.01350771 -0.1972378 -0.03578361  0.04348014 -0.41792942
#> s[3] -0.3057982 -0.05233582 -0.2495277  0.09318245 -0.08048530 -0.41472250
#> s[4] -0.1700720 -0.02146930  0.3660467 -0.24492897  0.06221162 -0.03184861
#> x[1]  1.0000000  0.00000000  0.0000000  0.00000000  1.00000000  0.00000000
#> x[2]  0.0000000  1.00000000  0.0000000  0.00000000  0.00000000  1.00000000
#> x[3]  0.0000000  0.00000000  1.0000000  0.00000000  0.00000000  0.00000000
#>            u[3]
#> s[1]  0.2189791
#> s[2] -0.1172880
#> s[3] -0.1764870
#> s[4]  0.2875431
#> x[1]  0.0000000
#> x[2]  0.0000000
#> x[3]  1.0000000
#> Left square root of noise covariance Sigma:
#>            u[1]       u[2] u[3]
#> u[1]  1.0000000  0.0000000    0
#> u[2] -0.1510780  1.0000000    0
#> u[3]  0.1920143 -0.1237398    1
# extract the corresponding free/deep parameters
th = extract_theta(model, tmpl)

# generate a sample with 50 observations
y = sim(model, n.obs = 50)$y

# conditional log likelihood
# the following statements return the same ll value!
ll(model, y, 'conditional')
#> [1] -4.115462
ll_theta(th, tmpl, y, 'conditional')
#> [1] -4.115462
fn = ll_FUN(tmpl, y, 'conditional')
fn(th)
#> [1] -4.115462

# concentrated conditional log likelihood
# the following statements return the same ll value!
ll(model, y, 'concentrated')
#> [1] -4.05168
ll_theta(th, tmpl, y, 'concentrated')
#> [1] -4.05168
fn = ll_FUN(tmpl, y, 'concentrated')
fn(th)
#> [1] -4.05168
# for this case, we may also compute the (analytic) gradient
gr = ll_FUN(tmpl, y, 'gr_concentrated')
gr(th)
#>  [1] -0.035510378  0.034818042 -0.084387036 -0.055252781 -0.168681470
#>  [6]  0.090886282 -0.107272572 -0.106298495  0.097460089  0.009128784
#> [11] -0.104673766  0.006148398  0.361105630  0.433500029 -0.145044069
#> [16] -0.071968044 -0.008647779  0.221878944 -0.295283813  0.051167458
#> [21] -0.020214612 -0.053694750 -0.001739417 -0.014929199  0.000000000
#> [26]  0.000000000  0.000000000  0.000000000  0.000000000  0.000000000

# log likelihood (via Kalman filter)
# the following statements return the same ll value!
ll(model, y, 'kf2')
#> [1] -4.124856
ll_theta(th, tmpl, y, 'kf2')
#> [1] -4.124856
ll(model, y, 'kf')
#> [1] -4.124856
ll_theta(th, tmpl, y, 'kf')
#> [1] -4.124856
fn = ll_FUN(tmpl, y, 'kf')
fn(th)
#> [1] -4.124856
```

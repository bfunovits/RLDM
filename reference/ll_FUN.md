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
#>             s[1]        s[2]        s[3]        s[4]        u[1]        u[2]
#> s[1]  0.00000000  0.00000000  0.00000000  1.00000000  0.15253617 -0.36242973
#> s[2] -0.02960837 -0.09734643  0.13501046  0.19710065 -0.33010318  0.36713901
#> s[3]  0.26414747  0.12324125  0.07361518 -0.29780813 -0.01250975  0.07081438
#> s[4]  0.33152503 -0.02196694 -0.02294610 -0.08258113  0.70768254  0.28681611
#> x[1]  1.00000000  0.00000000  0.00000000  0.00000000  1.00000000  0.00000000
#> x[2]  0.00000000  1.00000000  0.00000000  0.00000000  0.00000000  1.00000000
#> x[3]  0.00000000  0.00000000  1.00000000  0.00000000  0.00000000  0.00000000
#>            u[3]
#> s[1] -0.3437335
#> s[2] -0.2846045
#> s[3] -0.1803298
#> s[4]  0.4922631
#> x[1]  0.0000000
#> x[2]  0.0000000
#> x[3]  1.0000000
#> Left square root of noise covariance Sigma:
#>            u[1]       u[2] u[3]
#> u[1]  1.0000000  0.0000000    0
#> u[2] -0.4704373  1.0000000    0
#> u[3]  0.4350555 -0.5711638    1
# extract the corresponding free/deep parameters
th = extract_theta(model, tmpl)

# generate a sample with 50 observations
y = sim(model, n.obs = 50)$y

# conditional log likelihood
# the following statements return the same ll value!
ll(model, y, 'conditional')
#> [1] -4.23172
ll_theta(th, tmpl, y, 'conditional')
#> [1] -4.23172
fn = ll_FUN(tmpl, y, 'conditional')
fn(th)
#> [1] -4.23172

# concentrated conditional log likelihood
# the following statements return the same ll value!
ll(model, y, 'concentrated')
#> [1] -4.182474
ll_theta(th, tmpl, y, 'concentrated')
#> [1] -4.182474
fn = ll_FUN(tmpl, y, 'concentrated')
fn(th)
#> [1] -4.182474
# for this case, we may also compute the (analytic) gradient
gr = ll_FUN(tmpl, y, 'gr_concentrated')
gr(th)
#>  [1] -0.169894741 -0.096114707 -0.228725116 -0.332779561 -0.061868545
#>  [6] -0.074185344  0.001438396  0.064793316  0.077200476  0.394258857
#> [11]  0.161992344  0.093848540 -0.014857443 -0.050689809  0.490406828
#> [16]  0.032862825  0.062883330  0.136991842 -0.272113124 -0.040952224
#> [21]  0.318662848 -0.355185012  0.115213695 -0.057615072  0.000000000
#> [26]  0.000000000  0.000000000  0.000000000  0.000000000  0.000000000

# log likelihood (via Kalman filter)
# the following statements return the same ll value!
ll(model, y, 'kf2')
#> [1] -4.192763
ll_theta(th, tmpl, y, 'kf2')
#> [1] -4.192763
ll(model, y, 'kf')
#> [1] -4.192763
ll_theta(th, tmpl, y, 'kf')
#> [1] -4.192763
fn = ll_FUN(tmpl, y, 'kf')
fn(th)
#> [1] -4.192763
```

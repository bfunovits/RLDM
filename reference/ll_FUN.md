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
#>           s[1]       s[2]        s[3]       s[4]        u[1]        u[2]
#> s[1] 0.0000000  0.0000000  0.00000000  1.0000000  0.20857484 -0.41149993
#> s[2] 0.2293978 -0.1380675 -0.08606555 -0.2018322  0.04763087 -0.05484847
#> s[3] 0.3487320  0.3236556  0.42207815  0.2199682 -0.08944036  0.02261773
#> s[4] 0.1621189  0.1365278  0.05964505  0.6051999  0.17800942 -0.31638015
#> x[1] 1.0000000  0.0000000  0.00000000  0.0000000  1.00000000  0.00000000
#> x[2] 0.0000000  1.0000000  0.00000000  0.0000000  0.00000000  1.00000000
#> x[3] 0.0000000  0.0000000  1.00000000  0.0000000  0.00000000  0.00000000
#>             u[3]
#> s[1] -0.29330747
#> s[2]  0.07529302
#> s[3]  0.05620491
#> s[4]  0.17254724
#> x[1]  0.00000000
#> x[2]  0.00000000
#> x[3]  1.00000000
#> Left square root of noise covariance Sigma:
#>            u[1]        u[2] u[3]
#> u[1]  1.0000000  0.00000000    0
#> u[2] -0.2357979  1.00000000    0
#> u[3] -0.4381409 -0.04747029    1
# extract the corresponding free/deep parameters
th = extract_theta(model, tmpl)

# generate a sample with 50 observations
y = sim(model, n.obs = 50)$y

# conditional log likelihood
# the following statements return the same ll value!
ll(model, y, 'conditional')
#> [1] -4.499707
ll_theta(th, tmpl, y, 'conditional')
#> [1] -4.499707
fn = ll_FUN(tmpl, y, 'conditional')
fn(th)
#> [1] -4.499707

# concentrated conditional log likelihood
# the following statements return the same ll value!
ll(model, y, 'concentrated')
#> [1] -4.447215
ll_theta(th, tmpl, y, 'concentrated')
#> [1] -4.447215
fn = ll_FUN(tmpl, y, 'concentrated')
fn(th)
#> [1] -4.447215
# for this case, we may also compute the (analytic) gradient
gr = ll_FUN(tmpl, y, 'gr_concentrated')
gr(th)
#>  [1] -0.3277404991 -0.1407434796 -0.3685848644 -0.0215151691 -0.0030037120
#>  [6] -0.0591878729 -0.3351719111 -0.2173903412 -0.4199189542 -0.2918837948
#> [11] -0.1576236177 -0.4216633492  0.1924184086  0.2635383820  0.0021320010
#> [16] -0.0803958074  0.1807391011  0.1248992331  0.1225110598  0.0757156345
#> [21] -0.0009164495 -0.2747089633  0.0564316198 -0.0535579998  0.0000000000
#> [26]  0.0000000000  0.0000000000  0.0000000000  0.0000000000  0.0000000000

# log likelihood (via Kalman filter)
# the following statements return the same ll value!
ll(model, y, 'kf2')
#> [1] -4.477154
ll_theta(th, tmpl, y, 'kf2')
#> [1] -4.477154
ll(model, y, 'kf')
#> [1] -4.477154
ll_theta(th, tmpl, y, 'kf')
#> [1] -4.477154
fn = ll_FUN(tmpl, y, 'kf')
fn(th)
#> [1] -4.477154
```

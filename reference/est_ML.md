# Maximum Likelihood Estimation

**\[superseded\]** This is a naive implementation of Maximum Likelihood
Estimation. Rather use
[`ll()`](https://bfunovits.github.io/RLDM/reference/ll.md) and
[`ll_FUN()`](https://bfunovits.github.io/RLDM/reference/ll_FUN.md).

## Usage

``` r
est_ML(
  y,
  tmpl,
  th,
  which = c("concentrated", "conditional", "kf"),
  method = c("BFGS", "Nelder-Mead", "CG", "L-BFGS-B", "SANN", "Brent"),
  hessian = FALSE,
  control = list()
)
```

## Arguments

- y:

  sample, i.e. an \\(N,m)\\ dimensional matrix, or a "time series"
  object (i.e. `as.matrix(y)` should return an \\(N,m)\\-dimensional
  numeric matrix). Missing values (`NA`, `NaN` and `Inf`) are **not**
  supported.

- tmpl:

  a model template which describes the model class, see
  [`model structures()`](https://bfunovits.github.io/RLDM/reference/model_structures.md).
  Note that only the case of (non-empty, square) state space or ARMA
  models is implemented.

- th:

  Initial parameter estimate.

- which:

  (character string) determines which "likelihood" to be used, see also
  [`ll()`](https://bfunovits.github.io/RLDM/reference/ll.md). The option
  `"kf"` is only supported for state space models.

- method, hessian, control:

  are passed on to the optimization routine
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html).

## Value

A list with components

- model:

  The estimated model.

- th:

  The corresponding vector of *deep* parameters.

- ll:

  The log likelihood of the estimated model.

- which:

  The type of likelihood used.

- counts, convergence, message, hessian:

  as returned by [`stats::optim()`](https://rdrr.io/r/stats/optim.html).

## Note

- The optimization is computed with the general-purpose routine
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html).

- An initial estimate is **needed**.

- The procedure does **not** respect constraints like stability or
  minimum phase.

- The case of the *conditional, concentrated* likelihood is somewhat
  special. In this case the model template must have a particular
  structure: (1) The noise covariance is parametrized via the left
  cholesky factor. (2) The last \\m(m+1)/2\\ components of the parameter
  vector \\\theta\\ parametrize this left cholesky factor and the other
  components describe the system. (This implies that there is no
  overlap/dependency betweeen the "system parameters" and the "noise
  parameters".)

## Examples

``` r
# Generate a random model in echelon form model (m = 3)
tmpl = tmpl_stsp_echelon(nu = c(2,1,1))
model = r_model(template = tmpl, bpoles = 1, bzeroes = 1, sd = 0.25)
diag(model$sigma_L) = 1 # scale the diagonal entries of sigma_L
print(model)
#> state space model [3,3] with s = 4 states
#>              s[1]      s[2]       s[3]        s[4]       u[1]        u[2]
#> s[1]  0.000000000 0.0000000  0.0000000  1.00000000  0.2260673 -0.18269465
#> s[2] -0.300891527 0.1763168 -0.3097276 -0.10587253  0.0199123 -0.04753638
#> s[3] -0.009954261 0.2478604  0.6637246 -0.04959676 -0.3147068  0.13221617
#> s[4]  0.171745616 0.2860622 -0.0392293 -0.22370060  0.2564213  0.13755263
#> x[1]  1.000000000 0.0000000  0.0000000  0.00000000  1.0000000  0.00000000
#> x[2]  0.000000000 1.0000000  0.0000000  0.00000000  0.0000000  1.00000000
#> x[3]  0.000000000 0.0000000  1.0000000  0.00000000  0.0000000  0.00000000
#>             u[3]
#> s[1]  0.13742108
#> s[2] -0.16488559
#> s[3]  0.01435543
#> s[4] -0.70200263
#> x[1]  0.00000000
#> x[2]  0.00000000
#> x[3]  1.00000000
#> Left square root of noise covariance Sigma:
#>            u[1]       u[2] u[3]
#> u[1]  1.0000000  0.0000000    0
#> u[2] -0.1955948  1.0000000    0
#> u[3] -0.1660262 -0.1268121    1
# extract the corresponding free/deep parameters
th = extract_theta(model, tmpl)

# generate a sample with 500 observations
y = sim(model, n.obs = 500, n.burn_in = 100)$y

# We are cheating here and use the true model parameters
# as starting values for the optimization routine:

# estimate the model with the "exakt log likelihood"
out = est_ML(y, tmpl, th, which = 'kf')
KL_divergence(model, out$model)
#> [1] 0.01806206

# estimate the model with "conditional log likelihood"
out = est_ML(y, tmpl, th, which = 'conditional')
KL_divergence(model, out$model)
#> [1] 0.01786555

# estimate the model with "concentrated, conditional log likelihood"
out = est_ML(y, tmpl, th, which = 'concentrated')
KL_divergence(model, out$model)
#> [1] 0.01786693
```

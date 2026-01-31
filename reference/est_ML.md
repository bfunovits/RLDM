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
set.seed(123)
# Generate a random model in echelon form model (m = 3)
tmpl = tmpl_stsp_echelon(nu = c(2,1,1))
model = r_model(template = tmpl, bpoles = 1, bzeroes = 1, sd = 0.25)
diag(model$sigma_L) = 1 # scale the diagonal entries of sigma_L
print(model)
#> state space model [3,3] with s = 4 states
#>             s[1]       s[2]       s[3]        s[4]        u[1]       u[2]
#> s[1]  0.00000000 0.00000000  0.0000000  1.00000000  0.10019286  0.1244626
#> s[2] -0.14011891 0.01762710  0.1152291 -0.11141549  0.02767068 -0.4916543
#> s[3] -0.05754437 0.03232193 -0.3162653  0.30602045 -0.13896028  0.1753390
#> s[4]  0.38967708 0.42876625 -0.1717132  0.08995346  0.44672828 -0.1181979
#> x[1]  1.00000000 0.00000000  0.0000000  0.00000000  1.00000000  0.0000000
#> x[2]  0.00000000 1.00000000  0.0000000  0.00000000  0.00000000  1.0000000
#> x[3]  0.00000000 0.00000000  1.0000000  0.00000000  0.00000000  0.0000000
#>             u[3]
#> s[1] -0.26695593
#> s[2] -0.05449373
#> s[3] -0.25650111
#> s[4] -0.18222281
#> x[1]  0.00000000
#> x[2]  0.00000000
#> x[3]  1.00000000
#> Left square root of noise covariance Sigma:
#>            u[1]       u[2] u[3]
#> u[1]  1.0000000  0.0000000    0
#> u[2] -0.4216733  1.0000000    0
#> u[3]  0.2094468 -0.2845342    1
# extract the corresponding free/deep parameters
th = extract_theta(model, tmpl)

# generate a sample with 500 observations
y = sim(model, n.obs = 500, n.burn_in = 100)$y

# We are cheating here and use the true model parameters
# as starting values for the optimization routine:

# estimate the model with the "exakt log likelihood"
out = est_ML(y, tmpl, th, which = 'kf')
KL_divergence(model, out$model)
#> [1] 0.007795355

# estimate the model with "conditional log likelihood"
out = est_ML(y, tmpl, th, which = 'conditional')
KL_divergence(model, out$model)
#> [1] 0.007764377

# estimate the model with "concentrated, conditional log likelihood"
out = est_ML(y, tmpl, th, which = 'concentrated')
KL_divergence(model, out$model)
#> [1] 0.007765628
```

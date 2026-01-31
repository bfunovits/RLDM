# Gaussian log Likelihood of a State Space Model

These routines compute the log Likelihood for time invariant, linear
state space models of the form \$\$a\_{t+1} = A a_t + Bu_t\$\$ \$\$y_t =
C a_t + Du_t\$\$ with \\m\\-dimensional outputs \\y_t\\,
\\s\\-dimensional states \\a_t\\ and \\n\\-dimensional disturbances
\\u_t\\. The disturbances are white noise with a covariance matrix
\\\mathbf{E} u_t u_t'=\Sigma\\. Note that the disturbances and the
outputs may have *different* dimensions, however, only "wide" systems
with (\\m\leq n\\) are implemented.

The Gaussian log likelihood (for the case of Gaussian disturbances
\\u_t\sim N(0,\Sigma)\\ and \\a_1\sim N(a\_{1\|0},\Pi\_{1\|0})\\) here
is computed by the standard Kalman Filter or the square root Kalman
filter, see [`kf()`](https://bfunovits.github.io/RLDM/reference/kf.md).
The Kalman filter is a recursive scheme to compute the linear, least
squares predictions for \\a\_{t+1}\\ and \\y\_{t+1}\\ given the
observations \\y_t,\ldots,y_1\\ up to time \\t\\. These predictions are
notated with \\a\_{t+1\|t}\\ and \\y\_{t+1\|t}\\, the prediction error
for the output \\y\_{t+1}\\ is
\\\epsilon\_{t+1\|t}=(y\_{t+1}-y\_{t+1\|t})\\ and the corresponding
variances of the prediction errors are
\$\$\Pi\_{t+1\|t}=\mathbf{E}(a\_{t+1}-a\_{t+1\|t})
(a\_{t+1}-a\_{t+1\|t})',\$\$
\$\$\Sigma\_{t+1\|t}=\mathbf{E}(\epsilon\_{t+1\|t}
\epsilon\_{t+1\|t}').\$\$

The standard form of the Kalman filter is based on the parameter
matrices \\A,C\\, the variance of "state disturbances"
\\Q=\mathbf{E}(Bu_t (Bu_t)')=(B\Sigma B')\\, the variance of the
"measurement disturbances" \\R=\mathbf{E}(Du_t (Du_t)')=(D\Sigma D')\\
and the covariance \\S=\mathbf{E}(Bu_t(Du_t)')=(B\Sigma D')\\.
Furthermore we need the initial prediction \\a\_{1\|0}\\ and the
corresponding error variance \\\Pi\_{1\|0}\\.

For the square root form of the filter we need the "square roots"
\\\Pi\_{1\|0}^{1/2}\\ and \\\Sigma^{1/2}\\, i.e. matrices such that
\\\Pi\_{1\|0} = \Pi\_{1\|0}^{1/2} (\Pi\_{1\|0}^{1/2})'\\ and \\\Sigma =
\Sigma^{1/2}(\Sigma^{1/2})'\\. In addition, we define
\\H=(D',B')'\Sigma^{1/2}\\.

The (scaled) Gaussian log Likelihood of this model then may be expressed
as \$\$\frac{-1}{2N}\sum\_{t=1}^{N}\left(m\log(2\pi) +
\log\det\Sigma\_{t\|t-1} + (y_t - y\_{t\|t-1})' \Sigma\_{t\|t-1}^{-1}
(y_t - y\_{t\|t-1}) \right).\$\$

## Usage

``` r
ll_kf(model, y, method = c("kf", "kf2"), P1 = NULL, a1 = NULL, tol = 0)
```

## Arguments

- model:

  [`stspmod()`](https://bfunovits.github.io/RLDM/reference/stspmod.md)
  object, which represents the state space model.

- y:

  sample, i.e. an \\(N,m)\\ dimensional matrix, or a "time series"
  object (i.e. `as.matrix(y)` should return an \\(N,m)\\-dimensional
  numeric matrix). Missing values (`NA`, `NaN` and `Inf`) are **not**
  supported.

- method:

  Character string. If `method="kf"` then `ll_kf` calls the internal C++
  implementation ("standard form" of the Kalman filter) and for
  `method="kf2"` the "square root" form of the Kalman filter is used,
  i.e. the internal square root implementation is called. Up to
  numerical errors the outputs should not depend on the chosen method.

- P1:

  \\(s,s)\\ dimensional covariance matrix of the error of the initial
  state estimate, i.e. \\\Pi\_{1\|0}\\. If `NULL`, then the state
  covariance \\P = APA'+B\Sigma B'\\ is used. Note that this scheme
  assumes that the state space model is stable, i.e. that the state
  transition matrix \\A\\ is stable.

- a1:

  \\s\\ dimensional vector, which holds the initial estimate
  \\a\_{1\|0}\\ for the state at time \\t=1\\. If `a1=NULL`, then a zero
  vector is used.

- tol:

  (small) tolerance value (or zero). In order to speed up the
  computations, the algorithm(s) switch to a constant Kalman gain when
  there is no significant change in state error covariance. This
  behavior is controlled by the parameter `tol` and may be switched off
  by setting `tol=0`.

## Value

(double) The Gaussian log Likelihood of the model.

## Details

The core routines are internal C++ implementations (`ll_kf_cpp` and
`ll_kf2_cpp`) which are RcppArmadillo implementations of the standard
and the square root Kalman filter. The function `ll_kf` is a wrapper
function, which extracts the necessary parameters from an
[`stspmod()`](https://bfunovits.github.io/RLDM/reference/stspmod.md)
object, computes the initial covariance matrix `P1` and the initial
state estimate `a1` (if not provided) and then calls the internal C++
implementations.

Square root Kalman filter: For the square root \\\Pi\_{1\|0}^{1/2}\\ the
procedure first tries the Cholesky decomposition. If this fails (since
\\\Pi\_{1\|0}^{1/2}\\ is (close to) singular), then `ll_kf` tries to
compute a symmetric square root via the eigenvalue decomposition of
\\\Pi\_{1\|0}^{1/2}\\.

## Notes

The procedures only accept "wide" state space systems (\\m \leq n\\),
since for "tall" systems (\\m \> n\\) the variance of the prediction
errors (\\\Sigma\_{t+1\|t}\\) is singular for \\t\\ larger than some
threshold.

## See also

[`kf()`](https://bfunovits.github.io/RLDM/reference/kf.md) for
computation of the predictions \\a\_{t+1\|t}\\ and \\y\_{t+1\|t}\\.

## Examples

``` r
s = 4  # state dimension
m = 2  # number of outputs
n = m  # number of inputs (square case m=n)
n.obs = 100 # sample size

# generate a (stable) state space model (in innovation form)
tmpl = tmpl_stsp_full(m, n, s, sigma_L = "chol")
model = r_model(tmpl, bpoles = 1, sd = 0.5)
# generate a sample
data = sim(model, n.obs = n.obs, a1 = NA)

# compute Q, R, S and P1
sigma_L = model$sigma_L
sigma = tcrossprod(sigma_L)
R = model$sys$D %*% sigma %*% t(model$sys$D)
S = model$sys$B %*% sigma %*% t(model$sys$D)
Q = model$sys$B %*% sigma %*% t(model$sys$B)
P1 = lyapunov(model$sys$A, Q)

# compute H and square root of P1
H = rbind(model$sys$D, model$sys$B) %*% sigma_L
P1_R = chol(P1)

# compute logLikelihood (via Kalman Filter)
ll = ll_kf(model, data$y)

# compute logLikelihood (via square root Kalman Filter)
ll_test = ll_kf(model, data$y, method = 'kf2')
all.equal(ll, ll_test)
#> [1] TRUE

# Note: ll_kf_cpp and ll_kf2_cpp are internal C++ implementations
# called via .Call() by ll_kf()

# call the "full" kf routines
out = kf(model, data$y)
all.equal(ll, out$ll)
#> [1] TRUE
out = kf(model, data$y, method = 'kf2')
all.equal(ll, out$ll)
#> [1] TRUE
```

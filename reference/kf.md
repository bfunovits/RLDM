# Kalman Filter

These functions implement the "standard" Kalman filter and the "square
root" Kalman filter (also called "square root covariance filter") for
time invariant, linear state space systems without exogenous inputs, see
e.g. (Anderson and Moore 2005) .

## Usage

``` r
kf(model, y, method = c("kf", "kf2"), P1 = NULL, a1 = NULL)

kf_cpp(A, C, Q, R, S, y_t, P1, a1)

kf2_cpp(A, C, H_t, y_t, P1_R, a1)
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

  Character string. If `method="kf"` then `kf` calls `kf_cpp` ("standard
  form" of the Kalman filter) and for `method="kf2"` the "square root"
  form of the Kalman filter is used, i.e. `kf2_cpp` is called. Up to
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

- A:

  \\(s,s)\\ dimensional state transition matrix \\A\\.

- C:

  \\(m,s)\\ dimensional matrix \\C\\.

- Q, R, S:

  The variance, covariance matrices of the "state disturbances"
  (\\Bu_t\\) and the "measurement disturbances" (\\Du_t\\) as described
  above. These matrices must be of dimension \\(s,s)\\, \\(m,m)\\ and
  \\(s,m)\\ respectively.

- y_t:

  \\(m,N)\\ transposed data matrix `y_t = t(y)`.

- H_t:

  \\(n,s+m)\\ dimensional matrix. This parameter corresponds to the
  transpose \\H'\\ of \\H=(D',B')'\Sigma^{1/2}\\.

- P1_R:

  (right) square root of `P1`, i.e. `P1 = t(P1_R) \%*\% P1_R`.

## Value

List with components

- e:

  \\(N,m)\\ dimensional matrix with the standardized one-step ahead
  prediction errors. The \\t\\-th row of the matrix `e` corresponds to
  \$\$e_t = \Sigma\_{t\|t-1}^{-1/2}\epsilon\_{t\|t-1}.\$\$ If the model
  is correctly specified then these standardized residuals are white
  noise with a unit covariance matrix. So they may be used for validaton
  of the model.

- a:

  \\(N+1,s)\\ dimensional matrix with the estimated states. The \\t\\-th
  row of the matrix `a` corresponds to \\a\_{t\|t-1}\\. Given `y` and
  `a`, the one step ahead predictions \\y\_{t\|t-1}\\ may be computed
  with `yh = a \%*\% t(C)`.

- ll:

  (scaled) Gaussian log likelihood of the model
  \$\$-\frac{1}{2N}\sum\_{t=1}^{N}\left(m\log(2\pi) +
  \log\det\Sigma\_{t\|t-1} + (y_t - y\_{t\|t-1})' \Sigma\_{t\|t-1}^{-1}
  (y_t - y\_{t\|t-1}) \right).\$\$

- P1:

  \\(s,s)\\ dimensional covariance matrix of the error of the state
  prediction \\a\_{N+1\|N}\\, i.e. this matrix corresponds to
  \\\Pi\_{N+1\|N}\\.

## Details

The model considered is \$\$a\_{t+1} = A a_t + Bu_t\$\$ \$\$y_t = C
a_t + Du_t\$\$ with \\m\\-dimensional outputs \\y_t\\, \\s\\-dimensional
states \\a_t\\ and \\n\\-dimensional disturbances \\u_t\\. The
disturbances are white noise with a covariance matrix \\\mathbf{E}
u_tu_t'=\Sigma\\. Note that the disturbances and the outputs may have
*different* dimensions, however, only "wide" systems with (\\m\leq n\\)
are implemented.

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

The routines `kf_cpp`, `kf2_cpp` are RcppArmadillo implementations of
the standard form and of the square root form of the Kalman filter. The
wrapper function `kf` takes an
[`stspmod()`](https://bfunovits.github.io/RLDM/reference/stspmod.md)
object, which describes the state space model, and then calls the
approriate `RcppArmadillo` function.

Square root Kalman filter: For the square root \\\Pi\_{1\|0}^{1/2}\\ the
procedure first tries the Cholesky decomposition. If this fails (since
\\\Pi\_{1\|0}^{1/2}\\ is (close to) singular), then `ll_kf` tries to
compute a symmetric square root via the eigenvalue decomposition of
\\\Pi\_{1\|0}^{1/2}\\.

## Notes

The `RcppArmadillo` functions (`kf_cpp` and `kf2_cpp`) do not check the
input parameters, so these function must be used with some care.

The procedures only accept "wide" state space systems (\\m \leq n\\),
since for "tall" systems (\\m \> n\\) the variance of the prediction
errors (\\\Sigma\_{t+1\|t}\\) is singular for \\t\\ larger than some
threshold.

Up to now, there is no support for models with exogenous inputs.

## References

Anderson BDO, Moore JB (2005). *Optimal filtering*. Dover Publications
Inc., London. Originally published: Englewood Cliffs, Prentice-Hall
1979.

## See also

There exist a number of R packages which implement the Kalman filter,
e.g. KFAS. However, most of these packages do not allow for correlations
between the "state noise" and the "measurement noise".

If only the likelihood is needed, then one may use
[`ll_kf()`](https://bfunovits.github.io/RLDM/reference/ll_kf.md).

## Examples

``` r
s = 4  # state dimension
m = 2  # number of outputs
n = 3  # number of inputs,
n.obs = 100 # sample size

# generate a (stable) state space model
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

# call Kalman filter. Note y_t = t(y)!
out = kf_cpp(model$sys$A, model$sys$C, Q, R, S, t(data$y), P1, double(s))
# use the wrapper function
out_test = kf(model, data$y, method = 'kf')
all.equal(out, out_test)
#> [1] TRUE

# compute H and square root of P1
H = rbind(model$sys$D, model$sys$B) %*% sigma_L
P1_R = chol(P1)

# call square root Kalman filter. Note H_t = t(H) and y_t = t(y)!
out_test = kf2_cpp(model$sys$A, model$sys$C, t(H), t(data$y), P1_R, double(s))
all.equal(out, out_test)
#> [1] TRUE
# use the wrapper function
out_test = kf(model, data$y, method = 'kf2')
all.equal(out, out_test)
#> [1] TRUE

# The one step ahead predictions for y[t] may be computed by
yh = out$a %*% t(model$sys$C)
# and the (non scaled) prediction errors are
uh = data$y - out$a[1:n.obs,] %*% t(model$sys$C)
```

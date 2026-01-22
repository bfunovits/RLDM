# Log Likelihood Methods

Tools and methods for the computation of the (conditional or exact)
Gaussian log-likelihood of ARMA, RMFD, and state space models. For
functions which serve as input for optimizers like
[optim](https://rdrr.io/r/stats/optim.html), see
[ll_theta](https://bfunovits.github.io/RLDM/reference/ll_theta.md) and
[ll_FUN](https://bfunovits.github.io/RLDM/reference/ll_FUN.md) (where
the latter is a *function factory* which generates a closure that serves
as such an input).

## Usage

``` r
ll(obj, y, which, ...)

# S3 method for class 'armamod'
ll(obj, y, which = c("concentrated", "conditional"), skip = 0L, ...)

# S3 method for class 'stspmod'
ll(
  obj,
  y,
  which = c("concentrated", "conditional", "kf", "kf2"),
  skip = 0L,
  P1 = NULL,
  a1 = NULL,
  tol = 1e-08,
  ...
)
```

## Arguments

- obj:

  Object of type
  [`armamod()`](https://bfunovits.github.io/RLDM/reference/armamod.md),
  [`rmfdmod()`](https://bfunovits.github.io/RLDM/reference/rmfdmod.md),
  or
  [`stspmod()`](https://bfunovits.github.io/RLDM/reference/stspmod.md).

- y:

  Data sample given as \\(N,m)\\ dimensional matrix, or a "time series"
  object (in the sense that `as.matrix(y)` should return an
  \\(N,m)\\-dimensional numeric matrix). Missing values (`NA`, `NaN` and
  `Inf`) are **not** supported.

- which:

  (character) which likelihood to compute.

- ...:

  Not used.

- skip:

  (integer) skip initial observations. This parameter is only used, when
  the (concentrated) conditional likelihood is computed.

- P1:

  \\(s,s)\\ dimensional covariance matrix of the error of the initial
  state estimate. If `NULL` then the state covariance \\P=APA'+B\Sigma
  B'\\ is used. This parameter is only used, when the (exact) likelihood
  is computed via the Kalman Filter. See
  [`ll_kf()`](https://bfunovits.github.io/RLDM/reference/ll_kf.md) for
  more details.

- a1:

  \\s\\ dimensional vector, which holds the initial estimate for the
  state at time t=1. If `a1=NULL`, then a zero vector is used. This
  parameter is only used, when the (exact) likelihood is computed via
  the Kalman Filter. See
  [`ll_kf()`](https://bfunovits.github.io/RLDM/reference/ll_kf.md) for
  more details.

- tol:

  (small) tolerance value (or zero) used by the Kalman Filter routines,
  see [`ll_kf()`](https://bfunovits.github.io/RLDM/reference/ll_kf.md).

## Value

(double) the (scaled) log Likelihood of the model.

## Details

The procedure three choices ...

For an ARMA model \$\$a_0 y_t + a_1 y\_{t-1} + \cdots + a_p y\_{t-p} =
b_0 u_t + b_1 u\_{t-1} + \cdots + b_q u\_{t-q}\$\$ with Gaussian noise
\\u_t \sim N(0,\Sigma)\\ an *approximation* of the **scaled** log
likelihood is \$\$ll = -(1/2)(m \ln(2\pi) + \mathrm{tr}(S\Sigma^{-1}) +
\ln\det \Sigma + 2 \ln\det (a_0^{-1}b_0)\$\$ where \\S\\ denotes the
sample covariance of the residuals of the model
\$\$S=\frac{1}{N-s}\sum\_{t=s+1}^N e_t e_t'\$\$ The residuals are
computed from a sample \\y_t, t=1,\ldots,N\\ by solving the (inverse)
ARMA system \$\$b_0 e_t = -b_1 e\_{t-1} - \cdots - b_q e\_{t-q} + a_0
y_t + a_1 y\_{t-1} + \cdots + a_p y\_{t-p}\$\$ and setting all unknown
initial values \\y_t=0\\ and \\e_t=0\\ for \\t\leq 0\\ equal to zero.
See e.g.
[`solve_inverse_de()`](https://bfunovits.github.io/RLDM/reference/solve_de.md).
Note that the log Likelihood here is **scaled** by a factor \\1/(N-s)\\
and that the first \\s\\ observations are **skipped** when computing the
sample covariance matrix.

The log-likelihood may be easily maximized with respect to the noise
covariance matrix \\\Sigma\\. For given \\S\\, the optimal value for
\\\Sigma\\ is \\\Sigma=S\\. If we plug this maximizer into the
log-likelihood function, we obtain the "concentrated" log likelihood
function \$\$cll = -(1/2)(m \ln(2\pi) + m + \ln\det S + 2 \ln\det
(a_0^{-1}b_0)\$\$ which only depends on the sample `y` and the ARMA
parameter matrices \\a_i\\ and \\b_i\\.

For state space models the (approximate) log likelihood is computed
quite analogously.

## Note

To be precise, the functions returns \\1/(N-s)\\ times the (approximate)
log likelihood.

The above routines only handle the case of centered data, i.e. it is
assumed that the output process \\(y_t)\\ has mean zero!

The computation of the concentrated log likelihood assumes that the
model structure does **not** impose restrictions on the noise covariance
matrix.

## Examples

``` r
# Generate a random model in echelon form model (m = 3)
tmpl = tmpl_arma_echelon(nu = c(2,1,1))
model = r_model(template = tmpl, bpoles = 1, bzeroes = 1, sd = 0.25)
diag(model$sigma_L) = 1 # scale the diagonal entries of sigma_L
print(model)
#> ARMA model [3,3] with orders p = 2 and q = 2
#> AR polynomial a(z):
#>         z^0 [,1]  [,2]  [,3]    z^1 [,1]      [,2]       [,3]   z^2 [,1]
#> [1,]  1.00000000     0     0 -0.30971228 0.0000000  0.0000000 -0.2851629
#> [2,] -0.18030565     1     0 -0.07152422 0.2227471 -0.2489370  0.0000000
#> [3,]  0.02886177     0     1 -0.18531607 0.1409678  0.1947942  0.0000000
#>           [,2]       [,3]
#> [1,] 0.6455823 0.06594129
#> [2,] 0.0000000 0.00000000
#> [3,] 0.0000000 0.00000000
#> MA polynomial b(z):
#>         z^0 [,1]  [,2]  [,3]   z^1 [,1]        [,2]       [,3]    z^2 [,1]
#> [1,]  1.00000000     0     0 0.22986530 0.009034977  0.1316557 -0.08972547
#> [2,] -0.18030565     1     0 0.06354289 0.036216677  0.5051050  0.00000000
#> [3,]  0.02886177     0     1 0.20327238 0.068788081 -0.4723335  0.00000000
#>           [,2]       [,3]
#> [1,] 0.2853816 0.06148275
#> [2,] 0.0000000 0.00000000
#> [3,] 0.0000000 0.00000000
#> Left square root of noise covariance Sigma:
#>           u[1]      u[2] u[3]
#> u[1]  1.000000 0.0000000    0
#> u[2] -0.109147 1.0000000    0
#> u[3]  0.157453 0.1977473    1
# extract the corresponding free/deep parameters
th = extract_theta(model, tmpl)

# generate a sample with 50 observations
y = sim(model, n.obs = 50, n.burn_in = 100)$y

# conditional log likelihood
# the following statements return the same ll value!
ll(model, y, which = 'conditional', skip = 0)
#> [1] -4.405062
ll_theta(th, template= tmpl, y, which = 'conditional', skip = 0)
#> [1] -4.405062
llfun = ll_FUN(tmpl, y, which = 'conditional', skip = 0)
llfun(th)
#> [1] -4.405062

# concentrated, conditional log likelihood
# the following statements return the same ll value!
ll(model, y, which = 'concentrated', skip = 0)
#> [1] -4.365861
ll_theta(th, template= tmpl, y, which = 'concentrated', skip = 0)
#> [1] -4.365861
llfun = ll_FUN(tmpl, y, which = 'concentrated', skip = 0)
llfun(th)
#> [1] -4.365861
```

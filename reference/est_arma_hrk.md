# Hannan, Rissanen, Kavalieris estimation procedure

**\[deprecated\]** Estimate (V)ARMA models with the Hannan, Rissanen
Kavalieris procedure, see e.g. (Hannan and Rissanen 1982) and (Hannan et
al. 1986) .

## Usage

``` r
est_arma_hrk(
  y,
  e = NULL,
  tmpl,
  maxit = 1,
  tol = 0.001,
  trace = TRUE,
  p.max = NULL,
  ic = c("AIC", "BIC", "max"),
  mean_estimate = c("sample.mean", "intercept", "zero")
)
```

## Arguments

- y:

  sample, i.e. an \\(N,m)\\ dimensional matrix, or a "time series"
  object (i.e. `as.matrix(y)` should return an \\(N,m)\\-dimensional
  numeric matrix). Missing values (`NA`, `NaN` and `Inf`) are **not**
  supported.

- e:

  (initial) estimate of the disturbances \\u_t\\. If non `NULL` then `e`
  has to be an \\(N,m)\\ dimensional matrix, or a "time series" object
  (i.e an object which may be coerced to an \\(N,m)\\ dimensional matrix
  with `as.matrix(e)`).  
  If `NULL` then the procedure computes an estimate of the disturbances
  by fitting a "long" AR model to the data, see
  [`est_ar_ols()`](https://bfunovits.github.io/RLDM/reference/est_ar.md).  
  The matrix `e` may contain missing values (`NA`, `NaN` and `Inf`).
  Note that e.g. `est_ar_ols` returns residuals where the first \\p\\
  (here \\p\\ refers to the order of the fitted AR model) values are
  missing.

- tmpl:

  a model template, see
  [`model structures()`](https://bfunovits.github.io/RLDM/reference/model_structures.md).
  Note that only the case is implemented, where \\a_0=b_0\\ holds, the
  diagonal entries of \\a_0=b_0\\ are equal to one and all other fixed
  elements are equal to zero. Furthermore the square root `sigma_L` of
  the noise covariance matrix is asssumed to be a lower triangular
  matrix without any further restrictions.  
  The given template is coerced to a template of this kind. If the given
  template does not comply to these restrictions, then a warning message
  is issued.

- maxit:

  (integer) maximum number of iterations

- tol:

  (numeric) tolerance level

- trace:

  (boolean) if `trace=TRUE`, then some tracing information on the
  iterations is printed.

- p.max:

  (integer or `NULL`) Maximum order of the candidate AR models. For the
  default choice see below.

- ic:

  (character string) Which information criterion shall be used to find
  the optimal AR order. Note that `ic="max"` means that an AR(p) model
  with `p=p.max` is fitted. Default is `ic="AIC"`.

- mean_estimate:

  Character string giving the method used to estimate the mean \\\mu\\.
  Default is `mean_estimate = "sample.mean"`. See the details below.

## Value

List with components

- model:

  the estimated (V)ARMA model (i.e. an
  [`armamod()`](https://bfunovits.github.io/RLDM/reference/armamod.md)
  object).

- th:

  vector with the (free) parameters of the estimated (V)ARMA model.

- tmpl:

  the (coerced) model template.

- y.mean:

  estimate of the mean \\\mu\\.

- residuals:

  the residuals of the model, computed with
  [`solve_inverse_de()`](https://bfunovits.github.io/RLDM/reference/solve_de.md).

- sigma:

  the sample variance \\S\\ of the residuals, i.e. an estimate of the
  noise covariance matrix \\\Sigma\\.

- n.valid:

  number of "valid" observations, i.e. observations where all needed
  lagged values \\y\_{t-i}\\ and \\e\_{t-i}\\ are availiable. For an
  ARMA(p,q) model this implies that the number of valid observations is
  less than or equal to `n.obs -max(p,q)`.

- ll:

  Gaussian log likelihood: \$\$(-1/2)(m \ln(2\pi) + m + \ln\det(S))\$\$
  where \\S\\ denotes the sample variance of the residuals.

- iter:

  number of iterations.

- converged:

  (boolean) indicates whether the algorithm converged.

## Details

The main idea of the HRK procedure is as follows. If we have given
estimates, \\e_t\\ say, of the disturbances, then the ARMA model is
estimated from the equation \$\$y_t = -a^\*\_0(y_t+e_t) - a_1 y\_{t-1} -
\cdots - a_p y\_{t-p} + b_1 e\_{t-1} + \cdots + b_q e\_{t-q} +
v\_{t-1}\$\$ where \\a^\*\_0\\ is obtained from \\a_0 = b_0\\ by setting
all diagonal elements equal to zero. The entries in the parameter
matrices \\a_i\\ and \\b_i\\ are either treated as fixed (and equal to
zero) or "free". Now the above regression is estimated "componentwise",
i.e. for each component of \\y_t\\ the corresponding "free" parameters
are estimated by OLS.

Given the parameter estimates one computes new estimates for the
disturbances, by recursively solving the ARMA system, see
[`solve_inverse_de()`](https://bfunovits.github.io/RLDM/reference/solve_de.md).
The sample variance of these residuals is used as an estimate of the
noise covariance matrix \\\Sigma\\.

This procedure may be iterated: use the "new" estimates for the
disturbances to (re) estimate the ARMA parameters and to (re) estimate
the disturbances, ...

The parameters `maxit` and `tol` control this iterative scheme. The
iterations are stopped after at most `maxit` iterations or when there is
only a "small" change of the estimates. To be more precise, if `th`,
`th0` denote the vector of parameter estimates in the actual round and
the previous round, then the procedure stops if
`max(abs(th-th0)) <= tol`.

Note that in general there is no guarantee that this iterative scheme
converges or that the estimates are improved by iterating.

The user may supply his "own" (initial) estimates `e` of the
disturbances. If the parameter `e` is missing (or `NULL`) then the
procedure `est_arma_hrk` computes estimates of the disturbances by
fitting a "long" AR model to the data. To this end the procedure simply
calls
[`est_ar_ols()`](https://bfunovits.github.io/RLDM/reference/est_ar.md)
with the respective paramaters `p.max` (which controls the maximum
possible AR order), `ic` (which controls the information criterion used
to select the order of the AR model) and `mean_estimate` (which tells
`est_ar_ols` how to estimate the mean \\\mu\\). The default for the
maximum order `p.max` is \$\$\max(12, 10\log\_{10}(N), (N-1)/(m+1))\$\$

The procedure supports three options for the estimation of the mean
\\\mu = \mathbf{E} y_t\\. For `mean_estimate="zero"` the procedure sets
the (estimate of the) mean equal to zero. For
`mean_estimate="sample.mean"` the procedure simply uses the sample mean
of `y` as an estimate. Third option `mean_estimate="intercept"` uses an
intercept in the above regression(s) and computes the estimate of the
mean correspondingly. Note that this fails if the estimated AR
polynomial has a unit root, i.e. if \$\$\det \hat{a}(1) = 0.\$\$

There is no guarantee that the HRK algorithm returns a stable and
minimum phase ARMA model. In particular, if the estimated model is *not*
minimum phase then the recursive computation of the residuals often
yields useless results and correspondingly the cholesky decomposition of
the sample variance of the residuals (which is used as estimate of the
noise covariance matrix \\\Sigma\\) fails. In this case the procedure
stops with an error message.

## References

Hannan EJ, Rissanen J (1982). “Recursive estimation of mixed
autoregressive-moving average order.” *Biometrika*, **69**, 81–94.

Hannan EJ, Kavalieris L, Mackisack M (1986). “Recursive Estimation of
Linear Systems.” *Biometrika*, **73**(1), 119-133.

## Examples

``` r
# in order to get reproducible results
set.seed(4321)

# generate a random VARMA(p=2,q=1) model with m=2 outputs #####################
tmpl = tmpl_arma_pq(m = 2, n = 2, p = 2, q = 1)
model = r_model(template = tmpl, bpoles = 1, bzeroes = 1, sd = 0.25)
diag(model$sigma_L) = 1 # scale the diagonal entries of sigma_L
print(model)
#> ARMA model [2,2] with orders p = 2 and q = 1
#> AR polynomial a(z):
#>      z^0 [,1]  [,2]    z^1 [,1]      [,2]    z^2 [,1]        [,2]
#> [1,]        1     0 -0.10668935 0.1794017 -0.03208932 -0.07429186
#> [2,]        0     1 -0.05590295 0.2103614  0.40233680  0.04900116
#> MA polynomial b(z):
#>      z^0 [,1]  [,2]   z^1 [,1]        [,2]
#> [1,]        1     0  0.3101866 -0.01680908
#> [2,]        0     1 -0.1796745  0.08609177
#> Left square root of noise covariance Sigma:
#>          u[1] u[2]
#> u[1] 1.000000    0
#> u[2] 0.284866    1

# generate a sample with 200 observations
data = sim(model, n.obs = 200, n.burn_in = 100)

# estimate model with HRK
# note: we are cheating here and use the true disturbances!
out = est_arma_hrk(data$y, data$u, tmpl)
#> HRK estimation of ARMA model: m=2, n.obs=200, p=2, q=1
#> iter |th - th0|  n.val      MSE       ll 
#>    1      0.977    198    1.958   -2.754  
print(out$model)
#> ARMA model [2,2] with orders p = 2 and q = 1
#> AR polynomial a(z):
#>      z^0 [,1]  [,2]  z^1 [,1]        [,2]    z^2 [,1]        [,2]
#> [1,]        1     0 0.4774058 0.231267694 -0.08031496 -0.05099365
#> [2,]        0     1 0.1309504 0.005178953  0.33292233  0.06794187
#> MA polynomial b(z):
#>      z^0 [,1]  [,2]   z^1 [,1]        [,2]
#> [1,]        1     0 0.86715882 -0.05585317
#> [2,]        0     1 0.06883826 -0.24833582
#> Left square root of noise covariance Sigma:
#>           u[1]      u[2]
#> u[1] 0.9772386 0.0000000
#> u[2] 0.3424949 0.9410523
# ll() returns the same logLik value. However, we have to demean the data
all.equal(out$ll, ll(out$model, scale(data$y, center = out$y.mean, scale = FALSE),
                     'conditional', skip = 2))
#> [1] TRUE

# estimate the model with HRK
# use the residuals of a long AR model as estimates for the noise
out = est_arma_hrk(data$y, e = NULL, tmpl,
                   trace = TRUE, maxit = 10, mean_estimate = 'zero')
#> HRK estimation of ARMA model: m=2, n.obs=200, p=2, q=1
#> initial AR estimate of noise p.max=11, p=2, ll=-2.711726
#> iter |th - th0|  n.val      MSE       ll 
#>    1      0.939    197    1.865   -2.697  
#>    2      0.140    198    1.856   -2.692  
#>    3      0.060    198    1.856   -2.692  
#>    4      0.018    198    1.856   -2.692  
#>    5      0.007    198    1.856   -2.692  
#>    6      0.002    198    1.856   -2.692  
#>    7      0.002    198    1.856   -2.692  
#>    8      0.000    198    1.856   -2.692  
#> algorithm converged
print(out$model)
#> ARMA model [2,2] with orders p = 2 and q = 1
#> AR polynomial a(z):
#>      z^0 [,1]  [,2]   z^1 [,1]       [,2]   z^2 [,1]        [,2]
#> [1,]        1     0 0.20158389 0.13764630 0.05391753 -0.08392618
#> [2,]        0     1 0.06234842 0.09370431 0.37190127  0.05010730
#> MA polynomial b(z):
#>      z^0 [,1]  [,2]     z^1 [,1]       [,2]
#> [1,]        1     0  0.615023757 -0.1511000
#> [2,]        0     1 -0.001179974 -0.1500922
#> Left square root of noise covariance Sigma:
#>           u[1]      u[2]
#> u[1] 0.9230599 0.0000000
#> u[2] 0.3558209 0.9364313
# ll() returns the same logLik value. However, we have to demean the data
all.equal(out$ll, ll(out$model, scale(data$y, center = out$y.mean, scale = FALSE),
                     'conditional', skip = 2))
#> [1] TRUE

# Generate a random Model in echelon form model (m = 3) #######################
tmpl = tmpl_arma_echelon(nu = c(1,1,1))
model = r_model(template = tmpl, bpoles = 1, bzeroes = 1, sd = 0.25)
diag(model$sigma_L) = 1 # scale the diagonal entries of sigma_L
print(model)
#> ARMA model [3,3] with orders p = 1 and q = 1
#> AR polynomial a(z):
#>      z^0 [,1]  [,2]  [,3]   z^1 [,1]        [,2]       [,3]
#> [1,]        1     0     0  0.5893418 -0.00553193 -0.2113619
#> [2,]        0     1     0  0.1842932  0.06201503  0.1328532
#> [3,]        0     0     1 -0.1422159  0.23607602  0.3943483
#> MA polynomial b(z):
#>      z^0 [,1]  [,2]  [,3]    z^1 [,1]       [,2]      [,3]
#> [1,]        1     0     0 -0.07942568 -0.1992309 0.5002033
#> [2,]        0     1     0 -0.29348965  0.1566674 0.1483497
#> [3,]        0     0     1 -0.14925189  0.4185160 0.3153198
#> Left square root of noise covariance Sigma:
#>            u[1]      u[2] u[3]
#> u[1]  1.0000000 0.0000000    0
#> u[2]  0.2661327 1.0000000    0
#> u[3] -0.4719906 0.1258785    1

# generate a sample with 200 observations
data = sim(model, n.obs = 200, n.burn_in = 100)
# add mean value(s)
data$y = data$y + matrix(1:3, nrow = 200, ncol = 3, byrow = TRUE)

# estimate model with HRK
# note: we are cheating here and use the true disturbances!
out = est_arma_hrk(data$y, data$u, tmpl,
                   trace = FALSE, maxit = 1, mean_estimate = 'sample.mean')
print(out$y.mean)
#> [1] 0.9511495 1.9436019 2.8487408
print(out$model)
#> ARMA model [3,3] with orders p = 1 and q = 1
#> AR polynomial a(z):
#>      z^0 [,1]  [,2]  [,3]   z^1 [,1]        [,2]         [,3]
#> [1,]        1     0     0  0.6687489 -0.12600713 -0.154716557
#> [2,]        0     1     0  0.1944165 -0.04965131  0.057462781
#> [3,]        0     0     1 -0.3418481  0.50977514  0.007750021
#> MA polynomial b(z):
#>      z^0 [,1]  [,2]  [,3]    z^1 [,1]        [,2]        [,3]
#> [1,]        1     0     0  0.05047286 -0.31806802  0.52077685
#> [2,]        0     1     0 -0.21348959 -0.01193928  0.04368415
#> [3,]        0     0     1 -0.45130821  0.70739668 -0.11303414
#> Left square root of noise covariance Sigma:
#>            u[1]       u[2]    u[3]
#> u[1]  0.9644771 0.00000000 0.00000
#> u[2]  0.1932537 0.98300427 0.00000
#> u[3] -0.4325050 0.06192851 1.05045
# ll() returns the same logLik value. However, we have to demean the data
all.equal(out$ll, ll(out$model, scale(data$y, center = out$y.mean, scale = FALSE),
                     'conditional', skip = 1))
#> [1] TRUE

# estimate the model with HRK
# use the residuals of a long AR model as estimates for the noise
out = est_arma_hrk(data$y, e = NULL, tmpl,
                   maxit = 10, mean_estimate = 'intercept')
#> HRK estimation of ARMA model: m=3, n.obs=200, p=1, q=1
#> initial AR estimate of noise p.max=7, p=2, ll=-4.235502
#> iter |th - th0|  n.val      MSE       ll 
#>    1      1.051    197    3.200   -4.235  
#>    2      0.135    199    3.220   -4.249  
#>    3      0.081    199    3.223   -4.251  
#>    4      0.048    199    3.220   -4.250  
#>    5      0.035    199    3.222   -4.251  
#>    6      0.022    199    3.221   -4.250  
#>    7      0.015    199    3.222   -4.250  
#>    8      0.010    199    3.221   -4.250  
#>    9      0.007    199    3.222   -4.250  
#>   10      0.005    199    3.221   -4.250  
print(out$y.mean)
#> [1] 0.9537144 1.9386421 2.8536800
print(out$model)
#> ARMA model [3,3] with orders p = 1 and q = 1
#> AR polynomial a(z):
#>      z^0 [,1]  [,2]  [,3]   z^1 [,1]         [,2]       [,3]
#> [1,]        1     0     0  0.6740053 -0.261126479 -0.3217273
#> [2,]        0     1     0  0.1479565  0.008074719 -0.1286455
#> [3,]        0     0     1 -0.2984776  0.523115301  0.1643308
#> MA polynomial b(z):
#>      z^0 [,1]  [,2]  [,3]    z^1 [,1]        [,2]        [,3]
#> [1,]        1     0     0  0.05225048 -0.45689436  0.35600140
#> [2,]        0     1     0 -0.25603151  0.04533693 -0.14860902
#> [3,]        0     0     1 -0.41217354  0.72060291  0.04411931
#> Left square root of noise covariance Sigma:
#>            u[1]       u[2]     u[3]
#> u[1]  0.9625622 0.00000000 0.000000
#> u[2]  0.1902924 0.98177314 0.000000
#> u[3] -0.4307534 0.06340989 1.051195
# ll() returns the same logLik value. However, we have to demean the data
all.equal(out$ll, ll(out$model, scale(data$y, center = out$y.mean, scale = FALSE),
                     'conditional', skip = 1))
#> [1] TRUE

# We may also use this procedure to estimate AR models #####################
# where some coefficients are fixed = 0
a = dbind(d = 3, diag(2), array(NA_real_, dim = c(2,2,2)))
a[1,2,] = 0 # all coefficient matrices are lower triangular, i.e.
# y[2t] does not Granger cause y[1t]
tmpl = model2template(armamod(sys = lmfd(a=a),
                              sigma_L = matrix(NA_real_, nrow = 2, ncol = 2)),
                      sigma_L = 'chol')
model = r_model(template = tmpl, bpoles = 1, bzeroes = 1, sd = 0.25)
diag(model$sigma_L) = 1 # scale the diagonal entries of sigma_L
print(model)
#> ARMA model [2,2] with orders p = 2 and q = 0
#> AR polynomial a(z):
#>      z^0 [,1]  [,2]   z^1 [,1]      [,2]    z^2 [,1]       [,2]
#> [1,]        1     0 -0.1836216 0.0000000 -0.04334075 0.00000000
#> [2,]        0     1  0.2345054 0.2220447 -0.54760618 0.04449801
#> MA polynomial b(z):
#>      z^0 [,1]  [,2]
#> [1,]        1     0
#> [2,]        0     1
#> Left square root of noise covariance Sigma:
#>             u[1] u[2]
#> u[1]  1.00000000    0
#> u[2] -0.07017046    1

# generate a sample with 200 observations
data = sim(model, n.obs = 200, n.burn_in = 100)

# estimate model with HRK
out = est_arma_hrk(data$y, NULL, tmpl,
                   trace = FALSE, maxit = 1, mean_estimate = 'zero')
print(out$y.mean)
#> [1] 0 0
print(out$model)
#> ARMA model [2,2] with orders p = 2 and q = 0
#> AR polynomial a(z):
#>      z^0 [,1]  [,2]   z^1 [,1]    [,2]    z^2 [,1]       [,2]
#> [1,]        1     0 -0.1045799 0.00000 -0.04138278 0.00000000
#> [2,]        0     1  0.3202084 0.22883 -0.58485638 0.09887707
#> MA polynomial b(z):
#>      z^0 [,1]  [,2]
#> [1,]        1     0
#> [2,]        0     1
#> Left square root of noise covariance Sigma:
#>             u[1]      u[2]
#> u[1]  0.91829065 0.0000000
#> u[2] -0.04809834 0.9383795
# ll() returns the same logLik value. However, we have to demean the data
all.equal(out$ll, ll(out$model, scale(data$y, center = out$y.mean, scale = FALSE),
                     'conditional', skip = 2))
#> [1] TRUE

# reset the "seed"
set.seed(NULL)
```

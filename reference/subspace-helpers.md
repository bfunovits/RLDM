# Subspace Helper Methods

These procedure implement two subspace algorithms for the estimation of
state space models, the *AOKI* method, as described in (Aoki 1990) and
the *CCA* algorithm (see e.g. (Dahlen and Scherrer 2004) or (Bauer 2001)
). These subspace algorithms center on the weighted Hankel matrix
\$\$(R_f')^{-T} H\_{fp} R_p^{-1}\$\$ where the block *Hankel* matrix
\\H\_{fp}\\ is the covariance between the "past"
\\(y\_{t-1}',\cdots,y\_{t-p}')'\\ and the "future"
\\(y\_{t}',\cdots,y\_{t+f-1}')'\\ and \\R_f\\ and \\R_p\\ are the
cholesky factors of the covariance matrices of the "future" and the
"past" respectively. The singular values of this weighted Hankel matrix
are the canonical correlation coefficients between the past and the
future. Note that the implementation here always sets \\f = p+1\\.

## Usage

``` r
est_stsp_aoki(
  gamma,
  s.max,
  p,
  estorder = estorder_SVC,
  keep_models = FALSE,
  n.obs = NULL,
  ...
)

est_stsp_cca(
  gamma,
  s.max,
  p,
  estorder = estorder_SVC,
  keep_models = FALSE,
  n.obs = NULL,
  ...
)

est_stsp_cca_sample(
  y,
  s.max,
  p,
  estorder = estorder_SVC,
  keep_models = FALSE,
  mean_estimate = c("sample.mean", "zero"),
  ...
)
```

## Arguments

- gamma:

  \\(m,m,L+1)\\-dimensional array with the (sample) autocovariance
  function.

- s.max:

  (integer) maximum possible order.

- p:

  (integer) number of block columns of the Hankel matrix (size of the
  "past")

- estorder:

  function, to estimate the order of the system.

- keep_models:

  (boolean) should the function return a list with estimated system of
  order `0:s.max`?

- n.obs:

  sample size \\N\\.

- ...:

  additional parameters, passed on to the order estimation routine.

- y:

  \\(N,m)\\-dimensional matrix or an object, which may be coerced to a
  matrix with `as.matrix{y}`.

- mean_estimate:

  Character string giving the method used to estimate the mean \\\mu = E
  y_t\\. Default is to use the sample mean.

## Value

List with slots:

- model:

  a `stsp()` object, which represents the estimated state space model.

- models:

  either `NULL` (if `!keep_models`) or a list with the parameters of the
  estimated models with orders (`s=0:s.max+1`). This slot may e.g. be
  used to estimate the model order by some user defined model selection
  procedure.

- s:

  (integer) the estimate of the model order.

- info:

  list with information about the data and the design parameters of the
  estimation procedure.

- stats:

  ((s.max+1)-by-5)-dimensional matrix with statistics of the (estimated)
  state space models.

## Details

AOKIs method is a *realization* algorithm, i.e. it reconstructs the
underlying state space model from the (population) autocovariance
function. To this end a Riccati equation has to be solved, see
[`riccati()`](https://bfunovits.github.io/RLDM/reference/riccati.md). If
an estimated ACF is fed into this algorithm one obtains an estimate for
the state space model. However note that this may fail (in particular
the Riccati equation may have no positive definite solution) if the
estimate of the ACF is *not* positive definite, if the Hankel matrix is
too small or if state dimension is not correct.

The CCA method estimates the state space model by first constructing an
estimate of the states. Then the parameter matrices are estimated via
simple LS regressions. This procedure does not give the "true" model,
even in the case when the population ACF is used. However, the
"distance" between the true model and the estimated model converges to
zero, if the estimate of the ACF converges to the population ACF and the
size \\p\\ of the Hankel matrix converges to infinity.

There are two implementations of the CCA method:

1.  The routine `est_stsp_cca_sample()` operates directly on the
    supplied data.

2.  The routine `est_stsp_cca()` uses an (estimated) autocovariance
    function.

These algorithms may also be used as simple "model reduction
algorithms". If we want to approximate a high dimensional state space
model by a model of lower order, we may proceed as follows. First we
compute the ACF of the high dimensional model and then fed this ACF into
the subspace routines `est_stsp_cca` or `est_stsp_aoki`, however setting
the maximum order `s.max` to some value less than the true order.

**Order Estimation**

The order estimation is based on the Hankel singular values \\\sigma_s\\
and/or the log det values of the estimated noise covariance matrices
\\\ln\det \hat{\Sigma}\_s\\. Using only the Hankel singular values has
the advantage that only *one* model has to be estimated, whereas
otherwise estimates for *all* models with orders
\\s=0,\ldots,s\_{\max}\\ have to be computed.

In order to exploit this (small) advantage of singular values based
criteria the order estimation runs as follows: First the procedures
call  
` estorder(s.max, Hsv, n.par, m, n.obs, Hsize=c(f,p), ...)`  
Here `Hsv` is an \\pm\\ dimensional vector with the Hankel singular
values and `n.par` is an \\(s\_{\max}+1)\\ dimensional vector with the
respective number of parameters of the models with orders
\\s=0,\ldots,s\_{\max}\\. If this call returns an estimate of the order
then the procedures estimate a corresonding state space model.

If this call fails (i.e returns `NULL`) then the procedures estimate all
models with orders up to \\s\_{\max}\\ and the corresponding noise
covariance matrices. The order then is estimated by calling  
` estorder(s.max = s.max, Hsv, lndetSigma, n.par, m, n.obs, Hsize, ...)`  
where `lndetSigma` is the vector with the log det values of the
estimated noise covariance matrices (\\\ln\det \hat{\Sigma}\_s\\).

The package offers some predefined order selection procedures (see also
[subspace order
estimates](https://bfunovits.github.io/RLDM/reference/subspace-order-estimates.md)):

`estorder_max(s.max, ...)` simply returns the maximum order `s.max`
considered.

`estorder_rkH(s.max, Hsv, tol, ...)` estimates the order by an estimate
of the rank of the Hankel matrix.

`estorder_MOE(s.max, Hsv, ...)` estimates the order by searching for a
"gap" in the singular values.

`estorder_SVC(s.max, Hsv, n.par, n.obs, Hsize, penalty, ...)` implements
the so called *Singular Value Criteria*, see (Bauer 2001) : \$\$svc(s) =
\sigma\_{s+1}^2 + c(N)d(s)/N\$\$ Here \\\sigma_s\\ is the \\s\\-th
singular value of the weighted Hankel marix, \\N\\ is the sample size,
\\d(s) = 2ms\\ denotes the number of parameters for a state space model
with \\s\\ states (and \\m\\ outputs) and \\c(N)\\) is a "penalty"
(depending on the sample size).

The above order estimation procedures only use the Hankel singular
values, whereas the following procedure is based on the estimated noise
covariances.

`estorder_IVC(s.max, lndetSigma, n.par, n.obs, penalty, ...)` estimates
the order via an information criterion of the form \$\$ivc(s) =
\ln\det\hat\Sigma\_{s} + c(N)d(s)/N\$\$ where \\\hat\Sigma_s\\ is the
estimate of the noise covariace matrix obtained from a model with order
\\s\\, \\d(s)\\ denotes the number of parameters and \\c(N)\\ is a
"penalty" (depending on the sample size).

For both `estorder_SVC` and `estorder_IVC` the (optional) parameter
`penalty` controls the penalty term \\c(N)\\.

Note also that for `keep_models==TRUE` the estimation procedures compute
*all* models even in the case of a Hankel singular value based selection
criterion.

## References

Bauer D (2001). “Order estimation for subspace methods.” *Automatica*,
**37**(10), 1561 - 1573.
[doi:10.1016/S0005-1098(01)00118-2](https://doi.org/10.1016/S0005-1098%2801%2900118-2)
.

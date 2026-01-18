# Helper Functions for Order Estimation

These helper function are used in the subspace estimation routine
[subspace
methods](https://bfunovits.github.io/RLDM/reference/subspace-methods.md)
for the estimation of the order of a state space model. For a discussion
of order estimation in the context of subspace methods see e.g (Bauer
2001) .

## Usage

``` r
estorder_SVC(s.max, Hsv = NULL, n.par, n.obs, Hsize, penalty = "lnN", ...)

estorder_IVC(s.max, lndetSigma = NULL, n.par, n.obs, penalty = "BIC", ...)

estorder_max(s.max, ...)

estorder_rkH(s.max, Hsv = NULL, tol = sqrt(.Machine$double.eps), ...)

estorder_MOE(s.max, Hsv = NULL, ...)
```

## Arguments

- s.max:

  (integer) maximum order (maximum state space dimension).

- Hsv:

  vector with the Hankel singular values (must have at least s.max
  entries).

- n.par:

  (s.max+1) dimensional vector with the respective number of parameters.

- n.obs:

  (integer) sample size \\N\\ (or `Inf`).

- Hsize:

  two dimensional integer vector with the number of block rows and block
  columns of the Hankel matrix (`Hsize = c(f,p)`).

- penalty:

  determines the penalty term. See the details below.

- ...:

  optional additional parameters.

- lndetSigma:

  (s.max+1) dimensional vector with the logarithms of the determinants
  of the respective estimated noise covariance matrices.

- tol:

  (small) tolarenace used to determine the rank of the Hankel matrix.

## Value

Either NULL or a list with slots `$s` (the selected/estimated order) and
`$criterion` (an (`s.max+1`) dimensional vector).

## Details

`estorder_max` simply returns the maximum order `s.max` considered.

`estorder_rkH` estimates the order by an estimate of the rank of the
Hankel matrix. If the maximum singular value is smaller than or equal to
`.Machine$double.eps` then the estimate is set to `s=0`. Otherwise the
estimated order is equal to the number of singular values which are
larger than or equal to `tol` times the maximum singular value.

The function `estorder_MOE` searches for a "gap" in the singular values.
The order is set to maximum \\s\\ which satisfies \$\$\ln(\sigma_s) \>
(\ln(\sigma_1)+\ln(\sigma_m))/2\$\$ where \\\sigma_m\\ is the minimum,
non zero singular value. This scheme is also implemented in the N4SID
procedure of the system identification toolbox of MATLAB (Ljung, 1991).

The function `estorder_SVC` implements the so called *Singular Value
Criterion* \$\$svc(s) = \sigma\_{s+1}^2 + c(N)d(s)/N\$\$ (see e.g.
(Bauer 2001) ). Here \\\sigma_s\\ is the \\s\\-th singular value of the
weighted Hankel marix, \\N\\ is the sample size, \\d(s) = 2ms\\ denotes
the number of parameters for a state space model with \\s\\ states (and
\\m\\ outputs) and \\c(N)\\ is a "penalty" term. This term is
\\c(N)=\ln(N)\\ for `penalty = "lnN"` and \\c(N)=fp\ln(N)\\ for
`penalty = "fplnN"`. The estimate of the order is the minimizer of this
criterion.

`estorder_IVC` estimates the order via an information criterion of the
form \$\$ivc(s) = \ln\det\hat\Sigma\_{s} + c(N)d(s)/N\$\$ where
\\\hat\Sigma_s\\ is the estimate of the noise covariace matrix obtained
from a model with order \\s\\. Here the term \\c(N)\\ is chosen as
\\c(N)=2\\ for `penalty = "AIC"` and \\c(N)=\ln(N)\\ for
`penalty = "BIC"`.

Note also that for both routines `estorder_SVC` and `estorder_IVC` one
may also set `penalty` to an arbitrary numeric value! E.g. setting
`penalty=-1` would ensure that `estorder_SVC` always choses the maximum
possible order `s=s.max`.

## References

Bauer D (2001). “Order estimation for subspace methods.” *Automatica*,
**37**(10), 1561 - 1573.
[doi:10.1016/S0005-1098(01)00118-2](https://doi.org/10.1016/S0005-1098%2801%2900118-2)
.

# Compute the (concentrated) conditional log likelihood for a statespace system described by a model template.

This is an internal helper function, used by the function factory
[`ll_FUN`](https://bfunovits.github.io/RLDM/reference/ll_FUN.md). For a
more detailed documentation of the conditional log Likelihood, see
[`ll`](https://bfunovits.github.io/RLDM/reference/ll.md). The
conditional likelihood is computed for the initial state \\a_1\\ given
in the first column `a[,1]` of the matrix `a`.

## Usage

``` r
cll_theta_STSP_cpp(
  th,
  y,
  skip,
  concentrated,
  pi,
  H_pi,
  h_pi,
  L,
  H_L,
  h_L,
  a,
  u,
  dU
)
```

## Arguments

- th:

  \\(K)\\ dimensional vector of "deep" parameters.

- y:

  \\(m,N)\\ matrix with the observed outputs: \\(y_1,y_2,\ldots,y_N)\\.

- skip:

  (integer), skip the first residuals, when computing the sample
  covariance of the residuals.

- concentrated:

  (bool), if TRUE then the *concentrated*, conditional log Likelihood is
  computed

- pi:

  \\(m+s,m+s)\\ matrix, is overwritten with the system matrix \\\[A,B \|
  C,D\]\\.

- H_pi:

  \\(m+s)^2, K)\\ matrix.

- h_pi:

  \\((m+s)^2)\\-dimensional vector. Note that
  `vec(pi) = H_pi*th + h_pi`.

- L:

  \\(m,m)\\ matrix. If (concentrated==FALSE) then L is overwritten with
  the left square of the noise covariance matrix L corresponding to the
  deep parameters th. However, if (concentrated==TRUE) then L is
  overwritten with sample covariance matrix of the computed residuals!

- H_L:

  \\(m^2, K)\\ matrix.

- h_L:

  \\(m^2)\\-dimensional vector. Note that `vec(L) = H_L*th + h_L`.

- a:

  \\(s,N+1)\\ matrix. This matrix is overwritten with the (computed)
  states: \\(a_1,a_2,\ldots,a_N,a\_{N+1})\\. On input `a[,1]` must hold
  the initial state \\a_1\\.

- u:

  \\(m,N)\\ matrix. This matrix is overwritten with (computed)
  residuals: \\(u_1,u_2,\ldots,u_N)\\.

- dU:

  \\(mN,K)\\ matrix or \\(0,0)\\ matrix. This matrix is overwritten with
  the directional derivatives of the residuals. However, if the matrix
  is empty then no derivatives are computed.

## Value

(double) log Likelihood

## See also

[`outputs_ARMA_cpp`](https://bfunovits.github.io/RLDM/reference/outputs_ARMA_cpp.md),
[`residuals_ARMA_cpp`](https://bfunovits.github.io/RLDM/reference/residuals_ARMA_cpp.md),
[`cll_theta_ARMA_cpp`](https://bfunovits.github.io/RLDM/reference/cll_theta_ARMA_cpp.md),
[`outputs_STSP_cpp`](https://bfunovits.github.io/RLDM/reference/outputs_STSP_cpp.md),
[`residuals_STSP_cpp`](https://bfunovits.github.io/RLDM/reference/residuals_STSP_cpp.md),
`cll_theta_STSP_cpp` and
[`solve_de`](https://bfunovits.github.io/RLDM/reference/solve_de.md),
[`solve_inverse_de`](https://bfunovits.github.io/RLDM/reference/solve_de.md)
and [`ll`](https://bfunovits.github.io/RLDM/reference/ll.md).

# Compute the log likelihood for a statespace system described by a model template.

This is an internal helper function, used by the function factory
[`ll_FUN`](https://bfunovits.github.io/RLDM/reference/ll_FUN.md). For a
more detailed documentation of the log Likelihood, see
[`ll_kf`](https://bfunovits.github.io/RLDM/reference/ll_kf.md).

## Usage

``` r
ll_kf_theta_cpp(
  theta,
  y,
  SYS,
  H_SYS,
  h_SYS,
  sigma_L,
  H_sigma_L,
  h_sigma_L,
  VAR,
  P1,
  tol,
  err
)
```

## Arguments

- theta:

  \\(K)\\ dimensional vector of "deep" parameters.

- y:

  \\(m,N)\\ matrix with the observed outputs: \\(y_1,y_2,\ldots,y_N)\\.

- SYS:

  \\(m+s,m+s)\\ matrix, is overwritten with the system matrix \\\[A,B \|
  C,D\]\\.

- H_SYS:

  \\(m+s)^2, K)\\ matrix.

- h_SYS:

  \\((m+s)^2)\\-dimensional vector. Note that
  `vec(SYS) = H_SYS*theta + h_SYS`.

- sigma_L:

  \\(m,m)\\ matrix, is overwritten with the left square root of the
  noise covariance matrix.

- H_sigma_L:

  \\(m^2, K)\\ matrix.

- h_sigma_L:

  \\(m^2)\\-dimensional vector. Note that
  `vec(sigma_L) = H_sigma_L*theta + h_sigma_L`.

- VAR:

  \\(m+s,m+s)\\ matrix, is overwritten with the covariance matrix
  \\\[Q,S \| S',R\] = \[B \| C\] sigma_L sigma_L' \[B', C'\]\\

- P1:

  \\(s,s)\\ matrix, is overwritten with the initial state covariance
  matrix (computed via a Lyapunov equation).

- tol:

  (double) tolerance used by ll_kf_cpp.

- err:

  (double) return err, if the computation of P1 fails.

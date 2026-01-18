# Compute the (concentrated) conditional log likelihood for ARMA models described by a model template.

This internal helper function computes the (concentrated) conditional
log Likelihood of ARMA systems of the form \$\$a_0 y_t + a_1 y\_{t-1} +
\cdots + a_p y\_{t-p} = b_0 u_t + \cdots + b_q u\_{t-q}\$\$ The
conditional likelihood is computed for **zero** initial values
\\u_s=y_s=0\\ for \\s\leq 0\\.

## Usage

``` r
cll_theta_ARMA_cpp(
  th,
  y,
  skip,
  concentrated,
  ib0,
  H_b,
  h_b,
  B1,
  H_B,
  h_B,
  a0,
  A,
  H_A,
  h_A,
  L,
  H_L,
  h_L,
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

  (integer), omit the first "skip" residuals, when computing the
  likelihood.

- concentrated:

  (bool), if TRUE then the *concentrated*, conditional log Likelihood is
  computed

- ib0:

  \\(m, m)\\ matrix, is **overwritten** with the matrix \\b_0^{-1}a_0\\.

- H_b:

  \\(m^2, K)\\ matrix.

- h_b:

  \\((m^2)\\-dimensional vector. Note that `vec(b[0]) = H_b*th + h_b`.

- B1:

  \\(m, mq)\\ matrix, is **overwritten** with
  \\-b_0^{-1}(b_q,...,b_1)\\.

- H_B:

  \\((m^2)\*q, K)\\ matrix.

- h_B:

  \\((m^2)\*q)\\-dimensional vector. Note that
  `vec(-(b[q],...,b[1])) = H_B*th + h_B`.

- a0:

  \\(m, m)\\ matrix, is **overwritten** with \\a_0\\.

- A:

  \\(m, m(q+1))\\ matrix, is **overwritten** with
  \\b_0^{-1}(a_0,...,a_p\\.

- H_A:

  \\((m^2)\*(p+1), K)\\ matrix.

- h_A:

  \\((m^2)\*(p+1))\\-dimensional vector. Note that
  `vec((a[0],a[1],...,a[p])) = H_A*th + h_A`.

- L:

  \\(m,m)\\ matrix. If (concentrated==FALSE) then `L` is **overwritten**
  with the left square \\L\\ of the noise covariance matrix
  \\\Sigma=LL'\\ corresponding to the deep parameters th. However, if
  (concentrated==TRUE) then L is **overwritten** with sample covariance
  matrix of the computed residuals!

- H_L:

  \\(m^2, K)\\ matrix.

- h_L:

  \\(m^2)\\-dimensional vector. Note that `vec(L) = H_L*th + h_L`.

- u:

  \\(m,N)\\ matrix. This matrix is **overwritten** with (computed)
  residuals: \\(u_1,u_2,\ldots,u_N)\\.

- dU:

  \\(mN,(m^2)(p+q+2))\\ matrix or \\(0,0)\\ matrix. If non empty this
  matrix is **overwritten** with the directional derivatives of the
  residuals. However, if the matrix is empty then no derivatives are
  computed.

## Value

(double) log Likelihood

## Details

This function is mainly used by the function factory
[`ll_FUN`](https://bfunovits.github.io/RLDM/reference/ll_FUN.md). For a
more detailed documentation of the (concentrated) conditional log
Likelihood, see
[`ll`](https://bfunovits.github.io/RLDM/reference/ll.md).

The procedure first constructs the ARMA parameter matrices from the
given vector `th` of "deep" parameters.

- AR parameters `vec((a[0],a[1],...,a[p])) = h_A + H_A * th`.

- MA parameters `vec(b[0]) = h_b + H_b * th` and
  `vec(-(b[q],...,b[1])) = h_B + H_B * th`

- Left square root of noise covariance matrix \\\Sigma = LL'\\ and
  `vec(L) = h_L + H_L * th`.

The residuals (and their directional derivatives) are computed with
[`residuals_ARMA_cpp`](https://bfunovits.github.io/RLDM/reference/residuals_ARMA_cpp.md).

## Note

Use this procedure with care!

- The procedure does **not** check the input arguments.

- The procedure **overwrites** some of the input arguments

- The data matrices are organized columnwise (to avoid memory
  shuffling)!

- Note also the non standard representation of the coefficient matrices.

## See also

[`outputs_ARMA_cpp`](https://bfunovits.github.io/RLDM/reference/outputs_ARMA_cpp.md),
[`residuals_ARMA_cpp`](https://bfunovits.github.io/RLDM/reference/residuals_ARMA_cpp.md),
`cll_theta_ARMA_cpp`,
[`outputs_STSP_cpp`](https://bfunovits.github.io/RLDM/reference/outputs_STSP_cpp.md),
[`residuals_STSP_cpp`](https://bfunovits.github.io/RLDM/reference/residuals_STSP_cpp.md),
[`cll_theta_STSP_cpp`](https://bfunovits.github.io/RLDM/reference/cll_theta_STSP_cpp.md)
and
[`solve_de`](https://bfunovits.github.io/RLDM/reference/solve_de.md),
[`solve_inverse_de`](https://bfunovits.github.io/RLDM/reference/solve_de.md)
and [`ll`](https://bfunovits.github.io/RLDM/reference/ll.md).

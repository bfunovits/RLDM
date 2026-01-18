# Residuals of an ARMA system

This internal helper function computes the residuals and the directional
derivatives of the residuals of an ARMA system of the form \$\$a_0 y_t +
a_1 y\_{t-1} + \cdots + a_p y\_{t-p} = b_0 u_t + \cdots + b_q
u\_{t-q}\$\$

## Usage

``` r
residuals_ARMA_cpp(ib0, B1, A, t0, y, u, dU)
```

## Arguments

- ib0:

  \\(m, m)\\ matrix, **inverse** of the coefficient matrix \\b\[0\]\\.

- B1:

  \\(m, mq)\\ matrix, \\-b_0^{-1}(b_q,...,b_1)\\.

- A:

  \\(m, n(q+1))\\ matrix \\b_0^{-1}(a_0,...,a_p\\.

- t0:

  integer, start iteration at t = t0.

- y:

  \\(m, N)\\ matrix with the observed outputs \\(y_1,...,y_N\\.

- u:

  \\(m, N)\\ matrix. This matrix is **overwritten** with the computed
  residuals \\(u_1,...,u_N\\.

- dU:

  \\(mN, m^2(p+q+2))\\ matrix or an empty matrix. If non empty then this
  matrix is **overwritten** with the directional derivatives of the
  vectorized residuals. The \\j\\-th column of `dU` is the derivative of
  \\vec(u)\\ with respect to the \\j\\-th entry of
  \\\mathrm{vec}(a_0,a_1,\ldots,a_p,b_0,\ldots,b_q)\\

## Value

This RcppArmadillo routine returns `NULL` but **overwrites** the input
arguments `u` (and `dU`)!

## Details

Values \\y_t\\, \\u_t\\ for \\t\leq 0\\ are implicitly set to be zero.
However, by starting the iteration with some \\t_0\>1\\ we can enforce
non-zero initial values.

## Note

Use this procedure with care!

- The procedure does **not** check the input arguments. We require \\m =
  n \> 0\\, \\p,q \geq 0\\ and \\1 \leq t_0 \leq N\\.

- The procedure **overwrites** the input argument `u` (and `dU`).

- The data matrices are organized columnwise (to avoid memory
  shuffling)!

- Note also the non standard representation of the coefficient matrices.

## See also

[`outputs_ARMA_cpp`](https://bfunovits.github.io/RLDM/reference/outputs_ARMA_cpp.md),
`residuals_ARMA_cpp`,
[`cll_theta_ARMA_cpp`](https://bfunovits.github.io/RLDM/reference/cll_theta_ARMA_cpp.md),
[`outputs_STSP_cpp`](https://bfunovits.github.io/RLDM/reference/outputs_STSP_cpp.md),
[`residuals_STSP_cpp`](https://bfunovits.github.io/RLDM/reference/residuals_STSP_cpp.md),
[`cll_theta_STSP_cpp`](https://bfunovits.github.io/RLDM/reference/cll_theta_STSP_cpp.md)
and
[`solve_de`](https://bfunovits.github.io/RLDM/reference/solve_de.md),
[`solve_inverse_de`](https://bfunovits.github.io/RLDM/reference/solve_de.md)
and [`ll`](https://bfunovits.github.io/RLDM/reference/ll.md).

## Examples

``` r
# generate a random ARMA(2,1) model (3 outputs, 2 inputs)
p = 2
q = 1
m = 2
model = test_armamod(dim = c(m, m), degrees = c(p,q), digits = 2)

# prepare parameters for "outputs_ARMA_cpp"
A = unclass(model$sys$a)
a0 = A[,,1]
A1 = -A[,,(p+1):2]
dim(A1) = c(m, m*p)
A1 = solve(a0, A1)

B = unclass(model$sys$b)
dim(B) = c(m, m*(q+1))
B = solve(a0, B)

# generate random noise sequence (sample size N = 10)
n.obs = 10
u = matrix(rnorm(n.obs*m), nrow = m, ncol = n.obs)

# generate matrix for the outputs
y = matrix(0, nrow = m, ncol = n.obs)

# compute outputs
t0 = 2   # start iterations from t>=t0=2
outputs_ARMA_cpp(A1, B, t0, u, y)

# recompute the disturbances/residuals from the given outputs:
B = unclass(model$sys$b)
ib0 = B[,,1]
B1 = -B[,,(q+1):2]
dim(B1) = c(m, m*q)
B1 = solve(ib0, B1)

A = unclass(model$sys$a)
dim(A) = c(m, m*(p+1))
A = solve(ib0, A)

ib0 = solve(ib0)

uu = u + 0 # "deep copy" of the disturbances
uu[, t0:(n.obs)] = 0 # clear values for t >= t0
residuals_ARMA_cpp(ib0, B1, A, t0 = 2, y, uu, diag(0))
all.equal(u, uu) # check
#> [1] TRUE

# compute directional derivatives of residuals
dU = matrix(0, nrow = n.obs*m, ncol = (m^2)*(p+q+2))
residuals_ARMA_cpp(ib0, B1, A, t0 = 2, y, uu, dU)
```

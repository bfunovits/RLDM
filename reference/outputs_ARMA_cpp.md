# Outputs of an ARMA systems

This internal helper function computes the outputs of an ARMA system
\$\$a_0 y_t + a_1 y\_{t-1} + \cdots + a_p y\_{t-p} = b_0 u_t + \cdots +
b_q u\_{t-q}\$\$

## Usage

``` r
outputs_ARMA_cpp(A1, B, t0, u, y)
```

## Arguments

- A1:

  \\(m, mp)\\ matrix \\-a_0^{-1}(a_p,...,a_1)\\.

- B:

  \\(m, n(q+1))\\ matrix \\a_0^{-1}(b_0,...,b_q\\.

- t0:

  integer, start iteration at t = t0.

- u:

  \\(n, N)\\ matrix with the inputs \\(u_1,...,u_N\\.

- y:

  \\(m, N)\\ matrix with the outputs \\(y_1,...,y_N\\.

## Value

This RcppArmadillo routine returns `NULL` but **overwrites** the input
argument `y` with the computed outputs!

## Details

Values \\y_t\\, \\u_t\\ for \\t\leq 0\\ are implicitly set to be zero.
However, by starting the iteration with some \\t_0\>1\\ we can enforce
non-zero initial values.

## Note

Use this procedure with care!

- The procedure does **not** check the input arguments. We require \\m
  \> 0\\, \\p \geq 0\\, \\n(q+1) \geq 0\\ and \\1 \leq t_0 \leq N\\.

- The procedure **overwrites** the input argument `y`.

- The data matrices are organized columnwise (to avoid memory
  shuffling)!

- Note also the non standard representation of the coefficient matrices.

## See also

`outputs_ARMA_cpp`,
[`residuals_ARMA_cpp`](https://bfunovits.github.io/RLDM/reference/residuals_ARMA_cpp.md),
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
m = 3
n = 2
model = test_armamod(dim = c(m, n), degrees = c(p,q), digits = 2)
A = unclass(model$sys$a)
a0 = A[,,1]
A1 = -A[,,(p+1):2]
dim(A1) = c(m, m*p)
A1 = solve(a0, A1)
B = unclass(model$sys$b)
dim(B) = c(m, n*(q+1))
B = solve(a0, B)

# generate random noise sequence (sample size N = 10)
n.obs = 10
u = matrix(rnorm(n.obs*n), nrow = n, ncol = n.obs)
print(u)
#>            [,1]       [,2]       [,3]      [,4]       [,5]       [,6]
#> [1,]  0.4962119 -0.4722538  0.8577367 0.8148567 -0.4458320  0.5780798
#> [2,] -1.6218023  0.5796702 -0.4830497 0.9973271 -0.1576797 -1.0053253
#>             [,7]       [,8]       [,9]      [,10]
#> [1,] -1.68615924 -0.2180973 -0.8818466 -0.1040246
#> [2,]  0.04228322  0.6227171 -0.3819943 -1.4566517

# generate matrix for the outputs
y = matrix(0, nrow = m, ncol = n.obs)

# call outputs_ARMA_cpp()
outputs_ARMA_cpp(A1, B, t0 = 2, u, y) # start with t>=2
print(u)
#>            [,1]       [,2]       [,3]      [,4]       [,5]       [,6]
#> [1,]  0.4962119 -0.4722538  0.8577367 0.8148567 -0.4458320  0.5780798
#> [2,] -1.6218023  0.5796702 -0.4830497 0.9973271 -0.1576797 -1.0053253
#>             [,7]       [,8]       [,9]      [,10]
#> [1,] -1.68615924 -0.2180973 -0.8818466 -0.1040246
#> [2,]  0.04228322  0.6227171 -0.3819943 -1.4566517
print(y)  # y is overwritten with the computed outputs
#>      [,1]       [,2]       [,3]       [,4]        [,5]      [,6]     [,7]
#> [1,]    0 -0.1246525  1.1692922 -0.5552006 -0.69889674 -1.843248 9.744668
#> [2,]    0 -0.5444621 -0.4195822  2.5775542  0.08981804 -5.384511 3.330721
#> [3,]    0 -0.5493520  0.7354001  4.4236591 -5.00962349 -9.473436 5.621456
#>           [,8]      [,9]      [,10]
#> [1,] -6.957517  13.96052 -26.019929
#> [2,] 23.103350 -14.47345  -8.799821
#> [3,] 34.992335 -60.59104 -14.878826
```

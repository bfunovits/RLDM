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
#>            [,1]       [,2]       [,3]       [,4]       [,5]       [,6]
#> [1,] -0.6542588 -0.2331897 -0.7314605  0.9556143  1.2420865 -0.5361896
#> [2,]  3.1378073 -1.3051994 -0.7895210 -0.6484152 -0.2068317  0.5564057
#>            [,7]       [,8]       [,9]      [,10]
#> [1,] -0.2076044 -0.9593693 -0.0936314 -0.1847930
#> [2,]  0.7537292 -0.9031647  1.3370649 -0.1696349

# generate matrix for the outputs
y = matrix(0, nrow = m, ncol = n.obs)

# call outputs_ARMA_cpp()
outputs_ARMA_cpp(A1, B, t0 = 2, u, y) # start with t>=2
print(u)
#>            [,1]       [,2]       [,3]       [,4]       [,5]       [,6]
#> [1,] -0.6542588 -0.2331897 -0.7314605  0.9556143  1.2420865 -0.5361896
#> [2,]  3.1378073 -1.3051994 -0.7895210 -0.6484152 -0.2068317  0.5564057
#>            [,7]       [,8]       [,9]      [,10]
#> [1,] -0.2076044 -0.9593693 -0.0936314 -0.1847930
#> [2,]  0.7537292 -0.9031647  1.3370649 -0.1696349
print(y)  # y is overwritten with the computed outputs
#>      [,1]      [,2]      [,3]       [,4]      [,5]      [,6]      [,7]     [,8]
#> [1,]    0  2.518868 -8.510833   9.659450   8.20804  71.75834 105.29152 205.0337
#> [2,]    0 -6.256375  5.924761   5.239423  60.85647 101.53580 204.78371 142.0783
#> [3,]    0 -5.868227 -3.055842 -27.333396 -28.43466 -54.77651  29.65815 231.3036
#>            [,9]      [,10]
#> [1,]   78.20851  -383.7994
#> [2,] -215.49116 -1696.0254
#> [3,]  876.36735  1923.8205
```

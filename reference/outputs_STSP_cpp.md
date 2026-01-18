# Outputs of a statespace system

This internal helper function computes the outputs and states for a
statespace system of the form \$\$a\_{t+1} = A a_t + B u_t, \\ y_t = C
a_t + D u_t\$\$

## Usage

``` r
outputs_STSP_cpp(A, B, C, D, u, a, y)
```

## Arguments

- A:

  \\(s,s)\\ matrix.

- B:

  \\(s,n)\\ matrix.

- C:

  \\(m,s)\\ matrix.

- D:

  \\(m,n)\\ matrix.

- u:

  \\(n,N)\\ matrix with the inputs/disturbances:
  \\(u_1,u_2,\ldots,u_N)\\.

- a:

  \\(s,N+1)\\ matrix. This matrix is overwritten with the (computed)
  states: \\(a_1,a_2,\ldots,a_N,a\_{N+1})\\. On input `a[,1]` must hold
  the initial state \\a_1\\.

- y:

  \\(m,N)\\ matrix. This matrix is overwritten with (computed) outputs:
  \\(y_1,y_2,\ldots,y_N)\\.

## Value

This RcppArmadillo routine returns `NULL` but **overwrites** the input
arguments `a` and `u`!

## Note

Use this procedure with care!

- The procedure does **not** check the input arguments.

- The procedure **overwrites** the input arguments `a` and `u`.

- The data matrices are organized columnwise (to avoid memory
  shuffling)!

## See also

[`outputs_ARMA_cpp`](https://bfunovits.github.io/RLDM/reference/outputs_ARMA_cpp.md),
[`residuals_ARMA_cpp`](https://bfunovits.github.io/RLDM/reference/residuals_ARMA_cpp.md),
[`cll_theta_ARMA_cpp`](https://bfunovits.github.io/RLDM/reference/cll_theta_ARMA_cpp.md),
`outputs_STSP_cpp`,
[`residuals_STSP_cpp`](https://bfunovits.github.io/RLDM/reference/residuals_STSP_cpp.md),
[`cll_theta_STSP_cpp`](https://bfunovits.github.io/RLDM/reference/cll_theta_STSP_cpp.md)
and
[`solve_de`](https://bfunovits.github.io/RLDM/reference/solve_de.md),
[`solve_inverse_de`](https://bfunovits.github.io/RLDM/reference/solve_de.md)
and [`ll`](https://bfunovits.github.io/RLDM/reference/ll.md).

## Examples

``` r
# generate a random statespace model (3 outputs, 2 inputs and 4 states)
m = 3
n = 2
s = 4
model = test_stspmod(dim = c(m, n), s = s, digits = 2)

# generate random noise sequence (sample size N = 10)
n.obs = 10
u = matrix(rnorm(n.obs*n), nrow = n, ncol = n.obs)
print(u)
#>           [,1]       [,2]      [,3]         [,4]       [,5]       [,6]
#> [1,] 0.1621131  0.2754776 0.1546232 -0.241994688  1.3033774 -1.2047453
#> [2,] 0.2767673 -0.5245635 0.8115842  0.004326083 -0.2845295  0.4190233
#>            [,7]      [,8]      [,9]      [,10]
#> [1,] -0.6807778 0.4564316 0.5098132 -0.8951889
#> [2,] -0.2271627 1.4490524 0.1997966 -1.0671328

# generate matrix for the state sequence
a = matrix(0, nrow = s, ncol = n.obs+1)
a[,1] = rnorm(s) # random initial state a[1]
print(a)
#>            [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11]
#> [1,] -0.2910083    0    0    0    0    0    0    0    0     0     0
#> [2,] -0.4180777    0    0    0    0    0    0    0    0     0     0
#> [3,] -0.9646559    0    0    0    0    0    0    0    0     0     0
#> [4,]  0.3606955    0    0    0    0    0    0    0    0     0     0

# generate matrix for the outputs
y = matrix(0, nrow = m, ncol = n.obs)

# call outputs_STSP_cpp()
outputs_STSP_cpp(model$sys$A, model$sys$B, model$sys$C, model$sys$D, u, a, y)
print(u)
#>           [,1]       [,2]      [,3]         [,4]       [,5]       [,6]
#> [1,] 0.1621131  0.2754776 0.1546232 -0.241994688  1.3033774 -1.2047453
#> [2,] 0.2767673 -0.5245635 0.8115842  0.004326083 -0.2845295  0.4190233
#>            [,7]      [,8]      [,9]      [,10]
#> [1,] -0.6807778 0.4564316 0.5098132 -0.8951889
#> [2,] -0.2271627 1.4490524 0.1997966 -1.0671328
print(a)  # a is overwritten with the computed states
#>            [,1]       [,2]       [,3]       [,4]       [,5]       [,6]
#> [1,] -0.2910083  0.1470632 -2.2887687  0.9914409 -5.8038118   1.914746
#> [2,] -0.4180777 -0.8782638  0.4831436 -7.2554811  9.1149031 -29.334267
#> [3,] -0.9646559  0.4475396 -0.1555441 -0.2768109 -0.2716191   2.065672
#> [4,]  0.3606955 -0.8569591 -1.0539385 -1.2723679 -3.8823499  -4.454104
#>            [,7]      [,8]      [,9]      [,10]      [,11]
#> [1,] -19.063660  20.48665 -76.39810  103.28393 -299.50850
#> [2,]  36.521809 -95.20381 147.61934 -373.03877  653.31007
#> [3,]  -4.600986  15.56166 -31.20729   66.34195 -123.80344
#> [4,]  -9.358324 -13.67995 -16.20940  -40.31421  -39.08654
print(y)  # y is overwritten with the computed outputs
#>            [,1]       [,2]       [,3]      [,4]      [,5]       [,6]       [,7]
#> [1,] -0.1391880  1.6938865 -3.5804536  5.917635 -8.610496  15.436191 -36.811577
#> [2,]  0.7300401 -0.3661558 -0.7781827  1.092353 -1.795031  -2.244316  -2.058397
#> [3,] -1.3560614 -0.8337046 -3.0828712 -6.261713 -2.359993 -23.478763  -3.401920
#>            [,8]        [,9]       [,10]
#> [1,]  78.914011 -164.549864  336.015143
#> [2,]  -1.560201   -8.928909   -2.831807
#> [3,] -51.196342   -7.538429 -159.236736
```

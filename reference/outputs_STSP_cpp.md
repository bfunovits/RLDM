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
#>             [,1]       [,2]       [,3]       [,4]     [,5]       [,6]
#> [1,]  0.02257665 -0.6234985 -1.2689028 -1.6233913 1.788238  0.5762031
#> [2,] -1.73917082  0.8530448  0.4344063  0.2914045 1.657197 -0.1284107
#>             [,7]       [,8]       [,9]      [,10]
#> [1,] -0.09759244 -0.5895001 -1.2984312 -0.5606992
#> [2,] -0.22007338 -0.4481922 -0.5562581  0.4527295

# generate matrix for the state sequence
a = matrix(0, nrow = s, ncol = n.obs+1)
a[,1] = rnorm(s) # random initial state a[1]
print(a)
#>            [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11]
#> [1,] -0.9450664    0    0    0    0    0    0    0    0     0     0
#> [2,]  0.4197889    0    0    0    0    0    0    0    0     0     0
#> [3,] -0.3598256    0    0    0    0    0    0    0    0     0     0
#> [4,] -0.8936101    0    0    0    0    0    0    0    0     0     0

# generate matrix for the outputs
y = matrix(0, nrow = m, ncol = n.obs)

# call outputs_STSP_cpp()
outputs_STSP_cpp(model$sys$A, model$sys$B, model$sys$C, model$sys$D, u, a, y)
print(u)
#>             [,1]       [,2]       [,3]       [,4]     [,5]       [,6]
#> [1,]  0.02257665 -0.6234985 -1.2689028 -1.6233913 1.788238  0.5762031
#> [2,] -1.73917082  0.8530448  0.4344063  0.2914045 1.657197 -0.1284107
#>             [,7]       [,8]       [,9]      [,10]
#> [1,] -0.09759244 -0.5895001 -1.2984312 -0.5606992
#> [2,] -0.22007338 -0.4481922 -0.5562581  0.4527295
print(a)  # a is overwritten with the computed states
#>            [,1]      [,2]        [,3]       [,4]       [,5]      [,6]
#> [1,] -0.9450664  1.555291  6.09385706 -10.090029 -15.892886  86.11642
#> [2,]  0.4197889  2.359183  0.05989388 -22.516358  20.532756 108.67984
#> [3,] -0.3598256  1.702935 -4.48260356  -4.470646  41.449191 -31.19117
#> [4,] -0.8936101 -1.237574  1.37322571  -1.012969  -6.375628  18.55043
#>            [,7]      [,8]      [,9]      [,10]      [,11]
#> [1,]  -31.85077 -504.8281  957.9306  1791.3988  -8859.520
#> [2,] -289.42027 -264.5084 2263.5040 -1708.2699 -11974.277
#> [3,] -225.31993  545.9161  622.1753 -4458.9487   2724.218
#> [4,]   11.71331 -129.7429  124.3142   640.4324  -1777.915
print(y)  # y is overwritten with the computed outputs
#>            [,1]      [,2]      [,3]      [,4]       [,5]       [,6]       [,7]
#> [1,] -0.9100117 -2.374274  3.524024 -20.72283 -17.446762  150.92285  -94.13199
#> [2,] -2.9260985  3.220696  4.930470 -12.77119   3.299185   79.67083 -138.00689
#> [3,]  0.4386055 -1.981754 -8.882333  47.49607  13.879609 -315.65438  371.03509
#>           [,8]      [,9]      [,10]
#> [1,] -817.2299  1850.874  2458.5763
#> [2,] -302.0919  1330.514   -82.6597
#> [3,] 1469.2485 -4615.012 -2494.0035
```

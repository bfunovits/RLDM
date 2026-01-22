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
#>           [,1]       [,2]     [,3]       [,4]      [,5]      [,6]      [,7]
#> [1,] 0.2173447 -0.8353325 1.496910 -0.6515175 1.1346293 0.6763501 0.3363965
#> [2,] 0.4471227  0.1409061 1.362301 -0.5174458 0.6143891 1.5933323 0.5924614
#>              [,8]       [,9]      [,10]
#> [1,]  0.008539218  0.9694180 -1.5718134
#> [2,] -1.085621021 -0.3930895 -0.6153608

# generate matrix for the state sequence
a = matrix(0, nrow = s, ncol = n.obs+1)
a[,1] = rnorm(s) # random initial state a[1]
print(a)
#>            [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11]
#> [1,] -0.2412325    0    0    0    0    0    0    0    0     0     0
#> [2,]  0.2719050    0    0    0    0    0    0    0    0     0     0
#> [3,] -1.0406623    0    0    0    0    0    0    0    0     0     0
#> [4,]  0.5281270    0    0    0    0    0    0    0    0     0     0

# generate matrix for the outputs
y = matrix(0, nrow = m, ncol = n.obs)

# call outputs_STSP_cpp()
outputs_STSP_cpp(model$sys$A, model$sys$B, model$sys$C, model$sys$D, u, a, y)
print(u)
#>           [,1]       [,2]     [,3]       [,4]      [,5]      [,6]      [,7]
#> [1,] 0.2173447 -0.8353325 1.496910 -0.6515175 1.1346293 0.6763501 0.3363965
#> [2,] 0.4471227  0.1409061 1.362301 -0.5174458 0.6143891 1.5933323 0.5924614
#>              [,8]       [,9]      [,10]
#> [1,]  0.008539218  0.9694180 -1.5718134
#> [2,] -1.085621021 -0.3930895 -0.6153608
print(a)  # a is overwritten with the computed states
#>            [,1]      [,2]      [,3]      [,4]       [,5]      [,6]      [,7]
#> [1,] -0.2412325 1.9396248 -3.363503  1.199680  -5.903463 -6.793850 14.033443
#> [2,]  0.2719050 0.2963779 -6.528416  6.617881  -8.596936 18.492112 15.799060
#> [3,] -1.0406623 0.6957607 -2.229484  6.207147   2.824222  2.378189 11.417787
#> [4,]  0.5281270 1.3115528 -1.932456 11.033960 -10.124434 12.781377 -5.512745
#>           [,8]     [,9]      [,10]     [,11]
#> [1,] -10.38477 66.56542  -63.98627 107.75464
#> [2,] -14.51215 25.25327 -148.70864  96.55825
#> [3,] -46.86630 36.18867 -130.58331 200.13577
#> [4,] -17.59122 26.03063  -71.92477 168.91551
print(y)  # y is overwritten with the computed outputs
#>            [,1]       [,2]       [,3]      [,4]       [,5]      [,6]      [,7]
#> [1,] -0.3002033  1.2092360   1.070953  3.031778  -1.281007 -9.037767  4.143743
#> [2,]  0.1471419  0.4491558  -2.837089 13.045537  -3.369193 23.436127  6.933759
#> [3,]  0.9877319 -1.7203977 -11.166463 13.992650 -13.182688 52.098444 19.372472
#>           [,8]      [,9]      [,10]
#> [1,] -23.67842  60.03310  -41.93191
#> [2,] -50.17609  17.76221 -183.69595
#> [3,] -19.67283 -24.13063 -267.10533
```

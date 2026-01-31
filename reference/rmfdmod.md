# Constructor for RMFD Models

**\[experimental\]**

## Usage

``` r
rmfdmod(sys, sigma_L = NULL, names = NULL, label = NULL)
```

## Arguments

- sys:

  [`rationalmatrices::rmfd()`](https://bfunovits.github.io/rationalmatrices/reference/rmfd.html)
  object

- sigma_L:

  Left-factor of noise covariance, i.e. the covariance \\\sigma\\ is
  obtained as `sigma_L * t(sigma_L)`. If `sigma_L` is a vector of
  dimension \\n\\, where \\n\\ is the input dimension, only the diagonal
  elements are parametrized. If it is a vector of dimension \\n^2\\,
  then the elements of `sigma_L` are filled column by column.

- names:

  optional vector of character strings

- label:

  optional character string

## Value

Object of class `rmfdmod`.

## Details

A right-matrix fraction description (RMFD) plus parameterisation of
noise covariance. In (Hannan and Deistler 2012) , RMFDs are also called
dynamic adjustment forms. Internally, MFDs are lists with slots `sys`,
`sigma_L`, `names`, `label`. **Many of the generic functions which
construct derived objects like the autocovariance
[`autocov()`](https://bfunovits.github.io/RLDM/reference/autocov.md) are
not yet implemented for `rmfdmod` objects.**

## References

Hannan EJ, Deistler M (2012). *The Statistical Theory of Linear
Systems*, Classics in Applied Mathematics. SIAM, Philadelphia.
Originally published: John Wiley & Sons, New York, 1988.

## Examples

``` r
y = rmfdmod(sys = test_rmfd(dim = c(3,2), degrees = c(2,2)))
y
#> RMFD model [3,2] with orders p = 2 and q = 2
#> right factor polynomial c(z):
#>      z^0 [,1]  [,2]   z^1 [,1]       [,2]   z^2 [,1]       [,2]
#> [1,]        1     0 -0.6002596  1.5326106 -1.0264209  0.2568837
#> [2,]        0     1  2.1873330 -0.2357004 -0.7104066 -0.2466919
#> left factor polynomial d(z):
#>         z^0 [,1]       [,2]   z^1 [,1]        [,2]   z^2 [,1]       [,2]
#> [1,] -0.34754260 -0.7849045  0.9189966 -1.61788271  0.3011534 -0.8497043
#> [2,] -0.95161857 -1.6679419 -0.5753470 -0.05556197  0.1056762 -1.0241288
#> [3,] -0.04502772 -0.3802265  0.6079643  0.51940720 -0.6407060  0.1176466
#> Left square root of noise covariance Sigma:
#>      u[1] u[2]
#> u[1]    1    0
#> u[2]    0    1
```

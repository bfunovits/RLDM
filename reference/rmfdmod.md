# Constructor for RMFD Models

**\[experimental\]**

## Usage

``` r
rmfdmod(sys, sigma_L = NULL, names = NULL, label = NULL)
```

## Arguments

- sys:

  `rationalmatrices::rmfd()` object

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
#>      z^0 [,1]  [,2]   z^1 [,1]       [,2]   z^2 [,1]        [,2]
#> [1,]        1     0 -0.1014337  1.8955234 -0.6642255 -0.03049735
#> [2,]        0     1 -0.5479073 -0.8066148 -0.5150376 -1.47695345
#> left factor polynomial d(z):
#>        z^0 [,1]       [,2]  z^1 [,1]       [,2]   z^2 [,1]       [,2]
#> [1,] -0.4443564  1.1135145 0.4181386 -0.2823573  1.0879104 -0.4533885
#> [2,]  3.2097522 -1.5318151 2.1153912  1.3348598  0.1838209  0.2227293
#> [3,] -0.2778233 -0.3328608 0.2221793  1.0664345 -0.6965530 -0.1540822
#> Left square root of noise covariance Sigma:
#>      u[1] u[2]
#> u[1]    1    0
#> u[2]    0    1
```

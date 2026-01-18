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
#>      z^0 [,1]  [,2]   z^1 [,1]       [,2]  z^2 [,1]       [,2]
#> [1,]        1     0 -0.1082741 -0.1754903 1.7510323 -0.8972961
#> [2,]        0     1 -0.3803508  1.4467596 0.2942557  0.6696797
#> left factor polynomial d(z):
#>        z^0 [,1]       [,2]  z^1 [,1]       [,2]    z^2 [,1]       [,2]
#> [1,]  1.1192787  0.1722573 0.1290636 -1.1437321 -0.08693342  0.6529577
#> [2,] -1.3362424  0.5356509 0.6496108  0.6692810 -0.25264209 -0.2562246
#> [3,] -0.4130302 -0.3165201 1.3799802 -0.7006766  0.56466801  0.5767858
#> Left square root of noise covariance Sigma:
#>      u[1] u[2]
#> u[1]    1    0
#> u[2]    0    1
```

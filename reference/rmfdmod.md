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
#>      z^0 [,1]  [,2]   z^1 [,1]      [,2]  z^2 [,1]      [,2]
#> [1,]        1     0 -0.5349281 1.5632988 0.8315751 1.2387464
#> [2,]        0     1  2.0605885 0.7144591 1.5079436 0.1634034
#> left factor polynomial d(z):
#>        z^0 [,1]       [,2]   z^1 [,1]       [,2]   z^2 [,1]       [,2]
#> [1,]  0.8804698 -0.3063182 -0.2733697 -0.1739394 -1.9072101 -1.6174181
#> [2,] -0.7263570  0.9326638  1.2064370  1.3499435  0.4256843  0.3100937
#> [3,]  0.4788533  1.9304179  1.0867108 -0.6818069  1.0034773  0.7847326
#> Left square root of noise covariance Sigma:
#>      u[1] u[2]
#> u[1]    1    0
#> u[2]    0    1
```

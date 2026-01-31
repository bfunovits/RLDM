# Constructor for LMFD (ARMA) Models

A left-matrix fraction description (LMFD) plus parameterisation of noise
covariance.

## Usage

``` r
armamod(sys, sigma_L = NULL, names = NULL, label = NULL)
```

## Arguments

- sys:

  [`rationalmatrices::lmfd()`](https://bfunovits.github.io/rationalmatrices/reference/lmfd.html)
  or
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

Object of class `armamod`.

## Details

In Hannan, Deistler (2012, page 7), RMFDs are also called dynamic
adjustment forms. Internally, MFDs are lists with slots `sys`,
`sigma_L`, `names`, `label`.

## Examples

``` r
x = armamod(sys = lmfd(c(1, 0.5), 1), sigma_L = diag(1))
x
#> ARMA model [1,1] with orders p = 1 and q = 0
#> AR polynomial a(z):
#>      z^0 [,1] z^1 [,1]
#> [1,]        1      0.5
#> MA polynomial b(z):
#>      z^0 [,1]
#> [1,]        1
#> Left square root of noise covariance Sigma:
#>      u[1]
#> u[1]    1
```

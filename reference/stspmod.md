# Creator for stspmod class

Creator for stspmod class

## Usage

``` r
stspmod(sys, sigma_L = NULL, names = NULL, label = NULL)
```

## Arguments

- sys:

  [`rationalmatrices::stsp()`](https://bfunovits.github.io/rationalmatrices/reference/stsp.html)
  object

- sigma_L:

  noise covariance left

- names:

  optional vector of character strings

- label:

  optional chracter string

## Value

Object of class `stspmod`.

## Examples

``` r
x = stspmod(sys = test_stsp(dim = c(2,2), s = 2), sigma_L = diag(2))
x
#> state space model [2,2] with s = 2 states
#>            s[1]        s[2]      u[1]       u[2]
#> s[1] -0.2504141  0.24679921 0.2551648  0.5368560
#> s[2]  1.5495553 -0.73677154 0.2774468 -0.4604856
#> x[1] -1.0971396 -1.28000894 1.0000000  0.0000000
#> x[2]  0.9255112  0.07664366 0.0000000  1.0000000
#> Left square root of noise covariance Sigma:
#>      u[1] u[2]
#> u[1]    1    0
#> u[2]    0    1
```

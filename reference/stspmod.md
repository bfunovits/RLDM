# Creator for stspmod class

Creator for stspmod class

## Usage

``` r
stspmod(sys, sigma_L = NULL, names = NULL, label = NULL)
```

## Arguments

- sys:

  `rationalmatrices::stsp()` object

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
#>             s[1]       s[2]       u[1]       u[2]
#> s[1]  0.87012272 -0.2142707 -0.4749987 -0.1421839
#> s[2]  0.01924235 -0.8372097 -0.1486154 -0.1893236
#> x[1]  0.18011289 -0.5030431  1.0000000  0.0000000
#> x[2] -0.21975630 -0.8403386  0.0000000  1.0000000
#> Left square root of noise covariance Sigma:
#>      u[1] u[2]
#> u[1]    1    0
#> u[2]    0    1
```

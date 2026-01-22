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
#>            s[1]       s[2]      u[1]      u[2]
#> s[1] -0.1385598 -0.1071151 0.9055928 0.4116807
#> s[2]  1.1899153  0.7262129 0.7234385 0.9463023
#> x[1]  0.1320469 -1.6041407 1.0000000 0.0000000
#> x[2] -0.4050430 -0.8048852 0.0000000 1.0000000
#> Left square root of noise covariance Sigma:
#>      u[1] u[2]
#> u[1]    1    0
#> u[2]    0    1
```

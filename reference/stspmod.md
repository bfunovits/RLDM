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
#>            s[1]        s[2]       u[1]      u[2]
#> s[1] -0.1285609 -0.03155624  1.2988472 -1.325885
#> s[2]  1.1038513 -1.41231378 -0.2699241  1.311557
#> x[1]  1.0584341 -0.44596116  1.0000000  0.000000
#> x[2] -1.8095109  0.52668317  0.0000000  1.000000
#> Left square root of noise covariance Sigma:
#>      u[1] u[2]
#> u[1]    1    0
#> u[2]    0    1
```

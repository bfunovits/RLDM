# Coerce to State Space Model

The function
[`rationalmatrices::as.stsp.pseries()`](https://bfunovits.github.io/rationalmatrices/reference/as.stsp.html)
calls
[`rationalmatrices::pseries2stsp()`](https://bfunovits.github.io/rationalmatrices/reference/pseries2stsp.html)
with default parameters. Of course the
[`rationalmatrices::pseries()`](https://bfunovits.github.io/rationalmatrices/reference/pseries.html)
object must contain sufficiently many lags. NOT YET implemented

## Usage

``` r
as.stspmod(obj, ...)

# S3 method for class 'armamod'
as.stspmod(obj, ...)
```

## Arguments

- obj:

  object

- ...:

  optional additional parameters

- method:

  character string

## Value

object of class
[`stspmod()`](https://bfunovits.github.io/RLDM/reference/stspmod.md).

## Examples

``` r
# Convert an ARMA model to state space representation
arma_model = armamod(sys = lmfd(c(1, 0.5, 0.2), c(1, -0.3)), sigma_L = diag(1))
ss_model = as.stspmod(arma_model)
ss_model
#> state space model [1,1] with s = 2 states
#>      s[1] s[2] u[1]
#> s[1] -0.5 -0.2  0.2
#> s[2]  1.0  0.0 -0.8
#> x[1]  0.0  1.0  1.0
#> Left square root of noise covariance Sigma:
#>      u[1]
#> u[1]    1
```

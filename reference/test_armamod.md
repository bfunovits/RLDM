# Create Test ARMA model

This simple tool may be used to create a random ARMA model \$\$y_t + a_1
y\_{t-1} + \cdots + a_p y\_{t-p} = b_0 u_t + b_1 u\_{t-1} + \cdots + b_q
u\_{t-q}\$\$ with given order \\(p,q)\\.

## Usage

``` r
test_armamod(
  dim = c(1, 1),
  degrees = c(1, 1),
  b0 = NULL,
  sigma_L = NULL,
  digits = NULL,
  bpoles = NULL,
  bzeroes = NULL,
  n.trials = 100
)
```

## Arguments

- dim:

  integer vector `c(m,n)`.

- degrees:

  integer vector `c(p,q)`.

- b0:

  \\(m,n)\\ dimensional matrix (or `NULL`). See the details below.

- sigma_L:

  \\(n,n)\\ dimensional matrix (or `NULL`). See the details below.

- digits:

  integer, if non NULL then the randomly generated numbers are rounded
  to "digits" number of decimal places.

- bpoles:

  lower bound for the moduli of the poles of the corresponding transfer
  function (or NULL).

- bzeroes:

  lower bound for the moduli of the zeroes of the corresponding
  tranmsfer function (or NULL). This parameter is ignored for non-square
  matrices (m != n).

- n.trials:

  maximum number of trials.

## Value

[`armamod()`](https://bfunovits.github.io/RLDM/reference/armamod.md)
object, which represents the generated ARMA model.

## Details

We require \\m\>0\\ and \\p\geq 0\\.

The \\b_0\\ matrix defaults to a \\(m,n)\\-dimensional diagonal matrix
with ones on the diagonal (`diag(x=1, nrow = m, ncol = n)`). However,
one may also pass an arbitray (compatible) matrix to the procedure. This
matrix may contain `NA`'s, which then are replaced by random numbers.

The \\sigma_L\\ matrix defaults to a \\(n,n)\\-dimensional lower,
triangular matrix However, one may also pass an arbitrary (compatible)
\\sigma_L\\ matrix to the procedure.

The user may prescribe lower bounds for the moduli of the zeroes and/or
poles of the transfer function \$\$k(z) = a^{-1}(z) b(z).\$\$ In this
case the procedure simply generates (up to n.trials) random models until
a model is found which satisfies the constraint. The standard deviation
of the normal distribution, which is used to generate the random
entries, is decreased in each step. Of course this is a very crude
method and it may fail or need a very large number of randomly generated
matrices.

## Examples

``` r
### generate a random ARMA(1,1) model (with two outputs)
### we require that the model is stable and minimum phase
model = try(test_armamod(dim = c(2,2), degrees = c(1,1), digits = 2, bpoles = 1, bzeroes = 1))
if (!inherits(model, 'try-error')) {
   print(model)
   print(abs(poles(model$sys)))
   print(abs(zeroes(model$sys)))
}
#> ARMA model [2,2] with orders p = 1 and q = 1
#> AR polynomial a(z):
#>      z^0 [,1]  [,2] z^1 [,1]  [,2]
#> [1,]        1     0    -0.24  0.59
#> [2,]        0     1     0.12  0.55
#> MA polynomial b(z):
#>      z^0 [,1]  [,2] z^1 [,1]  [,2]
#> [1,]        1     0    -0.38 -0.44
#> [2,]        0     1     0.13 -0.31
#> Left square root of noise covariance Sigma:
#>      u[1] u[2]
#> u[1] 0.55 0.00
#> u[2] 1.43 0.79
#> [1] 1.584130 3.112729
#> [1] 2.390457 2.390457
```

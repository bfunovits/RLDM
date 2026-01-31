# Create Test state space Model

This simple tool may be used to create a random, state space model
\$\$a\_{t+1} = A a_t + B u_t \mbox{ and } y_t = C a_t + D u_t.\$\$

## Usage

``` r
test_stspmod(
  dim = c(1, 1),
  s = NULL,
  nu = NULL,
  D = NULL,
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

- s:

  state dimension (or NULL).

- nu:

  vector with the Kronecker indices (or `NULL`). Either the state space
  dimension `s` or the Kronecker indices `nu` must be non `NULL`. If
  both parameters are given, then the parameter `s` is ignored.

- D:

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

[`stspmod()`](https://bfunovits.github.io/RLDM/reference/stspmod.md)
object, which represents the generated model.

## Details

If the Kronecker indices (parameter `nu`) are given, then a state space
model in echelon canonical form is generated. This means that some of
the entries of the \\A,B,C\\ matrices are fixed to be one or zero and
the others are considerd as "free". See also
[`rationalmatrices::Kronecker-Indices()`](https://bfunovits.github.io/rationalmatrices/reference/Kronecker-Indices.html).
The entries of the \\A, B, C\\ matrices, which are not a priori fixed
are randomly generated.

If only the state dimension \\s\\ (parameter `s`) is given, then all
entries of the \\A, B, C\\ matrices are considered as "free".

The \\D\\ matrix defaults to a \\(m,n)\\-dimensional diagonal matrix
with ones on the diagonal (`diag(x=1, nrow = m, ncol = n)`). However,
one may also pass an arbitray (compatible) \\D\\ matrix to the
procedure. This matrix may contain `NA`'s, which then are replaced by
random numbers.

The \\sigma_L\\ matrix defaults to a \\(n,n)\\-dimensional lower,
triangular matrix However, one may also pass an arbitray (compatible)
\\sigma_L\\ matrix to the procedure.

The user may prescribe lower bounds for the moduli of the zeroes and/or
poles of the transfer function \$\$k(z) = C(I_m z{-1} - A)^{-1} B +
D.\$\$ In this case the procedure simply generates (up to n.trials)
random models until a model is found which satisfies the constraint. The
standard deviation of the normal distribution, which is used to generate
the random entries, is decreased in each step. Of course this is a very
crude method and it may fail or need a very large number of randomly
generated matrices.

Note also, that the generated model may be non-minimal.

## Examples

``` r
### random state space model with two outputs and state dimension s = 3
### The model is required to be stable and minimum phase
model = try(test_stspmod(dim = c(2,2), s = 3, digits = 2, bpoles = 1, bzeroes = 1))
if (!inherits(model, 'try-error')) {
   print(model)
   print(min(abs(poles(model$sys))) > 1)
   print(min(abs(zeroes(model$sys))) > 1)
   print(pseries2nu(pseries(model$sys, lag.max = 10))) # Kronecker indices
}
#> state space model [2,2] with s = 3 states
#>       s[1]  s[2]  s[3]  u[1]  u[2]
#> s[1] -0.11 -0.02  0.74  0.17  0.84
#> s[2]  0.15 -0.45 -0.40 -0.40 -0.68
#> s[3] -0.25 -0.28 -0.07  0.05 -0.34
#> x[1]  0.51 -0.53  1.11  1.00  0.00
#> x[2] -0.07 -0.01  0.49  0.00  1.00
#> Left square root of noise covariance Sigma:
#>       u[1] u[2]
#> u[1]  1.83 0.00
#> u[2] -0.37 2.16
#> [1] TRUE
#> [1] TRUE
#> [1] 2 1

### random state space model with three outputs and 2 inputs in echelon canonical form
### D is lower triangular (with ones on the diagonal)
### the model is required to stable (the transfer function has no poles within the unit circle)
model = try(test_stspmod(dim = c(3, 2), nu = c(2,3,0),
                         D = matrix(c(1,NA,NA,0,1,NA), nrow = 3, ncol = 2),
                         digits = 2, bpoles = 1))

if (!inherits(model, 'try-error')) {
   print(model)
   print(min(abs(poles(model$sys))) > 1)
   print(pseries2nu(pseries(model$sys, lag.max = 10))) # Kronecker indices
}
#> state space model [3,2] with s = 5 states
#>       s[1] s[2]  s[3]  s[4]  s[5]  u[1]  u[2]
#> s[1]  0.00 0.00  1.00  0.00  0.00  0.18  0.47
#> s[2]  0.00 0.00  0.00  1.00  0.00  0.11 -0.24
#> s[3] -0.87 0.18  0.56  0.50  0.00  0.70 -0.34
#> s[4]  0.00 0.00  0.00  0.00  1.00  0.05 -0.69
#> s[5] -0.18 0.39 -0.31 -0.22 -0.73  0.83  0.49
#> x[1]  1.00 0.00  0.00  0.00  0.00  1.00  0.00
#> x[2]  0.00 1.00  0.00  0.00  0.00 -0.14  1.00
#> x[3]  0.00 0.06  0.00  0.00  0.00 -0.08 -0.73
#> Left square root of noise covariance Sigma:
#>      u[1] u[2]
#> u[1] 2.82 0.00
#> u[2] 0.82 0.39
#> [1] TRUE
#> [1] 2 3 0
```

# Poles and Zeroes

Compute the poles and zeroes of VARMA and Statespace models. Note that
these models describe the corresponding processes as \$\$x_t = k(B)
u_t\$\$ where \\(u_t)\\ is a white noise process and \\k(B)\\ is a
rational filter (\\B\\ denotes the lag- or backward shift operator). The
poles and zeroes are the poles and zeroes of the rational transfer
function \\k(z)\\ of this filter.

## Usage

``` r
# S3 method for class 'armamod'
zeroes(x, tol = sqrt(.Machine$double.eps), print_message = TRUE, ...)

# S3 method for class 'armamod'
poles(x, tol = sqrt(.Machine$double.eps), print_message = TRUE, ...)

# S3 method for class 'rmfdmod'
zeroes(x, tol = sqrt(.Machine$double.eps), print_message = TRUE, ...)

# S3 method for class 'rmfdmod'
poles(x, tol = sqrt(.Machine$double.eps), print_message = TRUE, ...)

# S3 method for class 'stspmod'
zeroes(x, tol = sqrt(.Machine$double.eps), print_message = TRUE, ...)

# S3 method for class 'stspmod'
poles(x, tol = sqrt(.Machine$double.eps), print_message = TRUE, ...)
```

## Arguments

- x:

  an object which represents a VARMA, RMFD or statespace model (i.e. an
  [`armamod()`](https://bfunovits.github.io/RLDM/reference/armamod.md),
  [`rmfdmod()`](https://bfunovits.github.io/RLDM/reference/rmfdmod.md)
  or
  [`stspmod()`](https://bfunovits.github.io/RLDM/reference/stspmod.md)
  object).

- tol:

  Double. Default set to `sqrt(.Machine$double.eps)`. Required to decide
  on when a root is to be considered "at infinity".

- print_message:

  Boolean. Default set to TRUE. Prints a message if roots "at infinity "
  are discarded.

- ...:

  not used.

## Value

Vector of poles, respectively zeroes.

## See also

For more details we refer to the discussion about the computation of
poles and zeros of rational matrices in the companion package, see
rationalmatrices::poles and zeroes.

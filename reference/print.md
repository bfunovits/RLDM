# Print Methods

Print Methods

## Usage

``` r
# S3 method for class 'armamod'
print(
  x,
  digits = NULL,
  format = c("i|jz", "i|zj", "iz|j", "zi|j", "i|j|z", "character"),
  ...
)

# S3 method for class 'rmfdmod'
print(
  x,
  digits = NULL,
  format = c("i|jz", "i|zj", "iz|j", "zi|j", "i|j|z", "character"),
  ...
)

# S3 method for class 'stspmod'
print(x, digits = NULL, ...)

# S3 method for class 'impresp'
print(
  x,
  digits = NULL,
  format = c("i|jz", "i|zj", "iz|j", "zi|j", "i|j|z"),
  ...
)

# S3 method for class 'autocov'
print(
  x,
  digits = NULL,
  format = c("i|jz", "i|zj", "iz|j", "zi|j", "i|j|z"),
  ...
)

# S3 method for class 'fevardec'
print(
  x,
  digits = NULL,
  format = c("i|jz", "i|zj", "iz|j", "zi|j", "i|j|z"),
  ...
)

# S3 method for class 'freqresp'
print(
  x,
  digits = NULL,
  format = c("i|jz", "i|zj", "iz|j", "zi|j", "i|j|z"),
  ...
)

# S3 method for class 'spectrald'
print(
  x,
  digits = NULL,
  format = c("i|jz", "i|zj", "iz|j", "zi|j", "i|j|z"),
  ...
)
```

## Arguments

- x:

  RLDM object, i.e. a
  [`armamod()`](https://bfunovits.github.io/RLDM/reference/armamod.md),
  [`rmfdmod()`](https://bfunovits.github.io/RLDM/reference/rmfdmod.md),
  [`stspmod()`](https://bfunovits.github.io/RLDM/reference/stspmod.md),
  [`impresp()`](https://bfunovits.github.io/RLDM/reference/impresp.md),
  [`autocov()`](https://bfunovits.github.io/RLDM/reference/autocov.md),
  [`freqresp()`](https://bfunovits.github.io/RLDM/reference/freqresp.md),
  [`spectrum()`](https://rdrr.io/r/stats/spectrum.html) or
  [`fevardec()`](https://bfunovits.github.io/RLDM/reference/fevardec.md)
  object.

- digits:

  (integer) if non `NULL` then correspondingly rounded numbers are
  printed, see [`round()`](https://rdrr.io/r/base/Round.html).

- format:

  (character string) selects specific output formats. Note that `stsp()`
  and
  [`fevardec()`](https://bfunovits.github.io/RLDM/reference/fevardec.md)
  objects have no format option. The option `'character'` is only
  implemented for (V)ARMA models.

- ...:

  Further parameters are ignored.

## Value

`invisible(x)`

## Examples

``` r
# for VARMA models six different print formats are implemented ###################
m = armamod(test_lmfd(dim = c(2,2), degrees = c(1,1)), sigma_L = diag(2))
print(m, digits = 2, format = "i|jz")
#> ARMA model [2,2] with orders p = 1 and q = 1
#> AR polynomial a(z):
#>      z^0 [,1]  [,2] z^1 [,1]  [,2]
#> [1,]        1     0    -0.06  0.44
#> [2,]        0     1    -0.41  0.74
#> MA polynomial b(z):
#>      z^0 [,1]  [,2] z^1 [,1]  [,2]
#> [1,]     0.66  0.58     0.33  1.31
#> [2,]     0.08  0.85     1.67  0.67
#> Left square root of noise covariance Sigma:
#>      u[1] u[2]
#> u[1]    1    0
#> u[2]    0    1
```

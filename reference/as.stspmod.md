# Coerce to State Space Model

The function `as.stsp.pseries()` calls `pseries2stsp()` with default
parameters. Of course the `pseries()` object must contain sufficiently
many lags. NOT YET implemented

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

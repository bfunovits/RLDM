# Discrete Time Fourier Transform

Compute the Discrete Time Fourier Transform for data stored in a
3-dimensional array.

## Usage

``` r
dft_3D(a, n.f = dim(a)[3])
```

## Arguments

- a:

  \\(m,n,k)\\ dimensional (numeric) array

- n.f:

  (integer) number of frequencies

## Value

A
[`rationalmatrices::zvalues()`](https://bfunovits.github.io/rationalmatrices/reference/zvalues.html)
object.

## Examples

``` r
# Basic example
a <- array(1:12, dim = c(2, 2, 3))
result <- dft_3D(a)
result
#> ( 2 x 2 ) frequency response
#>      z[1] [,1]  [,2]    z[2] [,1]         [,2]    z[3] [,1]         [,2]
#> [1,]     15+0i 21+0i -6+3.464102i -6+3.464102i -6-3.464102i -6-3.464102i
#> [2,]     18+0i 24+0i -6+3.464102i -6+3.464102i -6-3.464102i -6-3.464102i
```

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

A `rationalmatrices::zvalues()` object.

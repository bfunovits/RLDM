# Forecast Error Variance Decomposition

Computes the Forecast Errors Variance Decomposition from a given
(orthogonalized) impulse response function.

## Usage

``` r
fevardec(obj, h.max = NULL, H = NULL)
```

## Arguments

- obj:

  [`impresp()`](https://bfunovits.github.io/RLDM/reference/impresp.md)
  object which represents the (orthogonalized) impulse response
  function.

- h.max:

  maximum forecast horizon. The default is one plus the number of lags
  of the `impresp` object.

- H:

  An (n x n) (non singular) matrix which renders the impulse response to
  a *orthogonalized* impulse response. The noise \\u_t\\ is transformed
  to \\H^{-1}u_t\\ and the impulse response coefficients (\\k_j
  \rightarrow k_j H\\) and the (left) square root of the noise
  covariance matrix (\\L \rightarrow H^{-1}L\\) are transformed
  correspondingly.  
  The default case `H=NULL` corresponds to the identity matrix (i.e. no
  transformation).  
  For `H='chol'`, the transformation matrix `H = t(chol(Sigma))` is
  obtained from the Choleski decomposition of the noise covariance
  \\\Sigma\\. For `H='eigen'` the symmetric square root of \\\Sigma\\
  (obtained from the the eigenvalue decomposition of \\\Sigma\\) is
  used. For `H='sigma_L'` the left square root of the noise covariance,
  which is stored in the object `obj`, is used. Other orthogonalization
  schemes may be obtained by setting \\H\\ to a suitable square root of
  \\\Sigma\\.  
  The procedure checks whether the transformation yields an
  orthogonalized impulse response. If not, an error is thrown.

## Value

`fevardec` object, i.e. a list with components

- vd:

  n-by-n-by-h.max array which contains the forecast error variance
  decomposition: `vd[i,j,h]` is the percentage of the variance of the
  h-step ahead forecast error of the i-th component due to the j-th
  orthogonalized shock.

- v:

  n-by-h.max matrix which contains the forecast error variances:
  `v[i,h]` is the variance of the h-step ahead forecast error for the
  i-th component.

- names:

  (m)-dimensional character vector

- label:

  character string or NULL

## Examples

``` r
model = stspmod(sys = stsp(A = c(0,0.2,1,-0.5), B = c(1,1,1,-1),
                           C = c(1,0,0,1)),
                sigma_L = t(chol(matrix(c(4,2,2,3),nrow=2))),
                names = c('y1','y2'), label = 'test model')
fevardec(impresp(model, lag.max=10), H = 'chol', h.max = 5) %>% print(digits = 2, format = 'iz|j')
#> test model: Forecast error variance decompositon [2,2] maximum horizon = 5
#>        u[1] u[2]
#> h=1 y1 1.00 0.00
#>     y2 0.33 0.67
#> h=2 y1 0.87 0.13
#>     y2 0.33 0.67
#> h=3 y1 0.78 0.22
#>     y2 0.29 0.71
#> h=4 y1 0.74 0.26
#>     y2 0.27 0.73
#> h=5 y1 0.72 0.28
#>     y2 0.26 0.74
fevardec(impresp(model, lag.max=4, H = 'eigen'))            %>% print(digits = 2, format = 'iz|j')
#> test model: Forecast error variance decompositon [2,2] maximum horizon = 5
#>        u[1] u[2]
#> h=1 y1 0.92 0.08
#>     y2 0.11 0.89
#> h=2 y1 0.66 0.34
#>     y2 0.36 0.64
#> h=3 y1 0.65 0.35
#>     y2 0.31 0.69
#> h=4 y1 0.62 0.38
#>     y2 0.30 0.70
#> h=5 y1 0.60 0.40
#>     y2 0.30 0.70
```

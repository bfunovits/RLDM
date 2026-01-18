# Impulse Response Function

Compute the (orthogonalized) impulse response function of a VARMA model
or a state space model. The impulse response coefficients are also
called *Power series parameters* of the system.

## Usage

``` r
impresp(obj, lag.max, H)

# S3 method for class 'armamod'
impresp(obj, lag.max = 12, H = NULL)

# S3 method for class 'rmfdmod'
impresp(obj, lag.max = 12, H = NULL)

# S3 method for class 'stspmod'
impresp(obj, lag.max = 12, H = NULL)

# S3 method for class 'impresp'
impresp(obj, lag.max = NULL, H = NULL)
```

## Arguments

- obj:

  [`armamod()`](https://bfunovits.github.io/RLDM/reference/armamod.md),
  [`stspmod()`](https://bfunovits.github.io/RLDM/reference/stspmod.md)
  object or `impresp()` object. The last case may be used to transform
  the impulse response function to a different orthogonalization scheme.

- lag.max:

  Maximum lag of the impulse response coefficients. This parameter is
  ignored in the case that `obj` is an impresp object.

- H:

  An (n x n) (non singular) matrix which specifies a transformation of
  the noise. The noise \\u_t\\ is transformed to \\H^{-1}u_t\\ and the
  impulse response coefficients (\\k_j \rightarrow k_j H\\) and the
  (left) square root of the noise covariance matrix (\\L \rightarrow
  H^{-1}L\\) are transformed correspondingly.  
  The default case `H=NULL` corresponds to the identity matrix (i.e. no
  transformation).  
  For `H='chol'`, the transformation matrix `H = t(chol(Sigma))` is
  determined from the Choleski decomposition of the noise covariance
  \\\Sigma\\. For `H='eigen'` the symmetric square root of \\\Sigma\\
  (obtained from the the eigenvalue decomposition of \\\Sigma\\) is
  used. For `H='sigma_L'` the left square root of the noise covariance,
  which is stored in the object `obj`, is used. In these cases one
  obtains an *orthogonalized* impulse response function. Other
  orthogonalization schemes may be obtained by setting \\H\\ to a
  suitable square root of \\\Sigma\\.

## Value

`impresp` object, i.e. a list with components

- irf:

  `pseries` object.

- sigma_L:

  (n,n)-dimensional matrix which contains left square of noise
  covariance matrix.

- names:

  (n)-dimensional character vector or NULL

- label:

  character string or NULL

## Details

The impulse response coefficients \\(k_j \\\|\\ j \geq 0)\\ define the
map between the noise and the output process. If the model is stable
then the stationary solution of the ARMA system, respectively state
space system, is given by \$\$ y_t = \sum\_{j \geq 0} k_j u\_{t-j}. \$\$
For a state space system the impulse response coefficients are \$\$k_0 =
D \mbox{ and }\$\$ \$\$k_j = CA^{j-1}B \mbox{ for }j \>0.\$\$ For an
ARMA model the coefficients are (recursively) computed by a comparison
of coefficients in the equation \$\$ (a_0 + a_1 z + \cdots + a_p
z^p)(k_0 + k_1 z + k_2 z^2 + \cdots ) = b_0 + b_1 z + \cdots + b_q z^q
\$\$

The S3 methods `impresp.*` compute the coefficients \\k_j\\ for \\j =
0,\cdots,N\\ and store the result, together with the left square root
(`sigma_L`) of the noise covariance \\\Sigma\\, in a **impresp** object.
`impresp` objects contain the complete information about the underlying
model, provided that the maximum lag \\N\\ is large enough. This means
that one may reconstruct the underlying model from an impulse response
object.

## References

Lütkepohl H (2005). *New Introduction to Multiple Time Series Analysis*.
Springer Berlin.

## Examples

``` r
# IRF from state space model ################################################
model = stspmod(stsp(A = c(0,0.2,1,-0.5), B = c(1,1,1,-1),
                     C = c(1,0,0,1)),
                sigma_L = matrix(c(4, 1, 1, 3), 2, 2),
                names = c('y1','y2'), label = 'test model')

# IRF
irf = impresp(model, lag.max=10)
irf
#> test model: Impulse response [2,2] with 10 lags
#>      lag=0 [,1]  [,2] lag=1 [,1]  [,2] lag=2 [,1]  [,2] lag=3 [,1]  [,2]
#> [1,]          1     0          1     1        1.0  -1.0      -0.30  0.70
#> [2,]          0     1          1    -1       -0.3   0.7       0.35 -0.55
#>      lag=4 [,1]   [,2] lag=5 [,1]    [,2] lag=6 [,1]     [,2] lag=7 [,1]
#> [1,]      0.350 -0.550    -0.2350  0.4150    0.18750 -0.31750  -0.140750
#> [2,]     -0.235  0.415     0.1875 -0.3175   -0.14075  0.24175   0.107875
#>           [,2] lag=8 [,1]       [,2]  lag=9 [,1]       [,2] lag=10 [,1]
#> [1,]  0.241750  0.1078750 -0.1843750 -0.08208750  0.1405375  0.06261875
#> [2,] -0.184375 -0.0820875  0.1405375  0.06261875 -0.1071438 -0.04772688
#>             [,2]
#> [1,] -0.10714375
#> [2,]  0.08167938

# Orthogonalized IRF: Cholesky
irf_chol = impresp(model, lag.max = 10, H = 'chol')
irf_chol
#> test model: Impulse response [2,2] with 10 lags
#>      lag=0 [,1]     [,2] lag=1 [,1]      [,2]  lag=2 [,1]      [,2]  lag=3 [,1]
#> [1,]   4.123106 0.000000   5.820855  2.667892  2.42535625 -2.667892 -0.04850713
#> [2,]   1.697749 2.667892   2.425356 -2.667892 -0.04850713  1.867524  0.50932481
#>           [,2] lag=4 [,1]      [,2] lag=5 [,1]       [,2] lag=6 [,1]       [,2]
#> [1,]  1.867524  0.5093248 -1.467341 -0.2643638  1.1071751  0.2340469 -0.8470557
#> [2,] -1.467341 -0.2643638  1.107175  0.2340469 -0.8470557 -0.1698962  0.6449629
#>      lag=7 [,1]       [,2]  lag=8 [,1]       [,2]  lag=9 [,1]       [,2]
#> [1,] -0.1698962  0.6449629  0.13175748 -0.4918926 -0.09985798  0.3749389
#> [2,]  0.1317575 -0.4918926 -0.09985798  0.3749389  0.07628049 -0.2858479
#>      lag=10 [,1]       [,2]
#> [1,]  0.07628049 -0.2858479
#> [2,] -0.05811184  0.2179117
print(irf_chol$sigma_L) # Sigma is (approximately equal to) the identity matrix
#>            [,1]      [,2]
#> [1,]  0.9701425 0.2425356
#> [2,] -0.2425356 0.9701425


# IRF from VARMA model ################################################
model = armamod(sys = test_lmfd(dim = c(2,2), degrees = c(2,1)))

irf = impresp(model)
print(irf, digits = 2, format = 'iz|j')
#> Orthogonalized impulse response [2,2] with 12 lags
#>                [,1]    [,2]
#>  lag=0 [1,]   -1.83    1.96
#>        [2,]   -1.14    0.94
#>  lag=1 [1,]    1.93   -4.66
#>        [2,]    1.67   -1.84
#>  lag=2 [1,]   -0.48    3.34
#>        [2,]   -4.90    5.54
#>  lag=3 [1,]   -0.10   -1.81
#>        [2,]    5.87  -12.12
#>  lag=4 [1,]    4.20   -1.19
#>        [2,]   -4.39   12.25
#>  lag=5 [1,]   -9.87   12.50
#>        [2,]    2.67  -11.50
#>  lag=6 [1,]   12.72  -22.70
#>        [2,]    7.61    3.57
#>  lag=7 [1,]  -16.17   31.52
#>        [2,]  -22.76   23.34
#>  lag=8 [1,]   11.59  -38.54
#>        [2,]   36.94  -55.10
#>  lag=9 [1,]    7.83   20.51
#>        [2,]  -53.07   92.97
#> lag=10 [1,]  -36.42   24.15
#>        [2,]   50.52 -126.93
#> lag=11 [1,]   80.61  -96.32
#>        [2,]  -12.67  107.13
#> lag=12 [1,] -127.68  206.26
#>        [2,]  -62.30  -14.68
```

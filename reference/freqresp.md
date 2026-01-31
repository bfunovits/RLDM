# Frequency Response Function

Compute the *frequency response function* (also called *transfer
function*) associated to a VARMA or state space model.

## Usage

``` r
freqresp(obj, n.f, ...)

# S3 method for class 'armamod'
freqresp(obj, n.f = 128, ...)

# S3 method for class 'stspmod'
freqresp(obj, n.f = 128, ...)

# S3 method for class 'impresp'
freqresp(obj, n.f = 128, ...)
```

## Arguments

- obj:

  [`armamod()`](https://bfunovits.github.io/RLDM/reference/armamod.md),
  [`stspmod()`](https://bfunovits.github.io/RLDM/reference/stspmod.md)
  or
  [`impresp()`](https://bfunovits.github.io/RLDM/reference/impresp.md)
  object. Note that for an impulse response object the result is only an
  approximation of the "true" frequency response due to the finite
  number of coefficients.

- n.f:

  number of frequencies.

- ...:

  not used.

## Value

`freqresp` object, i.e. a list with slots

- frr:

  [`rationalmatrices::zvalues()`](https://bfunovits.github.io/rationalmatrices/reference/zvalues.html)
  object.

- sigma_L:

  (n,n)-dimensional matrix which contains a left square root of the
  noise covariance matrix \\\Sigma\\.

- names:

  (m)-dimensional character vector or NULL. This optional slot stores
  the names for the components of the time series/process.

- label:

  character string or NULL.

## Details

The frequency response function (or transfer function) associated to an
ARMA or state space model is \$\$ K(\lambda) = \sum\_{j=0}^{\infty} k_j
e^{-i\lambda j} \$\$ where \\(k_j \\\|\\ j\geq 0)\\ is the *impulse
response* of the model. See also
[`impresp()`](https://bfunovits.github.io/RLDM/reference/impresp.md).

For an ARMA model the frequency response is equal to \$\$ K(\lambda) =
(a_0 + a_1 e^{-i\lambda} + \cdots + a_p e^{-i\lambda p})^{-1} (b_0 + b_1
e^{-i\lambda} + \cdots + b_q e^{-i\lambda q}) \$\$ and for a state space
model we have \$\$ K(\lambda) = C(e^{i\lambda}I_s - A)^{-1}B+D \$\$ Note
that \\K()\\ is the discrete-time Fourier transform (DTFT) of the
impulse response. If the impulse response is absolutely summable then
the coefficents \\k_j\\ may be reconstructed from the frequency response
via the inverse DTFT \$\$ k_j = \frac{1}{2\pi} \int\_{-\pi}^{\pi}
K(\lambda) e^{i\lambda j} d\lambda \$\$

The S3 methods `freqresp.*` evaluate the function on a grid of angular
frequencies \\\lambda_j = 2\pi j/N\\, \\j=0,\ldots,N-1\\ and store the
result (together with `sigma_L`) in a **freqresp** object.

## Examples

``` r
set.seed(3451) # set seed in order to get reproducible results

### generate random bivariate ARMA(1,1) model
# "bpoles = 1.1" implies that the poles have moduli larger than 1.1
# and therefore the impulse response coefficients decay with a rate (1.1)^k
arma_model = test_armamod(dim = c(2,2), degrees = c(1,1), bpoles = 1.1)
# frequency response
frr = freqresp(arma_model)
# compute the frequency response via the impulse response
irf = impresp(arma_model, lag.max = 100)
frr1 = freqresp(irf)
# since the impulse response quickly decays
# the "truncated" frequency response should be close to the true frequency response
all.equal(frr, frr1)
#> [1] TRUE
# create an equivalent state space model
stsp_model = as.stspmod(arma_model)
# of course the state space model has the same frequency response
# as the original ARMA model
frr1 = freqresp(stsp_model)
all.equal(frr, frr1)
#> [1] TRUE

# we can also reconstruct the impulse response from the
# frequency response, provided the frequency grid is "fine enough"
n.f = 2^6
frr = freqresp(arma_model, n.f = n.f)
# compute the impulse response via the inverse DTFT
K = unclass(frr$frr)
k1 = Re(apply(K, MARGIN = c(1,2), FUN = fft, inverse = TRUE)) / n.f
k1 = aperm(k1, c(2,3,1))
# impulse response
irf = impresp(arma_model, lag.max = n.f-1)
k = unclass(irf$irf)
# compare
all.equal(k, k1)
#> [1] TRUE

set.seed(NULL) # reset seed
# Create a simple model
model <- stspmod(sys = test_stsp(dim = c(2,2), s = 2), sigma_L = diag(2))

# Compute impresp
result <- impresp(model)
result
#> Orthogonalized impulse response [2,2] with 12 lags
#>      lag=0 [,1]  [,2] lag=1 [,1]      [,2] lag=2 [,1]       [,2] lag=3 [,1]
#> [1,]          1     0  1.0440746 3.3519438   1.055711 -0.4766323 -1.3543143
#> [2,]          0     1  0.5727021 0.2411002  -0.421677 -1.0533874 -0.2382033
#>            [,2] lag=4 [,1]     [,2] lag=5 [,1]       [,2] lag=6 [,1]       [,2]
#> [1,] -2.3984028 -0.1481255 1.584709  1.1408599  1.0889520 -0.4587138 -1.7966832
#> [2,]  0.3414096  0.4520817 0.657093 -0.0404582 -0.6001339 -0.3354913 -0.2146448
#>      lag=7 [,1]       [,2] lag=8 [,1]       [,2] lag=9 [,1]       [,2]
#> [1,] -0.6667991 0.04878797  0.6973745  1.3898064  0.1732484 -0.7392717
#> [2,]  0.2010347 0.58068867  0.1627305 -0.1238593 -0.2403245 -0.3946768
#>       lag=10 [,1]       [,2] lag=11 [,1]      [,2] lag=12 [,1]        [,2]
#> [1,] -0.636363972 -0.7212958   0.1845244 0.9457212  0.40791355  0.09091154
#> [2,] -0.006913935  0.2965374   0.1926782 0.1611623 -0.09172261 -0.31471666
```

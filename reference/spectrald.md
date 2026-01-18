# Spectral Density

Compute the spectral density of an ARMA process or a process defined by
a state space model.

## Usage

``` r
spectrald(obj, n.f, ...)

# S3 method for class 'armamod'
spectrald(obj, n.f = 128, ...)

# S3 method for class 'stspmod'
spectrald(obj, n.f = 128, ...)

# S3 method for class 'autocov'
spectrald(obj, n.f = 128, ...)

# S3 method for class 'impresp'
spectrald(obj, n.f = 128, ...)

# Default S3 method
spectrald(obj, n.f = NULL, demean = TRUE, ...)
```

## Arguments

- obj:

  [`armamod()`](https://bfunovits.github.io/RLDM/reference/armamod.md),
  [`stspmod()`](https://bfunovits.github.io/RLDM/reference/stspmod.md),
  [`autocov()`](https://bfunovits.github.io/RLDM/reference/autocov.md),
  [`impresp()`](https://bfunovits.github.io/RLDM/reference/impresp.md)
  object or a "time series" object, i.e. an object which may be coerced
  to a data matrix by `y = as.matrix(obj)`.

- n.f:

  number of frequencies.

- ...:

  not used.

- demean:

  (logical) should the data be demeaned before computing the
  periodogram?

## Value

`freqresp` object, i.e. a list with slots

- spd:

  `rationalmatrices::zvalues()` object.

- names:

  (m)-dimensional character vector or NULL. This optional slot stores
  the names for the components of the time series/process.

- label:

  character string or NULL.

- n.obs:

  (optional) integer or NULL.

## Details

The spectral density of a stationary process with an absolutely summable
autocovariance function \\(\gamma_j)\\ is given by \$\$ \Gamma(\lambda)
= \frac{1}{2\pi}\sum\_{j=-\infty}^{\infty} \gamma_j e^{-i\lambda j}.
\$\$

For an ARMA process, or process defined by a state space model the
spectral density is equal to \$\$ \Gamma(\lambda) = \frac{1}{2\pi}
K(\lambda) \Sigma K^\*(\lambda) \$\$ where \\\Sigma\\ is the noise
covariance, \\K()\\ is the frequency response of the model and
\\K^\*(\lambda)\\ is the Hermitean transpose of \\K(\lambda)\\. See also
[`autocov()`](https://bfunovits.github.io/RLDM/reference/autocov.md) and
[`freqresp()`](https://bfunovits.github.io/RLDM/reference/freqresp.md).

Note that \\\Gamma()\\ is (up to a factor \\2\pi\\) the discrete-time
Fourier transform (DTFT) of the autocovariance function and therefore
the ACF \\\gamma_j\\ may be reconstructed from the spectral density via
the inverse DTFT \$\$ \gamma_j = \int\_{-\pi}^{\pi} \Gamma(\lambda)
e^{i\lambda j} d\lambda \$\$

The S3 methods `spectrald.*` evaluate the spectral density function on a
grid of angular frequencies \\\lambda_j = 2\pi j/N\\, \\j=0,\ldots,N-1\\
and store the result in a **spectrald** object.

There are several possible ways to specify the process. One may provide
the ARMA (`armamod`) respectively the state space model (`stspmod`), the
autocovariance function (`autocov`) or the impulse response function
(`impresp`) which maps the noise to the outputs. Note however, that if
we have only given an `autocov` or `impresp` object then the computed
spectral density is only an approximation of the true spectral density
since only a finite number of covariances respectively impulse response
coefficients are given. The type of the autocovariance function
("covariances", "correlations" or "partial correlations") is irrelevenat
since the procedure alwayas uses the slot "`gamma`" which contains the
covariances.

The default method `spectrald.default` assumes that `obj` is a "time
series" object and tries to coerce this object to a data matrix via
`y = as.matrix(obj)`. In this case the procedure computes the
*periodogram* which is a simple estimate of the spectral density. The
periodgram may also be computed by the call
`spectrald(autocov(obj, max.lag = n.obs-1))`, i.e. by first computing
the sample auto covariance function and then computing the corresponding
spectral density.

Note that we use a different scaling than the
`stats::[spectrum][stats::spectrum]` routine.

## Examples

``` r
#' ### generate random 3-dimensional ARMA(1,1) model
# "bpoles = 1.1" implies that the poles have moduli larger than 1.1
# and therefore the impulse response coefficients decay with a rate (1.1)^k
arma_model = test_armamod(dim = c(3,3), degrees = c(1,1), bpoles = 1.1)

# spectral density
spd = spectrald(arma_model)

# compute the spectral density via the impulse response
spd1 = spectrald(impresp(arma_model, lag.max = 100))

# since the impulse response quickly decays
# the "truncated" spectral density should be close to the true one
all.equal(spd, spd1)
#> [1] "Component “spd”: Mean relative Mod difference: 3.54988e-05"

# compute the spectral density via the autocovariance function
spd1 = spectrald(autocov(arma_model, lag.max = 100))

# since the ACF quickly decays
# the "truncated" spectral density should be close to the true one
all.equal(spd, spd1)
#> [1] "Component “spd”: Mean relative Mod difference: 3.838035e-05"

# create an equivalent state space model
stsp_model = as.stspmod(arma_model)

# of course the state space model gives the same spectrum
# as the original ARMA model
spd1 = spectrald(stsp_model)
all.equal(spd, spd1)
#> [1] TRUE
```

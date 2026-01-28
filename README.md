# Rational Linear Dynamic Models (RLDM)

This `RLDM` (Rational Linear Dynamic Models) R package provides models for stationary processes with a rational spectral density and methods for their estimation.
We will refer to them as **rational models**.
It builds heavily on its sister R package `rationalmatrices`, see https://bfunovits.github.io/rationalmatrices/.

## Installation

You can install the latest version of the code using the `remotes` R package.

```
remotes::install_github("bfunovits/RLDM")
```

## Content

The package provides the following sets of functions whose documentation can be found in the reference page https://bfunovits.github.io/RLDM/reference of the website https://bfunovits.github.io/RLDM/ (created with https://pkgdown.r-lib.org/):

* Classes for the construction of rational models (consisting of an input covariance matrix and a rational matrix function from the `rationalmatrices` package):
    * VARMA models `armamod()`
    * State space models `stspmod()`
    * Right matrix fraction description (RMFD) models `rmfdmod()` (which is experimental)
    
* Templates for filling the linear parameters with deep parameters through an affine mapping. 
  Consists of 
    * a matrix $H$ where the number of rows is the number of linear parameters in a given model and the number of columns is the number of deep parameters in a given model 
    * a column vector $h$ of appropriate dimension

See `help("model structures")` and `help("local model structures")` for more details.    

* Generic functions to create objects which are derived from these rational models
    * The autocovariance sequence, see `autocov()`
    * Spectral density, see `spectrald()`
    * The transfer function/impulse response function (IRF), see `impresp()`
        * Forecast error variance decomposition, see `fevardec()`, for a given IRF
    * Frequency response (the transfer function evaluated on the unit circle), see `freqresp()`
    
* Several other generic functions which extend R's generic functions
    * `plot()`, `print()`, `str()`, `predict()`
    
* Some helpers for estimation methods: `solve_de()`, `solve_inverse_de()`, and more

* Moment estimation methods for 
    * AR models, see e.g. `est_ar()`
    * ARMA models, see the Hannan-Rissannen-Kavalieris algorithm in `est_arma_hrk3()`
    * state space models, see e.g. `est_stsp_cca()`
    
* Likelihood estimation methods
    * `ll()`
    * `ll_theta()` and `ll_FUN()` for the estimation of the deep parameters of a rational model
    * `ll_kf()`
    
* Some more tooling like
    * simulation in `sim()`
    * model comparison in `KL_divergence()`, `pm_test()`, `compare_estimates()`

## Getting Started

RLDM includes comprehensive documentation organized by user level:

1. **New to RLDM?** Start with `vignette("0_getting_started")`
   - Practical introduction with simple examples
   - AR models and multivariate VAR/VARMA systems
   - ~10-15 minute read with executable code

2. **Want a complete workflow example?** See `vignette("1_case_study")`
   - Blanchard-Quah economic data analysis
   - Comparing AR, state space, and ARMA models
   - Model selection, diagnostics, and forecasting

3. **Need technical details?** Consult `vignette("2_technical_reference")`
   - Complete class and method documentation
   - Mathematical foundations
   - Method selection guidance

## Installation Notes

The package depends on `rationalmatrices`. Install both:

```r
remotes::install_github("bfunovits/rationalmatrices")
remotes::install_github("bfunovits/RLDM")
```

**Vignette Building:** To build vignettes locally, you need `pandoc-citeproc` installed on your system:
- Ubuntu/Debian: `sudo apt-get install pandoc-citeproc`
- macOS: `brew install pandoc-citeproc`
- Windows: Included in pandoc installer

Pre-built vignettes are available on the [package website](https://bfunovits.github.io/RLDM/).

Then load and explore:

```r
library(RLDM)
browseVignettes("RLDM")  # View all vignettes (if built locally)
vignette("0_getting_started")  # Access specific vignette
```

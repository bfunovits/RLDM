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

## Usage

See the case study `vignette("d_casestudy2")` for a detailed example of how to use the package.

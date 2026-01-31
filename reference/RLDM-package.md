# RLDM: Rational Linear Dynamic Models

This package provides tools for stationary processes with rational
spectral density. It implements VARMA and state space models with
methods for estimation, simulation, prediction, and model comparison.

## Details

The package uses Rcpp/RcppArmadillo for performance-critical
computations including Kalman filtering and recursive least squares. It
depends on the sister package `rationalmatrices` for rational matrix
operations.

## Package Organization

R source files use a numeric prefix system organized by
purpose/workflow:

- **01\_**: Model representations and classes (`armamod`, `stspmod`,
  `rmfdmod`)

- **02\_**: Parameter templates (`tmpl_*` functions, `fill_template`,
  `extract_theta`)

- **03\_**: Derived properties (autocovariance, frequency response,
  impulse response, spectral density, poles)

- **04\_**: Time series operations (`solve_de`, `sim`,
  prediction/forecasting)

- **05\_**: Estimation methods (AR, ARMA, subspace, likelihood,
  recursive least squares)

- **06\_**: Visualization (plot methods for properties and predictions)

- **07\_**: Model comparison metrics and diagnostics

- **08\_**: Utilities (`print`, `str`, data documentation, package
  metadata)

## Getting Started

See vignettes for different user levels:

- [`vignette("0_getting_started")`](https://bfunovits.github.io/RLDM/articles/0_getting_started.md)
  for beginner-friendly introduction

- [`vignette("1_case_study")`](https://bfunovits.github.io/RLDM/articles/1_case_study.md)
  for practical end-to-end workflow

- [`vignette("2_technical_reference")`](https://bfunovits.github.io/RLDM/articles/2_technical_reference.md)
  for technical details and method selection

## Citation

When using RLDM in publications, please cite the relevant theory papers
referenced in the vignettes. Use `citation("RLDM")` for package
citation.

## See also

Useful links:

- <https://bfunovits.github.io/RLDM>

- <https://github.com/bfunovits/RLDM>

- Report bugs at <https://github.com/bfunovits/RLDM/issues>

## Author

Wolfgang Scherrer, Bernd Funovits Maintainer: <bernd.funovits@gmail.com>

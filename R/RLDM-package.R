#' RLDM: Rational Linear Dynamic Models
#'
#' @description
#' This package provides tools for stationary processes with rational spectral density.
#' It implements VARMA and state space models with methods for estimation,
#' simulation, prediction, and model comparison.
#'
#' @details
#' The package uses Rcpp/RcppArmadillo for performance-critical computations
#' including Kalman filtering and recursive least squares. It depends on the
#' sister package `rationalmatrices` for rational matrix operations.
#'
#' @section Package Organization:
#' R source files use a numeric prefix system organized by purpose/workflow:
#' - **01_**: Model representations and classes (`armamod`, `stspmod`, `rmfdmod`)
#' - **02_**: Parameter templates (`tmpl_*` functions, `fill_template`, `extract_theta`)
#' - **03_**: Derived properties (autocovariance, frequency response, impulse response, spectral density, poles)
#' - **04_**: Time series operations (`solve_de`, `sim`, prediction/forecasting)
#' - **05_**: Estimation methods (AR, ARMA, subspace, likelihood, recursive least squares)
#' - **06_**: Visualization (plot methods for properties and predictions)
#' - **07_**: Model comparison metrics and diagnostics
#' - **08_**: Utilities (`print`, `str`, data documentation, package metadata)
#'
#' @section Getting Started:
#' See vignettes for different user levels:
#' - `vignette("0_getting_started")` for beginner-friendly introduction
#' - `vignette("1_case_study")` for practical end-to-end workflow
#' - `vignette("2_technical_reference")` for technical details and method selection
#'
#' @section Citation:
#' When using RLDM in publications, please cite the relevant theory papers
#' referenced in the vignettes. Use `citation("RLDM")` for package citation.
#'
#' @author Wolfgang Scherrer, Bernd Funovits
#' Maintainer: <bernd.funovits@gmail.com>
"_PACKAGE"

#' @useDynLib RLDM, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @importFrom Rdpack reprompt
#' @import rationalmatrices
NULL

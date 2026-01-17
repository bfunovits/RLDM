# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

RLDM (Rational Linear Dynamic Models) is an R package for modeling stationary processes with rational spectral density. It uses Rcpp/RcppArmadillo for performance-critical computations including Kalman filtering and recursive least squares.

**Dependencies:** Requires the sister package `rationalmatrices` (installed from GitHub via `remotes::install_github("bfunovits/rationalmatrices")`).

## Build Commands

```r
# Development workflow
devtools::load_all()           # Load package in dev mode
devtools::document()           # Generate Roxygen docs
devtools::test()               # Run testthat tests
devtools::check()              # Full package check
devtools::build()              # Build source package

# Run a single test file
testthat::test_file("tests/testthat/test-templates.R")

# Rcpp workflow (after modifying C++ code)
Rcpp::compileAttributes()      # Regenerate RcppExports.cpp/R
devtools::load_all()           # Reload with new compiled code
```

## Architecture

### Model Classes (S3)

Three main model classes, all inheriting from `rldm`:
- **`armamod`** - VARMA models (left matrix fraction description)
- **`stspmod`** - State space models
- **`rmfdmod`** - Right matrix fraction description (experimental, many methods not yet implemented)

### Processing Pipeline

```
Model Construction (armamod/stspmod/rmfdmod)
    ↓
Parameter Templates (tmpl_* functions) - maps deep parameters to linear parameters
    ↓
Derived Objects (autocov, freqresp, impresp, spectrald)
    ↓
Estimation (moment-based or likelihood)
    ↓
Model Comparison & Diagnostics
```

### Key File Organization

R files use alphabetical prefixes for logical grouping and load-order control:
- `aa_*` - Core model classes (`armamod`, `stspmod`, `rmfdmod`)
- `ab_*` - Parameter templates and model generation
- `aca_*`, `acb_*`, etc. - Derived object methods:
  - `aca_` - autocov (autocovariance)
  - `acb_` - freqresp (frequency response)
  - `acc_` - impresp (impulse response)
  - `acd_` - spectrald (spectral density)
- `ad*` - Solving/simulation methods (`solve_de`, `sim`)
- `ae*` - Plot/print/predict methods
- `af*` - Estimation methods:
  - `afa_` - AR estimation
  - `afb_` - Likelihood-based estimation
  - `afc_` - Subspace/CCA methods

This alphabetical convention ensures files load in dependency order and groups related functionality together.

### C++ Code (`src/`)

- `kf.cpp` - Kalman filter implementation
- `solve_fwd_bwd_cll.cpp` - Forward-backward solving with conditional log-likelihood
- `rls_core.cpp` - Recursive least squares (exponential forgetting and windowed variants)
- `Makevars` - Links against OpenMP, LAPACK, BLAS

## Key Entry Points

- Model construction: `armamod()`, `stspmod()`, `rmfdmod()`
- Estimation: `est_ML()`, `est_ar()`, `est_arma_hrk3()`, `est_stsp_cca()`
- Analysis: `autocov()`, `spectrald()`, `impresp()`, `freqresp()`
- Solving: `solve_de()`, `solve_inverse_de()`
- Comparison: `KL_divergence()`, `pm_test()`, `compare_estimates()`

## Documentation

- Uses Roxygen2 with markdown support
- References managed via Rdpack
- Main case study: `vignette("d_casestudy2")`
- Online docs: https://bfunovits.github.io/RLDM/

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

## Testing

**Test Patterns:**
- Use `testthat::test_file("tests/testthat/test-*.R")` to run specific test files
- Use `devtools::test(filter = "pattern")` to run tests matching a pattern (e.g., `devtools::test(filter = "pfilter")`)
- Common assertions: `expect_no_error()`, `expect_true(inherits(obj, "class"))`, `expect_named()`, `expect_equal()`
- Set seeds with `set.seed()` for reproducible stochastic tests
- Convergence tests: compare particle filter results to Kalman filter baseline with increasing particle counts
- Validation bounds: Use reasonable thresholds (e.g., RMSE < 0.5, likelihood difference < 10) rather than strict convergence

**Example test structure:**
```r
test_that("function runs without errors", {
  set.seed(123)
  model <- r_model(tmpl_stsp_full(...))
  data <- sim(model, n.obs = 20)
  expect_no_error(result <- pfilter(model, data$y))
  expect_true(inherits(result, "pfilter"))
  expect_named(result, c("filtered_states", "predicted_states", ...))
})
```

## Benchmarking

- Use `microbenchmark` package for timing comparisons (`library(microbenchmark)`)
- Create benchmark scripts that save results to `.RData` files for later analysis
- Compare performance across parameter ranges (particle counts, noise levels, cross-covariance)
- Validate against theoretical expectations (e.g., optimal proposal should approach Kalman filter with increasing particles)

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

### Key File Organization (Numeric Functional System)

R files use a numeric prefix system organized by **purpose/workflow**:

| Prefix | Category | Content |
|--------|----------|---------|
| **01** | Representations | Model classes: `armamod()`, `stspmod()`, `rmfdmod()` |
| **02** | Templates | Parameter templates (`tmpl_*`), `fill_template()`, `extract_theta()` |
| **03** | Properties | Derived properties: autocov, frequency response, impulse response, spectral density, poles |
| **04** | Timeseries | Operations: `solve_de()`, `sim()`, prediction/forecasting |
| **05** | Estimation | All estimation methods: AR (OLS/YW/DLW), ARMA (HRK), subspace, likelihood, RLS |
| **06** | Visualization | Plot methods for properties and predictions |
| **07** | Comparison | Model comparison metrics and diagnostics |
| **08** | Utilities | Support: `print()`, `str()`, data documentation, package metadata |
| **99** | Auto-generated | RcppExports (DO NOT EDIT) |

**Key principle:** Files are organized by what the code **does** (operations), not what object type it operates on (S3 dispatch handles representation differences).

This numeric convention ensures files load in dependency order and groups related functionality together. See `R/README.md` for detailed organization documentation.

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

### Vignettes

RLDM has three main vignettes organized for different user levels:

1. **Getting Started** (`vignette("0_getting_started")`)
   - Beginner-friendly introduction
   - Simple univariate and bivariate examples
   - ~10-15 minute read with executable code
   - Best starting point for new users

2. **Case Study** (`vignette("1_case_study")`)
   - Practical end-to-end workflow
   - Blanchard-Quah economic data analysis
   - Comparing multiple estimation methods
   - Model diagnostics, prediction, and interpretation

3. **Technical Reference** (`vignette("2_technical_reference")`)
   - Comprehensive class and method documentation
   - Mathematical foundations
   - Method selection guidance
   - Reference material for advanced users

**Developer Documentation:**
- `inst/doc/technical_notes.Rmd` - Deep mathematical derivations and algorithm details
  - Durbin-Levinson-Whittle recursions
  - HRK Stage 3 algorithm specifications
  - ARMA ACF computation
  - Not built as vignette but included in package

### Documentation Tools

- Uses Roxygen2 with markdown support
- References managed via Rdpack
- Online docs: https://bfunovits.github.io/RLDM/

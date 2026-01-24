# Technology Stack

**Analysis Date:** 2026-01-24

## Languages

**Primary:**
- R (>= 2.10) - Core statistical computing language for the entire package
- C++11 - Performance-critical computations via Rcpp/RcppArmadillo

**Secondary:**
- R Markdown (.Rmd) - Documentation and vignettes
- Markdown (.md) - Project documentation

## Runtime

**Environment:**
- R runtime environment
- RStudio IDE (`.Rproj` file present)

**Package Manager:**
- R package system with `DESCRIPTION` file
- `remotes` for GitHub dependencies
- CRAN for standard R packages

## Frameworks

**Core:**
- Rcpp  - C++ integration for R
- RcppArmadillo - Linear algebra library for C++
- rationalmatrices - Sister package for rational matrix computations

**Testing:**
- testthat (>= 2.1.0) - Unit testing framework

**Documentation:**
- Roxygen2 (7.3.2) - Documentation generation with markdown support
- pkgdown - Package website generation
- Rdpack - Reference management for documentation

**Build/Dev:**
- devtools - Package development tools
- knitr - Dynamic report generation for vignettes
- rmarkdown - R Markdown processing

## Key Dependencies

**Critical:**
- rationalmatrices - Core dependency for rational matrix operations (GitHub: bfunovits/rationalmatrices)
- RcppArmadillo - Linear algebra in C++ code
- dplyr, purrr, tibble - Data manipulation and functional programming
- MASS - Statistical functions
- QZ - QZ decomposition for generalized eigenvalues

**Infrastructure:**
- LAPACK, BLAS - Linear algebra libraries (linked via `Makevars`)
- OpenMP - Parallel computing support (enabled in `Makevars`)

## Configuration

**Environment:**
- Configured via R package `DESCRIPTION` file at `/media/bernd/nvme/r_projects/acad_RLDM/DESCRIPTION`
- No `.env` files detected
- Build configuration in `src/Makevars` and `src/Makevars.win`

**Build:**
- `Makevars` - Unix build configuration with OpenMP and LAPACK/BLAS linking
- `Makevars.win` - Windows build configuration
- `.Rbuildignore` - Files excluded from package build

## Platform Requirements

**Development:**
- R (>= 2.10)
- Rcpp and RcppArmadillo development headers
- C++11 compatible compiler
- LAPACK and BLAS libraries
- OpenMP support for parallel computations

**Production:**
- R runtime environment
- Compatible with standard R package installation via `remotes::install_github()`

---

*Stack analysis: 2026-01-24*
# RLDM Package Cleanup & Streamlining

## What This Is

RLDM (Rational Linear Dynamic Models) is an R package for modeling stationary processes with rational spectral density, using Rcpp/RcppArmadillo for performance-critical computations. This project focuses on cleaning up the development workflow, ensuring documentation completeness, verifying build stability, and improving code organization for maintainability.

## Core Value

A clean, well-documented package that builds without errors and has a streamlined development workflow for the maintainer.

## Requirements

### Validated

(None yet — ship to validate)

### Active

- [ ] Repository cleanup: Move non-essential files from root to appropriate directories (benchmarks to inst/, logs to logs/)
- [ ] Documentation completeness: 100% Roxygen coverage for all exported functions with examples
- [ ] Build verification: `devtools::check()` passes with 0 errors and ≤1 warning
- [ ] Code organization: Review and clean up C++ code in src/, ensure R file numeric prefix consistency, organize S3 methods
- [ ] Vignette verification: Ensure vignettes build correctly during `devtools::check()`

### Out of Scope

- Test coverage improvements — focus on existing tests passing, not coverage metrics
- Data documentation — existing documentation is sufficient
- Performance optimization — maintain current performance levels
- New features — this is a maintenance/cleanup project only

## Context

Existing R package with numeric-prefixed R file organization (01_ representations, 02_ templates, etc.), C++ extensions for Kalman filtering and recursive least squares, and three main vignettes. Package has sister dependency `rationalmatrices`. Current root directory has clutter (multiple MD files, loose scripts) that needs organization.

## Constraints

- **Compatibility**: Must maintain backward compatibility with existing API
- **Dependencies**: Requires `rationalmatrices` package from GitHub
- **Workflow**: Focus on developer workflow, not end-user features
- **Time**: Prioritize cleanup that impacts development experience

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| Keep only essential MDs in root | Reduce clutter, follow R package conventions | — Pending |
| 100% Roxygen documentation | Ensure package is well-documented for future maintenance | — Pending |
| Allow ≤1 warning in check | Balance perfectionism with practical development | — Pending |
| Comprehensive code organization | Address C++, R, and S3 structure for maintainability | — Pending |

---
*Last updated: 2026-01-24 after initialization*
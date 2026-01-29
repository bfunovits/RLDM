# Phase 4: Documentation Completeness - Research

**Researched:** 2026-01-29
**Domain:** R package documentation with Roxygen2
**Confidence:** HIGH (based on existing patterns and R documentation standards)

## Summary

This research covers achieving 100% Roxygen documentation coverage for the RLDM R package. The package currently has good documentation coverage (23/24 R files documented), but needs systematic verification and completion to meet the phase requirements.

**Key findings:**
- The package uses Roxygen2 with markdown support (DESCRIPTION: `Roxygen: list(markdown = TRUE)`)
- Current documentation follows established R package patterns with mathematical notation
- Package-level documentation (`?RLDM`) is missing and needs to be created
- All exported functions appear to have documentation, but examples need verification
- The numeric file prefix system (01_, 02_, ...) provides natural grouping for @family tags

**Primary recommendation:** Use existing documentation patterns consistently, create package-level documentation, verify all exported functions have working examples, and use @inheritParams to reduce duplication.

## Standard Stack

The established tools for R package documentation:

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| Roxygen2 | 7.3.2 | Documentation generation | Standard for R packages, integrates with devtools |
| devtools | Latest | Package development workflow | Standard tool for R package development |
| Rdpack | 2.0+ | Bibliography management | Used for references in existing documentation |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| testthat | 3.0+ | Testing examples | For verifying examples work |
| roxygen2md | N/A | Markdown conversion | If converting existing Rd to markdown |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Roxygen2 | manual .Rd files | Roxygen2 is standard, reduces maintenance burden |
| devtools::document() | R CMD check --as-cran | devtools provides better development workflow |

**Installation:**
```bash
# Already in DESCRIPTION Imports/Suggests
# Roxygen2 configured in DESCRIPTION
```

## Architecture Patterns

### Recommended Project Structure
```
R/
├── 01_representations_classes.R    # Model classes
├── 02_templates.R                  # Parameter templates
├── 03_properties_*.R              # Derived properties
├── 04_timeseries_*.R              # Time series operations
├── 05_estimation_*.R              # Estimation methods
├── 06_visualization_*.R           # Plot methods
├── 07_comparison_metrics.R        # Model comparison
├── 08_utilities_*.R               # Support functions
├── RLDM-package.R                 # Package-level docs
└── RcppExports.R                  # Auto-generated (no docs)
```

### Pattern 1: Function Documentation Structure
**What:** Standard Roxygen sections for exported functions
**When to use:** All exported functions
**Example:**
```r
#' Brief Title
#'
#' Longer description with mathematical context (1-2 sentences).
#' Link to theory papers or vignettes when appropriate.
#'
#' @param param1 Description of first parameter
#' @param param2 Description of second parameter
#'
#' @return Object of class `classname` with description
#'
#' @export
#'
#' @examples
#' # Minimal proof-of-concept example
#' set.seed(123)  # Only for stochastic functions
#' result <- function_name(param1 = value1, param2 = value2)
#' result
```

### Pattern 2: @inheritParams Usage
**What:** Reuse parameter documentation across related functions
**When to use:** When multiple functions share parameters with identical meaning
**Example:**
```r
#' Function A
#'
#' @param x,y,z Shared parameters
#' @export
fun_a <- function(x, y, z) { }

#' Function B
#'
#' @inheritParams fun_a
#' @param w Additional parameter
#' @export
fun_b <- function(x, y, z, w) { }
```

### Pattern 3: Package-level Documentation
**What:** Documentation for the package itself (`?RLDM`)
**When to use:** Required for all R packages
**Example:**
```r
#' RLDM: Rational Linear Dynamic Models
#'
#' @description
#' This package provides tools for stationary processes with rational spectral density.
#'
#' @details
#' The package implements VARMA and state space models with methods for estimation,
#' simulation, prediction, and model comparison.
#'
#' @section Package Organization:
#' R source files use a numeric prefix system:
#' - 01_: Model representations and classes
#' - 02_: Parameter templates
#' - 03_: Derived properties (autocov, spectral density, etc.)
#' - 04_: Time series operations
#' - 05_: Estimation methods
#' - 06_: Visualization
#' - 07_: Model comparison
#' - 08_: Utilities
#'
#' @section Getting Started:
#' See vignettes:
#' - `vignette("0_getting_started")` for beginner introduction
#' - `vignette("1_case_study")` for practical workflow
#' - `vignette("2_technical_reference")` for technical details
#'
#' @author Wolfgang Scherrer, Bernd Funovits
#' Maintainer: <bernd.funovits@gmail.com>
#'
#' @docType package
#' @name RLDM
NULL
```

### Anti-Patterns to Avoid
- **Over-documenting examples:** Examples should be minimal proof-of-concept, not comprehensive tutorials
- **Missing examples:** All exported functions must have working examples
- **Long-running examples:** Examples should execute quickly (<5 seconds)
- **Error demonstrations:** Examples should only show successful usage
- **Inconsistent formatting:** Follow existing patterns for mathematical notation

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Documentation validation | Custom scripts | `devtools::check()` | Comprehensive validation including examples |
| Example testing | Manual testing | `testthat::test_examples()` | Automated verification examples work |
| Parameter documentation | Copy-paste text | `@inheritParams` | Reduces duplication, ensures consistency |
| Cross-references | Manual links | `@seealso`, `@family` | Automated linking in documentation |
| Bibliography management | Manual citations | `Rdpack` with `\insertCite` | Proper reference management |

**Key insight:** Roxygen2 and devtools provide comprehensive documentation tooling. Custom solutions introduce maintenance burden and miss edge cases.

## Common Pitfalls

### Pitfall 1: Missing Package-level Documentation
**What goes wrong:** `?RLDM` returns "No documentation for 'RLDM'"
**Why it happens:** No `RLDM-package.R` file or missing `@docType package`
**How to avoid:** Create `R/RLDM-package.R` with proper structure
**Warning signs:** Package check note about missing package documentation

### Pitfall 2: Broken Examples
**What goes wrong:** Examples fail during `R CMD check`
**Why it happens:** Missing dependencies, undefined variables, or runtime errors
**How to avoid:** Test examples with `devtools::run_examples()` or `testthat::test_examples()`
**Warning signs:** Check warnings about examples

### Pitfall 3: Inconsistent Parameter Documentation
**What goes wrong:** Same parameter documented differently across functions
**Why it happens:** Manual documentation without `@inheritParams`
**How to avoid:** Identify shared parameters, use `@inheritParams` consistently
**Warning signs:** Similar functions with different parameter descriptions

### Pitfall 4: Missing @export Tags
**What goes wrong:** Functions not exported despite being in NAMESPACE
**Why it happens:** Roxygen comment missing `@export` tag
**How to avoid:** Verify all exported functions have `@export` tags
**Warning signs:** `devtools::document()` warnings about missing exports

## Code Examples

Verified patterns from existing codebase:

### Function with Mathematical Documentation
```r
#' sigma_L Structure
#'
#' Create templates for the left square root \eqn{L} of the noise covariance
#' matrix \eqn{\Sigma = LL'}. This means that \eqn{L} is parametrized as
#' \deqn{\mbox{vec}(L) = h + H \theta}{vec(L) = h + H \theta} with a
#' (\eqn{k}-dimensional) parameter vector \eqn{\theta}.
#'
#' @param sigma_L numeric (n x n) matrix, where the free entries are coded with NAs
#' @param structure character string, determines the "structure" of sigma_L
#'
#' @return List with slots
#'   \itemize{
#'   \item `h` (\eqn{n^2}-dimensional vector),
#'   \item `H` (\eqn{(n^2, k)}-dimensional matrix) and
#'   \item `n.par` (integer) is the number of free/deep parameters (\eqn{=k}).
#'   }
#'
#' @export
#'
#' @examples
#' sigma_L = matrix(c(0, NA, 1, 0, 2, 3, NA, 1, 1), nrow = 3, ncol = 3)
#' sigma_L
#'
#' tmpl = tmpl_sigma_L(sigma_L, structure = 'as_given')
#' th = -(1:tmpl$n.par)
#' matrix(tmpl$h + tmpl$H %*% th, nrow = 3, ncol = 3)
```

### Function with @inheritParams
```r
#' Estimate AR model with OLS
#'
#' @inheritParams est_ar
#' @param mean_estimate character string, how to estimate the mean
#'
#' @export
est_ar_ols <- function(obj, p.max = NULL, penalty = -1,
                       mean_estimate = c("zero", "sample.mean", "intercept"),
                       trace = TRUE) {
  # Implementation
}
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Manual .Rd files | Roxygen2 in-source documentation | ~2010 | Reduced duplication, easier maintenance |
| Plain text documentation | Markdown support in Roxygen2 | Roxygen2 7.0.0 | Readable source, formatted output |
| No example testing | `testthat::test_examples()` | testthat 3.0.0 | Automated example verification |

**Deprecated/outdated:**
- `\code{}` and `\emph{}` in favor of markdown `` `code` `` and `*emphasis*`
- Manual bibliography management in favor of `Rdpack`

## Open Questions

Things that couldn't be fully resolved:

1. **Current documentation coverage percentage**
   - What we know: 23/24 R files have some documentation
   - What's unclear: Exact percentage of exported functions with complete documentation
   - Recommendation: Use `devtools::check()` to identify missing documentation

2. **Example completeness verification**
   - What we know: Some functions have examples
   - What's unclear: Whether all exported functions have working examples
   - Recommendation: Systematic check of all exported functions

3. **@family tags strategy**
   - What we know: Numeric file prefixes provide natural grouping
   - What's unclear: Optimal grouping strategy for functions
   - Recommendation: Group by file prefix (e.g., `@family 01_representations`)

## Sources

### Primary (HIGH confidence)
- Existing RLDM codebase - Current documentation patterns
- DESCRIPTION file - Roxygen2 configuration
- NAMESPACE file - Exported functions list

### Secondary (MEDIUM confidence)
- Roxygen2 documentation (training knowledge) - Standard practices
- R package development best practices (training knowledge)

### Tertiary (LOW confidence)
- WebSearch attempted but encountered API errors

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - Based on existing package configuration
- Architecture: HIGH - Based on existing patterns in codebase
- Pitfalls: MEDIUM - Based on common R package documentation issues

**Research date:** 2026-01-29
**Valid until:** 2026-02-28 (30 days for stable R documentation practices)

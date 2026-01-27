# Phase 2: Code Organization - Research

**Researched:** 2026-01-28
**Domain:** R package code organization, C++ standards, S3 method architecture
**Confidence:** HIGH (based on established R/C++ conventions and current codebase analysis)

## Summary

This research covers best practices for organizing R package code with C++ extensions, specifically for the RLDM package which uses Rcpp/RcppArmadillo, has a numeric prefix file organization system, and implements S3 object-oriented methods. The package currently has well-structured R code following a numeric prefix convention but needs C++ code cleanup, S3 method organization review, and dependency management optimization.

**Key findings:**
- C++ code needs Doxygen documentation and Google C++ Style Guide compliance
- R files follow numeric prefix convention consistently (01_*, 02_*, etc.)
- S3 methods are currently organized by function (all methods for a generic in one file)
- Dependencies in DESCRIPTION need version constraints and unused import review
- Build artifacts (*.o, *.so files) need cleanup from src/ directory

**Primary recommendation:** Implement Doxygen documentation for C++ functions, apply Google C++ Style Guide formatting, maintain current R file organization strategy, optimize S3 method registration, and add version constraints to dependencies while cleaning up unused imports.

## Standard Stack

The established tools and standards for R package code organization:

### C++ Code Standards
| Standard | Purpose | Why Standard |
|----------|---------|--------------|
| **Doxygen** | API documentation generation | Industry standard for C++ documentation, supports @param, @return, @brief tags |
| **Google C++ Style Guide** | Code formatting and naming conventions | Widely adopted, improves readability and consistency |
| **Rcpp/RcppArmadillo** | R-C++ integration | Standard for R packages with C++ code, provides seamless integration |
| **Roxygen2 for R** | R documentation generation | Standard for R package documentation, integrates with NAMESPACE |

### R Code Organization
| Pattern | Purpose | When to Use |
|---------|---------|-------------|
| **Numeric prefix system** | File organization by workflow/function | When files need logical grouping and load order control |
| **S3 method registration** | Object-oriented dispatch | For R's S3 object system, preferred over S4 for simplicity |
| **Roxygen @export** | Function export management | Automated NAMESPACE generation, better than manual management |

### Dependency Management
| Tool/Approach | Purpose | Why Standard |
|---------------|---------|--------------|
| **DESCRIPTION version constraints** | Dependency version specification | Ensures compatibility, prevents breakage |
| **Remotes field** | GitHub dependency specification | Standard for packages not on CRAN |
| **Import/Depends distinction** | Dependency type specification | Best practice for namespace management |

**Installation:** No installation needed - these are coding standards and organizational patterns.

## Architecture Patterns

### Recommended C++ File Structure
```
src/
├── kf.cpp                    # Kalman filter implementation
├── rls_core.cpp             # Recursive least squares
├── solve_fwd_bwd_cll.cpp    # Forward-backward solving
├── pf.cpp                   # Particle filter (if exists)
├── Makevars                 # Build configuration
└── Makevars.win            # Windows build configuration
```

### Pattern 1: Doxygen Documentation for C++ Functions
**What:** Use Doxygen-style comments for all C++ functions
**When to use:** All C++ functions in R packages
**Example:**
```cpp
/**
 * @brief Kalman filter implementation for state space models
 *
 * @param A State transition matrix (s x s)
 * @param C Observation matrix (m x s)
 * @param Q Process noise covariance (s x s)
 * @param R Observation noise covariance (m x m)
 * @param S Cross-covariance matrix (s x m)
 * @param y_t Observations (m x N)
 * @param P1 Initial state covariance (s x s)
 * @param a1 Initial state estimate (s x 1)
 * @return Rcpp::List containing filtered states, predictions, and log-likelihood
 *
 * @note Assumes column-major data organization for compatibility with R
 */
Rcpp::List kf_cpp(const arma::mat& A, const arma::mat& C,
                  const arma::mat& Q, const arma::mat& R, const arma::mat& S,
                  const arma::mat& y_t, const arma::mat& P1, const arma::colvec& a1) {
  // Implementation...
}
```

### Pattern 2: Google C++ Style Guide Compliance
**What:** Follow Google C++ Style Guide for formatting and naming
**When to use:** All C++ code in the package
**Key rules:**
- 2-space indentation (not tabs)
- 80-character line limit
- `snake_case` for functions and variables
- `CamelCase` for classes and structs
- `kConstantName` for constants
- Braces on same line for functions, new line for control statements

**Example:**
```cpp
// Good Google C++ style
void compute_kalman_gain(const arma::mat& covariance,
                         const arma::mat& observation_matrix,
                         arma::mat* kalman_gain) {
  if (covariance.n_rows == 0) {
    return;
  }

  *kalman_gain = covariance * trans(observation_matrix);
}

// Constants
const int kMaxIterations = 100;
const double kTolerance = 1e-6;
```

### Pattern 3: S3 Method Organization Strategies
**What:** Organize S3 methods either by class or by generic function
**When to use:** R packages with S3 object systems
**Options:**
1. **By generic function** (current approach): All `autocov.*` methods in one file
   - Pros: Easy to find all implementations of a generic
   - Cons: Files can become large with many methods

2. **By class**: All methods for `armamod` class in one file
   - Pros: All class behavior in one place
   - Cons: Hard to find specific generic implementations

**Current RLDM approach:** By generic function (recommended to maintain)

### Pattern 4: Roxygen S3 Method Registration
**What:** Use Roxygen tags for S3 method registration
**When to use:** All S3 methods in R packages
**Example:**
```r
#' @rdname autocov
#' @export
autocov.armamod = function(obj, type = c('covariance', 'correlation', 'partial'),
                           lag.max = 12, ...) {
  # Method implementation
}
```

### Anti-Patterns to Avoid
- **Mixed documentation styles:** Using Roxygen (`//'`) for C++ instead of Doxygen
- **Inconsistent naming:** Mixing naming conventions within C++ code
- **Unregistered S3 methods:** Methods not properly exported in NAMESPACE
- **Missing version constraints:** Dependencies without version bounds
- **Build artifacts in src/:** Leaving *.o, *.so files in source directory

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| C++ documentation | Custom comment system | Doxygen with @param, @return | Standard, toolable, widely understood |
| Code formatting | Manual formatting | Google C++ Style Guide | Consistent, automated tools available |
| S3 method registration | Manual NAMESPACE edits | Roxygen @export with method names | Automated, less error-prone |
| Dependency versioning | No version constraints | DESCRIPTION with version bounds | Prevents breakage, ensures compatibility |
| Build artifact cleanup | Manual deletion | .gitignore patterns, clean build process | Automated, prevents accidental commits |

**Key insight:** The R/C++ ecosystem has well-established standards for all aspects of package development. Following these standards reduces cognitive overhead and improves maintainability.

## Common Pitfalls

### Pitfall 1: Incomplete C++ Documentation
**What goes wrong:** C++ functions lack proper documentation, making maintenance difficult
**Why it happens:** Using Roxygen (`//'`) for C++ instead of Doxygen
**How to avoid:**
1. Use Doxygen-style comments (`/** */`) for all C++ functions
2. Include @param for all parameters, @return for return values
3. Add @brief for function summary, @note for important details
**Warning signs:** C++ functions with only brief comments or Roxygen-style comments

### Pitfall 2: Inconsistent Code Style
**What goes wrong:** Mixed naming conventions and formatting reduce readability
**Why it happens:** No enforced style guide, different contributors use different styles
**How to avoid:**
1. Adopt Google C++ Style Guide as standard
2. Use clang-format or similar tool for consistency
3. Establish clear naming conventions (snake_case, CamelCase, etc.)
**Warning signs:** Mixing `camelCase` and `snake_case`, inconsistent indentation

### Pitfall 3: S3 Method Registration Issues
**What goes wrong:** S3 methods not properly registered in NAMESPACE
**Why it happens:** Manual NAMESPACE editing errors or missing @export tags
**How to avoid:**
1. Always use `@export` with full method name (e.g., `autocov.armamod`)
2. Run `devtools::document()` to regenerate NAMESPACE
3. Verify NAMESPACE contains `S3method()` entries for all methods
**Warning signs:** "no applicable method" errors, methods not appearing in NAMESPACE

### Pitfall 4: Dependency Version Issues
**What goes wrong:** Package breaks when dependencies update
**Why it happens:** Missing version constraints in DESCRIPTION
**How to avoid:**
1. Add version constraints (e.g., `dplyr (>= 1.0.0)`)
2. Review imports for actually used functions
3. Remove unused imports to reduce dependency surface
**Warning signs:** Package checks fail after dependency updates, unused imports in NAMESPACE

### Pitfall 5: Build Artifact Pollution
**What goes wrong:** *.o, *.so files committed to repository or left in src/
**Why it happens:** Missing .gitignore patterns, incomplete cleanup
**How to avoid:**
1. Add `src/*.o`, `src/*.so`, `src/*.dll` to .gitignore
2. Clean build artifacts before commits
3. Use `devtools::clean_dll()` for Rcpp cleanup
**Warning signs:** Binary files in git status, large object files in repository

## Code Examples

Verified patterns from current codebase and standards:

### Current R File Organization (Good)
```
R/
├── 01_representations_classes.R    # Model classes: armamod(), stspmod(), rmfdmod()
├── 02_templates.R                  # Parameter templates (tmpl_*), fill_template()
├── 03_properties_autocov.R         # autocov() methods for all classes
├── 03_properties_frequency.R       # freqresp() methods
├── 03_properties_impulse.R         # impresp() methods
├── 03_properties_poles.R          # poles() methods
├── 03_properties_spectral.R       # spectrald() methods
├── 04_timeseries_predict.R        # predict() methods
├── 04_timeseries_simulate.R       # sim() methods
├── 04_timeseries_solve.R          # solve_de() methods
├── 05_estimation_ar.R             # AR estimation methods
├── 05_estimation_arma_hrk.R       # ARMA HRK estimation
├── 05_estimation_likelihood.R     # Likelihood estimation
├── 05_estimation_particle.R       # Particle filter estimation
├── 05_estimation_rls.R            # Recursive least squares
├── 05_estimation_subspace.R       # Subspace estimation
├── 06_visualization_plot.R        # plot() methods
├── 06_visualization_plot_prediction.R # plot_prediction()
├── 07_comparison_metrics.R        # Model comparison metrics
├── 08_utilities_data.R            # Data documentation
├── 08_utilities_print.R           # print() methods
├── 08_utilities_str.R             # str() methods
├── RcppExports.R                  # Auto-generated (DO NOT EDIT)
└── RLDM-package.R                 # Package documentation
```

### Current S3 Method Pattern (Good)
```r
# In 03_properties_autocov.R
#' @rdname autocov
#' @export
autocov.default = function(obj, type = c('covariance', 'correlation', 'partial'), ...) {
  # Default implementation for data objects
}

#' @rdname autocov
#' @export
autocov.armamod = function(obj, type = c('covariance', 'correlation', 'partial'), ...) {
  # Implementation for armamod objects
}

#' @rdname autocov
#' @export
autocov.stspmod = function(obj, type = c('covariance', 'correlation', 'partial'), ...) {
  # Implementation for stspmod objects
}
```

### Current C++ Documentation (Needs Improvement)
```cpp
// Current (Roxygen style, not ideal for C++)
//' @name kf
//' @rdname kf
//' @export
// [[Rcpp::export]]
Rcpp::List kf_cpp(const arma::mat& A, const arma::mat& C,
                  const arma::mat& Q, const arma::mat& R, const arma::mat& S,
                  const arma::mat& y_t, const arma::mat& P1, const arma::colvec& a1) {

// Recommended (Doxygen style)
/**
 * @brief Kalman filter implementation
 * @param A State transition matrix
 * @param C Observation matrix
 * @param Q Process noise covariance
 * @param R Observation noise covariance
 * @param S Cross-covariance matrix
 * @param y_t Observations (m x N)
 * @param P1 Initial state covariance
 * @param a1 Initial state estimate
 * @return Rcpp::List with filtered states, predictions, log-likelihood
 */
// [[Rcpp::export]]
Rcpp::List kf_cpp(const arma::mat& A, const arma::mat& C,
                  const arma::mat& Q, const arma::mat& R, const arma::mat& S,
                  const arma::mat& y_t, const arma::mat& P1, const arma::colvec& a1) {
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Manual C++ documentation | Doxygen with @param, @return | Industry standard | Automated documentation, better maintainability |
| Ad-hoc C++ style | Google C++ Style Guide | Widely adopted | Consistent, readable code |
| Manual NAMESPACE editing | Roxygen @export | R 3.0+ | Automated, less error-prone |
| No version constraints | DESCRIPTION version bounds | Best practice | Prevents breakage, ensures compatibility |
| Mixed R file organization | Numeric prefix system | RLDM-specific | Logical grouping, clear dependencies |

**Deprecated/outdated:**
- **Roxygen for C++:** Should use Doxygen for C++, Roxygen only for R
- **Inconsistent naming:** Should follow Google C++ Style Guide
- **Unconstrained dependencies:** Should specify version bounds
- **Build artifacts in src/:** Should be cleaned and gitignored

## Open Questions

Things that couldn't be fully resolved:

1. **C++ variable naming consistency**
   - What we know: Current C++ code uses mixed naming conventions
   - What's unclear: Whether to refactor all variable names or only new code
   - Recommendation: Apply Google C++ Style Guide to all C++ files for consistency

2. **S3 method inheritance patterns**
   - What we know: `rldm` is parent class for `armamod`, `stspmod`, `rmfdmod`
   - What's unclear: Whether all methods should have `rldm` default method
   - Recommendation: Keep current inheritance pattern, ensure `rldm` methods exist where appropriate

3. **Dependency version constraints**
   - What we know: Current DESCRIPTION has no version constraints
   - What's unclear: Appropriate version bounds for each dependency
   - Recommendation: Add conservative version constraints based on current usage

4. **Unused import cleanup**
   - What we know: Some imports may not be used
   - What's unclear: Which imports are actually used by package functions
   - Recommendation: Analyze namespace usage, remove truly unused imports

## Sources

### Primary (HIGH confidence)
- Current codebase analysis - Existing R/C++ code structure and patterns
- Google C++ Style Guide - Official style guide for C++ code formatting
- Doxygen documentation - Standard for C++ API documentation
- R Package Development (r-pkgs.org) - Standard R package conventions

### Secondary (MEDIUM confidence)
- Rcpp documentation - Best practices for R-C++ integration
- Roxygen2 documentation - S3 method registration patterns
- R community best practices - Dependency management and versioning

### Tertiary (LOW confidence)
- Web search attempted but encountered API errors - Would verify C++ standards with official documentation

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - Well-established C++/R standards
- Architecture: HIGH - Clear patterns from codebase analysis
- Pitfalls: HIGH - Common issues in R/C++ package development

**Research date:** 2026-01-28
**Valid until:** 2026-02-28 (30 days - stable conventions)
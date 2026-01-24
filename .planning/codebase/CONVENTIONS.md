# Coding Conventions

**Analysis Date:** 2026-01-24

## Naming Patterns

**Files:**
- Use numeric prefix system: `NN_category_subcategory.R` where NN is two-digit number
- Examples: `01_representations_classes.R`, `05_estimation_particle.R`, `03_properties_autocov.R`
- C++ files use descriptive names: `kf.cpp`, `pf.cpp`, `solve_fwd_bwd_cll.cpp`

**Functions:**
- Use snake_case for R functions: `armamod()`, `pfilter()`, `tmpl_stsp_full()`
- Use descriptive names that indicate purpose: `extract_theta()`, `fill_template()`, `sim()`
- C++ functions use snake_case with `_cpp` suffix: `kf_cpp()`, `pf_sir_cpp()`

**Variables:**
- Single lowercase letters for dimensions: `m` (outputs), `n` (inputs), `s` (states), `p`, `q` (AR/MA orders)
- Descriptive names for objects: `model`, `data`, `y`, `tmpl` (template), `th` (theta/parameters)
- Matrix variables: `A`, `B`, `C`, `D`, `Q`, `R`, `S` for state space matrices

**Types:**
- Model classes use descriptive suffixes: `armamod`, `stspmod`, `rmfdmod`
- Template functions prefixed with `tmpl_`: `tmpl_arma_pq()`, `tmpl_stsp_full()`
- Estimation functions prefixed with `est_`: `est_ar()`, `est_arma_hrk3()`

## Code Style

**Formatting:**
- 2-space indentation (standard R style)
- Line length appears flexible, typically under 80 characters
- Opening braces on same line, closing braces on own line
- Spaces around operators and after commas

**Linting:**
- No `.lintr` file detected
- No formal linting configuration found
- Code follows base R style conventions

## Import Organization

**Order:**
1. Roxygen documentation blocks
2. Function definitions with parameter validation
3. Core implementation logic
4. Return statements

**Path Aliases:**
- Not applicable for R package structure
- C++ uses `#include <RcppArmadillo.h>` and `#include <rationalmatrices_lyapunov.h>`

## Error Handling

**Patterns:**
- Use `stop()` with descriptive messages for invalid inputs
- Input validation at function start: `if (!inherits(sys, 'lmfd')) stop('"sys" must be an lmfd object')`
- Use `expect_no_error()` in tests rather than try-catch blocks
- C++ uses `stop()` from Rcpp for error handling

## Logging

**Framework:** Base R `message()`, `warning()`, `stop()`

**Patterns:**
- Use `warning()` for non-fatal issues (e.g., APF with cross-covariance)
- Use `message()` for informational output (minimal usage)
- No formal logging framework detected

## Comments

**When to Comment:**
- File headers with purpose and dependencies (see `CONTRIBUTING.md`)
- Complex mathematical operations
- Non-obvious implementation details
- C++ files have descriptive headers explaining purpose

**Roxygen Documentation:**
- Extensive use of Roxygen2 with markdown support
- Mathematical notation with LaTeX in `\deqn{}` and `\eqn{}`
- Examples sections with runnable code
- Cross-references with `@seealso` and `@rdname`
- Parameter documentation with type expectations

## Function Design

**Size:** Functions are moderate length, typically under 100 lines
**Parameters:** Clear parameter names with defaults where appropriate
**Return Values:** Consistent return structures (lists for complex objects)

## Module Design

**Exports:** All user-facing functions exported with `@export`
**Barrel Files:** Not used - each file contains related functionality grouped by numeric prefix
**S3 Methods:** Methods grouped by operation in single files (e.g., all `autocov()` methods in `03_properties_autocov.R`)

---

*Convention analysis: 2026-01-24*
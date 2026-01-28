# Phase 3: Build Verification - Research

**Researched:** 2026-01-28
**Domain:** R package quality assurance and build verification
**Confidence:** HIGH

## Summary

This research investigates how to implement comprehensive build verification for the RLDM R package. The package is a mature R package with C++ extensions (Rcpp/RcppArmadillo), S3 class system, vignettes, and dependencies including a sister package (`rationalmatrices`) installed from GitHub.

The standard approach for R package verification centers on `devtools::check()` (which internally uses `rcmdcheck::rcmdcheck()`), with additional validation of namespace management, installation from source, and vignette building. Key findings include:

- **Core verification tool:** `devtools::check()` provides comprehensive checking with configurable strictness
- **Namespace management:** Proper imports/exports are validated automatically by R CMD check
- **Build artifacts:** Clean environment management is critical for reproducible verification
- **Warning handling:** Standard practice allows certain informational warnings (non-ASCII, time verification) while requiring action on substantive warnings
- **Automated workflows:** GitHub Actions with r-lib/actions provide standardized checking patterns

**Primary recommendation:** Use `devtools::check()` with default arguments in a clean R session, clean build artifacts before/after, implement automatic fixes for common namespace issues, and log all results for review.

## Standard Stack

The established libraries/tools for R package verification:

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| devtools | Latest | Package development tools | Standard R package development toolkit, wraps rcmdcheck |
| rcmdcheck | Latest | Programmatic R CMD check | Underlying engine for devtools::check(), provides structured results |
| testthat | ≥2.1.0 | Unit testing framework | Already used in package, required for check() |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| Rcpp | ≥1.0.0 | C++ integration | Required for compiling C++ code during check |
| roxygen2 | 7.3.2 | Documentation generation | Already configured, needed for documentation checks |
| knitr | ≥1.30 | Vignette building | Required for vignette compilation during check |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| devtools::check() | rcmdcheck::rcmdcheck() | devtools provides higher-level interface, rcmdcheck gives more control |
| Manual checking | Automated scripts | Manual is error-prone, automation ensures consistency |

**Installation:**
```bash
# All dependencies are already in DESCRIPTION
# Sister package requires special handling:
remotes::install_github("bfunovits/rationalmatrices")
```

## Architecture Patterns

### Recommended Verification Structure
```
Verification Workflow:
1. Clean build artifacts
2. Install dependencies (including sister package)
3. Run devtools::check()
4. Parse results (errors/warnings/notes)
5. Attempt automatic fixes for common issues
6. Re-run check if fixes applied
7. Log all results
8. Clean build artifacts
```

### Pattern 1: Clean Environment Verification
**What:** Run checks in clean R session with fresh package build
**When to use:** Always for final verification to ensure reproducibility
**Example:**
```r
# Source: R package best practices
# Clean environment approach
system("R CMD build .")
system("R CMD check RLDM_*.tar.gz --as-cran")
```

### Pattern 2: Progressive Fixing with Logging
**What:** Collect all issues, attempt fixes, log results, repeat
**When to use:** When automatic fixes are possible (namespace issues, missing imports)
**Example:**
```r
# Source: rcmdcheck documentation
result <- rcmdcheck::rcmdcheck(args = "--no-manual")
if (length(result$errors) > 0) {
  log_errors(result$errors)
  # Attempt fixes...
  result <- rcmdcheck::rcmdcheck(args = "--no-manual")
}
```

### Anti-Patterns to Avoid
- **Checking in dirty environment:** Build artifacts from previous runs can hide issues
- **Stopping on first error:** Should collect all issues for comprehensive reporting
- **Ignoring namespace warnings:** Missing imports/exports cause runtime failures
- **Not cleaning between runs:** Can lead to false positives/negatives

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Parsing check output | Custom regex parsing | rcmdcheck::parse_check() | Handles all R versions, output formats, edge cases |
| Dependency installation | Manual install scripts | r-lib/actions/setup-r-dependencies | Handles system dependencies, version conflicts |
| Clean environment | Manual file deletion | devtools::clean_dll(), remove.packages() | Ensures complete cleanup, handles platform differences |
| Warning classification | Custom heuristics | rcmdcheck result structure | Properly categorizes errors/warnings/notes |

**Key insight:** R package checking has well-established patterns and edge cases (platform differences, R version variations, dependency resolution) that custom solutions consistently get wrong.

## Common Pitfalls

### Pitfall 1: Dirty Build Environment
**What goes wrong:** Previous compilation artifacts cause false passes or strange errors
**Why it happens:** R retains .o, .so files between sessions, C++ code not recompiled
**How to avoid:** Always clean src/*.o, src/*.so, and use devtools::clean_dll()
**Warning signs:** "object file not found" errors, inconsistent check results

### Pitfall 2: Missing Sister Package
**What goes wrong:** `rationalmatrices` not installed, causing check failures
**Why it happens:** Remotes: field in DESCRIPTION not automatically handled by all tools
**How to avoid:** Explicit installation step before checking
**Warning signs:** "Package 'rationalmatrices' not available" errors

### Pitfall 3: Namespace Mismatches
**What goes wrong:** Functions exported but not imported, or vice versa
**Why it happens:** Roxygen2 @export/@importFrom tags not properly maintained
**How to avoid:** Run devtools::document() before checking, verify NAMESPACE
**Warning signs:** "no visible binding for global variable" warnings

### Pitfall 4: Vignette Build Failures
**What goes wrong:** Vignettes fail to knit during check
**Why it happens:** Missing vignette dependencies, caching issues
**How to avoid:** Test vignette building separately, ensure all Suggest: packages installed
**Warning signs:** "Error: processing vignette failed" during check

## Code Examples

Verified patterns from official sources:

### Running Comprehensive Check
```r
# Source: devtools documentation
# Run check with all standard tests
check_result <- devtools::check(
  document = TRUE,      # Update documentation first
  build_args = NULL,    # Default build arguments
  manual = TRUE,        # Include manual checks
  cran = FALSE,         # Don't use --as-cran (per decision)
  remote = FALSE,
  incoming = FALSE
)

# Access results
errors <- check_result$errors
warnings <- check_result$warnings
notes <- check_result$notes
```

### Clean Build Environment
```r
# Source: R package development best practices
# Clean C++ compilation artifacts
clean_files <- c(
  list.files("src", pattern = "\\.(o|so|dll)$", full.names = TRUE),
  "src/RcppExports.cpp",
  "src/RcppExports.o"
)
unlink(clean_files)

# Clean R build artifacts
unlink("RLDM.Rcheck", recursive = TRUE)
unlink(list.files(pattern = "RLDM_.*\\.tar\\.gz$"))
```

### Installing Dependencies
```r
# Source: r-lib/actions patterns
# Install sister package from GitHub
if (!requireNamespace("rationalmatrices", quietly = TRUE)) {
  remotes::install_github("bfunovits/rationalmatrices")
}

# Install package dependencies
devtools::install_deps(dependencies = TRUE)
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Manual R CMD check | devtools::check() wrapper | devtools v1.0 (2015) | Standardized, programmatic access |
| Custom check scripts | rcmdcheck package | rcmdcheck v1.0 (2016) | Structured results, better parsing |
| Local checking only | GitHub Actions integration | r-lib/actions v1 (2019) | Automated CI/CD, multi-platform |

**Deprecated/outdated:**
- **Custom Makefiles:** Use standard R CMD build/check infrastructure
- **Manual vignette building:** Let devtools::check() handle it automatically
- **Platform-specific checks:** Use rcmdcheck which handles platform differences

## Open Questions

1. **Warning threshold configuration**
   - What we know: devtools::check() has `cran` argument for strictness
   - What's unclear: Exact configuration for "≤1 warning" allowance
   - Recommendation: Use `cran = FALSE`, manually filter acceptable warnings

2. **Clean session vs current session**
   - What we know: Clean session ensures reproducibility
   - What's unclear: Performance impact vs reliability tradeoff
   - Recommendation: Use clean R session for verification, current for development

3. **Automatic fix scope**
   - What we know: Can fix missing imports, documentation issues
   - What's unclear: Which fixes are safe to automate vs require human review
   - Recommendation: Automate namespace fixes only, log others for review

## Sources

### Primary (HIGH confidence)
- R Packages book (r-pkgs.org) - R CMD check documentation and best practices
- devtools.r-lib.org - devtools::check() function documentation
- rcmdcheck.r-lib.org - rcmdcheck package capabilities and usage

### Secondary (MEDIUM confidence)
- R Extensions Manual - Official R package building and installation guidelines
- r-lib/actions repository - GitHub Actions workflows for R package checking

### Tertiary (LOW confidence)
- Community patterns for R package verification (needs validation against current RLDM implementation)

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - Well-documented R package development ecosystem
- Architecture: HIGH - Established patterns from r-lib and CRAN
- Pitfalls: MEDIUM - Based on documented common issues, needs RLDM-specific validation

**Research date:** 2026-01-28
**Valid until:** 2026-02-28 (30 days - stable R ecosystem)
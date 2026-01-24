# Technology Stack

**Project:** RLDM Package Cleanup & Streamlining
**Domain:** R package maintenance, documentation, and build verification
**Researched:** 2026-01-24
**Overall confidence:** MEDIUM (WebSearch verification unavailable, based on training data and package analysis)

## Recommended Stack

### Core Development Framework

| Technology | Version | Purpose | Why Recommended |
|------------|---------|---------|-----------------|
| devtools | ≥2.4.5 | Package development workflow | Standard tool for R package development; provides `load_all()`, `document()`, `test()`, `check()`, `build()` |
| usethis | ≥2.2.2 | Package setup and automation | Automates repetitive tasks, creates consistent project structure, integrates with devtools |
| roxygen2 | ≥7.3.0 | Documentation generation | Standard for R package documentation with markdown support; integrates with devtools |

### Testing & Quality Assurance

| Technology | Version | Purpose | Why Recommended |
|------------|---------|---------|-----------------|
| testthat | ≥3.2.0 | Unit testing framework | Modern testing framework with snapshot tests, parallel execution, and improved reporting |
| covr | ≥3.6.4 | Test coverage analysis | Integrates with testthat to measure code coverage; identifies undocumented code |
| lintr | ≥3.1.0 | Static code analysis | Enforces coding style, catches common errors before testing |

### Documentation & Website

| Technology | Version | Purpose | Why Recommended |
|------------|---------|---------|-----------------|
| pkgdown | ≥2.0.7 | Package website generation | Creates professional documentation websites from Roxygen comments and vignettes |
| knitr | ≥1.45 | Dynamic report generation | Required for vignettes and documentation examples |
| rmarkdown | ≥2.25 | R Markdown processing | Supports multiple output formats for documentation |

### Build & Verification

| Technology | Version | Purpose | Why Recommended |
|------------|---------|---------|-----------------|
| R CMD check | (built-in) | Package validation | Gold standard for package quality; devtools::check() wraps this |
| rhub | ≥1.1.3 | Multi-platform checking | Tests package on CRAN's build servers across platforms |
| revdepcheck | ≥1.0.0 | Reverse dependency checking | Checks impact of changes on dependent packages |

### C++ Integration (Existing)

| Technology | Version | Purpose | Why Recommended |
|------------|---------|---------|-----------------|
| Rcpp | ≥1.0.11 | C++ integration for R | Standard for high-performance C++ code in R packages |
| RcppArmadillo | ≥0.12.6.6.0 | Linear algebra in C++ | Efficient matrix operations for Kalman filtering and RLS |
| RcppParallel | ≥5.1.7 | Parallel computing | Optional for performance improvements in particle filters |

## Supporting Libraries

| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| withr | ≥2.5.2 | Temporary environment settings | For tests that need temporary directory, environment variables, or options |
| vdiffr | ≥1.0.6 | Visual regression testing | When package has plotting functions that need visual verification |
| spelling | ≥2.2.2 | Spell checking documentation | For comprehensive documentation quality checks |
| goodpractice | ≥1.0.4 | Code quality analysis | For automated code review and best practices enforcement |

## Development Tools

| Tool | Purpose | Notes |
|------|---------|-------|
| RStudio IDE | Integrated development environment | Provides package development pane, integrated debugging, and project management |
| git | Version control | Essential for collaboration and change tracking |
| GitHub Actions | Continuous integration | Automated testing and checking on push/pull requests |
| renv | ≥1.0.3 | Project-specific dependencies | For reproducible development environments |

## Installation

```r
# Core development dependencies
install.packages(c("devtools", "usethis", "roxygen2", "testthat", "covr", "lintr"))

# Documentation and website
install.packages(c("pkgdown", "knitr", "rmarkdown"))

# Build verification
install.packages(c("rhub", "revdepcheck"))

# Supporting libraries
install.packages(c("withr", "vdiffr", "spelling", "goodpractice"))
```

## Alternatives Considered

| Recommended | Alternative | When to Use Alternative |
|-------------|-------------|-------------------------|
| testthat (≥3.2.0) | RUnit, tinytest | When working with legacy codebases that already use these frameworks |
| pkgdown | staticdocs, bookdown | staticdocs is deprecated; bookdown for book-length documentation |
| GitHub Actions | Travis CI, CircleCI | When already invested in other CI platforms |
| RStudio IDE | VS Code with R extension | When preferring lightweight editors or cross-language development |

## What NOT to Use

| Avoid | Why | Use Instead |
|-------|-----|-------------|
| `R CMD build` directly | Manual process, error-prone | `devtools::build()` with proper configuration |
| Manual documentation | Inconsistent, hard to maintain | roxygen2 with markdown support |
| Ad-hoc testing | Unreliable, not automated | testthat with systematic test organization |
| Global package installation | Version conflicts, reproducibility issues | renv for project-specific environments |

## Stack Patterns by Variant

**If focusing on CRAN submission:**
- Use `rhub::check_for_cran()` for comprehensive platform testing
- Use `devtools::release()` for submission workflow
- Because CRAN has strict requirements across platforms

**If focusing on internal/team use:**
- Use GitHub Actions for CI/CD
- Use pkgdown for internal documentation site
- Because automation and documentation are key for team collaboration

**If package has C++ code:**
- Use `Rcpp::compileAttributes()` after C++ changes
- Use `devtools::load_all(recompile = TRUE)` for testing
- Because C++ code requires compilation and proper linking

## Version Compatibility

| Package | Compatible With | Notes |
|---------|-----------------|-------|
| testthat ≥3.2.0 | R ≥4.1.0 | Third edition introduces breaking changes |
| pkgdown ≥2.0.0 | bootstrap 5 | New template system, different configuration |
| Rcpp ≥1.0.11 | C++11 standard | Requires C++11 compatible compiler |
| devtools ≥2.4.5 | usethis ≥2.2.0 | Tight integration for package workflows |

## Sources

**Note:** WebSearch verification unavailable due to API issues. Recommendations based on:
- Analysis of existing RLDM package structure and dependencies
- Training data on R package development ecosystem (last updated Jan 2025)
- Current DESCRIPTION file shows roxygen2 7.3.2, testthat ≥2.1.0
- Package follows standard R package conventions

**Confidence levels:**
- HIGH: devtools, roxygen2, testthat (standard tools, confirmed in DESCRIPTION)
- MEDIUM: usethis, pkgdown, covr (commonly used, but versions need verification)
- LOW: rhub, revdepcheck, specific version recommendations (need current CRAN verification)

**Areas needing verification:**
- Current CRAN versions of all recommended packages (2026)
- GitHub Actions workflows for R packages
- Best practices for C++ code organization in R packages
- Integration of particle filter benchmarking with package checks

---
*Stack research for: R package maintenance, documentation, and build verification*
*Researched: 2026-01-24*
*Confidence: MEDIUM (limited verification possible)*
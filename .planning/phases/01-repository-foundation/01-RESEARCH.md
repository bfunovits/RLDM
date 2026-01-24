# Phase 1: Repository Foundation - Research

**Researched:** 2026-01-24
**Domain:** R package repository structure and organization
**Confidence:** HIGH (based on established R package conventions and current repository analysis)

## Summary

This research covers best practices for organizing an R package repository, specifically focusing on the RLDM package which contains C++ code (Rcpp/RcppArmadillo), benchmark scripts, log files, and development artifacts. The package follows standard R package structure but has accumulated development files in the root directory that need reorganization.

**Key findings:**
- R packages have well-defined standard directory structures (DESCRIPTION, R/, src/, man/, etc.)
- Root directory should contain only essential package files, not development artifacts
- Benchmark scripts belong in `inst/benchmarks/` or `inst/scripts/`
- Log files should be in a separate `logs/` directory (excluded from package build)
- `.gitignore` needs comprehensive patterns for R package development with C++

**Primary recommendation:** Create a clean separation between package source files (in standard locations) and development artifacts (in dedicated directories), following R package conventions while maintaining the existing numeric file prefix system for R code organization.

## Standard Stack

The established structure for R packages with C++ extensions:

### Core R Package Structure
| Directory/File | Purpose | Why Standard |
|----------------|---------|--------------|
| `DESCRIPTION` | Package metadata, dependencies, license | Required for all R packages |
| `NAMESPACE` | Function exports/imports | Required for proper namespace management |
| `R/` | R source code files | Standard location for R code |
| `man/` | Documentation files (auto-generated) | Standard location for help files |
| `src/` | C/C++ source code | Required for packages with compiled code |
| `inst/` | Additional package files installed with package | For benchmarks, examples, data files |
| `tests/` | Test files (testthat) | Standard location for package tests |
| `vignettes/` | Package vignettes | Standard location for long-form documentation |

### Development Support Structure
| Directory | Purpose | When to Use |
|-----------|---------|-------------|
| `logs/` | Development log files | For session logs, debugging output |
| `.planning/` | Project planning documents | For GSD methodology tracking |
| `data-raw/` | Raw data processing scripts | For data preparation (excluded from build) |
| `figure/` | Generated figures | For documentation images, plots |

### Essential Root Files
| File | Required? | Purpose |
|------|-----------|---------|
| `README.md` | Recommended | Package overview, installation instructions |
| `LICENSE` or `LICENSE.md` | Required for CRAN | Package license (GPL-3 in this case) |
| `CLAUDE.md` | Project-specific | Claude Code instructions for this project |
| `.gitignore` | Required | Git ignore patterns |
| `.Rbuildignore` | Required | Files to exclude from package build |

**Installation:** No installation needed - this is about repository structure.

## Architecture Patterns

### Recommended Project Structure
```
RLDM/
├── DESCRIPTION                    # Package metadata
├── NAMESPACE                     # Exports/imports
├── README.md                     # Package overview
├── LICENSE.md                    # License file (GPL-3)
├── CLAUDE.md                     # Claude instructions
├── .gitignore                    # Git ignore patterns
├── .Rbuildignore                 # Build ignore patterns
├── R/                            # R source code (01_*, 02_*, etc.)
├── src/                          # C++ code (kf.cpp, rls_core.cpp, etc.)
├── man/                          # Documentation (auto-generated)
├── tests/                        # Testthat tests
├── vignettes/                    # Package vignettes
├── inst/                         # Installed files
│   ├── benchmarks/               # Benchmark scripts
│   ├── examples/                 # Example scripts (existing)
│   ├── extdata/                  # External data (existing)
│   └── include/                  # C++ headers (existing)
├── logs/                         # Development logs
│   ├── claude/                   # Claude session logs
│   └── smc/                      # SMC development logs
└── .planning/                    # Planning documents (existing)
```

### Pattern 1: R Package Root File Management
**What:** Keep root directory clean with only essential package files
**When to use:** All R package development
**Example:**
```bash
# Good root files:
DESCRIPTION
NAMESPACE
README.md
LICENSE.md
CLAUDE.md
.gitignore
.Rbuildignore

# Files to move elsewhere:
benchmark_*.R → inst/benchmarks/
*debug*.R → logs/debug/
*test*.R (not testthat) → logs/test_scripts/
*_LOG*.md → logs/smc/
log_*.md → logs/claude/
```

### Pattern 2: Development Artifact Organization
**What:** Separate development artifacts from package source
**When to use:** When package has debugging scripts, benchmarks, logs
**Example:**
```bash
# Create organized structure:
mkdir -p inst/benchmarks
mkdir -p logs/{claude,smc,debug,test_scripts}

# Move files:
mv benchmark_*.R inst/benchmarks/
mv *debug*.R logs/debug/
mv *test*.R logs/test_scripts/  # except tests/testthat/
mv *LOG*.md logs/smc/
mv log_*.md logs/claude/
```

### Anti-Patterns to Avoid
- **Root directory clutter:** Having development scripts in root makes package structure unclear
- **Mixed test locations:** Having some tests in `tests/` and others in root causes confusion
- **Missing LICENSE file:** CRAN requires explicit license file, not just DESCRIPTION entry
- **Incomplete .gitignore:** Missing patterns for R package build artifacts

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Git ignore patterns | Custom .gitignore | Standard R package .gitignore templates | Covers all R/C++ build artifacts, IDE files |
| License file | Write custom license | Use standard GPL-3 template | Ensures proper licensing, CRAN compliance |
| Benchmark organization | Ad-hoc benchmark locations | `inst/benchmarks/` convention | Standard location, installed with package |
| Log file management | Scattered log files | Dedicated `logs/` directory | Clear separation, easy cleanup |

**Key insight:** R package ecosystem has established conventions for all common development tasks. Deviating causes confusion for contributors and tools.

## Common Pitfalls

### Pitfall 1: Broken Path References
**What goes wrong:** Moving files breaks relative paths in scripts
**Why it happens:** Scripts use `../` or assume current directory
**How to avoid:**
1. Check scripts for hardcoded paths before moving
2. Use `here::here()` or `rprojroot::find_root()` for path resolution
3. Update any path references after moving files
**Warning signs:** Scripts failing with "file not found" errors after move

### Pitfall 2: Incomplete .gitignore
**What goes wrong:** Build artifacts, temporary files committed to git
**Why it happens:** Missing patterns for R/C++ development
**How to avoid:** Use comprehensive .gitignore covering:
- R build artifacts (`*.tar.gz`, `*.Rcheck/`)
- C++ build artifacts (`*.o`, `*.so`, `*.dll`)
- IDE files (`.Rproj.user`, `.DS_Store`)
- Session files (`.Rhistory`, `.Rdata`)
**Warning signs:** Binary files in git diff, large commits of build artifacts

### Pitfall 3: License File Issues
**What goes wrong:** Package fails CRAN checks due to missing license file
**Why it happens:** LICENSE only specified in DESCRIPTION, not as file
**How to avoid:** Create `LICENSE.md` or `LICENSE` file with full license text
**Warning signs:** `devtools::check()` warns about missing license file

### Pitfall 4: Test File Confusion
**What goes wrong:** Mixing `testthat` tests with development test scripts
**Why it happens:** Both use "test" prefix but serve different purposes
**How to avoid:**
- `tests/testthat/test-*.R` for package tests (run with `devtools::test()`)
- `logs/test_scripts/` for development/debugging tests
**Warning signs:** Test scripts in root that aren't part of testthat suite

## Code Examples

Verified patterns from R package conventions:

### Standard .gitignore for R Package with C++
```gitignore
# R package build artifacts
*.tar.gz
*.Rcheck/
RLDM.Rcheck/
.Rproj.user
.Rhistory
.Ruserdata
.RData
.Rapp.history

# C++ build artifacts
src/*.o
src/*.so
src/*.dll
src/*.dylib

# Documentation
docs/
pkgdown/

# IDE and OS files
.DS_Store
Thumbs.db
*.swp
*~

# Package-specific exclusions
logs/
figure/
data-raw/
inst/doc/
*.pdf
*.html
*.bak
todo.md
```

### Creating License File for GPL-3
```bash
# Create LICENSE.md from DESCRIPTION license
echo "## License" > LICENSE.md
echo "" >> LICENSE.md
echo "This package is licensed under the GPL-3 license." >> LICENSE.md
echo "See https://www.gnu.org/licenses/gpl-3.0.html for details." >> LICENSE.md
echo "" >> LICENSE.md
echo "Copyright (c) 2026 Wolfgang Scherrer, Bernd Funovits" >> LICENSE.md
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Ad-hoc repository structure | Standard R package layout | R 2.10+ | Consistent across all packages |
| Manual namespace management | roxygen2 with @export | R 3.0+ | Automated, less error-prone |
| Custom test organization | testthat in tests/testthat/ | 2012+ | Standardized testing framework |
| Mixed development/source files | Clear separation with inst/, logs/ | Best practice | Cleaner, more maintainable |

**Deprecated/outdated:**
- **Root-level test scripts:** Should be in `tests/testthat/` or `logs/test_scripts/`
- **Scattered log files:** Should be consolidated in `logs/` directory
- **Missing license file:** CRAN now requires explicit license file

## Open Questions

Things that couldn't be fully resolved:

1. **Benchmark script dependencies**
   - What we know: `benchmark_particle_filters.R` depends on RLDM package
   - What's unclear: Whether it has hardcoded paths to other scripts/files
   - Recommendation: Check script for `source()` calls or file reads before moving

2. **Existing inst/ directory usage**
   - What we know: `inst/examples/`, `inst/extdata/`, `inst/include/` exist
   - What's unclear: Whether these are actively used by package functions
   - Recommendation: Preserve existing inst/ structure, add benchmarks subdirectory

3. **.Rbuildignore vs .gitignore scope**
   - What we know: `.Rbuildignore` excludes from package build, `.gitignore` from git
   - What's unclear: Optimal division between the two files
   - Recommendation: Keep build-specific exclusions in `.Rbuildignore`, development artifacts in `.gitignore`

## Sources

### Primary (HIGH confidence)
- R Package Development (r-pkgs.org) - Standard package structure conventions
- Current repository analysis - Existing file patterns and organization
- R documentation - DESCRIPTION, NAMESPACE, directory structure requirements

### Secondary (MEDIUM confidence)
- R community best practices - Common patterns for R package development
- CRAN submission guidelines - License and file requirements

### Tertiary (LOW confidence)
- Web search attempted but encountered API errors - Would verify with official R documentation

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - Well-established R package conventions
- Architecture: HIGH - Clear patterns from repository analysis
- Pitfalls: HIGH - Common issues in R package development

**Research date:** 2026-01-24
**Valid until:** 2026-02-24 (30 days - stable conventions)

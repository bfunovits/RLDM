# Codebase Structure

**Analysis Date:** 2026-01-24

## Directory Layout

```
/media/bernd/nvme/r_projects/acad_RLDM/
├── R/                          # R source code (organized by numeric prefix)
│   ├── 01_representations_classes.R    # Model classes
│   ├── 02_templates.R                  # Parameter templates
│   ├── 03_properties_*.R              # Derived properties
│   ├── 04_timeseries_*.R              # Time series operations
│   ├── 05_estimation_*.R              # Estimation methods
│   ├── 06_visualization_*.R           # Plot methods
│   ├── 07_comparison_metrics.R        # Model comparison
│   ├── 08_utilities_*.R               # Utilities and data
│   ├── RcppExports.R                  # Auto-generated C++ bindings
│   └── RLDM-package.R                 # Package metadata
├── src/                        # C++ source code
│   ├── kf.cpp                  # Kalman filter implementation
│   ├── pf.cpp                  # Particle filter implementation
│   ├── rls_core.cpp            # Recursive least squares
│   ├── solve_fwd_bwd_cll.cpp   # Forward-backward solving
│   ├── RcppExports.cpp         # Auto-generated C++ exports
│   └── Makevars                # Build configuration
├── man/                        # Documentation (Rd files)
├── tests/                      # Test suite
│   └── testthat/               # testthat tests
├── vignettes/                  # Package vignettes
│   ├── 0_getting_started.Rmd   # Beginner tutorial
│   ├── 1_case_study.Rmd        # Practical workflow
│   └── 2_technical_reference.Rmd # Advanced reference
├── inst/                       # Package installation files
│   ├── examples/               # Example code
│   ├── extdata/                # External data
│   └── include/                # C++ headers
├── data/                       # Package datasets
├── data-raw/                   # Raw data processing scripts
├── docs/                       # Package website (pkgdown)
├── figure/                     # Generated figures
├── .github/                    # GitHub workflows
├── .planning/                  # Planning documents
└── [root config files]         # DESCRIPTION, NAMESPACE, etc.
```

## Directory Purposes

**R/:**
- Purpose: All R source code organized by workflow
- Contains: Function definitions, S3 methods, exports
- Key files: All `0[1-8]_*.R` files following numeric convention

**src/:**
- Purpose: Performance-critical C++ implementations
- Contains: Kalman filter, particle filter, RLS algorithms
- Key files: `kf.cpp`, `pf.cpp`, `rls_core.cpp`

**tests/:**
- Purpose: Unit and integration tests
- Contains: testthat test files
- Key files: `test-pfilter.R`, `test-est_ar.R`, `test-templates.R`

**vignettes/:**
- Purpose: User documentation and tutorials
- Contains: RMarkdown vignettes at different user levels
- Key files: `0_getting_started.Rmd`, `1_case_study.Rmd`, `2_technical_reference.Rmd`

**inst/:**
- Purpose: Files installed with package
- Contains: Examples, external data, C++ headers
- Key files: `inst/examples/`, `inst/include/`

**data/:**
- Purpose: Package datasets
- Contains: `.rda` files with example data
- Key files: `BQdata.rda`, `RSdata.rda` (Blanchard-Quah and Rudebusch-Svensson data)

## Key File Locations

**Entry Points:**
- `/media/bernd/nvme/r_projects/acad_RLDM/R/RLDM-package.R`: Package metadata and imports
- `/media/bernd/nvme/r_projects/acad_RLDM/R/01_representations_classes.R`: Model constructors

**Configuration:**
- `/media/bernd/nvme/r_projects/acad_RLDM/DESCRIPTION`: Package metadata and dependencies
- `/media/bernd/nvme/r_projects/acad_RLDM/NAMESPACE`: Exported functions and S3 methods
- `/media/bernd/nvme/r_projects/acad_RLDM/src/Makevars`: C++ build configuration

**Core Logic:**
- `/media/bernd/nvme/r_projects/acad_RLDM/R/05_estimation_*.R`: All estimation methods
- `/media/bernd/nvme/r_projects/acad_RLDM/src/kf.cpp`: Kalman filter implementation
- `/media/bernd/nvme/r_projects/acad_RLDM/src/pf.cpp`: Particle filter implementation

**Testing:**
- `/media/bernd/nvme/r_projects/acad_RLDM/tests/testthat/test-pfilter.R`: Particle filter tests
- `/media/bernd/nvme/r_projects/acad_RLDM/tests/testthat/test-est_ar.R`: AR estimation tests

## Naming Conventions

**Files:**
- Pattern: `NN_category_description.R` where NN is 01-08 prefix
- Example: `05_estimation_particle.R`, `03_properties_autocov.R`

**Directories:**
- Pattern: Lowercase with descriptive names
- Example: `tests/`, `vignettes/`, `data-raw/`

**Functions:**
- Pattern: `snake_case` for R functions, `CamelCase` for classes
- Example: `est_arma_hrk()`, `armamod` class

**C++ Files:**
- Pattern: Descriptive names with `.cpp` extension
- Example: `kf.cpp`, `pf.cpp`, `rls_core.cpp`

## Where to Add New Code

**New Estimation Method:**
- Primary code: `/media/bernd/nvme/r_projects/acad_RLDM/R/05_estimation_NEWMETHOD.R`
- Tests: `/media/bernd/nvme/r_projects/acad_RLDM/tests/testthat/test-NEWMETHOD.R`
- Follow existing patterns in `05_estimation_arma_hrk.R` or `05_estimation_particle.R`

**New Model Property:**
- Primary code: `/media/bernd/nvme/r_projects/acad_RLDM/R/03_properties_NEWPROPERTY.R`
- Tests: Add to existing property test files
- Follow patterns in `03_properties_autocov.R`

**New C++ Algorithm:**
- Implementation: `/media/bernd/nvme/r_projects/acad_RLDM/src/NEWALGORITHM.cpp`
- R binding: Add to existing R files or create new `05_estimation_*.R`
- Run `Rcpp::compileAttributes()` after adding `[[Rcpp::export]]`

**New Visualization:**
- Primary code: `/media/bernd/nvme/r_projects/acad_RLDM/R/06_visualization_NEWPLOT.R`
- Follow S3 method patterns in `06_visualization_plot.R`

**Utilities:**
- Shared helpers: `/media/bernd/nvme/r_projects/acad_RLDM/R/08_utilities_NEW.R`
- Data documentation: `/media/bernd/nvme/r_projects/acad_RLDM/R/08_utilities_data.R`

## Special Directories

**.planning/:**
- Purpose: Planning and analysis documents
- Contains: Codebase mapping outputs
- Generated: Yes (by Claude agent)
- Committed: Yes

**.github/:**
- Purpose: GitHub Actions workflows
- Contains: CI/CD configuration
- Generated: Partially
- Committed: Yes

**data-raw/:**
- Purpose: Scripts to process raw data into package datasets
- Contains: Data processing R scripts
- Generated: No
- Committed: Yes

**figure/:**
- Purpose: Generated figures for documentation
- Contains: Plot outputs
- Generated: Yes
- Committed: Yes

---

*Structure analysis: 2026-01-24*
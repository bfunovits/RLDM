# Architecture

**Analysis Date:** 2026-01-24

## Pattern Overview

**Overall:** S3 object-oriented system with functional programming patterns

**Key Characteristics:**
- S3 class hierarchy with `rldm` as base class
- Three main model representations: `armamod`, `stspmod`, `rmfdmod`
- Performance-critical operations implemented in C++ via Rcpp/RcppArmadillo
- Template-based parameter estimation system
- Numeric prefix organization for R files by workflow purpose

## Layers

**Model Layer:**
- Purpose: Define model classes and constructors
- Location: `/media/bernd/nvme/r_projects/acad_RLDM/R/01_representations_classes.R`
- Contains: `armamod()`, `stspmod()`, `rmfdmod()` constructors
- Depends on: `rationalmatrices` package for LMFD/RMFD objects
- Used by: All estimation and analysis functions

**Template Layer:**
- Purpose: Parameter templates for estimation
- Location: `/media/bernd/nvme/r_projects/acad_RLDM/R/02_templates.R`
- Contains: `tmpl_*` functions, `fill_template()`, `extract_theta()`
- Depends on: Model layer
- Used by: Estimation functions

**Properties Layer:**
- Purpose: Compute derived model properties
- Location: `/media/bernd/nvme/r_projects/acad_RLDM/R/03_properties_*.R`
- Contains: `autocov()`, `freqresp()`, `impresp()`, `spectrald()`, `poles()`
- Depends on: Model layer
- Used by: Analysis and visualization functions

**Time Series Layer:**
- Purpose: Time series operations
- Location: `/media/bernd/nvme/r_projects/acad_RLDM/R/04_timeseries_*.R`
- Contains: `solve_de()`, `sim()`, prediction functions
- Depends on: Model layer, C++ implementations
- Used by: Estimation and simulation workflows

**Estimation Layer:**
- Purpose: Model estimation methods
- Location: `/media/bernd/nvme/r_projects/acad_RLDM/R/05_estimation_*.R`
- Contains: AR (OLS/YW/DLW), ARMA (HRK), subspace, likelihood, RLS, particle filters
- Depends on: Template layer, C++ implementations
- Used by: End-user estimation workflows

**C++ Performance Layer:**
- Purpose: High-performance numerical computations
- Location: `/media/bernd/nvme/r_projects/acad_RLDM/src/`
- Contains: Kalman filter, particle filter, recursive least squares, forward-backward solving
- Depends on: Rcpp, RcppArmadillo, LAPACK/BLAS
- Used by: R estimation and time series functions

## Data Flow

**Model Estimation Workflow:**

1. Model construction: `armamod()` or `stspmod()`
2. Template creation: `tmpl_*()` functions
3. Parameter estimation: `est_*()` functions (call C++ implementations)
4. Model validation: `compare_estimates()`, diagnostic functions
5. Analysis: `autocov()`, `spectrald()`, `impresp()`

**State Space Filtering Flow:**

1. Model: `stspmod()` object with A, B, C, D, Q, R, S matrices
2. Data: Time series matrix `y`
3. Filter: `pfilter()` or `kf()` (dispatches to C++ implementations)
4. Output: Filtered states, predictions, likelihood

**State Management:**
- Models are immutable S3 objects
- Estimation returns new model objects with fitted parameters
- No global state; pure functional approach

## Key Abstractions

**Model Classes:**
- Purpose: Represent different model parameterizations
- Examples: `/media/bernd/nvme/r_projects/acad_RLDM/R/01_representations_classes.R`
- Pattern: S3 classes inheriting from `rldm` base class

**Parameter Templates:**
- Purpose: Map deep parameters to model matrices
- Examples: `tmpl_stsp_full()`, `tmpl_arma_echelon()`
- Pattern: Linear mapping `vec(SYS) = H_SYS * theta + h_SYS`

**Filter Objects:**
- Purpose: Represent filtering results
- Examples: `pfilter()` output class
- Pattern: List with standardized components (filtered_states, predicted_states, etc.)

## Entry Points

**Package Loading:**
- Location: `/media/bernd/nvme/r_projects/acad_RLDM/R/RLDM-package.R`
- Triggers: `library(RLDM)` or `devtools::load_all()`
- Responsibilities: Register C++ functions, import dependencies

**Model Construction:**
- Location: `/media/bernd/nvme/r_projects/acad_RLDM/R/01_representations_classes.R`
- Triggers: User calls `armamod()`, `stspmod()`, `rmfdmod()`
- Responsibilities: Validate inputs, create S3 objects

**Estimation:**
- Location: `/media/bernd/nvme/r_projects/acad_RLDM/R/05_estimation_*.R`
- Triggers: User calls `est_*()` functions
- Responsibilities: Parameter estimation, model fitting

## Error Handling

**Strategy:** Defensive programming with explicit validation

**Patterns:**
- Input validation in constructors and public functions
- `stop()` with informative messages for invalid inputs
- `tryCatch()` for numerical errors in C++ calls
- Return values include convergence status for iterative methods

## Cross-Cutting Concerns

**Logging:** Minimal; uses R's warning/error system
**Validation:** Input validation in all public functions
**Authentication:** Not applicable (local R package)

---

*Architecture analysis: 2026-01-24*
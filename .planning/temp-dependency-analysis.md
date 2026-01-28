# Dependency Usage Analysis for RLDM Package

## Current DESCRIPTION Dependencies

### Depends
- R (>= 2.10)
- rationalmatrices

### Imports
- dplyr
- MASS
- purrr
- Rcpp
- Rdpack
- tibble
- QZ

### Suggests
- DiagrammeR
- kableExtra
- knitr
- lubridate
- mvbutils
- printr
- rmarkdown
- testthat (>= 2.1.0)
- xts

### LinkingTo
- Rcpp
- RcppArmadillo
- rationalmatrices

### Remotes
- bfunovits/rationalmatrices

## Usage Analysis

### ✅ Used Dependencies (Confirmed Usage)

#### 1. rationalmatrices (Depends)
- **Usage**: Extensive usage throughout codebase
- **Files**: Multiple R files reference rationalmatrices functions
- **Functions**: `lmfd()`, `rmfd()`, `stsp()`, `zvalues()`, `pseries()`, `nu2basis()`, `plot_3D()`, etc.
- **Status**: Essential dependency, correctly in Depends

#### 2. Rcpp (Imports + LinkingTo)
- **Usage**: For `evalCpp()` and `sourceCpp()`
- **Files**: `R/RLDM-package.R`, `R/05_estimation_rls.R`
- **Functions**: `evalCpp`, `sourceCpp`
- **Status**: Required for C++ integration

#### 3. dplyr (Imports)
- **Usage**: `bind_rows()` function
- **Files**: `R/05_estimation_arma_hrk.R` (lines 306, 367, 527, 576)
- **Status**: Used

#### 4. purrr (Imports)
- **Usage**: `map_dbl()` function
- **Files**: `R/05_estimation_arma_hrk.R` (lines 171, 299, 301, 362, 520, 523, 571)
- **Status**: Used

#### 5. tibble (Imports)
- **Usage**: `tibble()` function
- **Files**: `R/05_estimation_arma_hrk.R` (lines 175, 181, 307, 368, 528, 577)
- **Status**: Used

#### 6. MASS (Imports)
- **Usage**: `ginv()` function
- **Files**: `R/04_timeseries_solve.R` (lines 518, 792), `R/05_estimation_rls.R` (line 76)
- **Status**: Used

#### 7. QZ (Imports)
- **Usage**: `qz.dgges()` and `qz.dtgsen()` functions
- **Files**: `R/05_estimation_subspace.R` (lines 1514, 1533)
- **Status**: Used

#### 8. Rdpack (Imports)
- **Usage**: Documentation macros (`RdMacros: Rdpack` in DESCRIPTION)
- **Files**: Imported in `R/RLDM-package.R` but `reprompt()` not directly called
- **Status**: Required for documentation building (RdMacros field)

#### 9. magrittr (Imports via importFrom)
- **Usage**: `%>%` operator re-exported
- **Files**: `R/05_estimation_rls.R` imports it, used throughout codebase with `%>%`
- **Status**: Used

### ✅ Base Packages (Not in Imports but used)
- **stats**: `lsfit()`, `rnorm()`, `acf()`, `optim()`, etc.
- **graphics**: Plotting functions
- **grDevices**: Color functions

### 🔍 Suggests Dependencies (Optional)

#### 1. testthat (>= 2.1.0)
- **Usage**: Testing framework
- **Files**: Test files in `tests/testthat/`
- **Status**: Used for testing

#### 2. knitr
- **Usage**: Vignette building
- **Files**: All vignettes (`vignettes/*.Rmd`)
- **Status**: Used for vignettes

#### 3. rmarkdown
- **Usage**: Vignette building
- **Files**: Likely used with knitr
- **Status**: Presumed used for vignettes

#### 4. Other Suggests (Unconfirmed Usage)
- **DiagrammeR**, **kableExtra**, **lubridate**, **mvbutils**, **printr**, **xts**
- **Status**: Need to check vignettes/tests for usage

## Findings

1. **All current Imports are actually used** in the codebase
2. **rationalmatrices** is correctly in Depends (essential dependency)
3. **Rdpack** is needed for documentation macros even though `reprompt()` not directly called
4. **No unused imports found** in current DESCRIPTION
5. **Version constraints missing** for most dependencies (need to add in Task 2)
6. **R version constraint** is very old (>= 2.10) - should update to >= 3.6.0

## Recommendations for Task 2

1. Add version constraints to all dependencies
2. Update R version constraint from >= 2.10 to >= 3.6.0
3. Keep all current imports (all are used)
4. Add conservative version bounds based on current usage patterns
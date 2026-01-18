# R/ Folder Organization

## Numeric Functional System

RLDM uses a numeric prefix system organized by **purpose/workflow**, not
object type. This organization scheme emphasizes that mathematical
operations (estimation, properties, visualization) are conceptually
independent of model representation.

### Core Philosophy

- **Representations** (lmfd/armamod, stsp/stspmod, rmfd/rmfdmod) are
  just different encodings of the same rational process
- **Operations** (estimation, properties, prediction) work with all
  representations via S3 dispatch
- **Templates** map deep (free) parameters to model parameters using
  affine transformations
- Users think in terms of workflow (“I want to estimate”) not
  representation

### File Organization by Numeric Prefix

| Prefix | Category        | Purpose                                                                         | Key Files                      |
|--------|-----------------|---------------------------------------------------------------------------------|--------------------------------|
| **01** | Representations | Model class definitions and representation-specific methods                     | `01_representations_classes.R` |
| **02** | Templates       | Parameter templates for each model type (ARMA, RMFD, STSP, local, structural)   | `02_templates.R`               |
| **03** | Properties      | Derived properties (ACF, spectral, impulse response, frequency response, poles) | `03_properties_*.R`            |
| **04** | Timeseries      | Time series operations (solve, simulate, predict/forecast)                      | `04_timeseries_*.R`            |
| **05** | Estimation      | All estimation methods (AR, ARMA, subspace, likelihood, RLS)                    | `05_estimation_*.R`            |
| **06** | Visualization   | Plot methods for properties and predictions                                     | `06_visualization_*.R`         |
| **07** | Comparison      | Model comparison, diagnostics, and evaluation metrics                           | `07_comparison_metrics.R`      |
| **08** | Utilities       | Supporting functions (print, str, data, package docs)                           | `08_utilities_*.R`             |
| **99** | Auto-generated  | Rcpp bindings (DO NOT EDIT MANUALLY)                                            | `RcppExports.R`                |

### Load Order

Files load alphabetically. Dependencies must come before dependents:

1.  **01** - Core classes (needed by everything)
2.  **02** - Templates (needed by estimation and model creation)
3.  **03-04** - Properties and operations (use classes)
4.  **05-07** - Higher-level functionality (use estimation/properties)
5.  **08** - Utilities (minimal dependencies)
6.  **99** - Auto-generated (loaded last)

### File Naming Convention

`NN_category_subcategory.R` where:

- **NN** = Two-digit numeric prefix
- **category** = Functional area (single word: representations,
  templates, properties, etc.)
- **subcategory** = Specific aspect (method name, model type, or
  descriptor)

**Examples:** - `02_templates.R` - All template-related functions -
`03_properties_autocov.R` - Autocovariance computation -
`04_timeseries_solve.R` - Difference equation solving -
`05_estimation_ar.R` - AR model estimation - `06_visualization_plot.R` -
General plotting methods

### S3 Method Organization

**Key principle:** Methods are grouped by operation (noun), not by
class.

**Example:** `03_properties_autocov.R` contains: -
[`autocov()`](https://bfunovits.github.io/RLDM/reference/autocov.md)
generic -
[`autocov.armamod()`](https://bfunovits.github.io/RLDM/reference/autocov.md) -
via Yule-Walker equations -
[`autocov.stspmod()`](https://bfunovits.github.io/RLDM/reference/autocov.md) -
via Lyapunov equation - `autocov.rmfdmod()` - experimental -
[`autocov.default()`](https://bfunovits.github.io/RLDM/reference/autocov.md) -
from time series data

This emphasizes: **Operation is primary, representation is
implementation detail**.

### Roxygen2 Integration

- Each file has a header comment block (see example below)
- Cross-references via `@rdname` for related functions
- `@include` directives specify file loading order for documentation
- Roxygen tags group documentation by functional area
  (`@name model_structures`, etc.)

### Template Header Format

Each file should begin with a structured header comment:

``` r
# ====================================================================
# Category Name
#
# Purpose: [Clear, single-sentence description]
#   - Function 1 (brief role)
#   - Function 2 (brief role)
#
# Key Operations: [What this file does in the workflow]
# Representations: [Which model types, if relevant]
# Dependencies: [Which files this imports/depends on]
# Related Files: [Cross-references to related functionality]
# ====================================================================
```

### Transition from Alphabetical System

The previous alphabetical prefix system (aa\_, ab\_, aca\_, etc.) has
been reorganized:

**Old → New:** - `aa_*` → `01_representations_*` - `ab_*` →
`02_templates` - `aca_*, acb_*, ...` → `03_properties_*` -
`ada_*, adb_*` → `03_properties_poles`, `04_timeseries_*` -
`aeb_*, aec_*` → `08_utilities_*` - `aeaa_*` → `06_visualization_*` -
`afa*` → `05_estimation_*` - `ag_*` → `07_comparison_*`

### Verification Commands

After reorganization, verify everything still works:

``` r
# Load package in development mode
devtools::load_all()

# Generate documentation (respects new file organization)
devtools::document()

# Check for any warnings/errors
devtools::check()

# Run test suite
devtools::test()
```

### Guidelines for Adding New Files

1.  **Choose the right prefix** based on what the code does
2.  **Use descriptive names** - subcategory should clarify function
3.  **Add header comment** explaining purpose and dependencies
4.  **Check load order** - ensure dependencies load first (automatic via
    naming)
5.  **Use @rdname** for related functions in documentation
6.  **Add to tests** if implementing new functionality

### Notes

- All 23 R files (excluding RcppExports.R) fit within categories 01-08
- Load order is automatically handled via alphabetical sorting of
  filenames
- No subdirectories are used (R package requirement)
- Documentation is unified across files via Roxygen2 tags

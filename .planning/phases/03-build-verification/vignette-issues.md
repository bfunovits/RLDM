# Vignette Build Failure Analysis

**Date:** 2026-01-28
**Source:** `check-initial.log` from Plan 03-01 Task 2
**Check command:** `devtools::check(document = TRUE, manual = TRUE, cran = FALSE)`

## Summary

All three vignettes fail to build:
1. **0_getting_started.Rmd** - Code error (dimension mismatch)
2. **1_case_study.Rmd** - Code error (coercion to polm object)
3. **2_technical_reference.Rmd** - System dependency (pandoc-citeproc missing)

## Detailed Analysis

### 1. 0_getting_started.Rmd

**Error:**
```
Error in `pred_ss$yhat[plot_range[-1], 1]`:
! incorrect number of dimensions
```

**Location:** Lines 308-322 (unnamed-chunk-15)

**Root cause:** Code error - dimension mismatch when accessing array elements. The code expects `pred_ss$yhat` to have certain dimensions but it doesn't match.

**Category:** Code issue (internal)

**Files affected:** `vignettes/0_getting_started.Rmd`

**Recommended action:** Fix the dimension access in the vignette code. Likely need to check structure of `pred_ss$yhat` and adjust indexing.

### 2. 1_case_study.Rmd

**Error:**
```
Error in `lmfd()`:
! could not coerce "a" to a "polm" object!
```

**Location:** Lines 305-311 (unnamed-chunk-10)

**Root cause:** Code error - `rationalmatrices::lmfd()` function failing to convert an object to polynomial matrix (polm) format. The object `models$SSECF$sys` likely has incorrect structure.

**Category:** Code issue (internal)

**Files affected:** `vignettes/1_case_study.Rmd`

**Recommended action:** Debug the `models$SSECF$sys` object structure. May need to fix data preparation or use different conversion method.

### 3. 2_technical_reference.Rmd

**Error:**
```
Error running filter pandoc-citeproc:
Could not find executable pandoc-citeproc
Error: processing vignette '2_technical_reference.Rmd' failed with diagnostics:
pandoc document conversion failed with error 83
```

**Location:** During pandoc document conversion

**Root cause:** System dependency missing - `pandoc-citeproc` executable not found on system.

**Category:** System dependency (external)

**Files affected:** `vignettes/2_technical_reference.Rmd`

**Recommended action:** Install `pandoc-citeproc` system package. On Ubuntu/Debian: `sudo apt-get install pandoc-citeproc`. This is a system-level dependency, not a code issue.

## Impact on BUILD-01 (0 errors)

1. **System dependency (2_technical_reference.Rmd):** External issue, not blocking BUILD-01 if we treat it as system configuration issue.
2. **Code errors (0_getting_started.Rmd, 1_case_study.Rmd):** Internal issues that would prevent BUILD-01 (0 errors) if vignettes are required to build.

## Phase 3 Context Decisions

From 03-CONTEXT.md:
- "Vignettes must build without errors"
- "Claude's discretion: Specific warning patterns to ignore (non‑ASCII, time verification)"

## Recommendations for Plan 03-02

### Option A: Fix code issues, document system dependency
1. **Fix 0_getting_started.Rmd dimension error**
   - Analyze `pred_ss$yhat` structure
   - Correct indexing or data preparation
2. **Fix 1_case_study.Rmd coercion error**
   - Debug `models$SSECF$sys` object
   - Ensure proper conversion to polm format
3. **Document pandoc-citeproc as system requirement**
   - Add to package documentation
   - Note in INSTALL instructions

### Option B: Skip vignette building for BUILD-01
1. Use `--ignore-vignettes` flag for check
2. Document vignette issues as known limitations
3. Address in later phase (Phase 4: Documentation)

### Recommended: Option A
Since Phase 3 CONTEXT states "Vignettes must build without errors", we should attempt to fix the code issues. The system dependency (pandoc-citeproc) can be documented as external requirement.

## Action Plan for Plan 03-02

1. **Task:** Fix dimension error in `0_getting_started.Rmd`
   - File: `vignettes/0_getting_started.Rmd` lines 308-322
   - Action: Debug `pred_ss$yhat` structure and fix indexing

2. **Task:** Fix coercion error in `1_case_study.Rmd`
   - File: `vignettes/1_case_study.Rmd` lines 305-311
   - Action: Debug `models$SSECF$sys` and fix `lmfd()` conversion

3. **Task:** Document system dependency
   - File: `INSTALL` or `README.md`
   - Action: Add note about `pandoc-citeproc` requirement

4. **Task:** Verify fixes
   - Run `devtools::check(vignettes = TRUE)` after fixes
   - Ensure only system dependency error remains

## Files to Modify in Plan 03-02

- `vignettes/0_getting_started.Rmd`
- `vignettes/1_case_study.Rmd`
- Package documentation files (INSTALL/README)

## Success Criteria for Plan 03-02

- [ ] `0_getting_started.Rmd` builds without dimension errors
- [ ] `1_case_study.Rmd` builds without coercion errors
- [ ] System dependency documented
- [ ] `devtools::check(vignettes = TRUE)` shows only pandoc-citeproc error (external)
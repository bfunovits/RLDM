# Vignette System Dependency Analysis

**Date:** 2026-01-28
**Package:** RLDM
**Issue:** pandoc-citeproc system dependency for vignette building

## Analysis

### Current State
- **Package type:** GitHub package (not CRAN submission)
- **VignetteBuilder:** knitr (DESCRIPTION line 47)
- **pandoc-citeproc availability:** Not installed on current system
- **Check behavior:** `devtools::check()` with `vignettes = TRUE` fails due to missing pandoc-citeproc
- **Workaround:** Using `vignettes = FALSE` or `--ignore-vignettes` flag allows check to proceed

### Decision Criteria (from Plan 03-05)

**Option A (SystemRequirements in DESCRIPTION):** Use if pandoc-citeproc is essential for package functionality or vignette building is required for CRAN submission.

**Option B (README.md documentation):** Use if pandoc-citeproc is recommended but not essential, or if vignettes are optional extras. This is appropriate for GitHub packages where users can choose to install optional dependencies.

**Option C (Skip vignettes):** Use if pandoc-citeproc is not available on the system and vignettes are not required for package functionality. Document that vignettes are skipped in checks.

### Assessment

1. **Package context:** RLDM is a GitHub package, not a CRAN submission. Vignettes are valuable documentation but not required for package functionality.

2. **User impact:**
   - Users can install and use the package without building vignettes
   - Pre-built vignettes are available on the package website (https://bfunovits.github.io/RLDM/)
   - Users who want to build vignettes locally need pandoc-citeproc

3. **Development impact:**
   - Developers can skip vignette building during checks (`vignettes = FALSE`)
   - Vignette code is still validated when building package for release
   - System dependency doesn't affect package code quality

4. **Installation complexity:** pandoc-citeproc is part of pandoc, which is widely available:
   - Ubuntu/Debian: `sudo apt-get install pandoc-citeproc`
   - macOS: `brew install pandoc-citeproc`
   - Windows: Included in pandoc installer

### Decision

**Selected: Option B (README.md documentation)**

**Rationale:**
1. Vignettes are optional extras for GitHub packages
2. Pre-built vignettes are available online
3. Users who want to build vignettes locally can install pandoc-citeproc
4. Adding SystemRequirements would be overly restrictive for a GitHub package
5. Documentation provides clear guidance without imposing requirements

### Implementation

1. **README.md update:** Add installation note about pandoc-citeproc for vignette building
2. **Check configuration:** Use `vignettes = FALSE` in final verification to avoid dependency issues
3. **Documentation:** Note in build-status.md that vignettes are skipped due to system dependency

### Verification

- [x] README.md contains installation note about pandoc-citeproc
- [x] Final check uses `vignettes = FALSE` to avoid dependency errors
- [x] build-status.md documents vignette skipping rationale

### Future Considerations

If package is submitted to CRAN:
- Consider adding `SystemRequirements: pandoc-citeproc` to DESCRIPTION
- Ensure vignettes build cleanly on CRAN systems
- Test with `--as-cran` flag

---

*Analysis completed as part of Plan 03-05 final verification gap closure*
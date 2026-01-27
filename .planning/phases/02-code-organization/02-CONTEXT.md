# Phase 2: Code Organization - Context

**Gathered:** 2026-01-28
**Status:** Ready for planning

<domain>
## Phase Boundary

Review and organize C++ and R code for maintainability. Clean up src/ C++ code, ensure R file numeric prefix consistency, organize S3 methods properly, and review DESCRIPTION dependencies. Creating new functionality or changing API contracts is out of scope—this phase focuses on organizing existing code.

</domain>

<decisions>
## Implementation Decisions

### C++ Code Standards
- Doxygen-style documentation for C++ functions (full API documentation with @param, @return, etc.)
- Adopt Google C++ Style Guide for consistency
- Prioritize both code structure and variable naming cleanup (minimize cognitive overhead)
- No additional unit tests for C++ functions beyond existing R-level tests

### R File Consistency
- Files not following numeric prefix convention: handle case-by-case, preference for leaving as exceptions when justified
- Hybrid category approach for files with multiple purposes (create appropriate categorization)
- No README in R/ directory; add information about numeric prefix system to parent directory README for GitHub users

### S3 Method Organization
- Claude's discretion for optimal organization strategy
- Claude's discretion for S3 registration approach (explicit vs Roxygen)
- Claude's discretion for method inheritance patterns (rldm parent class handling)
- Claude's discretion for method registry documentation needs

### Dependency Review
- Keep dependencies used anywhere in package (including internal functions)
- Rationalmatrices dependency stays as Remotes (GitHub-only installation)
- Claude's discretion on version constraints and unused import cleanup

### Claude's Discretion
- Determine appropriate adjustments to numeric prefix categories
- Choose optimal S3 method organization strategy (by class vs by function)
- Decide on S3 registration approach and inheritance patterns
- Assess need for method registry documentation
- Review and potentially add version constraints for dependencies
- Identify and clean up unused imports if appropriate

</decisions>

<specifics>
## Specific Ideas

- Package quality is critical: RLDM is foundational for other packages (svarmawhf) and depends on rationalmatrices
- "Little cognitive overhead for humans" - prioritize clean, understandable code structure
- Information about numeric prefix system should be in parent directory README (not R/README.md) for GitHub users

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope.

</deferred>

---

*Phase: 02-code-organization*
*Context gathered: 2026-01-28*
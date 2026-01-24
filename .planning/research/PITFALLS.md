# Domain Pitfalls

**Domain:** R package maintenance, documentation, and build verification
**Project:** RLDM Package Cleanup & Streamlining
**Researched:** 2026-01-24

## Critical Pitfalls

Mistakes that cause rewrites or major issues.

### Pitfall 1: Breaking Backward Compatibility
**What goes wrong:** Changes to function interfaces break existing user code
**Why it happens:** Refactoring without considering existing usage
**Consequences:** Users can't upgrade, issue reports, loss of trust
**Prevention:**
- Use `lifecycle` package to deprecate functions properly
- Maintain exported API stability
- Test with `revdepcheck` before major changes
**Detection:** `revdepcheck::revdep_check()` shows broken dependencies

### Pitfall 2: C++ ABI Incompatibility
**What goes wrong:** C++ code compiled with different standards breaks on user systems
**Why it happens:** Inconsistent compiler flags or C++ standard settings
**Consequences:** Installation failures, cryptic linking errors
**Prevention:**
- Use consistent `CXX_STD = CXX11` in `src/Makevars`
- Test on multiple platforms with `rhub`
- Document system requirements clearly
**Detection:** `rhub::check()` on multiple platforms reveals ABI issues

### Pitfall 3: Documentation Drift
**What goes wrong:** Documentation doesn't match actual function behavior
**Why it happens:** Code changes without updating documentation
**Consequences:** User confusion, incorrect usage, support burden
**Prevention:**
- Documentation-first development
- Regular `devtools::document()` runs
- Examples in documentation must run without error
**Detection:** `devtools::check()` warns about mismatched arguments

### Pitfall 4: Namespace Conflicts
**What goes wrong:** Functions mask or are masked by other packages
**Why it happens:** Common function names, improper imports
**Consequences:** Unpredictable behavior, hard-to-debug issues
**Prevention:**
- Use `@importFrom` selectively in roxygen2
- Prefix internal functions with `.`
- Use `pkg::fun()` syntax in code
**Detection:** `devtools::check()` reports namespace issues

## Moderate Pitfalls

Mistakes that cause delays or technical debt.

### Pitfall 1: Test Fragility
**What goes wrong:** Tests fail randomly or on specific platforms
**Why it happens:** Non-deterministic tests, platform-specific assumptions
**Consequences:** False CI failures, ignored test suites
**Prevention:**
- Use `set.seed()` for stochastic tests
- Avoid testing exact floating-point equality
- Use `testthat::skip_on_*()` for platform-specific skips
**Detection:** Intermittent CI failures, platform-specific test failures

### Pitfall 2: Vignette Build Failures
**What goes wrong:** Vignettes don't build during `devtools::check()`
**Why it happens:** Missing dependencies, long runtime, side effects
**Consequences:** CRAN rejection, incomplete documentation
**Prevention:**
- All vignette dependencies in `Suggests`
- Use `eval = FALSE` for slow code chunks
- Test vignette building locally
**Detection:** `devtools::check()` fails with vignette errors

### Pitfall 3: Over-Engineering Cleanup
**What goes wrong:** Spending too much time on perfect organization vs. functional cleanup
**Why it happens:** Perfectionism, unclear priorities
**Consequences:** Project stalls, diminishing returns
**Prevention:**
- Focus on `devtools::check()` passing first
- Incremental improvements
- Document what "good enough" looks like
**Detection:** Spending days on minor organizational details

### Pitfall 4: Ignoring Warnings
**What goes wrong:** Treating `devtools::check()` warnings as acceptable
**Why it happens:** "It still works" mentality, time pressure
**Consequences:** Hidden bugs, future breakage, CRAN issues
**Prevention:**
- Address all warnings systematically
- Document intentional warnings with `# nocov` or comments
- Treat warnings as errors in CI
**Detection:** Growing warning count in `devtools::check()` output

## Minor Pitfalls

Mistakes that cause annoyance but are fixable.

### Pitfall 1: Inconsistent Code Style
**What goes wrong:** Hard to read, maintain, and review code
**Why it happens:** Multiple contributors, no style enforcement
**Prevention:** Use `lintr` with consistent `.lintr` configuration

### Pitfall 2: Poor Commit Messages
**What goes wrong:** Hard to understand history, blame, and changes
**Why it happens:** Quick commits without clear messaging
**Prevention:** Follow conventional commit messages, explain "why" not just "what"

### Pitfall 3: Missing Examples
**What goes wrong:** Users can't understand how to use functions
**Why it happens:** Documentation written as afterthought
**Prevention:** Require `@examples` section for all exported functions

### Pitfall 4: Unorganized Root Directory
**What goes wrong:** Cluttered, unprofessional appearance
**Why it happens:** Temporary files accumulate
**Prevention:** Regular cleanup, proper `.Rbuildignore`

## Phase-Specific Warnings

| Phase Topic | Likely Pitfall | Mitigation |
|-------------|---------------|------------|
| Repository cleanup | Breaking file references | Update all `source()`, `load()`, `readRDS()` calls after moving files |
| Documentation audit | Missing `@export` tags | Use `roxygen2::roxygenise()` to see what gets documented |
| C++ code review | Breaking Rcpp interfaces | Test after each change with `Rcpp::compileAttributes()` |
| Test organization | Breaking existing tests | Run `devtools::test()` after moving test files |
| CI/CD setup | Complex GitHub Actions | Start with basic R-CMD-check workflow, expand gradually |
| Benchmark integration | Side effects in package | Keep benchmarks separate from package code |

## RLDM-Specific Pitfalls

### Particle Filter Integration
**Risk:** Particle filter benchmarks may have side effects or complex dependencies
**Mitigation:** Keep benchmarking code separate from package core; use `inst/benchmarks/`

### rationalmatrices Dependency
**Risk:** GitHub dependency may cause installation issues
**Mitigation:** Document installation clearly; consider conditional installation in tests

### Numeric Prefix System
**Risk:** May confuse new contributors unfamiliar with the system
**Mitigation:** Document the system in `R/README.md`; maintain consistency

### S3 Method Proliferation
**Risk:** Many S3 methods may be poorly organized
**Mitigation:** Audit all S3 generics and methods; group by generic

## Prevention Checklist

Before any major change:

- [ ] Run `devtools::check()` to establish baseline
- [ ] Run `devtools::test()` to ensure tests pass
- [ ] Check `revdepcheck` if package has dependencies
- [ ] Test on multiple platforms if possible
- [ ] Document the change and its rationale
- [ ] Update examples and vignettes if affected

After changes:

- [ ] Run `devtools::check()` again
- [ ] Verify examples still work
- [ ] Ensure vignettes build
- [ ] Update documentation if needed

## Recovery Strategies

When a pitfall occurs:

1. **Immediate rollback** - Use git to revert problematic changes
2. **Isolate the issue** - Create minimal reproducible example
3. **Fix systematically** - Address root cause, not just symptoms
4. **Add prevention** - Update processes to avoid recurrence
5. **Document the lesson** - Add to this pitfalls document

## Sources

**Based on analysis of:**
- Common R package maintenance issues
- CRAN package rejection reasons
- Experience with package refactoring
- RLDM package specific characteristics

**Confidence:** HIGH (pitfalls are well-documented in R community)
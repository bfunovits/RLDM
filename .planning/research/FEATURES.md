# Feature Landscape

**Domain:** R package maintenance, documentation, and build verification
**Project:** RLDM Package Cleanup & Streamlining
**Researched:** 2026-01-24

## Table Stakes

Features users expect. Missing = package feels incomplete or unprofessional.

| Feature | Why Expected | Complexity | Notes |
|---------|--------------|------------|-------|
| 100% Roxygen documentation | CRAN requires documentation for all exported functions; users expect help pages | Medium | Need to check all exported functions have proper @export tags and documentation |
| Passing `devtools::check()` | Basic quality gate; warnings/errors indicate problems | High | Current package may have warnings; goal is 0 errors, ≤1 warning |
| Working vignettes | Users learn from examples; vignettes must build during check | Medium | Three existing vignettes need verification |
| Consistent code style | Readability and maintainability; follows community standards | Low | Use lintr to enforce style; existing code may need formatting |
| Proper namespace | All imports/exports correctly declared in NAMESPACE | Medium | Roxygen should handle this automatically |
| Version control integration | Collaboration and change tracking | Low | Already using git; need proper .gitignore and workflow |

## Differentiators

Features that set well-maintained packages apart. Not expected, but valued.

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| Automated CI/CD | Catches issues early; ensures package works across environments | Medium | GitHub Actions for R packages |
| Code coverage reporting | Shows test completeness; identifies untested code | Medium | covr integration with CI |
| Professional documentation website | Easy navigation; searchable documentation | Medium | pkgdown with custom theme |
| Benchmark suite | Performance tracking over time | High | Existing particle filter benchmarks need integration |
| Reverse dependency checks | Ensures changes don't break dependent code | High | revdepcheck for major changes |
| Spell-checked documentation | Professional polish; catches typos | Low | spelling package integration |

## Anti-Features

Features to explicitly NOT build. Common mistakes in package maintenance.

| Anti-Feature | Why Avoid | What to Do Instead |
|--------------|-----------|-------------------|
| Manual NAMESPACE editing | Error-prone; conflicts with roxygen2 | Let roxygen2 manage NAMESPACE automatically |
| Global environment pollution | Causes side effects; hard to debug | Use proper function scoping; avoid `.GlobalEnv` assignments |
| Hard-coded paths | Not portable; breaks on other systems | Use `system.file()`, `tempdir()`, or user configuration |
| Excessive dependencies | Increases installation friction; more failure points | Minimize Imports; use Suggests for optional features |
| Monolithic functions | Hard to test and maintain | Break into smaller, focused functions with clear responsibilities |
| Ignoring warnings in check | Warnings often indicate real problems | Address warnings systematically; document intentional warnings |

## Feature Dependencies

```
Infrastructure Setup → Code Organization → Documentation → Testing → Build Verification
      ↓                      ↓                  ↓            ↓            ↓
   git setup           R file cleanup      Roxygen      testthat     devtools::check()
   .gitignore          C++ organization   examples      covr         rhub
   project structure   S3 method org      vignettes     lintr        revdepcheck
```

```
Phase 1: Foundation
├── Repository cleanup (root directory organization)
├── Basic documentation audit
└── Initial check() baseline

Phase 2: Code Quality
├── R code organization (numeric prefix consistency)
├── C++ code review
├── S3 method organization
└── lintr enforcement

Phase 3: Documentation
├── 100% Roxygen coverage
├── Working examples
├── Vignette verification
└── pkgdown website

Phase 4: Verification
├── Comprehensive test suite
├── CI/CD setup
├── Benchmark integration
└── Final check() validation
```

## MVP Recommendation

For initial cleanup (MVP), prioritize:

1. **Repository cleanup** - Move non-essential files from root to appropriate directories
2. **Documentation completeness** - 100% Roxygen coverage for exported functions
3. **Build verification** - `devtools::check()` passes with 0 errors, ≤1 warning

Defer to post-MVP:
- **Automated CI/CD** - Can be added after basic cleanup is complete
- **Code coverage metrics** - Focus on tests passing first, then measure coverage
- **Professional website** - Basic pkgdown is sufficient initially

## Critical Path

The most important sequence for success:

1. **Establish baseline** - Run `devtools::check()` to identify current issues
2. **Fix critical errors** - Address any check() errors first (blockers)
3. **Documentation gap analysis** - Identify undocumented exported functions
4. **Organize code structure** - Ensure R file numeric prefix consistency
5. **Verify vignettes** - Ensure they build correctly

## Risk Areas

| Feature | Risk Level | Mitigation |
|---------|------------|------------|
| C++ code organization | High | Review `src/` directory structure; check `Makevars` configuration |
| S3 method consistency | Medium | Audit all S3 generics and methods; ensure proper registration |
| Vignette dependencies | Medium | Check all vignette dependencies are in Suggests |
| Particle filter benchmarks | High | Separate benchmarking from package core; ensure no side effects |

## Sources

**Note:** Based on analysis of existing RLDM package and standard R package conventions:

- Current package structure analysis
- DESCRIPTION file dependencies
- R file organization (numeric prefix system)
- Existing vignettes and documentation
- Common R package maintenance patterns

**Confidence:** MEDIUM (standard practices, but specific to this package's needs)
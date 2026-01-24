# Roadmap: RLDM Package Cleanup & Streamlining

## Overview

This roadmap guides the cleanup and streamlining of the RLDM R package, focusing on repository organization, code quality, documentation completeness, build verification, and professional website deployment. The journey starts with foundational repository cleanup, progresses through code organization and build verification, completes comprehensive documentation, and finishes with a professional documentation website.

## Phases

**Phase Numbering:**
- Integer phases (1, 2, 3): Planned milestone work
- Decimal phases (2.1, 2.2): Urgent insertions (marked with INSERTED)

Decimal phases appear between their surrounding integers in numeric order.

- [ ] **Phase 1: Repository Foundation** - Clean up root directory, establish baseline structure
- [ ] **Phase 2: Code Organization** - Review and organize C++ and R code for maintainability
- [ ] **Phase 3: Build Verification** - Ensure package builds cleanly and passes checks
- [ ] **Phase 4: Documentation Completeness** - Achieve 100% Roxygen coverage with working examples
- [ ] **Phase 5: Website Deployment** - Create professional documentation website with pkgdown

## Phase Details

### Phase 1: Repository Foundation
**Goal**: Clean, organized repository with proper structure
**Depends on**: Nothing (first phase)
**Requirements**: REPO-01, REPO-02, REPO-03
**Success Criteria** (what must be TRUE):
  1. Root directory contains only essential files (README.md, LICENSE.md, CLAUDE.md, .gitignore)
  2. Benchmark scripts are moved to inst/benchmarks directory
  3. Log files are moved to logs/ directory
  4. Proper .gitignore for R package development is established
**Plans**: 5 plans (3 original + 2 gap closure)

Plans:
- [ ] 01-01-PLAN.md — Foundation Setup: Create directory structure, LICENSE.md, update .gitignore/.Rbuildignore
- [ ] 01-02-PLAN.md — File Organization: Move benchmark scripts, log files, debug/test scripts to organized locations
- [ ] 01-03-PLAN.md — Verification & Cleanup: Final root cleanup and human verification of structure
- [ ] 01-04-PLAN.md — Gap Closure: Human decisions and initial cleanup (gap closure)
- [ ] 01-05-PLAN.md — Gap Closure: Complete root cleanup and final verification (gap closure)

### Phase 2: Code Organization
**Goal**: Well-organized, maintainable code structure
**Depends on**: Phase 1
**Requirements**: CODE-01, CODE-02, CODE-03, CODE-04
**Success Criteria** (what must be TRUE):
  1. C++ code in src/ is reviewed and cleaned up
  2. R files follow numeric prefix consistency (01_, 02_, etc.)
  3. S3 methods are properly organized and registered
  4. DESCRIPTION file dependencies are reviewed and updated
**Plans**: TBD

Plans:
- [ ] 02-01: TBD

### Phase 3: Build Verification
**Goal**: Package builds cleanly and passes checks
**Depends on**: Phase 2
**Requirements**: BUILD-01, BUILD-02, BUILD-03
**Success Criteria** (what must be TRUE):
  1. `devtools::check()` passes with 0 errors and ≤1 warning
  2. Namespace imports/exports are correctly declared
  3. Package installs cleanly from source without errors
**Plans**: TBD

Plans:
- [ ] 03-01: TBD

### Phase 4: Documentation Completeness
**Goal**: Comprehensive, usable documentation
**Depends on**: Phase 3
**Requirements**: DOCS-01, DOCS-02, DOCS-03
**Success Criteria** (what must be TRUE):
  1. All exported functions have Roxygen documentation
  2. All exported functions have working examples
  3. Package-level documentation (?RLDM) exists and is accurate
**Plans**: TBD

Plans:
- [ ] 04-01: TBD

### Phase 5: Website Deployment
**Goal**: Professional documentation website
**Depends on**: Phase 4
**Requirements**: WEB-01, WEB-02
**Success Criteria** (what must be TRUE):
  1. pkgdown website builds locally without errors
  2. Website can be deployed to GitHub Pages or similar service
**Plans**: TBD

Plans:
- [ ] 05-01: TBD

## Progress

**Execution Order:**
Phases execute in numeric order: 1 → 2 → 3 → 4 → 5

| Phase | Plans Complete | Status | Completed |
|-------|----------------|--------|-----------|
| 1. Repository Foundation | 0/5 | Not started | - |
| 2. Code Organization | 0/1 | Not started | - |
| 3. Build Verification | 0/1 | Not started | - |
| 4. Documentation Completeness | 0/1 | Not started | - |
| 5. Website Deployment | 0/1 | Not started | - |

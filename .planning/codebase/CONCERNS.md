# Codebase Concerns

**Analysis Date:** 2026-01-24

## Tech Debt

**Particle Filter APF Implementation:**
- Issue: Auxiliary Particle Filter (APF) returns biased likelihood estimates for linear Gaussian models with cross-covariance (S ≠ 0)
- Files: `/media/bernd/nvme/r_projects/acad_RLDM/src/pf.cpp` (lines 220-653), `/media/bernd/nvme/r_projects/acad_RLDM/R/05_estimation_particle.R` (lines 164-166)
- Impact: APF method produces near-zero likelihood values when cross-covariance is present, making it unreliable for models with correlated state and observation noise
- Fix approach: Implement correct conditional distribution in APF weight calculation or deprecate APF for linear Gaussian models in favor of optimal proposal

**Debug Code in Production Files:**
- Issue: Multiple "FOR DEBUGGING" comments and debug code blocks remain in C++ source files
- Files: `/media/bernd/nvme/r_projects/acad_RLDM/src/solve_fwd_bwd_cll.cpp` (lines 472, 481, 490, 505, 512, 526, 699, 724)
- Impact: Code clutter, potential performance impact from unused debug statements
- Fix approach: Remove debug comments and conditional debug code blocks, implement proper logging system if needed

**Large Source Files:**
- Issue: Several R files exceed 1000 lines, making maintenance difficult
- Files: `/media/bernd/nvme/r_projects/acad_RLDM/R/02_templates.R` (1635 lines), `/media/bernd/nvme/r_projects/acad_RLDM/R/05_estimation_subspace.R` (1562 lines), `/media/bernd/nvme/r_projects/acad_RLDM/R/05_estimation_likelihood.R` (1211 lines)
- Impact: Reduced code readability, increased cognitive load for maintainers, harder to locate specific functionality
- Fix approach: Refactor into smaller, focused modules based on functionality

**Compiled Object Files in Source Control:**
- Issue: Multiple `.o` (object) files and compiled libraries (`.so`, `.dll`) in `src/` directory
- Files: `/media/bernd/nvme/r_projects/acad_RLDM/src/*.o`, `/media/bernd/nvme/r_projects/acad_RLDM/src/RLDM.so`, `/media/bernd/nvme/r_projects/acad_RLDM/src/RLDM.dll`
- Impact: Unnecessary repository bloat, potential version conflicts, violates standard R package structure
- Fix approach: Add to `.gitignore`, clean compiled artifacts, rely on build system

## Known Bugs

**APF Weight Calculation Bug:**
- Symptoms: APF returns biased likelihood estimates (near-zero) for linear Gaussian models with cross-covariance
- Files: `/media/bernd/nvme/r_projects/acad_RLDM/src/pf.cpp` (APF implementation)
- Trigger: Using `method = "apf"` with models where `norm(S, "F") > 1e-10`
- Workaround: Use `method = "optimal"` for linear Gaussian models with cross-correlation

**Debug Files in Project Root:**
- Symptoms: Multiple `debug_*.R` files clutter project root directory
- Files: `/media/bernd/nvme/r_projects/acad_RLDM/debug_*.R` (10+ files)
- Trigger: Development debugging left in production environment
- Workaround: Move to `inst/debug/` or `tests/debug/` directory

## Security Considerations

**Input Validation:**
- Risk: Potential buffer overflows or memory issues in C++ code with malformed inputs
- Files: `/media/bernd/nvme/r_projects/acad_RLDM/src/pf.cpp`, `/media/bernd/nvme/r_projects/acad_RLDM/src/kf.cpp`
- Current mitigation: Basic input validation in R wrapper functions, size checks in C++
- Recommendations: Add comprehensive bounds checking, validate matrix dimensions before operations, implement exception handling for edge cases

**External Dependency Management:**
- Risk: Dependency on GitHub package (`bfunovits/rationalmatrices`) without version pinning
- Files: `/media/bernd/nvme/r_projects/acad_RLDM/DESCRIPTION` (line 41)
- Current mitigation: `Remotes` field specifies GitHub source
- Recommendations: Consider CRAN submission of dependency or version pinning

## Performance Bottlenecks

**Particle Filter Memory Usage:**
- Problem: Stores full particle trajectories (s × N_particles × (N+1) array) by default
- Files: `/media/bernd/nvme/r_projects/acad_RLDM/src/pf.cpp` (lines 65, 250)
- Cause: `particle_trajectories` cube allocation for all time steps
- Improvement path: Add option to disable trajectory storage, implement checkpointing for long sequences

**C++ Loop Optimization:**
- Problem: Nested loops in particle filter implementations may be suboptimal
- Files: `/media/bernd/nvme/r_projects/acad_RLDM/src/pf.cpp` (main filter loops)
- Cause: Particle-by-particle operations instead of vectorized matrix operations
- Improvement path: Explore Armadillo vectorization, parallelization with OpenMP

**Singular Matrix Regularization:**
- Problem: Adds εI regularization to covariance matrices (1e-10 * ||matrix||_F)
- Files: `/media/bernd/nvme/r_projects/acad_RLDM/src/pf.cpp` (lines 80-90, 264-274)
- Cause: Cholesky decomposition fails for singular or near-singular matrices
- Improvement path: Consider rank-revealing decompositions (SVD, LDL') for better numerical stability

## Fragile Areas

**RMFD Model Implementation:**
- Files: `/media/bernd/nvme/r_projects/acad_RLDM/R/01_representations_classes.R` (rmfdmod class)
- Why fragile: Many methods not yet implemented, experimental status noted in documentation
- Safe modification: Limited to bug fixes until full implementation
- Test coverage: Minimal - primarily placeholder implementation

**Particle Filter Numerical Stability:**
- Files: `/media/bernd/nvme/r_projects/acad_RLDM/src/pf.cpp` (weight calculations, log-sum-exp)
- Why fragile: Weight degeneracy, numerical underflow in likelihood calculations
- Safe modification: Maintain log-domain computations, preserve regularization parameters
- Test coverage: Good for linear Gaussian, limited for nonlinear cases

**Template System Complexity:**
- Files: `/media/bernd/nvme/r_projects/acad_RLDM/R/02_templates.R` (1635 lines)
- Why fragile: Complex parameter mapping between deep and linear parameters
- Safe modification: Add tests before changes, maintain backward compatibility
- Test coverage: Moderate - template tests exist but may not cover all edge cases

## Scaling Limits

**Particle Filter Computational Complexity:**
- Current capacity: O(N × N_particles × s³) for matrix operations
- Limit: Memory (particle trajectories) and time (sequential resampling)
- Scaling path: Implement parallel resampling, reduce state dimension when possible, use Rao-Blackwellization

**Large State Dimension Models:**
- Current capacity: Limited by O(s³) matrix operations in Kalman/predictive steps
- Limit: Cubic scaling with state dimension s
- Scaling path: Exploit sparse structure, use iterative solvers, approximate methods

## Dependencies at Risk

**rationalmatrices GitHub Dependency:**
- Risk: Development package not on CRAN, breaking changes possible
- Impact: Build failures if dependency interface changes
- Migration plan: Monitor for CRAN submission, consider vendoring critical components

**RcppArmadillo Version Compatibility:**
- Risk: Armadillo linear algebra library updates may change numerical results
- Impact: Subtle numerical differences in filter outputs
- Migration plan: Version pinning, comprehensive regression testing

## Missing Critical Features

**Nonlinear/Non-Gaussian Model Support:**
- Problem: Particle filter implementation assumes linear Gaussian in C++ core
- Files: `/media/bernd/nvme/r_projects/acad_RLDM/src/pf.cpp` (hardcoded linear transitions)
- Blocks: Cannot handle true nonlinear state space models without code modification
- Priority: High for particle filter utility

**Parameter Estimation with Particle Filters:**
- Problem: No particle MCMC or SMC² implementation for parameter estimation
- Blocks: Cannot estimate model parameters using particle filters
- Priority: Medium - extends applicability beyond filtering

**GPU Acceleration:**
- Problem: No GPU support for particle filter computations
- Blocks: Performance limited to CPU, cannot leverage modern hardware
- Priority: Low - nice-to-have enhancement

## Test Coverage Gaps

**Nonlinear Model Tests:**
- What's not tested: Particle filter performance on true nonlinear/non-Gaussian models
- Files: `/media/bernd/nvme/r_projects/acad_RLDM/tests/testthat/test-pfilter.R`
- Risk: Implementation may fail for intended use cases
- Priority: High

**Numerical Edge Cases:**
- What's not tested: Extreme parameter values (near-singular matrices, large/small variances)
- Files: All estimation and filter test files
- Risk: Numerical instability in production use
- Priority: Medium

**Memory and Performance Tests:**
- What's not tested: Scaling behavior with large particle counts or long time series
- Files: No performance regression tests
- Risk: Performance degradation unnoticed in development
- Priority: Medium

**Cross-Platform Compatibility:**
- What's not tested: Windows/macOS/Linux numerical consistency
- Files: No platform-specific testing
- Risk: Different numerical results across platforms
- Priority: Low

---

*Concerns audit: 2026-01-24*
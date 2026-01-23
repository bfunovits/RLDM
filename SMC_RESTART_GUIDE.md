# Particle Filter Extension - Restart Guide

## Quick Start for Resuming Work

### 1. Load and Test Current State
```r
# Load the package
devtools::load_all()

# Run basic particle filter test
source("test_apf_fix.R")

# Run test suite (check current failures)
devtools::test(filter = "pfilter")
```

### 2. Critical Bug: Filtered State Consistency
**Problem**: `filtered_states[t,] ≠ particles[,,t] %*% weights` (9 failing tests)

**Debugging steps**:
1. Add debug output to `src/pf.cpp`:
   ```cpp
   // After computing filtered_states.col(t+1)
   Rcpp::Rcout << "t=" << t
               << " filtered=" << filtered_states.col(t+1).t()
               << " weighted_avg=" << (particles * weights).t() << std::endl;
   ```

2. Check weight propagation:
   - Are weights reset to uniform after resampling?
   - Do filtered_states use current weights or previous weights?

3. Verify particle trajectory storage:
   - `particle_trajectories.slice(t+1)` should store current `particles`
   - Check indexing: t=0 is initial, t=1 after first update

**Files to examine**:
- `src/pf.cpp`: Lines 180-185 (SIR) and 373-376 (APF) - filtered state computation
- `tests/testthat/test-pfilter.R`: Lines 129-135 - consistency test

### 3. Implement Optimal Proposal (High Priority)
**Goal**: Dramatically improve performance for linear Gaussian models

**Implementation steps**:
1. Add new function `pf_optimal_cpp()` to `src/pf.cpp`:
   ```cpp
   // Optimal proposal: p(x_t | x_{t-1}, y_t) = N(μ_t, Σ_t)
   // μ_t = A x_{t-1} + K_t (y_t - C A x_{t-1})
   // Σ_t = Q - K_t C Q
   // K_t = Q C' (C Q C' + R)^{-1}
   ```

2. Update R interface in `R/05_estimation_particle.R`:
   - Add `method = "optimal"` option to `pfilter()`
   - Add `filter_type = "optimal"` option to `ll_pfilter()`

3. Test performance improvement:
   - Compare with Kalman filter for linear Gaussian models
   - Expect near-zero weight degeneracy and accurate likelihood

### 4. Create Nonlinear Examples
**Purpose**: Demonstrate particle filter value over Kalman filter

**Example ideas**:
1. **Nonlinear state transition**:
   ```r
   # x_t = sin(x_{t-1}) + noise
   transition <- function(x) sin(x) + rnorm(length(x), 0, sigma_x)
   ```

2. **Non-Gaussian observation**:
   ```r
   # y_t ~ Poisson(λ = exp(C x_t))
   observation <- function(x) rpois(1, exp(C %*% x))
   ```

3. **Stochastic volatility model**:
   ```r
   # x_t = α + β x_{t-1} + σ_v v_t
   # y_t = exp(x_t/2) ε_t
   ```

### 5. Performance Optimization
**Potential improvements**:
1. **Vectorize weight computations**:
   ```cpp
   // Current: Loop over particles
   for (int i = 0; i < N_particles; i++) {
     y_pred = C * new_particles.col(i);
     // ...
   }

   // Possible: Compute all predictions at once
   arma::mat all_y_pred = C * new_particles;
   ```

2. **Parallelization**:
   - Use OpenMP for particle propagation
   - Use RcppParallel for weight computations

3. **Memory optimization**:
   - Add option to not store full particle trajectories
   - Store only summary statistics if not needed

### 6. Testing Strategy
**Immediate tests needed**:
1. **Bug fix verification**:
   ```r
   # After fixing filtered state consistency
   testthat::test_file("tests/testthat/test-pfilter.R")
   ```

2. **Optimal proposal validation**:
   ```r
   # Compare with Kalman filter
   source("test_optimal_proposal.R")
   ```

3. **Nonlinear model tests**:
   ```r
   # Create new test file
   testthat::test_file("tests/testthat/test-pfilter-nonlinear.R")
   ```

### 7. Documentation Updates
**Required updates**:
1. **Function documentation**:
   - Add optimal proposal method to `?pfilter`
   - Add nonlinear examples to `?pfilter` examples

2. **Vignette creation**:
   ```r
   # Create vignette demonstrating:
   # 1. Linear Gaussian baseline (comparison with Kalman)
   # 2. Nonlinear examples (particle filter superiority)
   # 3. Performance guidelines
   ```

3. **Package NEWS**:
   - Document particle filter addition
   - Note performance characteristics and limitations

### 8. Useful Commands
```r
# Build and test
Rcpp::compileAttributes()
devtools::document()
devtools::load_all()
devtools::check()

# Debug specific test
testthat::test_file("tests/testthat/test-pfilter.R", filter = "returns consistent")

# Profile performance
Rprof("pf_profile.out")
pf_result <- pfilter(model, y, N_particles = 10000, method = "sir")
Rprof(NULL)
summaryRprof("pf_profile.out")
```

### 9. Quick Reference
**Current working files**:
- `src/pf.cpp` - Core C++ implementation (bug location)
- `R/05_estimation_particle.R` - R interface
- `tests/testthat/test-pfilter.R` - Test suite (failing tests)
- `SMC_STATUS.md` - Detailed status report
- `SMC_EXTENSION_LOG.md` - Development log

**Key functions to examine**:
- `pf_sir_cpp()` lines 180-185: Filtered state computation (bug suspect)
- `pf_apf_cpp()` lines 373-376: Filtered state computation in APF
- `pfilter.stspmod()`: R interface parameter handling

**Test files for validation**:
- `test_apf_fix.R`: Current working test
- `test_convergence.R`: Convergence with particle count
- `test_pfilter_dev.R`: Comprehensive development test

---

## Immediate Next Actions (Recommended Order)

1. **Fix filtered state consistency bug** (blocking)
2. **Implement optimal proposal** (major performance improvement)
3. **Add nonlinear examples** (demonstrate value)
4. **Optimize performance** (vectorization, parallelization)
5. **Update documentation** (vignette, examples, NEWS)

**Estimated effort**:
- Bug fix: 2-4 hours
- Optimal proposal: 4-8 hours
- Nonlinear examples: 4-8 hours
- Optimization: 8-16 hours
- Documentation: 4-8 hours

**Success criteria**:
1. All tests passing
2. Optimal proposal matching Kalman filter performance for linear Gaussian
3. Clear nonlinear examples showing particle filter superiority
4. Comprehensive documentation and performance guidelines
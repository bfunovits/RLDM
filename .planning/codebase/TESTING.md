# Testing Patterns

**Analysis Date:** 2026-01-24

## Test Framework

**Runner:**
- testthat (>= 2.1.0)
- Config: Standard testthat structure with `tests/testthat/` directory

**Assertion Library:**
- testthat assertions: `expect_no_error()`, `expect_true()`, `expect_equal()`, `expect_named()`, `expect_type()`

**Run Commands:**
```bash
devtools::test()               # Run all tests
devtools::test(filter = "pattern")  # Run tests matching pattern
testthat::test_file("tests/testthat/test-*.R")  # Run specific test file
```

## Test File Organization

**Location:**
- Co-located in `tests/testthat/` directory
- One test file per functional area

**Naming:**
- `test-{functionality}.R` pattern
- Examples: `test-pfilter.R`, `test-templates.R`, `test-est_ar.R`

**Structure:**
```
tests/
├── testthat.R              # Test runner
└── testthat/
    ├── test-pfilter.R      # Particle filter tests
    ├── test-templates.R    # Template function tests
    ├── test-est_ar.R       # AR estimation tests
    └── _snaps/             # Snapshot test directory
```

## Test Structure

**Suite Organization:**
```r
test_that("descriptive test name", {
  # Setup
  set.seed(123)
  m <- 1
  s <- 2
  n <- s

  # Execution
  tmpl <- tmpl_stsp_full(m, n, s, sigma_L = "chol")
  model <- r_model(tmpl, bpoles = 1, sd = 0.5)
  data <- sim(model, n.obs = 20)

  # Assertions
  expect_no_error({
    pf_sir <- pfilter(model, data$y, N_particles = 100, method = "sir")
  })
  expect_true(inherits(pf_sir, "pfilter"))
  expect_named(pf_sir, c("filtered_states", "predicted_states", ...))
})
```

**Patterns:**
- **Setup pattern:** Use `set.seed()` for reproducibility, define test parameters
- **Execution pattern:** Create models, generate data, call functions
- **Assertion pattern:** Chain assertions with descriptive messages
- **Teardown pattern:** No explicit teardown needed (R garbage collection)

## Mocking

**Framework:** Not heavily used - tests focus on actual implementations

**Patterns:**
- Minimal mocking observed
- Tests use actual model generation with `r_model()` and `sim()`
- Randomness controlled with `set.seed()`

**What to Mock:**
- External dependencies not mocked (rationalmatrices package)
- Random number generation controlled via seeds

**What NOT to Mock:**
- Core mathematical operations
- Model creation and simulation

## Fixtures and Factories

**Test Data:**
```r
# Common test setup pattern
set.seed(123)
m <- 1
s <- 2
n <- s
tmpl <- tmpl_stsp_full(m, n, s, sigma_L = "chol")
model <- r_model(tmpl, bpoles = 1, sd = 0.5)
data <- sim(model, n.obs = 20)
y <- data$y
```

**Location:**
- Created inline within each test
- No separate fixture files
- Repeated patterns across tests

## Coverage

**Requirements:** No formal coverage requirements detected

**View Coverage:**
```r
# Not configured in package
# Would require covr package installation
```

## Test Types

**Unit Tests:**
- Scope: Individual functions and methods
- Approach: Test input validation, output structure, basic functionality
- Examples: `test-pfilter.R` tests particle filter methods

**Integration Tests:**
- Scope: Multiple functions working together
- Approach: Test template filling/extraction, model properties
- Examples: `test-templates.R` tests template round-trip consistency

**E2E Tests:**
- Framework: Not used - focused on unit and integration tests
- Scope: Complete workflows not explicitly tested

## Common Patterns

**Async Testing:**
```r
# Not applicable - no async operations
```

**Error Testing:**
```r
# Test for expected warnings
expect_warning(
  pfilter(model, y, N_particles = 100, method = "apf"),
  "APF method may produce biased likelihood estimates"
)

# Test for no warnings
expect_no_warning(pfilter(model, y, N_particles = 100, method = "sir"))
```

**Convergence Testing:**
```r
# Test convergence with increasing particle counts
particle_counts <- c(100, 500, 1000)
rmse_optimal <- numeric(length(particle_counts))

for (i in seq_along(particle_counts)) {
  N <- particle_counts[i]
  pf_opt <- pfilter(model, y, N_particles = N, method = "optimal")
  pf_states <- pf_opt$filtered_states[1:50, ]
  kf_states <- kf_result$a[1:50, ]
  rmse_optimal[i] <- sqrt(mean((pf_states - kf_states)^2))

  expect_true(is.finite(pf_opt$log_likelihood))
  expect_true(all(is.finite(pf_opt$effective_sample_size)))
}

# RMSE should generally decrease with more particles
expect_true(all(rmse_optimal < 1.0))
```

**Parameter Range Testing:**
```r
# Test across parameter ranges
for (sd in c(0.1, 0.5, 1.0, 2.0)) {
  tmpl <- tmpl_stsp_full(m, n, s, sigma_L = "chol")
  model <- r_model(tmpl, bpoles = 1, sd = sd)
  data <- sim(model, n.obs = 20)
  y <- data$y

  # All methods should run without errors
  expect_no_error(pfilter(model, y, N_particles = 200, method = "sir"))
  expect_no_error(pfilter(model, y, N_particles = 200, method = "apf"))
  expect_no_error(pfilter(model, y, N_particles = 200, method = "optimal"))
}
```

**Edge Case Testing:**
```r
# Test zero-dimensional state space
m <- 1
s <- 0  # zero-dimensional state space
n <- m
tmpl <- tmpl_stsp_full(m, n, s, sigma_L = "chol")
model <- r_model(tmpl, bpoles = 1, sd = 0.5)
data <- sim(model, n.obs = 5)
y <- data$y

# Should work with zero states
expect_no_error({
  pf <- pfilter(model, y, N_particles = 50, method = "sir")
})
expect_equal(dim(pf$filtered_states), c(6, 0))
expect_equal(dim(pf$particles), c(0, 50, 6))
```

---

*Testing analysis: 2026-01-24*
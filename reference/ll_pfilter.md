# Particle Filter Approximation of Log-Likelihood

Approximates the log-likelihood of a state space model using particle
filters. This function is useful for nonlinear/non-Gaussian models where
the exact Kalman filter likelihood is not available.

## Usage

``` r
ll_pfilter(
  model,
  y,
  N_particles = 1000,
  filter_type = c("sir", "apf", "optimal"),
  resampling = c("systematic", "multinomial", "stratified"),
  ess_threshold = 0.5,
  N_runs = 10,
  ...
)

# S3 method for class 'stspmod'
ll_pfilter(
  model,
  y,
  N_particles = 1000,
  filter_type = c("sir", "apf", "optimal"),
  resampling = c("systematic", "multinomial", "stratified"),
  ess_threshold = 0.5,
  N_runs = 10,
  ...
)
```

## Arguments

- model:

  [`stspmod()`](https://bfunovits.github.io/RLDM/reference/stspmod.md)
  object, which represents the state space model. For nonlinear models,
  additional parameters may be required.

- y:

  sample, i.e. an \\(N,m)\\ dimensional matrix, or a "time series"
  object (i.e. `as.matrix(y)` should return an \\(N,m)\\-dimensional
  numeric matrix). Missing values (`NA`, `NaN` and `Inf`) are **not**
  supported.

- N_particles:

  Number of particles to use (default: 1000).

- filter_type:

  Type of particle filter to use for likelihood approximation. Options:
  `"sir"` (default), `"apf"`.

- resampling:

  Resampling method: `"multinomial"`, `"systematic"` (default), or
  `"stratified"`.

- ess_threshold:

  Effective sample size threshold for triggering resampling (default:
  0.5). Resampling occurs when ESS \< ess_threshold \* N_particles.

- N_runs:

  Number of independent particle filter runs to average over (default:
  10). Averaging reduces the variance of the likelihood estimator.

- ...:

  Additional arguments passed to filter implementations.

## Value

Particle filter approximation of the log-likelihood (scaled by \\1/N\\).

## See also

[`ll_kf()`](https://bfunovits.github.io/RLDM/reference/ll_kf.md) for
exact Kalman filter likelihood (linear Gaussian models),
[`pf()`](https://rdrr.io/r/stats/Fdist.html) for particle filtering.

## Examples

``` r
# Linear Gaussian example: compare particle filter likelihood with Kalman filter
set.seed(123)
s = 2
m = 1
n = m
n.obs = 100

tmpl = tmpl_stsp_full(m, n, s, sigma_L = "chol")
model = r_model(tmpl, bpoles = 1, sd = 0.5)
data = sim(model, n.obs = n.obs)

# Kalman filter likelihood (exact)
ll_exact = ll_kf(model, data$y)

# Particle filter likelihood approximation
ll_approx = ll_pfilter(model, data$y, N_particles = 1000, N_runs = 5)

# Compare (should be close for linear Gaussian model with enough particles)
cat("Exact (Kalman):", ll_exact, "\n")
#> Exact (Kalman): -0.2533242 
cat("Approx (PF):", ll_approx, "\n")
#> Approx (PF): -33.54036 
cat("Difference:", ll_approx - ll_exact, "\n")
#> Difference: -33.28704 
```

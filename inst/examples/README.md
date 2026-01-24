# Particle Filter Examples for Nonlinear/Non-Gaussian Models

This directory contains example scripts demonstrating particle filters for nonlinear and non-Gaussian state space models. These examples show the superiority of particle filters over Kalman filter approximations for such models.

## Files

### Core Helper Functions
- `helper_pf.R` - Generic particle filter implementation in R with systematic, multinomial, and stratified resampling.

### Example Scripts

1. **`nonlinear_sin.R`** - Nonlinear state transition
   - Model: `x_t = sin(x_{t-1}) + η_t`, `y_t = x_t + ε_t`
   - Comparison: Particle filter vs Extended Kalman Filter (EKF)
   - Demonstrates handling of nonlinear state transitions

2. **`nonlinear_poisson.R`** - Non-Gaussian observation
   - Model: `x_t = A x_{t-1} + η_t`, `y_t ~ Poisson(λ = exp(C x_t))`
   - Comparison: Particle filter vs Gaussian approximation (EKF)
   - Demonstrates handling of Poisson observations

3. **`stochastic_volatility.R`** - Stochastic volatility model
   - Model: `x_t = α + β x_{t-1} + σ_v v_t`, `y_t = exp(x_t/2) ε_t`
   - Comparison: Particle filter vs log-squared transformation approximation
   - Classic financial time series model

## Usage

Each example is self-contained and can be run directly:

```r
# Run nonlinear sin example
source("inst/examples/nonlinear_sin.R")

# Run Poisson observation example
source("inst/examples/nonlinear_poisson.R")

# Run stochastic volatility example
source("inst/examples/stochastic_volatility.R")
```

## Key Features Demonstrated

- **Nonlinear state transitions**: Particle filters propagate particles through nonlinear dynamics without linearization
- **Non-Gaussian observations**: Particle filters compute exact observation likelihoods, not Gaussian approximations
- **Effective Sample Size (ESS) monitoring**: Adaptive resampling based on ESS threshold
- **Multiple resampling methods**: Systematic, multinomial, and stratified resampling
- **Diagnostic plots**: States, ESS, weights, likelihood contributions

## Comparison with Kalman Filter

Each example compares particle filter performance with an approximate Kalman filter:
1. **Extended Kalman Filter (EKF)**: Linearizes nonlinear dynamics (sin example)
2. **Gaussian approximation**: Approximates non-Gaussian observations as Gaussian (Poisson example)
3. **Log-squared transformation**: Transforms multiplicative noise to additive Gaussian (volatility example)

Particle filters generally provide more accurate state estimates and likelihood approximations for nonlinear/non-Gaussian models, at the cost of increased computational complexity.

## Integration with RLDM Package

These examples use standalone R implementations for educational clarity. The RLDM package's production particle filter implementation is in C++ (`src/pf.cpp`) and supports linear Gaussian models with optimal proposal, SIR, and APF algorithms.

## References

- Gordon, N. J., Salmond, D. J., & Smith, A. F. M. (1993). Novel approach to nonlinear/non-Gaussian Bayesian state estimation
- Pitt, M. K., & Shephard, N. (1999). Filtering via simulation: Auxiliary particle filters
- Doucet, A., Godsill, S., & Andrieu, C. (2000). On sequential Monte Carlo sampling methods for Bayesian filtering
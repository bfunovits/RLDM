# Technical Reference: RLDM Classes and Methods

This document is a **reference guide** for RLDM classes, mathematical
foundations, and estimation methods. For practical examples, see the
“Getting Started” and “Case Study” vignettes.

------------------------------------------------------------------------

## Preliminaries

### Notation

The following notation is used throughout RLDM:

| Symbol       | Meaning                                                                         |
|--------------|---------------------------------------------------------------------------------|
| $m$          | Dimension of the output process $\left( y_{t} \right)$                          |
| $n$          | Dimension of the noise process $\left( u_{t} \right)$, typically $n = m$        |
| $s$          | Dimension of the state process $\left( a_{t} \right)$ (state space models only) |
| $N$          | Sample size, denoted `n.obs` in code                                            |
| $\mathbb{E}$ | Expectation operator                                                            |
| $\gamma_{k}$ | Autocovariance at lag $k$                                                       |
| $\Sigma$     | Noise covariance matrix ${\mathbb{E}}\left( u_{t}u_{t}\prime \right)$           |

### Sign Convention

RLDM uses a **non-standard sign convention** for AR coefficients to
maintain consistency with matrix fraction descriptions.

**Standard form** (often seen in time series textbooks):
$$y_{t} = a_{1}y_{t - 1} + a_{2}y_{t - 2} + \cdots + a_{p}y_{t - p} + u_{t}$$

**RLDM form** (matrix fraction description):
$$y_{t} + a_{1}y_{t - 1} + a_{2}y_{t - 2} + \cdots + a_{p}y_{t - p} = u_{t}$$

Notice the **opposite signs** on the AR coefficients. This ensures
consistency when working with rational matrix fractions where:
$$a_{0}y_{t} + a_{1}y_{t - 1} + \cdots + a_{p}y_{t - p} = b_{0}u_{t} + b_{1}u_{t - 1} + \cdots + b_{q}u_{t - q}$$

with $a_{0} = I_{m}$ (identity matrix).

------------------------------------------------------------------------

## Model Classes

### Overview: Model Representations

RLDM supports three equivalent ways to represent rational linear dynamic
models:

#### 1. Left Matrix Fraction Description (LMFD): `armamod`

The ARMA/VARMA representation: $$a(z)y_{t} = b(z)u_{t}$$

where $a(z)$ and $b(z)$ are polynomial matrices in the lag operator $z$.

#### 2. State Space Form: `stspmod`

The canonical state space representation: $$\begin{aligned}
s_{t + 1} & {= As_{t} + Bu_{t}} \\
y_{t} & {= Cs_{t} + Du_{t}}
\end{aligned}$$

where $s_{t}$ is the state vector of dimension $s$.

#### 3. Right Matrix Fraction Description (RMFD): `rmfdmod`

An alternative ARMA representation (experimental): $$\begin{aligned}
w_{t} & {= c(z)u_{t}} \\
y_{t} & {= d(z)w_{t}}
\end{aligned}$$

------------------------------------------------------------------------

### VARMA/ARMA Models: `armamod` Class

The `armamod` class represents **VARMA** (multivariate ARMA) processes
in left matrix fraction form.

#### Structure

An `armamod` object is an S3 list with the following slots:

| Slot      | Type            | Description                                                |
|-----------|-----------------|------------------------------------------------------------|
| `sys`     | `lmfd [m,n]`    | Left matrix fraction: filter $k(z) = a(z)^{- 1}b(z)$       |
| `sigma_L` | `double [n,n]`  | Left square root of noise covariance ($\Sigma = LL\prime$) |
| `names`   | `character [m]` | Component names of $y_{t}$                                 |
| `label`   | `character`     | Descriptive label for the process                          |

**Class attribute:** `c("armamod", "rldm")`

#### Constructor

``` r
armamod(sys, sigma_L, names = NULL, label = NULL)
```

**Example:**

``` r
# Create a simple AR(1) model: y_t = 0.8*y_{t-1} + u_t
# In RLDM convention: y_t - 0.8*y_{t-1} = u_t
a_coef <- matrix(-0.8)  # negative sign for standard convention
b_coef <- matrix(1)
sys <- lmfd(a_coef, b_coef)

model <- armamod(sys, sigma_L = matrix(1), names = "Output", label = "AR(1)")
model
#> AR(1): ARMA model [1,1] with orders p = 0 and q = 0
#> AR polynomial a(z):
#>      z^0 [,1]
#> [1,]     -0.8
#> MA polynomial b(z):
#>      z^0 [,1]
#> [1,]        1
#> Left square root of noise covariance Sigma:
#>      u[1]
#> u[1]    1
```

#### Available Methods

``` r
methods(class = 'armamod')
#>  [1] as.stspmod autocov    freqresp   impresp    ll         poles     
#>  [7] predict    print      sim        spectrald  str        zeroes    
#> see '?methods' for accessing help and source code
```

Key methods: -
**[`autocov()`](https://bfunovits.github.io/RLDM/reference/autocov.md)** -
Compute autocovariance function -
**[`impresp()`](https://bfunovits.github.io/RLDM/reference/impresp.md)** -
Impulse response functions -
**[`spectrald()`](https://bfunovits.github.io/RLDM/reference/spectrald.md)** -
Spectral density -
**[`freqresp()`](https://bfunovits.github.io/RLDM/reference/freqresp.md)** -
Frequency response -
**[`predict()`](https://bfunovits.github.io/RLDM/reference/predict.md)** -
Forecasting -
**[`solve_de()`](https://bfunovits.github.io/RLDM/reference/solve_de.md)** -
Simulate from the model -
**[`print()`](https://rdrr.io/r/base/print.html) /
[`plot()`](https://rdrr.io/r/graphics/plot.default.html)** -
Visualization

------------------------------------------------------------------------

### State Space Models: `stspmod` Class

The `stspmod` class represents processes in canonical **state space**
form.

State space models offer several advantages over other representations:

- **Minimality**: The state dimension $s$ represents the essential
  “memory” of the system, providing the most compact representation.
- **Balancing**: Realizations can be balanced to improve numerical
  stability (`balance` parameter in template functions).
- **Latent state interpretation**: Hidden states can represent
  meaningful unobserved components (trends, cycles, factors).
- **Flexibility**: Can model complex dynamics with relatively few
  parameters through appropriate state dimension selection.

#### Structure

An `stspmod` object is an S3 list with slots:

| Slot      | Type            | Description                     |
|-----------|-----------------|---------------------------------|
| `sys`     | `stsp [m,n]`    | State space realization         |
| `sigma_L` | `double [n,n]`  | Left factor of noise covariance |
| `names`   | `character [m]` | Component names                 |
| `label`   | `character`     | Descriptive label               |

**Class attribute:** `c("stspmod", "rldm")`

#### Constructor

``` r
stspmod(sys, sigma_L, names = NULL, label = NULL)
```

**Example:**

``` r
# Create a simple state space model
# s_{t+1} = 0.8*s_t + u_t
# y_t = s_t

A_matrix <- matrix(0.8)
B_matrix <- matrix(1)
C_matrix <- matrix(1)
D_matrix <- matrix(0)

sys <- stsp(A_matrix, B_matrix, C_matrix, D_matrix)
model_ss <- stspmod(sys, sigma_L = matrix(1))
model_ss
#> state space model [1,1] with s = 1 states
#>      s[1] u[1]
#> s[1]  0.8    1
#> x[1]  1.0    0
#> Left square root of noise covariance Sigma:
#>      u[1]
#> u[1]    1
```

#### Available Methods

``` r
methods(class = 'stspmod')
#>  [1] autocov    freqresp   impresp    ll_pfilter ll         pfilter   
#>  [7] poles      predict    print      sim        spectrald  str       
#> [13] zeroes    
#> see '?methods' for accessing help and source code
```

Similar to `armamod`, state space models support: - Autocovariance,
impulse responses, spectral density computation - Forecasting and
simulation - Conversion to impulse response representation

#### State Space Parameterizations

RLDM supports several standard state space canonical forms through
template functions:

- **[`tmpl_DDLC()`](https://bfunovits.github.io/RLDM/reference/local_model_structures.md)** -
  Diagonal-Direct-Lead-Coefficient form
- **[`tmpl_stsp_echelon()`](https://bfunovits.github.io/RLDM/reference/model_structures.md)** -
  Echelon canonical form
- **`tmpl_stsp_innovation()`** - Innovation form

------------------------------------------------------------------------

### Right Matrix Fraction Description: `rmfdmod` Class

The `rmfdmod` class represents processes in right matrix fraction form
(experimental, limited methods):

``` r
rmfdmod(sys, sigma_L, names = NULL, label = NULL)
```

Most estimation and analysis methods are not yet implemented for this
class.

------------------------------------------------------------------------

## Derived Objects

These classes represent computed properties of process models, not model
specifications themselves.

### Autocovariance: `autocov` Class

Stores autocovariances, autocorrelations, or partial autocorrelations.

#### Structure

| Slot    | Type                    | Description                               |
|---------|-------------------------|-------------------------------------------|
| `acf`   | array                   | Correlations/covariances at each lag      |
| `type`  | `character`             | “covariance”, “correlation”, or “partial” |
| `gamma` | `array [m,m,lag.max+1]` | 3D array of autocovariances               |
| `names` | `character [m]`         | Variable names                            |
| `label` | `character`             | Descriptive label                         |

#### Computation

``` r
# From data
autocov.default(obj, lag.max = 20, type = "correlation")

# From model
autocov.armamod(obj, lag.max = 20, type = "correlation")
autocov.stspmod(obj, lag.max = 20, type = "correlation")
```

**Example:**

``` r
# Generate AR model and compute ACF
y <- sim(model, n.obs = 500)
acov <- autocov(model, lag.max = 15, type = "correlation")
plot(acov)
```

![Time series plot of an AR(1) model showing simulated process output
over 500
observations](2_technical_reference_files/figure-html/unnamed-chunk-9-1.png)

### Impulse Response: `impresp` Class

Stores impulse response coefficients: $k_{j}$ where
$y_{t} = \sum_{j \geq 0}k_{j}u_{t - j}$.

#### Structure

| Slot      | Type                 | Description             |
|-----------|----------------------|-------------------------|
| `irf`     | `array [m,n,n.lags]` | IRF coefficients        |
| `sigma_L` | `double [n,n]`       | Noise covariance factor |
| `names`   | `character [m]`      | Output variable names   |
| `label`   | `character`          | Descriptive label       |

#### Computation

``` r
impresp.armamod(obj, lag.max = 40, H = NULL)
impresp.stspmod(obj, lag.max = 40, H = NULL)
```

**Example:**

``` r
irf <- impresp(model, lag.max = 20)
plot(irf, main = "Impulse Response")
```

![Four-panel impulse response plot showing how shocks in state space
model propagate to outputs over 40
lags](2_technical_reference_files/figure-html/unnamed-chunk-11-1.png)

The `H` parameter allows custom orthogonalization (default: Cholesky).

### Spectral Density: `spectrald` Class

Frequency-domain representation:
$\Gamma(\lambda) = K(\lambda)\Sigma K^{*}(\lambda)$ where $K(\lambda)$
is the frequency response.

#### Structure

| Slot    | Type              | Description                        |
|---------|-------------------|------------------------------------|
| `spd`   | `array [m,m,n.f]` | Spectral density on frequency grid |
| `names` | `character [m]`   | Variable names                     |
| `label` | `character`       | Descriptive label                  |

#### Computation

``` r
spectrald(obj, n.f = 128, ...)
```

**Example:**

``` r
spec <- spectrald(model, n.f = 256)
plot(spec)
```

![Four-panel spectral density plot showing power across frequencies for
outputs of the state space
model](2_technical_reference_files/figure-html/unnamed-chunk-13-1.png)

### Frequency Response: `freqresp` Class

The frequency transfer function
$K(\lambda) = \sum_{j}k_{j}e^{- i\lambda j}$.

``` r
freqresp.armamod(obj, n.f = 256)
```

### Forecast Error Variance Decomposition: `fevardec` Class

Decomposes forecast error variance into contributions from each shock.

| Slot    | Type            | Description             |
|---------|-----------------|-------------------------|
| `vd`    | `array [m,m,h]` | Variance decomposition  |
| `v`     | `matrix [m,h]`  | Forecast error variance |
| `names` | `character [m]` | Variable names          |

``` r
fevardec(obj, h.max = 40)
```

------------------------------------------------------------------------

## Estimation Methods

### AR Model Estimation

#### Yule-Walker Method: `est_ar_yw()`

**Theory:** Solves the Yule-Walker equations using Cholesky
decomposition:
$$\gamma_{k} = a_{1}\gamma_{k - 1} + \cdots + a_{p}\gamma_{k - p}$$

**Advantages:** - Statistically efficient (minimum variance) - Produces
all AR orders simultaneously - Automatic log-determinant sequence for
model selection

**Example:**

``` r
est_ar_yw(gamma, p.max = 10, penalty = -1)
```

#### Durbin-Levinson-Whittle Method: `est_ar_dlw()`

**Theory:** Recursive algorithm for AR parameter and partial
autocorrelation computation.

**Advantages:** - Numerically stable - Computes partial ACF - Efficient
recursions

**Example:**

``` r
est_ar_dlw(gamma, p.max = 10, penalty = -1)
```

#### Wrapper Function: `est_ar()`

**Primary interface for AR estimation with automatic order selection:**

``` r
est_ar(data, p.max = NULL, method = "auto",
       ic = "AIC", mean_estimate = "sample.mean")
```

**Parameters:** - `ic`: Information criterion (“AIC”, “BIC”, “HQ”) -
`mean_estimate`: How to estimate mean (“sample.mean”, “zero”, or a
vector) - `method`: “yw” (Yule-Walker), “dlw” (Durbin-Levinson-Whittle),
or “auto”

**Returns:** List with components: - `model`: Estimated `armamod`
object - `p`: Selected order - `stats`: Criterion values for each
order - `ll`: Log-likelihood values

------------------------------------------------------------------------

### ARMA Model Estimation

#### Hannan-Rissanen-Kavalieris (HRK): `est_arma_hrk()`

Three-stage procedure for VARMA estimation without nonlinear
optimization:

``` r
est_arma_hrk(data, tmpl = NULL, p = NULL, q = NULL,
             mean_estimate = "zero")
```

**Stage 1:** Estimate long AR model **Stage 2:** Compute initial ARMA
estimates **Stage 3:** Refine using feasible GLS

**Advantages:** - Provides consistent initial estimates - No nonlinear
optimization required - Suitable for maximum likelihood refinement

**When to use:** - Initial estimates for multivariate ARMA - Moderate to
high-dimensional systems - When you want to avoid local optima

------------------------------------------------------------------------

### State Space Estimation

#### CCA Method: `est_stsp_ss(method = "cca")`

**Canonical Correlation Analysis** for determining state space order and
initial estimates.

``` r
est_stsp_ss(data, method = "cca", s = NULL,
            keep_models = FALSE, ...)
```

**Advantages:** - Data-driven order selection - Suitable when no prior
knowledge of system structure - Good numerical properties

**State Dimension Selection:** When `s = NULL` (default), CCA
automatically selects state dimension using singular value criteria. For
manual selection, consider: - **Singular Value Criterion (SVC)**: Look
for a gap in singular values of the weighted Hankel matrix -
**Information Criteria**: AIC/BIC computed for different state
dimensions - **Cross-validation**: Out-of-sample prediction
performance - **Interpretability**: Choose dimension that matches
expected number of latent factors

In practice, start with automatic selection via
[`est_stsp_ss()`](https://bfunovits.github.io/RLDM/reference/subspace-methods.md)
with `s = NULL`, then inspect the singular values plot (if available) or
compare information criteria across dimensions.

#### Subspace Methods: `est_stsp_ss(method = "cca"/"ho"/"moesp")`

Multiple subspace identification methods available through matrix
fraction estimation.

------------------------------------------------------------------------

### Maximum Likelihood Estimation

#### General ML Framework: `ll_FUN()`

Create likelihood functions for structured parameter templates:

``` r
llfun <- ll_FUN(tmpl, data, skip = 0, which = "concentrated")
```

Then optimize using standard R optimizers:

``` r
out <- optim(theta0, llfun, method = "BFGS", control = list(fnscale = -1))
```

#### Parameter Templates

Templates map deep model parameters to linear parameters for
optimization:

- **[`tmpl_DDLC()`](https://bfunovits.github.io/RLDM/reference/local_model_structures.md)** -
  Diagonal-Direct-Lead-Coefficient
- **[`tmpl_arma_echelon()`](https://bfunovits.github.io/RLDM/reference/model_structures.md)** -
  ARMA in echelon form
- **[`tmpl_stsp_echelon()`](https://bfunovits.github.io/RLDM/reference/model_structures.md)** -
  State space in echelon form

**Workflow:**

``` r
# Create template from initial estimate
tmpl <- tmpl_DDLC(model_initial, sigma_L = 'identity')

# Extract starting parameter vector
theta0 <- extract_theta(model_initial, tmpl)

# Create and optimize likelihood
llfun <- ll_FUN(tmpl, data, which = "concentrated")
opt <- optim(theta0, llfun, method = "BFGS", control = list(fnscale = -1, maxit = 500))

# Reconstruct model from optimized parameters
model_ml <- fill_template(opt$par, tmpl)
```

### Particle Filtering

#### Overview

Particle filters (Sequential Monte Carlo methods) extend Kalman
filtering to **nonlinear and/or non-Gaussian state space models**. While
the Kalman filter provides optimal closed-form solutions for linear
Gaussian models, particle filters use Monte Carlo approximation to
handle:

- **Nonlinear state transitions**:
  $x_{t} = f\left( x_{t - 1},u_{t} \right)$
- **Non-Gaussian observations**: $y_{t} \sim p\left( y|x_{t} \right)$
  with arbitrary distribution
- **Complex noise structures**

The RLDM implementation includes three particle filter variants:

1.  **SIR (Bootstrap) filter**: Basic sampling importance resampling
2.  **APF (Auxiliary Particle Filter)**: Two-stage resampling with
    lookahead
3.  **Optimal proposal**: For linear Gaussian models, approaches Kalman
    filter performance

#### Basic Usage: Linear Gaussian Comparison

For linear Gaussian models, particle filters approximate the Kalman
filter solution:

``` r
library(RLDM)

# Create a linear Gaussian state space model
m <- 1; s <- 2; n <- s
tmpl <- tmpl_stsp_full(m, n, s, sigma_L = "chol")
model <- r_model(tmpl, bpoles = 1, sd = 0.5)
data <- sim(model, n.obs = 100)
y <- data$y

# Kalman filter (optimal for linear Gaussian)
kf_result <- kf(model, y)

# Particle filters
pf_sir <- pfilter(model, y, N_particles = 1000, method = "sir")
pf_apf <- pfilter(model, y, N_particles = 1000, method = "apf")
pf_opt <- pfilter(model, y, N_particles = 1000, method = "optimal")

# Compare filtered states
plot(kf_result$a[,1], type = "l", col = "black", lwd = 2,
     main = "Filtered State Estimates")
lines(pf_sir$filtered_states[,1], col = "blue", lty = 2)
lines(pf_apf$filtered_states[,1], col = "red", lty = 2)
lines(pf_opt$filtered_states[,1], col = "green", lty = 2)
legend("topright", legend = c("Kalman", "SIR", "APF", "Optimal"),
       col = c("black", "blue", "red", "green"), lty = c(1, 2, 2, 2), lwd = c(2, 1, 1, 1))
```

#### Likelihood Approximation

Particle filters can approximate the log-likelihood for
nonlinear/non-Gaussian models where exact Kalman filter likelihood is
unavailable:

``` r
# Kalman filter likelihood (exact for linear Gaussian)
ll_kf <- ll_kf(model, y)

# Particle filter likelihood approximation (averaged over multiple runs)
ll_pf <- ll_pfilter(model, y, N_particles = 1000, filter_type = "apf", N_runs = 10)

cat("Kalman filter LL:", ll_kf, "\n")
cat("Particle filter LL:", ll_pf, "\n")
cat("Difference:", ll_pf - ll_kf, "\n")
```

#### Nonlinear and Non-Gaussian Examples

RLDM includes example scripts demonstrating particle filter superiority
for nonlinear/non-Gaussian models:

1.  **`inst/examples/nonlinear_sin.R`** - Nonlinear state transition:
    $x_{t} = \sin\left( x_{t - 1} \right) + \eta_{t}$
2.  **`inst/examples/nonlinear_poisson.R`** - Non-Gaussian observation:
    $y_{t} \sim \text{Poisson}\left( \lambda = \exp\left( Cx_{t} \right) \right)$
3.  **`inst/examples/stochastic_volatility.R`** - Stochastic volatility
    model: $y_{t} = \exp\left( x_{t}/2 \right)\varepsilon_{t}$

Each example compares particle filter performance with approximate
Kalman filter variants (Extended Kalman Filter, Gaussian approximation,
log-squared transformation).

#### Diagnostics and Visualization

The
[`plot.pfilter()`](https://bfunovits.github.io/RLDM/reference/plot.pfilter.md)
method provides several diagnostic plots:

``` r
pf <- pfilter(model, y, N_particles = 500, method = "sir")

# Filtered state estimates
plot(pf, type = "states")

# Effective Sample Size (ESS) over time
plot(pf, type = "ess")

# Particle weight distribution at final time
plot(pf, type = "weights")

# Log-likelihood contributions per time step
plot(pf, type = "likelihood")
```

#### Performance Considerations

- **Particle count**: More particles improve accuracy but increase
  computation (typical: 1000-10000)
- **Resampling methods**: Systematic (default), multinomial, or
  stratified
- **ESS threshold**: Adaptive resampling when ESS \< threshold ×
  N_particles (default: 0.5)
- **Optimal proposal**: Use for linear Gaussian models to minimize
  weight degeneracy

#### References

- Gordon, N. J., Salmond, D. J., & Smith, A. F. M. (1993). Novel
  approach to nonlinear/non-Gaussian Bayesian state estimation
- Pitt, M. K., & Shephard, N. (1999). Filtering via simulation:
  Auxiliary particle filters
- Doucet, A., Godsill, S., & Andrieu, C. (2000). On sequential Monte
  Carlo sampling methods for Bayesian filtering

------------------------------------------------------------------------

## Solution and Simulation

### Solving Difference Equations

#### Forward Solution: `solve_de()`

Simulate forward from initial conditions:

$$y_{t} = k_{0}u_{t} + k_{1}u_{t - 1} + k_{2}u_{t - 2} + \cdots$$

``` r
y <- solve_de(sys, u, y_init = NULL)
```

#### Inverse Solution: `solve_inverse_de()`

Extract residuals given data (solve backwards):

$$u_{t} = a(z)^{- 1}b(z)^{- 1}y_{t}$$

``` r
u <- solve_inverse_de(sys, y)$u
```

------------------------------------------------------------------------

## Mathematical Foundations

### Stationary ARMA Processes

An ARMA process is described by: $$a(z)y_{t} = b(z)u_{t}$$

where: - $a(z) = I_{m} + a_{1}z + \cdots + a_{p}z^{p}$ (AR polynomial) -
$b(z) = b_{0} + b_{1}z + \cdots + b_{q}z^{q}$ (MA polynomial) -
$\left( u_{t} \right)$ is white noise with
${\mathbb{E}}\left( u_{t}u_{t}\prime \right) = \Sigma$

The process is **stationary** if all roots of
$\det\left( a(z) \right) = 0$ lie outside the unit circle.

The stationary solution is:
$$y_{t} = \sum\limits_{j \geq 0}k_{j}u_{t - j}$$

where $k_{j}$ satisfy: $a(z)k(z) = b(z)$.

### Autocovariance Function

For stationary ARMA processes, the autocovariance function
$\gamma_{k} = {\mathbb{E}}\left( y_{t + k}y_{t}\prime \right)$ satisfies
the **generalized Yule-Walker equations**:

$$a_{0}\gamma_{k} + a_{1}\gamma_{k - 1} + \cdots + a_{p}\gamma_{k - p} = \begin{cases}
{b_{0}\Sigma k_{0}\prime + \cdots + b_{q}\Sigma k_{q - k}\prime} & {0 \leq k \leq q} \\
0 & {k > q}
\end{cases}$$

### Spectral Factorization

The spectral density can be written as:
$$\Gamma(\lambda) = \frac{1}{2\pi}K(\lambda)\Sigma K^{*}(\lambda)$$

where $K(\lambda)$ is the frequency response and $K^{*}$ denotes
conjugate transpose.

This enables: - Frequency domain analysis - Filter design - Whitening
transformations

### Kronecker Indices

For VARMA models in canonical form, **Kronecker indices**
$\left( n_{1},\ldots,n_{m} \right)$ define the minimal orders required
for each output: - Total state dimension: $s = n_{1} + \cdots + n_{m}$ -
They specify the echelon structure

------------------------------------------------------------------------

## Method Selection Guide

### Choosing Between AR, ARMA, and State Space

| Method          | When to Use                     | Advantages                     | Disadvantages               |
|-----------------|---------------------------------|--------------------------------|-----------------------------|
| **AR**          | Baseline, simple systems        | Fast, stable, interpretable    | May need high order         |
| **ARMA**        | Parsimonious fits               | Fewer parameters               | Parameter estimation harder |
| **State Space** | Complex systems, latent factors | Flexible, interpretable states | Requires order selection    |

### Choosing State Space Parameterizations

RLDM supports several state space parameterizations, each with different
strengths:

| Parameterization                            | When to Use                                                    | Key Features                                                | Estimation Method                                                                                                                               |
|---------------------------------------------|----------------------------------------------------------------|-------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------|
| **CCA (Canonical Correlation Analysis)**    | Initial exploration, data-driven order selection               | Automatic state dimension selection, good initial estimates | `est_stsp_ss(method = "cca")`                                                                                                                   |
| **DDLC (Diagonal-Direct-Lead-Coefficient)** | Maximum likelihood refinement, stable optimization             | Numerically stable, orthogonal to equivalence classes       | ML with [`tmpl_DDLC()`](https://bfunovits.github.io/RLDM/reference/local_model_structures.md) + [`optim()`](https://rdrr.io/r/stats/optim.html) |
| **Echelon Form**                            | Parsimonious canonical representation, structural restrictions | Minimal parameters via Kronecker indices, identifiable      | [`tmpl_stsp_echelon()`](https://bfunovits.github.io/RLDM/reference/model_structures.md) + ML                                                    |

**Guidelines:** 1. **Start with CCA** for initial model discovery and
automatic order selection 2. **Refine with DDLC** for maximum likelihood
estimation with good numerical properties 3. **Use echelon form** when
you need canonical, identifiable representations with structural
restrictions 4. **Convert between forms** using
[`pseries2stsp()`](https://bfunovits.github.io/rationalmatrices/reference/pseries2stsp.html)
and template functions as needed

### Model Comparison Workflow

1.  **Start with AR** using
    [`est_ar()`](https://bfunovits.github.io/RLDM/reference/est_ar.md)
2.  **Try state space** with
    [`est_stsp_ss()`](https://bfunovits.github.io/RLDM/reference/subspace-methods.md)
    for dimensionality reduction
3.  **Compare with ARMA** using HRK initial estimates + ML refinement
4.  **Refine best model** with parameter templates and optimization
5.  **Validate** using residual diagnostics, ACF, and spectral plots

### Estimation Method Selection

| Task                     | Recommended Method                                                                        |
|--------------------------|-------------------------------------------------------------------------------------------|
| Order selection          | AIC via [`est_ar()`](https://bfunovits.github.io/RLDM/reference/est_ar.md)                |
| Quick state space        | CCA via [`est_stsp_ss()`](https://bfunovits.github.io/RLDM/reference/subspace-methods.md) |
| Multivariate ARMA        | HRK via [`est_arma_hrk()`](https://bfunovits.github.io/RLDM/reference/est_arma_hrk.md)    |
| Maximum likelihood       | Parameter templates + [`optim()`](https://rdrr.io/r/stats/optim.html)                     |
| Highly structured models | Custom templates with `tmpl_*()`                                                          |

------------------------------------------------------------------------

## References

(Scherrer and Deistler 2019)

Scherrer, Wolfgang, and Manfred Deistler. 2019. “Chapter 6 - Vector
Autoregressive Moving Average Models.” In *Conceptual Econometrics Using
r*, edited by Hrishikesh D. Vinod and C. R. Rao, 41:145–91. Handbook of
Statistics. Elsevier.
https://doi.org/<https://doi.org/10.1016/bs.host.2019.01.004>.

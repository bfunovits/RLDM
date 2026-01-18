# Getting Started with RLDM

## What is RLDM?

**RLDM** (Rational Linear Dynamic Models) is an R package for modeling
and analyzing **stationary time series** using linear systems theory. It
implements powerful methods for:

- **Autoregressive (AR) models** - Simple, interpretable baseline models
- **ARMA (VARMA) models** - Flexible models combining AR and Moving
  Average terms
- **State Space models** - Compact representations of complex dynamics
- **Impulse response analysis** - Understanding how shocks propagate
  through your system
- **Spectral analysis** - Analyzing behavior at different frequencies
- **Model comparison** - Selecting the best model for your data

RLDM is particularly useful when you have **multivariate time series
data** and want to go beyond simple AR(1) models to capture more complex
dependencies.

## Installation

RLDM depends on the companion package `rationalmatrices` for matrix
fraction computations. Install both from GitHub:

``` r
# Install dependencies first
remotes::install_github("bfunovits/rationalmatrices")

# Then install RLDM
remotes::install_github("bfunovits/RLDM")
```

Then load the package:

``` r
library(RLDM)
```

## Simple AR Model: Univariate Example

Let’s start with a simple univariate time series. We’ll generate
synthetic data from an AR(2) process and estimate the model.

### Generate AR(2) Data

``` r
# Simulate an AR(1) process: y_t = 0.8*y_{t-1} + u_t
#
# In RLDM convention: y_t - 0.8*y_{t-1} = u_t
a_coef <- matrix(-0.8)  # AR coefficient with sign convention
b_coef <- matrix(1)     # MA term (just white noise)
sys <- lmfd(a_coef, b_coef)

# Create ARMA model with noise variance = 1
model_true <- armamod(sys, sigma_L = matrix(1), label = "AR(1)")

# Simulate 200 observations
y_list <- sim(model_true, n.obs = 200)
y <- y_list$y

# Plot the simulated data
plot(y, main = "Simulated AR(1) Process", ylab = "y_t", type = 'l')
```

![Time series plot of simulated AR(1) process showing 200 observations
over time with fluctuations around
zero](0_getting_started_files/figure-html/unnamed-chunk-3-1.png)

### Estimate the Model

Now let’s estimate an AR model from this data using the Yule-Walker
method. RLDM uses information criteria (AIC) to automatically select the
order:

``` r
# Estimate AR model with automatic order selection
result <- est_ar(y, p.max = 10)

# Display basic information
cat("Selected AR order:", result$p, "\n")
#> Selected AR order: 0
cat("Number of parameters:", result$p, "\n")
#> Number of parameters: 0
cat("Information Criterion:", result$stats[result$p + 1, "ic"], "\n\n")
#> Information Criterion: 0.3897879

# View the estimated model
result$model
#> ARMA model [1,1] with orders p = 0 and q = 0
#> AR polynomial a(z):
#>      z^0 [,1]
#> [1,]        1
#> MA polynomial b(z):
#>      z^0 [,1]
#> [1,]        1
#> Left square root of noise covariance Sigma:
#>          u[1]
#> u[1] 1.215182
```

### Model Diagnostics

Check if the model adequately captures the dynamics:

``` r
# Compute residuals
residuals <- solve_inverse_de(result$model$sys, y)$u

# Plot residuals
par(mfrow = c(2, 1))
plot(residuals, main = "Residuals", type = 'l', ylab = "Residuals")
abline(h = 0, col = 'red', lty = 2)

# ACF of residuals (should show no significant correlations)
acf(residuals, main = "ACF of Residuals")
```

![Two-panel diagnostic plot: upper panel shows residuals from AR(1)
model estimation over time with a red reference line at zero; lower
panel shows autocorrelation function of
residuals](0_getting_started_files/figure-html/unnamed-chunk-5-1.png)

``` r
par(mfrow = c(1, 1))
```

### Make Predictions

Forecast future values:

``` r
# Predict 10 steps ahead
n.ahead <- 10
pred <- predict(result$model, y, h = 1, n.ahead = n.ahead)

# Combine data and predictions
n_obs <- length(y)
plot_range <- max(1, n_obs - 40):n_obs
plot(plot_range, y[plot_range], type = 'l',
     xlab = "Time", ylab = "Value",
     main = "Data and Predictions")

# Add predictions
pred_time <- (n_obs+1):(n_obs+n.ahead)
lines(pred_time, pred$yhat.ahead[, 1], col = 'blue', lwd = 2)

# Add confidence bands (approximate, assuming normality)
# Extract the variance for each step ahead (the variance increases with horizon)
variance_ahead <- as.vector(pred$sigmahat.ahead[1, 1, ])
se <- sqrt(variance_ahead)
lines(pred_time, pred$yhat.ahead[, 1] + 1.96*se, col = 'blue', lty = 2)
lines(pred_time, pred$yhat.ahead[, 1] - 1.96*se, col = 'blue', lty = 2)

legend("topleft", c("Data", "Forecast", "95% CI"),
       col = c("black", "blue", "blue"), lty = c(1, 1, 2))
```

![Forecast plot showing recent data (black line) from AR(1) model with
10-step ahead predictions (blue line) and 95% confidence bands (blue
dashed
lines)](0_getting_started_files/figure-html/unnamed-chunk-6-1.png)

## Multivariate Example: Bivariate VAR

Now let’s work with a **multivariate system** - a 2-dimensional VAR(1)
model.

### Generate Bivariate Data

``` r
# Define a 2D VAR(1) model
#
# [ y1_t ]   [ 0.8  0.1 ] [ y1_{t-1} ]   [ u1_t ]
# [ y2_t ] = [ 0.2  0.7 ] [ y2_{t-1} ] + [ u2_t ]

# AR matrices (note sign convention)
A1 <- matrix(c(-0.8, -0.2, -0.1, -0.7), nrow = 2, byrow = TRUE)
B0 <- diag(2)

# Create VARMA system
sys_var <- lmfd(A1, B0)

# Noise covariance (correlation between u1 and u2)
rho <- 0.5
Sigma_L <- matrix(c(1, rho, 0, sqrt(1 - rho^2)), nrow = 2)

# Create model and simulate
model_var <- armamod(sys_var, sigma_L = Sigma_L,
                     names = c("Series 1", "Series 2"))

y_var_list <- sim(model_var, n.obs = 300)
y_var <- y_var_list$y

# Plot both series
plot.ts(y_var, main = "Bivariate VAR(1) Process")
```

![Time series plot of two correlated series from a bivariate VAR(1)
process with 300 observations
each](0_getting_started_files/figure-html/unnamed-chunk-7-1.png)

### Estimate the Multivariate Model

``` r
# Estimate VAR model
result_var <- est_ar(y_var, p.max = 5)

cat("Selected VAR order:", result_var$p, "\n")
#> Selected VAR order: 0
cat("Total parameters:", result_var$p * 2^2, "\n\n")
#> Total parameters: 0

# View estimated model
result_var$model
#> ARMA model [2,2] with orders p = 0 and q = 0
#> AR polynomial a(z):
#>      z^0 [,1]  [,2]
#> [1,]        1     0
#> [2,]        0     1
#> MA polynomial b(z):
#>      z^0 [,1]  [,2]
#> [1,]        1     0
#> [2,]        0     1
#> Left square root of noise covariance Sigma:
#>           u[1]     u[2]
#> u[1] 1.1841356 0.000000
#> u[2] 0.1760054 1.273754
```

### Analyze the System

**Impulse Response Function**: How does the system respond to shocks?

``` r
# Compute impulse responses
irf <- impresp(result_var$model, lag.max = 20)

# Plot impulse responses
plot(irf,
     main = "Impulse Response Functions",
     legend = c("shock to Series 1", "shock to Series 2"))
```

![Four-panel impulse response functions plot showing how shocks to each
series in a bivariate VAR(1) system propagate over 20
lags](0_getting_started_files/figure-html/unnamed-chunk-9-1.png)

**Spectral Density**: Understand frequency domain behavior:

``` r
# Compute and plot spectral density
spec <- spectrald(result_var$model, n.f = 128)

plot(spec,
     main = "Spectral Density",
     legend = c("Series 1", "Series 2"))
```

![Spectral density plot showing power across frequencies for two series
in the estimated VAR(1)
model](0_getting_started_files/figure-html/unnamed-chunk-10-1.png)

**Autocovariance Structure**:

``` r
# Compute autocovariances
acov <- autocov(result_var$model, lag.max = 10, type = "correlation")

# Plot correlations
plot(acov, main = "Sample and Model ACF")
```

![Four-panel autocorrelation function plot comparing sample ACF to
model-implied ACF for the bivariate VAR(1)
system](0_getting_started_files/figure-html/unnamed-chunk-11-1.png)

## Where to Go Next

### For Real-World Applications

Check out the **Case Study** vignette
([`vignette("1_case_study")`](https://bfunovits.github.io/RLDM/articles/1_case_study.md))
which walks through: - Economic data analysis (GDP growth and
unemployment) - Comparing multiple estimation methods - Selecting the
best model - Forecasting and model diagnostics

### For More Details

The **Technical Reference** vignette
([`vignette("2_technical_reference")`](https://bfunovits.github.io/RLDM/articles/2_technical_reference.md))
provides: - Detailed explanation of each model class - Mathematical
foundations - Guidance on method selection - Parameter templates and
estimation algorithms

### Key Functions Reference

| Task                     | Function                                                                                                                                                              |
|--------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Estimate AR models**   | [`est_ar()`](https://bfunovits.github.io/RLDM/reference/est_ar.md)                                                                                                    |
| **Estimate ARMA models** | [`est_arma_hrk()`](https://bfunovits.github.io/RLDM/reference/est_arma_hrk.md)                                                                                        |
| **Estimate state space** | [`est_stsp_cca()`](https://bfunovits.github.io/RLDM/reference/subspace-helpers.md), [`est_stsp_ss()`](https://bfunovits.github.io/RLDM/reference/subspace-methods.md) |
| **Generate forecasts**   | [`predict()`](https://bfunovits.github.io/RLDM/reference/predict.md)                                                                                                  |
| **Impulse responses**    | [`impresp()`](https://bfunovits.github.io/RLDM/reference/impresp.md)                                                                                                  |
| **Spectral analysis**    | [`spectrald()`](https://bfunovits.github.io/RLDM/reference/spectrald.md)                                                                                              |
| **Model diagnostics**    | [`autocov()`](https://bfunovits.github.io/RLDM/reference/autocov.md), [`pm_test()`](https://bfunovits.github.io/RLDM/reference/pm_test.md)                            |

## Quick Tips

1.  **Start simple**: Begin with
    [`est_ar()`](https://bfunovits.github.io/RLDM/reference/est_ar.md)
    to understand your data before moving to more complex models
2.  **Check diagnostics**: Always plot residuals and ACF to verify your
    model
3.  **Use AIC**: The automatic order selection in
    [`est_ar()`](https://bfunovits.github.io/RLDM/reference/est_ar.md)
    and
    [`est_arma_hrk()`](https://bfunovits.github.io/RLDM/reference/est_arma_hrk.md)
    works well in most cases
4.  **Compare models**: Use
    [`compare_estimates()`](https://bfunovits.github.io/RLDM/reference/compare_estimates.md)
    to systematically compare methods
5.  **Multivariate**: RLDM really shines with multivariate data where
    relationships matter

## References

(Scherrer and Deistler 2019)

Scherrer, Wolfgang, and Manfred Deistler. 2019. “Chapter 6 - Vector
Autoregressive Moving Average Models.” In *Conceptual Econometrics Using
r*, edited by Hrishikesh D. Vinod and C. R. Rao, 41:145–91. Handbook of
Statistics. Elsevier.
https://doi.org/<https://doi.org/10.1016/bs.host.2019.01.004>.

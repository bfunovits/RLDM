# Diagnostics for Particle Filter Results

Plotting and diagnostic functions for particle filter output.

## Usage

``` r
# S3 method for class 'pfilter'
plot(x, type = c("states", "ess", "weights", "likelihood"), ...)
```

## Arguments

- x:

  Particle filter result object from
  [`pf()`](https://rdrr.io/r/stats/Fdist.html).

- type:

  Type of diagnostic plot: `"states"` (filtered states with credibility
  intervals), `"ess"` (effective sample size over time), `"weights"`
  (histogram of final weights), `"likelihood"` (log-likelihood
  contributions over time).

- ...:

  Additional arguments passed to plotting functions.

## Value

Plot (invisibly returns `x`).

## Examples

``` r
# See examples in pf()
```

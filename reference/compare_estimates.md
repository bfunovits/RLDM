# Compare Estimated Models

This utility function computes a number of statistics which may be used
to compare/evaluate a set of estimated models. In particular the
function returns the log Likelihood (`ll`), the Akaike Information
Criterion (`AIC`), the Bayes Information Criterion (`BIC`), the Final
Prediction Error (`FPE`) and the p-values of a Portmanteau test for
serial correlation of the residuals.

## Usage

``` r
compare_estimates(estimates, y, n.lags = NULL)
```

## Arguments

- estimates:

  (named) list of estimates. Each slot should contain a list with slots
  `$model` (the estimated model) and `$n.par` the corresponding number
  of (free) parameters of the model (class).

- y:

  (N-by-m)-dimensional matrix with the observed data (or an object which
  may be coerced to a matrix with `as.matrix(y)`).

- n.lags:

  number of lags for the Portmantaeu test for serial correlation of the
  residuals, see also
  [`pm_test()`](https://bfunovits.github.io/RLDM/reference/pm_test.md).

## Value

Matrix with the computed statistics for the estimated models. This
matrix has attributes `m`, `n.obs` and `n.lags`.

## Details

The concentrated, conditional (scaled) log Likelihood \$\$ll = -(1/2)(m
\ln(2\pi) + m + \ln\det S + 2 \ln\det (k_0)\$\$ is computed with
`ll(model, y, skip = 0, concentrated = TRUE)`, see
[`ll()`](https://bfunovits.github.io/RLDM/reference/ll.md). Here \\S\\
denotes the sample covariance of the residuals of the model.

The information criteria are \$\$AIC = -2 ll + (2/N) \kappa\$\$ \$\$BIC
= -2 ll + (\ln(N)/N) \kappa\$\$ where \\\kappa\\ denotes the respective
number of free parameters.

The Final Prediction Error is \$\$FPE =
\det(S)\frac{N+\kappa}{N-\kappa}\$\$

For the portmanteau test, see
[`pm_test()`](https://bfunovits.github.io/RLDM/reference/pm_test.md). If
the number of lags is not specified then the procedure choses a default
value based on the sample size. (This values is attached to the output
as attribute).

Note that the procedure (re)evaulates these measures, even if the
estimates contain this information (e.g. the residuals or the log
likelihood may be stored in the corresponding list). The reason is to
have a common data set and a common evaluation procedure for all
estimates.

Typically the data `y` is the "estimation data set", i.e. the data which
has been used to estimate the models. However, one may also pass a new
data set to the procedure.

## Examples

``` r
# Generate some data
set.seed(123)
model = stspmod(sys = stsp(A = matrix(c(0.5, 0.2, 0, 0.3), 2, 2),
                           B = matrix(c(1, 0.5), 2, 1),
                           C = matrix(c(1, 0), 1, 2),
                           D = matrix(1, 1, 1)),
                sigma_L = matrix(1, 1, 1))
y = sim(model, n.obs = 100)$y

# Create two different estimates (simulated for example)
estimate1 = list(model = model, n.par = 4)
estimate2 = list(model = stspmod(sys = stsp(A = matrix(c(0.4, 0.1, 0, 0.35), 2, 2),
                                           B = matrix(c(1.1, 0.4), 2, 1),
                                           C = matrix(c(0.9, 0), 1, 2),
                                           D = matrix(1, 1, 1)),
                                    sigma_L = matrix(1.2, 1, 1)),
                 n.par = 4)

# Compare estimates
estimates = list("Estimate 1" = estimate1, "Estimate 2" = estimate2)
comparison = compare_estimates(estimates, y, n.lags = 5)
comparison
#>            #par        ll      AIC      BIC       FPE    PM test
#> Estimate 1    4 -1.327622 2.735244 2.839451 0.9024965 0.03727574
#> Estimate 2    4 -1.319477 2.718955 2.823161 0.8879145 0.08825713
#> attr(,"m")
#> [1] 100
#> attr(,"n.obs")
#> [1] 100
#> attr(,"n.lags")
#> [1] 5
```

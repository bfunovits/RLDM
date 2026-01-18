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

# (Multivariate) Recursive Least Squares

This function implements the Recursive Least Squares (RLS) algorithm
with exponentially down-weighted past of fixed window size. It allows
for multivariate regressions, i.e. for multiple left-hand-side variables
`y`. However, the algorithm is formulated such that the regressors are
the same for each individual component of the LHS variable, i.e. \$\$Y =
XB + U\$\$ where

- `Y` and `U` are \\(T \times n\_{y})\\ dimensional,

- `X` is \\(T \times n\_{x})\\ dimensional,

- `B` is \\(n\_{x} \times n\_{y})\\ dimensional,

In a VAR(p) context where the endogenous variable `y` is
`n`-dimensional, this corresponds to dimensions \\(T \times n)\\, \\(T
\times n p)\\, \\( n p \times n)\\ for `Y`, `X`, `B` respectively. The
fact that we use the same regressors for each component of the LHS
simplifies the updating formula in RLS a bit.

## Usage

``` r
rls_exp_cpp(Y, X, r, n_init, allow_neg = TRUE, debug_flag = FALSE)

rls_window_cpp(Y, X, ws, allow_neg = TRUE, debug_flag = FALSE)
```

## Arguments

- Y, X:

  Matrices of doubles. Left-hand-side and right-hand-side variables in
  possibly multivariate regression. The number of rows corresponds to
  the number of observations

- r:

  Matrix of doubles of dimension \\(n\_{r} \times 1)\\. Contains
  forgetting factors. Note that each observation has the same weight
  (i.e. different components )

- n_init:

  Integer. Number of observations used for initial estimate of
  regression coefficients. Default is set to twice the number of
  regressors in `X`. Note that forecasts and forecast errors are only
  produced AFTER the initial estimate.

- allow_neg:

  Boolean. Default set to true. If false, negative forecasts are not
  allowed and set to zero.

- debug_flag:

  Boolean. Default set to false. If true, output is printed to the
  console.

- ws:

  Integer. Fixed number of observation to be used to calculate
  estimates. Number of observations used for initial estimate of
  regression coefficients.

## Value

List containing

- `y_pred`: Cube of doubles of dimension \\(n\_{obs} \times n_y \times
  n\_{r})\\. Contains (one-step-ahead) predictions for all forgetting
  factors. In case of windowed RLS, the third dimension (containing
  forgetting factors) is discarded (such that the output is a matrix).

- `fe_honest:` Cube of doubles of dimension \\(n\_{obs} \times n_y
  \times n\_{r})\\. Contains (one-step-ahead) forecast errors for all
  forgetting factors. In case of windowed RLS, the third dimension
  (containing forgetting factors) is discarded (such that the output is
  a matrix).

## Details

Inputs are not checked, so this should be done within R before calling
this function.

The main reference is Chapter 4 in Young (2012).

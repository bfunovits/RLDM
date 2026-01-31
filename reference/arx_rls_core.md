# RLS function

This function implements the Recursive Least Squares (RLS) algorithm
with exponentially down-weighted past.

## Usage

``` r
arx_rls_core(
  y,
  X,
  r = matrix(c(0.9, 0.95, 0.975, 0.9875, 0.99375, 0.996875, 0.9984375, 1), ncol = 1),
  n_init = NULL,
  start_of_eval = NULL,
  end_of_train = NULL,
  enhance_conv = TRUE
)
```

## Arguments

- y:

  Vector of doubles. Response variable in regression. There must not be
  NAs in this vector.

- X:

  Matrix of doubles. Dimension = (length(y) x maximal number of
  regressors). NAs are not handled separately (i.e. they must be checked
  before this function is called).

- r:

  Matrix of doubles (a column vector of dimension \\(
  \text{forgetting\\factors} \times 1)\\ containing the forgetting
  factors.

- n_init:

  Integer. Number of observations used for initial estimate of beta. If
  no value provided, at least 21 observations (3 weeks) or 3 times the
  number of regressors is used.

- start_of_eval:

  Integer. Starting value of evaluation period for honest prediction
  error (if we start too early, there are bad initial estimations
  involved)

- end_of_train:

  Integer. Index specifying the end of the training set. (Afterwards the
  forecasting period starts.) Important for calculating the honest
  prediction error (otherwise it would not be out-of-sample...).
  comment_bf: If one wants to use arx_rls_core() outside the automated
  forecasting framework, set `end_of_train` to `length(y)`.

- enhance_conv:

  Boolean. Indicates whether the convergence enhancing factor as
  described in Young (2011) page 55.

## Value

List containing

- `y_pred`: vector of predictions (vector of same length as input y)

- `fev_honest`: double. honest prediction error (calculated from
  `end_of_train - start_of_eval + 1` observations)

- `forgetting`: double. forgetting factor which produced the minimal
  honest prediction error

## Examples

``` r
# Generate simple time series data
set.seed(123)
n <- 100
x1 <- rnorm(n)
x2 <- rnorm(n)
X <- cbind(x1, x2)
y <- 2*x1 + 3*x2 + rnorm(n, sd = 0.5)
y <- matrix(y, ncol = 1)  # Convert to matrix

# Run RLS with default forgetting factors
result <- arx_rls_core(y, X, n_init = 20, start_of_eval = 25, end_of_train = 80)
str(result)
#> List of 3
#>  $ y_pred    : num [1:100] NA NA NA NA NA NA NA NA NA NA ...
#>  $ fev_honest: num 0.187
#>  $ forgetting: num 1
```

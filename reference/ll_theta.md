# Log-likelihood Given Deep Parameters

**\[superseded\]** See
[`ll_FUN()`](https://bfunovits.github.io/RLDM/reference/ll_FUN.md). The
template `template` is filled with the deep parameters in `th`.
Subsequently, the S3 method
[`ll()`](https://bfunovits.github.io/RLDM/reference/ll.md) is called for
the class provided in the template and the value of the **scaled**
log-likelihood function is returned, see
[`ll()`](https://bfunovits.github.io/RLDM/reference/ll.md).

## Usage

``` r
ll_theta(th, template, y, which, ...)
```

## Arguments

- th:

  Vector of deep parameter

- template:

  A model template, see [model
  structures](https://bfunovits.github.io/RLDM/reference/model_structures.md).

- y:

  Data sample given as \\(N,m)\\ dimensional matrix, or a "time series"
  object (in the sense that `as.matrix(y)` should return an
  \\(N,m)\\-dimensional numeric matrix). Missing values (`NA`, `NaN` and
  `Inf`) are **not** supported.

- which:

  (character) which likelihood to compute.

- ...:

  Not used.

## Value

Value of log-likelihood for a given deep/free parameter vector `th` and
a model structure defined via `template`. Note that this function simply
calls `ll(fill_template(th, template), y, which, ...)`.

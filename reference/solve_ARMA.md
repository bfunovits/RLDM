# Solve ARMA system

Compute the outputs of ARMA(p, q) systems of the form \$\$y_t = a_1
y\_{t-1} + ... + a_p y\_{t-p} + b_0 u_t + \cdots + b_q u\_{t-q}\$\$

## Usage

``` r
solve_ARMA_R(a, b, u, y, t0)
```

## Arguments

- a:

  \\(m, mp)\\ matrix \\(a_p,...,a_1)\\.

- b:

  \\(m, n(q+1))\\ matrix \\(b_0,...,b_q\\.

- u:

  \\(n, N)\\ matrix with the inputs \\(u_1,...,u_N\\.

- y:

  \\(m, N)\\ matrix with the outputs \\(y_1,...,y_N\\.

- t0:

  integer, start iteration at t=t0.

## Value

The R implementation `solve_ARMA_R` returns the matrix `y` with the
computed outputs. The RcppArmadillo implementation returns `NULL` but
**overwrites** the input argument `y`!

## Details

Values \\y_t\\, \\u_t\\ for \\t\leq 0\\ are implicitly set to be zero.
However, if we start the iteration with some \\t_0\>1\\ we can enforce
non-zero initial values.

The routines are used internally and hence do **not** check their
arguments. We require \\m \> 0\\, \\p \geq 0\\, \\n \geq 0\\, \\(q+1)
\geq 0\\ and \\1 \leq t_0 \leq N\\. Note also that the RcppArmadillo
implementation **overwrites** the input argument `y`. Use this procedure
with care!

Note the non standard arguments: The order of the AR coefficients is
reversed. The data matrices are organized column-wise (to avoid memory
shuffling)!

## Examples

``` r
# Basic example
result <- solve_ARMA_R()
#> Error in solve_ARMA_R(): argument "y" is missing, with no default
result
#> Error: object 'result' not found
```

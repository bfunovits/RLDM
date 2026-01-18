# Solve RMFD system for given inputs

Compute the outputs \\y_t\\ of RMFD(p, q) systems of the form \$\$y_t =
d_0 v_t + d_1 v\_{t-1} + \cdots + d_q v\_{t-q}\$\$ where \$\$v_t + c_1
v\_{t-1} + \cdots + c_p v\_{t-p} = u_t\$\$ for given inputs \\u_t\\
which are contained column-wis in the data matrix `data_input`, see
below.

## Usage

``` r
solve_RMFD_R(polm_c, polm_d, data_input, t0 = 1)
```

## Arguments

- polm_c, polm_d:

  polm objects. Describe jointly RMFD. \\c(z)\\ is square, has the
  identity as zero-lag coefficient, and its coefficients will be
  reversed in the procedure: \\(c_p,...,c_1)\\. \\d(z)\\ might be tall,
  and its zero-lag coefficient matrix is in general free.

- data_input:

  \\(n, N)\\ matrix with the inputs \\(u_1,...,u_N\\.

- t0:

  integer, start iteration at t=t0.

## Value

The R implementation `solve_RMFD_R()` returns the matrix `y` with the
computed outputs. \$\$y_t = d(z) c(z)^{-1} u_t\$\$ The RcppArmadillo
implementation
[`solve_rmfd_cpp()`](https://bfunovits.github.io/RLDM/reference/solve_rmfd_cpp.md)
returns `NULL` but **overwrites** the input argument. Note that the
RcppArmadillo implementation has a different user interface (it is
intended for internal use only).

## Details

Values \\y_t\\, \\u_t\\ for \\t\leq 0\\ are implicitely set to be zero.
However, if we start the iteration with some \\t_0\>1\\ we can enforce
non-zero initial values.

The routines are used internally and hence do **not** check their
arguments. We require the number of outputs to be positive \\m \> 0\\,
the number of inputs to be non-negative \\n \geq 0\\, the degree of
\\c(z)\\ %\>% to be non-negative \\p \geq 0\\, the degree of \\d(z)\\ is
unrestricted, i.e. \\(q+1) \geq 0\\, and for the starting value and the
sample size \\1 \leq t_0 \leq N\\ holds. Note also that the
RcppArmadillo implementation **overwrites** the input argument `y`. Use
this procedure with care!

Note the non standard arguments: The polynomial matrices \\c(z)\\ and
\\d(z)\\ are saved as "wide" matrices. The order of the coefficients in
\\c(z)\\ is reversed and there is no \\c_0\\ coefficient (because it is
required to be the identity matrix). The order of the coefficients in
\\d(z)\\ is as usual, \\d_0\\ is available too. The data matrices are
organized columnwise (to avoid memory shuffling)!

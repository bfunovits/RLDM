# Simulating Output from an RMFD Model (or obtain residuals)

This RcppArmadillo function calculates for given inputs `data_in` of
dimension \\(n \times nobs)\\, where \\nobs\\ is the sample size, the
outputs of dimension \\m\\. Note that data matrices are "transposed" in
the sense that every column corresponds to one observation because of
memory management. `data_out` is thus of dimension \\(m x nobs)\\. This
function is intended for internal use and thus arguments are not
checked.

## Usage

``` r
solve_rmfd_cpp(poly_inv, poly_fwd, data_in, data_out, t0)
```

## Arguments

- poly_inv:

  Matrix of dimension \\(n \times n p)\\, representing a square matrix
  polynomial \\c(z)\\ with \\c_0\\ equal to the identity matrix (and
  therefore not stored). The coefficients need to be in **reverse**
  direction, i.e. \\(c_p, ... , c_1)\\, where \\p\\ denotes the degree
  of \\c(z)\\.

- poly_fwd:

  Matrix of dimensions \\(m \times n(q+1))\\, representing a (possibly
  tall) matrix polynomial \\d(z)\\ of dimension \\(m \times n)\\, where
  \\m \geq n\\. The coefficient are stored "as usual" and including
  \\d_0\\, i.e. \\(d_0, d_1, ... , d\_{q-1}, d\_{q})\\, where \\q\\
  denotes the degree of \\d(z)\\.

- data_in:

  Matrix of dimension \\(n \times n_obs)\\, i.e. \\(u_1, ..., u_T)\\.
  Inputs to the RMFD system.

- data_out:

  Matrix of dimension \\(m \times n_obs)\\, i.e. \\(y_1, ..., y_T)\\.
  Outputs of the RMFD system. Initially zero and will be
  **overwritten**.

- t0:

  Integer. Time index from which we should start calculating a solution.
  Usually equal to 1.

## Value

`data_out` is overwritten with the outputs of the RMFD system.

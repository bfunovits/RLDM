# Solve (linear) Difference Equations

The procedure `solve_de()` solves the difference equations associated to
(V)ARMA models \$\$a_0 y_t + a_1 y\_{t-1} + \cdots + a_p y\_{t-p} = b_0
u_t + b_1 u\_{t-1} + ... b_1 u\_{t-q}\$\$ or state space models
\$\$a\_{t+1} = A a_t + B u_t \mbox{ and } y_t = C a_t + D u_t.\$\$

## Usage

``` r
solve_de(sys, u, ...)

# S3 method for class 'stsp'
solve_de(sys, u, a1 = NULL, ...)

# S3 method for class 'lmfd'
solve_de(sys, u, u0 = NULL, y0 = NULL, ...)

# S3 method for class 'rmfd'
solve_de(sys, u, u0 = NULL, y0 = NULL, ...)

solve_inverse_de(sys, y, ...)

# S3 method for class 'stsp'
solve_inverse_de(sys, y, a1 = NULL, ...)

# S3 method for class 'lmfd'
solve_inverse_de(sys, y, u0 = NULL, y0 = NULL, ...)

# S3 method for class 'rmfd'
solve_inverse_de(sys, y, u0 = NULL, y0 = NULL, ...)
```

## Arguments

- sys:

  [`rationalmatrices::lmfd()`](https://bfunovits.github.io/rationalmatrices/reference/lmfd.html)
  or
  [`rationalmatrices::stsp()`](https://bfunovits.github.io/rationalmatrices/reference/stsp.html)
  object which describes the difference equation.

- u:

  \\(N,n)\\ matrix with the noise (\\u_t\\, \\t=1,...,N\\).

- ...:

  not used.

- a1:

  \\m\\ dimensional vector with the initial state \\a_1\\. If `a1=NULL`
  then a zero vector is used.

- u0:

  \\(h,n)\\ dimensional matrix with starting values for the disturbances
  \\(u\_{1-h}, \ldots, u\_{-1}, u_0)\\. Note that the last row
  corresponds to \\u_0\\, the last but one row to \\u\_{-1}\\ and so on.
  If \\h\>q\\ then only the last \\q\\ rows of `u0` are used. In the
  case \\h\<q\\ the "missing" initial values are set to zero vectors.  
  The default value `u0=NULL` sets all initial values \\u_t\\, \\t \leq
  0\\ equal to zero vectors.

- y0:

  \\(h,m)\\ dimensional matrix with starting values for the outputs
  \\(y\_{1-h}, \ldots, y\_{-1}, y_0)\\. This (optional) parameter is
  interpreted analogously to `u0`.

- y:

  \\(N,m)\\ matrix with the outputs (\\y_t\\, \\t=1,...,N\\).

## Value

List with slots

- y:

  \\(N,n)\\ matrix with the (computed) outputs.

- u:

  \\(N,n)\\ matrix with the (computed) noise.

- a:

  \\(N+1,n)\\ matrix with the (computed) states (\\a_t\\,
  \\t=1,...,N+1\\). Note that this matrix has (\\N+1\\) rows! This slot
  is only present for state space models.

## Details

`solve_de()` computes the outputs \\y_t\\ for \\t=1,\ldots,N\\ for given
disturbances \\u_t\\ \\t=1,\ldots,N\\. The starting values (\\u_t\\ and
\\y_t\\ for \\t\leq 0\\ for VARMA models and \\a_1\\ for state space
models) may be given as optional arguments. The default is to use zero
vectors.

For the reverse direction, i.e. to reconstruct the disturbances if the
outputs are given, the function `solve_inverse_de` may be used. In this
case the system must be square and the matrix \\D\\ respectively \\b_0\\
must be invertible.

These functions are mainly intended for internal use and hence only some
basic checks on the input parameters are performed.

## Examples

``` r
### generate a random ARMA(2,1) model (with two outputs) #########
model = test_armamod(dim = c(2,2), degrees = c(2,1),
                     digits = 2, bpoles = 1, bzeroes = 1)

# generate random noise sequence (sample size N = 100)
u = matrix(rnorm(100*2), nrow = 100, ncol = 2)

# generate random initial values
u0 = matrix(rnorm(2), nrow = 1, ncol = 2) # u[0]
y0 = matrix(rnorm(2), nrow = 1, ncol = 2) # y[0]

# compute outputs "y[t]"
# note that y0 has only one row, thus y[-1] is set to zero!
data = solve_de(model$sys, u = u, y0 = y0, u0 = u0)

# we can reconstruct the noise "u" from given outputs "y"
data1 = solve_inverse_de(model$sys, y = data$y, u0 = u0, y0 = y0)
all.equal(data$u, data1$u)
#> [1] TRUE

### generate a random state space model (3 outputs and 4 states) ##
model = test_stspmod(dim = c(3,3), s = 4,
                     digits = 2, bpoles = 1, bzeroes = 1)

# generate random noise sequence (sample size N = 100)
u = matrix(rnorm(100*3), nrow = 100, ncol = 3)

# generate random initial state a[1]
a1 = rnorm(4)

# compute outputs "y[t]"
data = solve_de(model$sys, u = u, a1 = a1)

# we can reconstruct the noise "u" from given outputs "y"
data1 = solve_inverse_de(model$sys, y = data$y, a1 = data$a[1,])
all.equal(data$u, data1$u)
#> [1] TRUE

# Basic example
result <- solve_inverse_de()
#> Error in solve_inverse_de(): argument "sys" is missing, with no default
result
#> Error: object 'result' not found
```

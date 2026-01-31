# Model Structures

These tools define and implement **model structures** where the **model
parameters** are represented by an affine function of some **free**
(**deep**) parameters. As an example consider multivariate ARMA models.
The AR coefficients \\a_k\\, the MA coefficients \\b_k\\ and the (left)
square root of the noise covariance matrix, \\L\\ say, are vectorized
and stacked into a (long) parameter vector as \$\$\pi =
(\mbox{vec}(a_1)',\ldots,\mbox{vec}(a_p)',
\mbox{vec}(b_1)',\ldots,\mbox{vec}(b_q)',\mbox{vec}(L)')'\$\$ This
parameter vector then is written as \$\$\pi = h + H\theta\$\$ where
\\\theta\\ represents the **free** parameters. Of course the matrix
\\H\\ is assumed to have full column rank. This parameterization scheme
is quite flexible. In particular, ARMA and state space models in
**echelon form** may be represented by this scheme.  
Templates and the related tools are mainly used for estimation and for
the generation of (random) models for simulations and testing.

## Usage

``` r
model2template(
  model,
  sigma_L = c("as_given", "chol", "symm", "identity", "full_normalized")
)

tmpl_arma_pq(
  m,
  n,
  p,
  q,
  sigma_L = c("chol", "symm", "identity", "full_normalized")
)

tmpl_arma_echelon(
  nu,
  n = length(nu),
  sigma_L = c("chol", "symm", "identity", "full_normalized")
)

tmpl_rmfd_pq(
  m,
  n,
  p,
  q,
  sigma_L = c("chol", "symm", "identity", "full_normalized")
)

tmpl_rmfd_echelon(
  nu,
  m = length(nu),
  sigma_L = c("chol", "symm", "identity", "full_normalized")
)

tmpl_stsp_full(
  m,
  n,
  s,
  sigma_L = c("chol", "symm", "identity", "full_normalized")
)

tmpl_stsp_ar(m, p, sigma_L = c("chol", "symm", "identity", "full_normalized"))

tmpl_stsp_echelon(
  nu,
  n = length(nu),
  sigma_L = c("chol", "symm", "identity", "full_normalized")
)
```

## Arguments

- model:

  [`armamod()`](https://bfunovits.github.io/RLDM/reference/armamod.md),
  [`rmfdmod()`](https://bfunovits.github.io/RLDM/reference/rmfdmod.md)
  or
  [`stspmod()`](https://bfunovits.github.io/RLDM/reference/stspmod.md)
  object.

- sigma_L:

  (character string) determines the form of the (left) square root of
  the noise covariance \\\Sigma\\. The choice `"chol"` gives a lower
  triangular matrix, `"symm"` gives a symmetric matrix and `"identity"`
  corresponds to a fixed (identity) matrix. The procedure
  `model2template` has an additional option `"as_given"` which means
  that the structure of the square root `sigma_L` is completely
  determined by the `sigma_L` slot of the given model.

- m:

  output dimension

- n:

  input dimension (= number of shocks = dimension of the noise process)

- p:

  order of the AR polynomial for ARMA models, respectively the order of
  the right factor polynomial \\c(z)\\ in an RMFD model.

- q:

  order of the MA polynomial for ARMA models, respectively the order of
  the left factor polynomial \\d(z)\\ in an RMFD model.

- nu:

  vector of Kronecker indices. For ARMA models the Kronecker indices
  describe the basis rows and for RMFD models the basis columns of the
  Hankel matrix of the impulse response coefficients.

- s:

  state dimension for state space models.

## Value

The functions `model2template` and `tmpl_***` return a model template.

## Details

The functions `model2template` and `tmpl_***` generate model "templates"
which represent certain model structures, where the model parameters are
affine functions of some **free**, respectively **deep**, parameters.

The template contains information about the model structure in the
following slots

1.  `h, H` represent the vector \\h\\ and the marix \\H\\ as described
    above. (The vector \\\pi\\ of stacked model parameters is
    represented as \\\pi = h +H \theta\\ where \\\theta\\ is a vector of
    *deep* parameters.)

2.  `class=["armamod"|"rmfdmod"|"stspmod"]`: determines whether the
    template parametrizes ARMA, RMFD or state space models.

3.  `order`: an integer vector which contains the dimensions and orders
    of the model. For ARMA and RMFD models `order = c(m,n,p,q)` and for
    state space models `order = c(m,n,s)`.

4.  `n.par`: the number of *free* parameters, i.e. the dimension of the
    vector \\\theta\\.

5.  `nu`: This optional slot contains the Kronecker indices \\\nu\\.

**model2template:**

The function `model2template()` takes an
[`armamod()`](https://bfunovits.github.io/RLDM/reference/armamod.md),
[`rmfdmod()`](https://bfunovits.github.io/RLDM/reference/rmfdmod.md), or
[`stspmod()`](https://bfunovits.github.io/RLDM/reference/stspmod.md)
object where the free parameters are coded as `NA`'s, `NaN`'s or `Inf`'s
and constructs a corresponding model template.

For the parametrization of the (left) square root, \\L\\ say, of the
noise covariance \\\Sigma = LL'\\ the following choices are possible: In
the case `sigma_L = "as_given"` the slot `model$sigma_L` of the given
`model` is used to construct the template: `NA` entries are considered
as free and all other entries as fixed. For the choice
`sigma_L = "chol"` first all entries of `model$sigma_L` above the
diagonal are set to zero and then the template is constructed as above.
In the case `sigma_L = "symm"` the matrix `model$sigma_L` is first
replaced by a symmetric one and then the template is constructed
(according to the `NA`'s) such that the square root `L=sigma_L` is
always symmetric. The choice `sigma_L = "identity"` sets the matrix
`L = sigma_L` to the identity matrix. Finally, the choice
`sigma_L = "full_normalized"` sets the diagonal elements equal to ones
and all other elements to NAs in `L = sigma_L`.

**tmpl\_***:*\*

The functions `tmpl_***` implement the following model structures:

- tmpl_arma_pq:

  ARMA models
  ([`armamod()`](https://bfunovits.github.io/RLDM/reference/armamod.md))
  with prescribed orders \\(p,q\\).

- tmpl_arma_echelon:

  ARMA models
  ([`armamod()`](https://bfunovits.github.io/RLDM/reference/armamod.md))
  in echelon form, with given Kronecker indices \\\nu\\.

- tmpl_rmfd_pq:

  RMFD models
  ([`rmfdmod()`](https://bfunovits.github.io/RLDM/reference/rmfdmod.md))
  with prescribed orders \\(p,q\\).

- tmpl_rmfd_echelon:

  RMFD models
  ([`rmfdmod()`](https://bfunovits.github.io/RLDM/reference/rmfdmod.md))
  in echelon form, with Kronecker indices \\\nu\\. Note that for RMFD
  models the Kronecker indices refer to the basis of the *column space*
  of the Hankel matrix of the impulse response coefficients.

- tmpl_stsp_full:

  Fully parametrized state space models
  ([`stspmod()`](https://bfunovits.github.io/RLDM/reference/stspmod.md))
  with given state space dimension \\s\\, i.e. a models where each entry
  in the matrices \\A,B,C\\ is considered non-fixed.

- tmpl_stsp_echelon:

  State space models
  ([`stspmod()`](https://bfunovits.github.io/RLDM/reference/stspmod.md))
  in echelon form, with given Kronecker indices \\\nu\\.

- tmpl_state space_ar:

  State space model representations
  ([`stspmod()`](https://bfunovits.github.io/RLDM/reference/stspmod.md))
  of AR models with given order \\p\\. Here only the "square" case
  \\m=n\\ is implemented.

For these model structures the impulse response (transfer function) is
scaled such that the \\(m,n)\\-dimensional lag zero coefficient, \\k_0\\
say, is of the form

- \\m=n\\:

  \\k_0\\ is the \\m\\-dimensional identity matrix.

- \\m\<n\\:

  The first \\m\\ columns of \\k_0\\ form the \\m\\-dimensional identity
  matrix and the remaining columns are zero.

- \\m\>n\\:

  The first \\n\\ rows of \\k_0\\ form the \\n\\-dimensional identity
  matrix and the remaining rows are *free*.

For the parametrization of the (left) square root \\L\\ of the noise
covariance \\\Sigma = LL'\\ the following choices are possible: For
`sigma_L="chol"` the matrix \\L\\ is lower triangular (all entries on
and below the main diagonal are considered as free and the entries above
the diagonal are zero). For `sigma_L="symm"` the matrix \\L\\ is
symmetric (all entries on and below the main diagonal are considered as
free and the entries above the diagonal are such that \\L=L'\\ holds).
For `sigma_L="identity"` the matrix \\L\\ is *fixed* to the
\\n\\-dimensional identity matrix. For `sigma_L="full_normalized"` the
diagonal elements of the matrix \\L\\ are *fixed* to ones and all other
elements are free.

## See also

[`r_model()`](https://bfunovits.github.io/RLDM/reference/r_model.md),
[`fill_template()`](https://bfunovits.github.io/RLDM/reference/fill_template.md),
[`ll()`](https://bfunovits.github.io/RLDM/reference/ll.md),
[`ll_theta()`](https://bfunovits.github.io/RLDM/reference/ll_theta.md)
and other estimation procedures

## Examples

``` r
# ######################################################
# construct a template from a model
# ######################################################

# Let us consider scalar ARMA(5,1) models
# for quarterly data with a strong seasonal component.
# In order to have parsimonious models we want a[2]=a[3]=0:
model = armamod(lmfd(a = c(1,NA,0,0,NA,NA), b = c(1,NA)))
tmpl = model2template(model)

# Let's see how the "free" parameters are mapped to the model parameters
print(cbind(tmpl$h, tmpl$H))
#>       [,1] [,2] [,3] [,4] [,5]
#>  [1,]    1    0    0    0    0
#>  [2,]    0    1    0    0    0
#>  [3,]    0    0    0    0    0
#>  [4,]    0    0    0    0    0
#>  [5,]    0    0    1    0    0
#>  [6,]    0    0    0    1    0
#>  [7,]    1    0    0    0    0
#>  [8,]    0    0    0    0    1
#>  [9,]    1    0    0    0    0
th = -(1:tmpl$n.par)
fill_template(th, tmpl)
#> ARMA model [1,1] with orders p = 5 and q = 1
#> AR polynomial a(z):
#>      z^0 [,1] z^1 [,1] z^2 [,1] z^3 [,1] z^4 [,1] z^5 [,1]
#> [1,]        1       -1        0        0       -2       -3
#> MA polynomial b(z):
#>      z^0 [,1] z^1 [,1]
#> [1,]        1       -4
#> Left square root of noise covariance Sigma:
#>      u[1]
#> u[1]    1

# Generate a random model with this structure
th0 = rnorm(tmpl$n.par, sd = 0.1)
model = fill_template(th0, tmpl)

# Extract the "free" parameters from the model
th = extract_theta(model, tmpl)
all.equal(th, th0)
#> [1] TRUE

# This model structure fixes sigma_L = 1.
# If we change sigma_L = 2 then the model does not fit to the template.
model$sigma_L = 2
# the default choice on_error = 'ignore', tells
# extract_theta to ignore this misfit:
th = extract_theta(model, tmpl, on_error = 'ignore')
# with on_error = 'warn' we get a warning and
# with on_error = 'stop' would throw an error.
th = extract_theta(model, tmpl, on_error = 'warn')
#> Warning: model does not match template. max(abs(res)) = 1
# We may also "ignore" sigma_L
th = extract_theta(model, tmpl, on_error = 'stop', ignore_sigma_L=TRUE)

# If the orders/class of template and model does not fit
if (FALSE) { # \dontrun{
model = armamod(lmfd(a = c(1,1), b = c(1,1)))
extract_theta(model, tmpl)
model = stspmod(stsp(D = 1))
extract_theta(model, tmpl)
} # }

# ######################################################
# the parameter "sigma_L"
# ######################################################

# consider a state space model (with 1 state) for a 3-dimensional process
model = stspmod(stsp(A = 1, B = c(1,0,0), C = c(1,1,1), D = diag(3)))

# We may specify an arbitrary structure for the left square root (L = sigma_L)
# of the noise covariance Sigma. Any NA entry is considered as a "free" parameter.
L = matrix(c(0, NA, 1, 0, 2, 3, NA, 1, 1), nrow = 3, ncol = 3)
L
#>      [,1] [,2] [,3]
#> [1,]    0    0   NA
#> [2,]   NA    2    1
#> [3,]    1    3    1
# L has 2 NA entries and thus we get a model structure with 2 free parameters.
model$sigma_L = L

tmpl = model2template(model, sigma_L = 'as_given')
th = -(1:tmpl$n.par)
fill_template(th, tmpl)
#> state space model [3,3] with s = 1 states
#>      s[1] u[1] u[2] u[3]
#> s[1]    1    1    0    0
#> x[1]    1    1    0    0
#> x[2]    1    0    1    0
#> x[3]    1    0    0    1
#> Left square root of noise covariance Sigma:
#>      u[1] u[2] u[3]
#> u[1]    0    0   -2
#> u[2]   -1    2    1
#> u[3]    1    3    1

# The choice sigma_L = 'chol' forces L to be lower triangular.
# In the case considered here, we get a model structure with 1 free parameter.
tmpl = model2template(model, sigma_L = 'chol')
th = -(1:tmpl$n.par)
fill_template(th, tmpl)
#> state space model [3,3] with s = 1 states
#>      s[1] u[1] u[2] u[3]
#> s[1]    1    1    0    0
#> x[1]    1    1    0    0
#> x[2]    1    0    1    0
#> x[3]    1    0    0    1
#> Left square root of noise covariance Sigma:
#>      u[1] u[2] u[3]
#> u[1]    0    0    0
#> u[2]   -1    2    0
#> u[3]    1    3    1

# The choice sigma_L = 'symm' forces L = sigma_L to be symmetric.
# In the case considered here we thus get a model structure with 2 free parameters.
tmpl = model2template(model, sigma_L = 'symm')
th = -(1:tmpl$n.par)
fill_template(th, tmpl)
#> state space model [3,3] with s = 1 states
#>      s[1] u[1] u[2] u[3]
#> s[1]    1    1    0    0
#> x[1]    1    1    0    0
#> x[2]    1    0    1    0
#> x[3]    1    0    0    1
#> Left square root of noise covariance Sigma:
#>      u[1] u[2] u[3]
#> u[1]    0   -1   -2
#> u[2]   -1    2    2
#> u[3]   -2    2    1

# The choice sigma_L = 'identity' set L equal to the identity matrix,
# i.e. sigma_L is fixed.
tmpl = model2template(model, sigma_L = 'identity')
th = numeric(0)
fill_template(th, tmpl)
#> state space model [3,3] with s = 1 states
#>      s[1] u[1] u[2] u[3]
#> s[1]    1    1    0    0
#> x[1]    1    1    0    0
#> x[2]    1    0    1    0
#> x[3]    1    0    0    1
#> Left square root of noise covariance Sigma:
#>      u[1] u[2] u[3]
#> u[1]    1    0    0
#> u[2]    0    1    0
#> u[3]    0    0    1
tmpl$n.par # there are no free parameters: tmpl$n.par = 0
#> [1] 0

# The choice sigma_L = 'full_normalized' sets the diagonal elements of L equal to ones,
# and leaves all other elements free.
tmpl = model2template(model, sigma_L = 'full_normalized')
th = -(1:tmpl$n.par)
fill_template(th, tmpl)
#> state space model [3,3] with s = 1 states
#>      s[1] u[1] u[2] u[3]
#> s[1]    1    1    0    0
#> x[1]    1    1    0    0
#> x[2]    1    0    1    0
#> x[3]    1    0    0    1
#> Left square root of noise covariance Sigma:
#>      u[1] u[2] u[3]
#> u[1]    1   -3   -5
#> u[2]   -1    1   -6
#> u[3]   -2   -4    1


# ######################################################
# ARMA(p,q) models
# ######################################################

m = 2 # output dimension
p = 1 # AR order
q = 1 # MA order

# model structure with lower triangular sigma_L
tmpl = tmpl_arma_pq(m, n = m, p, q, sigma_L = "chol")
th = rnorm(tmpl$n.par)
th = -(1:tmpl$n.par)
fill_template(th, tmpl)
#> ARMA model [2,2] with orders p = 1 and q = 1
#> AR polynomial a(z):
#>      z^0 [,1]  [,2] z^1 [,1]  [,2]
#> [1,]        1     0       -1    -3
#> [2,]        0     1       -2    -4
#> MA polynomial b(z):
#>      z^0 [,1]  [,2] z^1 [,1]  [,2]
#> [1,]        1     0       -5    -7
#> [2,]        0     1       -6    -8
#> Left square root of noise covariance Sigma:
#>      u[1] u[2]
#> u[1]   -9    0
#> u[2]  -10  -11

# model structure with symmetric sigma_L
tmpl = tmpl_arma_pq(m, n = m, p, q, sigma_L = "symm")
fill_template(th, tmpl)
#> ARMA model [2,2] with orders p = 1 and q = 1
#> AR polynomial a(z):
#>      z^0 [,1]  [,2] z^1 [,1]  [,2]
#> [1,]        1     0       -1    -3
#> [2,]        0     1       -2    -4
#> MA polynomial b(z):
#>      z^0 [,1]  [,2] z^1 [,1]  [,2]
#> [1,]        1     0       -5    -7
#> [2,]        0     1       -6    -8
#> Left square root of noise covariance Sigma:
#>      u[1] u[2]
#> u[1]   -9  -10
#> u[2]  -10  -11

# model structure with sigma_L = I
tmpl = tmpl_arma_pq(m, n = m, p, q, sigma_L = "identity")
# here the number of free paramaters is of course (by 3) smaller
# than for the above model structures!
fill_template(th[1:(length(th)-3)], tmpl)
#> ARMA model [2,2] with orders p = 1 and q = 1
#> AR polynomial a(z):
#>      z^0 [,1]  [,2] z^1 [,1]  [,2]
#> [1,]        1     0       -1    -3
#> [2,]        0     1       -2    -4
#> MA polynomial b(z):
#>      z^0 [,1]  [,2] z^1 [,1]  [,2]
#> [1,]        1     0       -5    -7
#> [2,]        0     1       -6    -8
#> Left square root of noise covariance Sigma:
#>      u[1] u[2]
#> u[1]    1    0
#> u[2]    0    1

# Create a template
tmpl <- tmpl_arma_echelon()
#> Error in tmpl_arma_echelon(): argument "nu" is missing, with no default
tmpl
#> $h
#>  [1] 1 0 0 1 0 0 0 0 1 0 0 1 0 0 0 0 1 0 0 1
#> 
#> $H
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
#>  [1,]    0    0    0    0    0    0    0    0
#>  [2,]    0    0    0    0    0    0    0    0
#>  [3,]    0    0    0    0    0    0    0    0
#>  [4,]    0    0    0    0    0    0    0    0
#>  [5,]    1    0    0    0    0    0    0    0
#>  [6,]    0    1    0    0    0    0    0    0
#>  [7,]    0    0    1    0    0    0    0    0
#>  [8,]    0    0    0    1    0    0    0    0
#>  [9,]    0    0    0    0    0    0    0    0
#> [10,]    0    0    0    0    0    0    0    0
#> [11,]    0    0    0    0    0    0    0    0
#> [12,]    0    0    0    0    0    0    0    0
#> [13,]    0    0    0    0    1    0    0    0
#> [14,]    0    0    0    0    0    1    0    0
#> [15,]    0    0    0    0    0    0    1    0
#> [16,]    0    0    0    0    0    0    0    1
#> [17,]    0    0    0    0    0    0    0    0
#> [18,]    0    0    0    0    0    0    0    0
#> [19,]    0    0    0    0    0    0    0    0
#> [20,]    0    0    0    0    0    0    0    0
#> 
#> $class
#> [1] "armamod"
#> 
#> $order
#> [1] 2 2 1 1
#> 
#> $n.par
#> [1] 8
#> 

# Use the template with fill_template()
# filled <- fill_template(tmpl, theta = rnorm(tmpl$n.par))

# ######################################################
# RMFD(p,q) models y[t] = d(z) c(z)^-1 e[t]
# ######################################################

m = 2 # output dimension
p = 1 # order of c(z)
q = 1 # order of d(z)

# model structure with lower triangular sigma_L
tmpl = tmpl_rmfd_pq(m, n = m, p, q, sigma_L = "chol")
th = rnorm(tmpl$n.par)
th = -(1:tmpl$n.par)
fill_template(th, tmpl)
#> RMFD model [2,2] with orders p = 1 and q = 1
#> right factor polynomial c(z):
#>      z^0 [,1]  [,2] z^1 [,1]  [,2]
#> [1,]        1     0       -1    -5
#> [2,]        0     1       -2    -6
#> left factor polynomial d(z):
#>      z^0 [,1]  [,2] z^1 [,1]  [,2]
#> [1,]        1     0       -3    -7
#> [2,]        0     1       -4    -8
#> Left square root of noise covariance Sigma:
#>      u[1] u[2]
#> u[1]   -9    0
#> u[2]  -10  -11

# model structure with symmetric sigma_L
tmpl = tmpl_rmfd_pq(m, n = m, p, q, sigma_L = "symm")
fill_template(th, tmpl)
#> RMFD model [2,2] with orders p = 1 and q = 1
#> right factor polynomial c(z):
#>      z^0 [,1]  [,2] z^1 [,1]  [,2]
#> [1,]        1     0       -1    -5
#> [2,]        0     1       -2    -6
#> left factor polynomial d(z):
#>      z^0 [,1]  [,2] z^1 [,1]  [,2]
#> [1,]        1     0       -3    -7
#> [2,]        0     1       -4    -8
#> Left square root of noise covariance Sigma:
#>      u[1] u[2]
#> u[1]   -9  -10
#> u[2]  -10  -11

# model structure with sigma_L = I
tmpl = tmpl_rmfd_pq(m, n = m, p, q, sigma_L = "identity")
# here the number of free paramaters is of course (by 3) smaller
# than for the above model structures!
fill_template(th[1:(length(th)-3)], tmpl)
#> RMFD model [2,2] with orders p = 1 and q = 1
#> right factor polynomial c(z):
#>      z^0 [,1]  [,2] z^1 [,1]  [,2]
#> [1,]        1     0       -1    -5
#> [2,]        0     1       -2    -6
#> left factor polynomial d(z):
#>      z^0 [,1]  [,2] z^1 [,1]  [,2]
#> [1,]        1     0       -3    -7
#> [2,]        0     1       -4    -8
#> Left square root of noise covariance Sigma:
#>      u[1] u[2]
#> u[1]    1    0
#> u[2]    0    1

# Create a template
tmpl <- tmpl_rmfd_echelon()
#> Error in tmpl_rmfd_echelon(): argument "nu" is missing, with no default
tmpl
#> $h
#>  [1] 1 0 0 0 1 0 0 0 0 1 0 0 0 1 0 0 1 0 0 1
#> 
#> $H
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
#>  [1,]    0    0    0    0    0    0    0    0
#>  [2,]    0    0    0    0    0    0    0    0
#>  [3,]    1    0    0    0    0    0    0    0
#>  [4,]    0    1    0    0    0    0    0    0
#>  [5,]    0    0    0    0    0    0    0    0
#>  [6,]    0    0    0    0    0    0    0    0
#>  [7,]    0    0    1    0    0    0    0    0
#>  [8,]    0    0    0    1    0    0    0    0
#>  [9,]    0    0    0    0    0    0    0    0
#> [10,]    0    0    0    0    0    0    0    0
#> [11,]    0    0    0    0    1    0    0    0
#> [12,]    0    0    0    0    0    1    0    0
#> [13,]    0    0    0    0    0    0    0    0
#> [14,]    0    0    0    0    0    0    0    0
#> [15,]    0    0    0    0    0    0    1    0
#> [16,]    0    0    0    0    0    0    0    1
#> [17,]    0    0    0    0    0    0    0    0
#> [18,]    0    0    0    0    0    0    0    0
#> [19,]    0    0    0    0    0    0    0    0
#> [20,]    0    0    0    0    0    0    0    0
#> 
#> $class
#> [1] "rmfdmod"
#> 
#> $order
#> [1] 2 2 1 1
#> 
#> $n.par
#> [1] 8
#> 

# Use the template with fill_template()
# filled <- fill_template(tmpl, theta = rnorm(tmpl$n.par))
# Create a template
tmpl <- tmpl_stsp_full()
#> Error in tmpl_stsp_full(): argument "m" is missing, with no default
tmpl
#> $h
#>  [1] 1 0 0 0 1 0 0 0 0 1 0 0 0 1 0 0 1 0 0 1
#> 
#> $H
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
#>  [1,]    0    0    0    0    0    0    0    0
#>  [2,]    0    0    0    0    0    0    0    0
#>  [3,]    1    0    0    0    0    0    0    0
#>  [4,]    0    1    0    0    0    0    0    0
#>  [5,]    0    0    0    0    0    0    0    0
#>  [6,]    0    0    0    0    0    0    0    0
#>  [7,]    0    0    1    0    0    0    0    0
#>  [8,]    0    0    0    1    0    0    0    0
#>  [9,]    0    0    0    0    0    0    0    0
#> [10,]    0    0    0    0    0    0    0    0
#> [11,]    0    0    0    0    1    0    0    0
#> [12,]    0    0    0    0    0    1    0    0
#> [13,]    0    0    0    0    0    0    0    0
#> [14,]    0    0    0    0    0    0    0    0
#> [15,]    0    0    0    0    0    0    1    0
#> [16,]    0    0    0    0    0    0    0    1
#> [17,]    0    0    0    0    0    0    0    0
#> [18,]    0    0    0    0    0    0    0    0
#> [19,]    0    0    0    0    0    0    0    0
#> [20,]    0    0    0    0    0    0    0    0
#> 
#> $class
#> [1] "rmfdmod"
#> 
#> $order
#> [1] 2 2 1 1
#> 
#> $n.par
#> [1] 8
#> 

# Use the template with fill_template()
# filled <- fill_template(tmpl, theta = rnorm(tmpl$n.par))
# Create a template
tmpl <- tmpl_stsp_ar()
#> Error in tmpl_stsp_ar(): argument "m" is missing, with no default
tmpl
#> $h
#>  [1] 1 0 0 0 1 0 0 0 0 1 0 0 0 1 0 0 1 0 0 1
#> 
#> $H
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
#>  [1,]    0    0    0    0    0    0    0    0
#>  [2,]    0    0    0    0    0    0    0    0
#>  [3,]    1    0    0    0    0    0    0    0
#>  [4,]    0    1    0    0    0    0    0    0
#>  [5,]    0    0    0    0    0    0    0    0
#>  [6,]    0    0    0    0    0    0    0    0
#>  [7,]    0    0    1    0    0    0    0    0
#>  [8,]    0    0    0    1    0    0    0    0
#>  [9,]    0    0    0    0    0    0    0    0
#> [10,]    0    0    0    0    0    0    0    0
#> [11,]    0    0    0    0    1    0    0    0
#> [12,]    0    0    0    0    0    1    0    0
#> [13,]    0    0    0    0    0    0    0    0
#> [14,]    0    0    0    0    0    0    0    0
#> [15,]    0    0    0    0    0    0    1    0
#> [16,]    0    0    0    0    0    0    0    1
#> [17,]    0    0    0    0    0    0    0    0
#> [18,]    0    0    0    0    0    0    0    0
#> [19,]    0    0    0    0    0    0    0    0
#> [20,]    0    0    0    0    0    0    0    0
#> 
#> $class
#> [1] "rmfdmod"
#> 
#> $order
#> [1] 2 2 1 1
#> 
#> $n.par
#> [1] 8
#> 

# Use the template with fill_template()
# filled <- fill_template(tmpl, theta = rnorm(tmpl$n.par))

# ######################################################
# state space models in echelon form
# ######################################################
nu = c(3,2,4)   # Kronecker indices
m = length(nu)  # number of outputs/inputs
tmpl = tmpl_stsp_echelon(nu = nu)

# generate a random vector of parameters.
# Note that "tmpl$n.par" contains the number free parameters.
th = rnorm(tmpl$n.par)

# generate a model according to this structure with the parameters th
model = fill_template(th, tmpl)
print(model)
#> state space model [3,3] with s = 9 states
#>             s[1]        s[2]       s[3]        s[4]     s[5]       s[6]
#> s[1]  0.00000000  0.00000000  0.0000000  1.00000000 0.000000  0.0000000
#> s[2]  0.00000000  0.00000000  0.0000000  0.00000000 1.000000  0.0000000
#> s[3]  0.00000000  0.00000000  0.0000000  0.00000000 0.000000  1.0000000
#> s[4]  0.00000000  0.00000000  0.0000000  0.00000000 0.000000  0.0000000
#> s[5] -0.56187636  1.59850877  0.6307541 -0.52111732 1.300199 -0.1331510
#> s[6]  0.00000000  0.00000000  0.0000000  0.00000000 0.000000  0.0000000
#> s[7] -0.34391723 -0.08856511 -0.1136399 -0.48987045 2.293079 -1.7565274
#> s[8]  0.00000000  0.00000000  0.0000000  0.00000000 0.000000  0.0000000
#> s[9]  0.09049665  1.08079950 -1.5329020  0.04715443 1.547581 -0.3887799
#> x[1]  1.00000000  0.00000000  0.0000000  0.00000000 0.000000  0.0000000
#> x[2]  0.00000000  1.00000000  0.0000000  0.00000000 0.000000  0.0000000
#> x[3]  0.00000000  0.00000000  1.0000000  0.00000000 0.000000  0.0000000
#>            s[7]       s[8]     s[9]        u[1]       u[2]        u[3]
#> s[1] 0.00000000  0.0000000 0.000000 -0.44655722  0.1181445 -0.57746800
#> s[2] 0.00000000  0.0000000 0.000000  0.17480270  0.1340386  2.00248273
#> s[3] 0.00000000  0.0000000 0.000000  0.07455118  0.2210195  0.06670087
#> s[4] 1.00000000  0.0000000 0.000000  0.42816676  1.6408462  1.86685184
#> s[5] 0.08920722  0.0000000 0.000000  0.02467498 -0.2190504 -1.35090269
#> s[6] 0.00000000  1.0000000 0.000000 -1.66747510  0.1680654  0.02098359
#> s[7] 0.84501300  0.6843094 0.000000  0.73649596  1.1683839  1.24991457
#> s[8] 0.00000000  0.0000000 1.000000  0.38602657  1.0541810 -0.71524219
#> s[9] 0.96252797 -1.3952743 0.849643 -0.26565163  1.1452631 -0.75268897
#> x[1] 0.00000000  0.0000000 0.000000  1.00000000  0.0000000  0.00000000
#> x[2] 0.00000000  0.0000000 0.000000  0.00000000  1.0000000  0.00000000
#> x[3] 0.00000000  0.0000000 0.000000  0.00000000  0.0000000  1.00000000
#> Left square root of noise covariance Sigma:
#>            u[1]       u[2]      u[3]
#> u[1] -0.9385387  0.0000000 0.0000000
#> u[2] -1.0525133  0.3311792 0.0000000
#> u[3] -0.4371595 -2.0142105 0.2119804

# we can extract the free parameters from this given model
all.equal(th, extract_theta(model, tmpl, on_error = 'stop'))
#> [1] TRUE

# check the impulse response
k = impresp(model, lag.max = 2*sum(nu) + 1)

# the lag zero coeffcient k[0] is equal to the identity
all.equal(unclass(k$irf)[,,1], diag(m))
#> [1] TRUE

# check the Kronecker indices
all.equal(rationalmatrices::pseries2nu(k$irf), nu)
#> [1] TRUE
```

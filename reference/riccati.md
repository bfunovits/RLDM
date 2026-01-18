# Solve a discrete time, algebraic Riccati equation

This function solves the discrete time, algebraic Riccati equation \$\$X
= AXA' + (M -AXC')(G-CXC')^{-1}(M-AXC')'\$\$ where A is a square
(s-by-s) matrix, M and C' are matrices of dimension (s-by-m) and G is a
square (positive definite) matrix of dimension (m-by-m). Given certain
regularity conditions (see the discussion below) the solution X computed
by `riccati` is a positive definite matrix (of size (s-by-s)).

## Usage

``` r
riccati(A, M, C, G, only.X = TRUE)
```

## Arguments

- A:

  (s-by-s) matrix

- M:

  (s-by-m) matrix

- C:

  (m-by-s) matrix

- G:

  (m-by-m) matrix

- only.X:

  boolean

## Value

if (only.X) is TRUE then `riccati` just returns the solution X,
otherwise a list with slots

- X:

  the solution of the Ricatti equation

- B:

  the matrix \\B = (M- AXC')\Sigma^{-1}\\

- sigma:

  the matrix \\\Sigma=G-CXC'\\

- lambda:

  the (2m) vector of eigenvalues of the associated generalized
  eigenvalue problem. The first \\m\\ entries are the eigenvalues of
  \\A-BC\\.

## Details

Here (within this package) this function is mainly used to construct a
state space model if we have given the autocovariance function of a
process \\(y_t)\\ with a rational spectral density. The ACF is
represented by four matrices (A,M,C,G) as \\\gamma(0)=G\\ and
\\\gamma(k) = CA^{k-1}M\\ for \\k\>0\\. This *realization problem* is
related to the so called *spectral factorization problem*. If the matrix
\\A\\ is stable, the pair \\(A,C)\\ is controllable, the pair \\(A,M)\\
is controllable (i.e. (A,M,C,G) is a "minimal" realization of the ACF)
and if the spectral density is positive definite (i.e. has no zeros)
then the Riccati equation has a (unique) solution \\X\\ which is
positive definite and where the matrix \\A - (M-AXC')(G-CXC')^{-1}C\\ is
stable. The process \\(y_t)\\ then has a state space representation with
parameters \\(A,B,C,D=I)\\, where \\B=(M-AXC')\Sigma^{-1}\\ and
\\\Sigma=(G-CXC')\\ is the innovation covariance.

The solution is computed via an (2m-by-2m) dimensional generalized
eigenvalue problem, which in turn is solved with a QZ decomposition. See
[`QZ::qz.dgges()`](https://rdrr.io/pkg/QZ/man/bb_qz_dgges.html) and
[`QZ::qz.dtgsen()`](https://rdrr.io/pkg/QZ/man/bc_qz_dtgsen.html). The
eigenvalues of modulus less than one are the eigenvalues of the matrix
\\A-BC\\.

In addition `riccati` is also used to compute the stochastically
balanced realization of a state space model, see
`rationalmatrices::grammians()` and `rationalmatrices::balance()`.

Note that this function is mainly used as a utility function and
therefore no checks on the given input parameters are performed.

## Examples

``` r
# create a "random" state space model, which satisfies the
# stability and the (strict) miniphase assumption
m = 2 # number of outputs
s = 4 # number of states

model = r_model(template = tmpl_stsp_full(m, m, s),
                bpoles = 1, bzeroes = 1, sd = 0.25)
# scale sigma
model$sigma_L = model$sigma_L / sqrt(sum(diag(model$sigma_L %*% t(model$sigma_L))))

# extract the model parameter matrices
A = model$sys$A
B = model$sys$B
C = model$sys$C
sigma = model$sigma_L %*% t(model$sigma_L)

# compute the variance of the state P = A P A' + B sigma B'
P = lyapunov(A, B %*% sigma %*% t(B))

# variance of the output y[t]: G = C P C' + sigma
G = C %*% P %*% t(C) + sigma
# covariance between s[t+1] and y[t]: M = A P C' + B sigma
M = A %*% P %*% t(C) + B %*% sigma

# check that P solves the Riccati equation P = APA' + (M -APC')(G-CPC')^{-1}(M-APC')'
all.equal(P,
          A %*% P %*% t(A) +
            (M - A %*% P %*% t(C)) %*% solve(G - C%*% P %*% t(C), t(M - A %*% P %*% t(C))))
#> [1] TRUE

# compute P from the Riccati equation: P = APA' + (M -APC')(G-CPC')^{-1}(M-APC')'
out = riccati(A, M, C, G, only.X=FALSE)

# check the solution
all.equal(P, out$X)
#> [1] TRUE
all.equal(B, out$B)
#> [1] TRUE
all.equal(sigma, out$sigma)
#> [1] TRUE

# eigenvalues of (A-BC) ( <=> reciprocals of the zeroes of the system)
lambda = eigen(A - B %*% C, only.values=TRUE)$values
all.equal(sort(lambda), sort(out$lambda[1:s]))
#> [1] TRUE
```

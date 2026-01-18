# sigma_L Structure

Create templates for the left square root \\L\\ of the noise covariance
matrix \\\Sigma = LL'\\. This means that \\L\\ is parametrized as
\$\$\mbox{vec}(L) = h + H \theta\$\$ with a (\\k\\-dimensional)
parameter vector \\\theta\\.

## Usage

``` r
tmpl_sigma_L(
  sigma_L,
  structure = c("as_given", "chol", "symm", "identity", "full_normalized")
)
```

## Arguments

- sigma_L:

  numeric (n x n) matrix, where the free entries are coded with NAs

- structure:

  character string, determines the "structure" of sigma_L, see the
  examples.

## Value

List with slots

- `h` (\\n^2\\-dimensional vector),

- `H` (\\(n^2, k)\\-dimensional matrix, where \\k\\ denotes the number
  of free/deep parameters) and

- `n.par` (integer) is the number of free/deep parameters (\\=k\\).

## Details

The parameter `structure` has the following meaning

- as_given:

  Use the given parameter `sigma_L` to construct the template: `NA`
  entries are considered as free and all other entries as fixed.

- chol:

  Set all entries of `sigma_L` above the diagonal to zero and then
  proceed as above.

- symm:

  First make `sigma_L` symmetric (`sigma_L = (sigma_L + t(sigma_L))/2`)
  and then use this `sigma_L` as template. However, `h`, `H` are
  constructed such that \\h + H\theta\\ gives a symmetric matrix! Note
  that NAs overwrite fixed values, see examples.

- identity:

  Use the identity matrix as template. In this case there are no free
  parameters, i.e. \\\theta\\ is an empty vector (vector with zero
  length).

- full_normalized:

  Ones on the diagonal, otherwise all parameters are free.

## Examples

``` r
sigma_L = matrix(c(0, NA, 1, 0, 2, 3, NA, 1, 1), nrow = 3, ncol = 3)
sigma_L
#>      [,1] [,2] [,3]
#> [1,]    0    0   NA
#> [2,]   NA    2    1
#> [3,]    1    3    1

tmpl = tmpl_sigma_L(sigma_L, structure = 'as_given')
th = -(1:tmpl$n.par)
matrix(tmpl$h + tmpl$H %*% th, nrow = 3, ncol = 3)
#>      [,1] [,2] [,3]
#> [1,]    0    0   -2
#> [2,]   -1    2    1
#> [3,]    1    3    1

tmpl = tmpl_sigma_L(sigma_L, structure = 'chol')
th = -(1:tmpl$n.par)
matrix(tmpl$h + tmpl$H %*% th, nrow = 3, ncol = 3)
#>      [,1] [,2] [,3]
#> [1,]    0    0    0
#> [2,]   -1    2    0
#> [3,]    1    3    1

tmpl = tmpl_sigma_L(sigma_L, structure = 'symm')
th = -(1:tmpl$n.par)
matrix(tmpl$h + tmpl$H %*% th, nrow = 3, ncol = 3)
#>      [,1] [,2] [,3]
#> [1,]    0   -1   -2
#> [2,]   -1    2    2
#> [3,]   -2    2    1

tmpl = tmpl_sigma_L(sigma_L, structure = 'identity')
tmpl$n.par # = 0
#> [1] 0
matrix(tmpl$h, nrow = 3, ncol = 3)
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    1    0
#> [3,]    0    0    1

tmpl = tmpl_sigma_L(sigma_L, structure = 'full_normalized')
th = -(1:tmpl$n.par)
matrix(tmpl$h + tmpl$H %*% th, nrow = 3, ncol = 3)
#>      [,1] [,2] [,3]
#> [1,]    1   -3   -5
#> [2,]   -1    1   -6
#> [3,]   -2   -4    1
```

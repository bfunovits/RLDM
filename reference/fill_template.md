# Connect Deep Parameters with a Model

`fill_template` fills a given `template`, i.e. in essence an affine
mapping from the free parameters to the linear parameters of a model
class specified in the template, with the given free (deep) parameters
`th` and returns a model (e.g. an
[`armamod()`](https://bfunovits.github.io/RLDM/reference/armamod.md) or
a [`stspmod()`](https://bfunovits.github.io/RLDM/reference/stspmod.md)).
The procedure can be used to generate random models, see e.g.
[`r_model()`](https://bfunovits.github.io/RLDM/reference/r_model.md), or
for *M-estimates*, i.e. for estimation procedures where the estimate is
obtained by minimizing (maximizing) a criterion (e.g. ML or GMM
estimation).  
  
The "inverse function" `extract_theta` extracts the free/deep parameters
from the linear parameters of a model, by using the information provided
in the template. To this end the procedure first constructs the vector
\\\pi\\ of the stacked (linear) model parameters and then determines the
deep parameters \\\theta\\ as the least squares solution of the equation
system \$\$(\pi - h) = H\theta\$\$ If the residuals are non zero, then
the model does not (exactly) fit to the model structure. The threshold
`tol` is used to decide whether the model fits to the template or not.
The parameter `on_error` decides what to do in the case of a
"significant" misfit. For `"ignore"` the procedure ignores the misfit,
for `"warn"` the procedure issues a warning, and for `"stop"` the
procedure stops with an appropriate error message.  
In many cases the noise covariance is not explicitly parametrized, since
\\\Sigma\\ is estimated by another route. This may be accomplished by
fixing sigma_L to the identity matrix, with the option
`sigma_L = "identity"` for the `tmpl_***` functions. In order to
**extract the system parameters** (e.g. the AR and MA parameters for an
ARMA model) from a model where `sigma_L` is not equal to the identity,
one may use the option `ignore_sigma_L = TRUE`. This ignores a possible
mismatch for `sigma_L` but still checks whether the system parameters
are in accordance to the model structure.

## Usage

``` r
fill_template(th, template)

extract_theta(
  model,
  template,
  tol = sqrt(.Machine$double.eps),
  on_error = c("ignore", "warn", "stop"),
  ignore_sigma_L = FALSE
)
```

## Arguments

- th:

  Vector containing free (deep) parameters.

- template:

  A template like listed in
  [`model structures()`](https://bfunovits.github.io/RLDM/reference/model_structures.md),
  or a template explicitly specified by the user with
  [`model2template()`](https://bfunovits.github.io/RLDM/reference/model_structures.md).
  Essentially, a template is an affine mapping parametrised as a vector
  `h` and matrix `H` which connect the free (deep) parameters to the
  linear parameters of the model.

- model:

  A model object, i.e. an
  [`armamod()`](https://bfunovits.github.io/RLDM/reference/armamod.md),
  [`stspmod()`](https://bfunovits.github.io/RLDM/reference/stspmod.md),
  or
  [`rmfdmod()`](https://bfunovits.github.io/RLDM/reference/rmfdmod.md)
  object, from which deep parameters should be extracted.

- tol:

  In `extract_theta`, small double specifying the tolerance for an
  acceptable distance between the linear parameters and `H` times the
  deep parameters. Default is st to `sqrt(.Machine$double.eps)`.

- on_error:

  In `extract_theta`, character string with possible choices `ignore`,
  `warn`, `stop`. Specifies what should happen when the distance between
  the linear parmameters and `H` times the deep parameters, is larger
  than the specified `tol`. Default is `ignore`

- ignore_sigma_L:

  Boolean, default set to `FALSE`. If TRUE, the linear and free
  parameters pertaining the left square root `sigma_L` of the error
  covariance matrix are ignored. See also
  [`tmpl_sigma_L()`](https://bfunovits.github.io/RLDM/reference/tmpl_sigma_L.md)
  and
  [`model structures()`](https://bfunovits.github.io/RLDM/reference/model_structures.md)
  for more detail.

## Value

`fill_template` returns a model object, i.e. an
[`armamod()`](https://bfunovits.github.io/RLDM/reference/armamod.md),
[`stspmod()`](https://bfunovits.github.io/RLDM/reference/stspmod.md), or
[`rmfdmod()`](https://bfunovits.github.io/RLDM/reference/rmfdmod.md)
object, according to the class of the template (and the given parameters
`th`). The function `extract_theta` returns the vector of *free*
parameters for a given model and template.

## Connection to Likelihood Estimation

These functions are important for likelihood estimation where the
following instances of functionality are necessary.

- When an initial estimate is given through a model (together with a
  template), one may use `extract_theta` to extract the deep parameters.
  This vector of initial free/deep parameter values needs to be supplied
  to an optimizer.

- The optimized deep parameter values need to be filled into the model
  by using the structure provided by a template. This is done with
  `fill_template`

## See also

[`model structures()`](https://bfunovits.github.io/RLDM/reference/model_structures.md),
[`ll()`](https://bfunovits.github.io/RLDM/reference/ll.md),
[`ll_theta()`](https://bfunovits.github.io/RLDM/reference/ll_theta.md)
and [`ll_FUN()`](https://bfunovits.github.io/RLDM/reference/ll_FUN.md).

## Examples

``` r
# Extract deep parameter from ARMA object with ARMA(p,q) template ##########
(armamod_obj = test_armamod(dim = c(2,2), degree = c(3,1)))
#> ARMA model [2,2] with orders p = 3 and q = 1
#> AR polynomial a(z):
#>      z^0 [,1]  [,2]  z^1 [,1]      [,2]    z^2 [,1]       [,2]  z^3 [,1]
#> [1,]        1     0 0.4788278 0.8033552 -0.02115351 -0.3974357 1.5255242
#> [2,]        0     1 0.7681575 0.6468269  1.35827126  1.2689082 0.4007101
#>           [,2]
#> [1,]  1.021274
#> [2,] -1.538554
#> MA polynomial b(z):
#>      z^0 [,1]  [,2]   z^1 [,1]       [,2]
#> [1,]        1     0  1.7449898  0.8355711
#> [2,]        0     1 -0.4705249 -0.8206135
#> Left square root of noise covariance Sigma:
#>            u[1]      u[2]
#> u[1]  0.3213127 0.0000000
#> u[2] -0.8014435 0.6287823
(tmpl_obj = tmpl_arma_pq(2, 2, 3, 1))
#> $h
#>  [1] 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0
#> 
#> $H
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
#>  [1,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#>  [2,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#>  [3,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#>  [4,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#>  [5,]    1    0    0    0    0    0    0    0    0     0     0     0     0
#>  [6,]    0    1    0    0    0    0    0    0    0     0     0     0     0
#>  [7,]    0    0    1    0    0    0    0    0    0     0     0     0     0
#>  [8,]    0    0    0    1    0    0    0    0    0     0     0     0     0
#>  [9,]    0    0    0    0    1    0    0    0    0     0     0     0     0
#> [10,]    0    0    0    0    0    1    0    0    0     0     0     0     0
#> [11,]    0    0    0    0    0    0    1    0    0     0     0     0     0
#> [12,]    0    0    0    0    0    0    0    1    0     0     0     0     0
#> [13,]    0    0    0    0    0    0    0    0    1     0     0     0     0
#> [14,]    0    0    0    0    0    0    0    0    0     1     0     0     0
#> [15,]    0    0    0    0    0    0    0    0    0     0     1     0     0
#> [16,]    0    0    0    0    0    0    0    0    0     0     0     1     0
#> [17,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [18,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [19,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [20,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [21,]    0    0    0    0    0    0    0    0    0     0     0     0     1
#> [22,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [23,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [24,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [25,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [26,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [27,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [28,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#>       [,14] [,15] [,16] [,17] [,18] [,19]
#>  [1,]     0     0     0     0     0     0
#>  [2,]     0     0     0     0     0     0
#>  [3,]     0     0     0     0     0     0
#>  [4,]     0     0     0     0     0     0
#>  [5,]     0     0     0     0     0     0
#>  [6,]     0     0     0     0     0     0
#>  [7,]     0     0     0     0     0     0
#>  [8,]     0     0     0     0     0     0
#>  [9,]     0     0     0     0     0     0
#> [10,]     0     0     0     0     0     0
#> [11,]     0     0     0     0     0     0
#> [12,]     0     0     0     0     0     0
#> [13,]     0     0     0     0     0     0
#> [14,]     0     0     0     0     0     0
#> [15,]     0     0     0     0     0     0
#> [16,]     0     0     0     0     0     0
#> [17,]     0     0     0     0     0     0
#> [18,]     0     0     0     0     0     0
#> [19,]     0     0     0     0     0     0
#> [20,]     0     0     0     0     0     0
#> [21,]     0     0     0     0     0     0
#> [22,]     1     0     0     0     0     0
#> [23,]     0     1     0     0     0     0
#> [24,]     0     0     1     0     0     0
#> [25,]     0     0     0     1     0     0
#> [26,]     0     0     0     0     1     0
#> [27,]     0     0     0     0     0     0
#> [28,]     0     0     0     0     0     1
#> 
#> $class
#> [1] "armamod"
#> 
#> $order
#> [1] 2 2 3 1
#> 
#> $n.par
#> [1] 19
#> 

extract_theta(armamod_obj, tmpl_obj)
#>  [1]  0.47882776  0.76815751  0.80335524  0.64682691 -0.02115351  1.35827126
#>  [7] -0.39743566  1.26890825  1.52552422  0.40071009  1.02127370 -1.53855420
#> [13]  1.74498984 -0.47052488  0.83557110 -0.82061355  0.32131274 -0.80144349
#> [19]  0.62878228


# Fill template with deep parameters #################
(tmpl_obj = tmpl_arma_echelon(nu = c(3,2)))
#> $h
#>  [1] 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
#> 
#> $H
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
#>  [1,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#>  [2,]    1    0    0    0    0    0    0    0    0     0     0     0     0
#>  [3,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#>  [4,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#>  [5,]    0    1    0    0    0    0    0    0    0     0     0     0     0
#>  [6,]    0    0    1    0    0    0    0    0    0     0     0     0     0
#>  [7,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#>  [8,]    0    0    0    1    0    0    0    0    0     0     0     0     0
#>  [9,]    0    0    0    0    1    0    0    0    0     0     0     0     0
#> [10,]    0    0    0    0    0    1    0    0    0     0     0     0     0
#> [11,]    0    0    0    0    0    0    1    0    0     0     0     0     0
#> [12,]    0    0    0    0    0    0    0    1    0     0     0     0     0
#> [13,]    0    0    0    0    0    0    0    0    1     0     0     0     0
#> [14,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [15,]    0    0    0    0    0    0    0    0    0     1     0     0     0
#> [16,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [17,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [18,]    1    0    0    0    0    0    0    0    0     0     0     0     0
#> [19,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [20,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [21,]    0    0    0    0    0    0    0    0    0     0     1     0     0
#> [22,]    0    0    0    0    0    0    0    0    0     0     0     1     0
#> [23,]    0    0    0    0    0    0    0    0    0     0     0     0     1
#> [24,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [25,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [26,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [27,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [28,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [29,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [30,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [31,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [32,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [33,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [34,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [35,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [36,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#>       [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23]
#>  [1,]     0     0     0     0     0     0     0     0     0     0
#>  [2,]     0     0     0     0     0     0     0     0     0     0
#>  [3,]     0     0     0     0     0     0     0     0     0     0
#>  [4,]     0     0     0     0     0     0     0     0     0     0
#>  [5,]     0     0     0     0     0     0     0     0     0     0
#>  [6,]     0     0     0     0     0     0     0     0     0     0
#>  [7,]     0     0     0     0     0     0     0     0     0     0
#>  [8,]     0     0     0     0     0     0     0     0     0     0
#>  [9,]     0     0     0     0     0     0     0     0     0     0
#> [10,]     0     0     0     0     0     0     0     0     0     0
#> [11,]     0     0     0     0     0     0     0     0     0     0
#> [12,]     0     0     0     0     0     0     0     0     0     0
#> [13,]     0     0     0     0     0     0     0     0     0     0
#> [14,]     0     0     0     0     0     0     0     0     0     0
#> [15,]     0     0     0     0     0     0     0     0     0     0
#> [16,]     0     0     0     0     0     0     0     0     0     0
#> [17,]     0     0     0     0     0     0     0     0     0     0
#> [18,]     0     0     0     0     0     0     0     0     0     0
#> [19,]     0     0     0     0     0     0     0     0     0     0
#> [20,]     0     0     0     0     0     0     0     0     0     0
#> [21,]     0     0     0     0     0     0     0     0     0     0
#> [22,]     0     0     0     0     0     0     0     0     0     0
#> [23,]     0     0     0     0     0     0     0     0     0     0
#> [24,]     1     0     0     0     0     0     0     0     0     0
#> [25,]     0     1     0     0     0     0     0     0     0     0
#> [26,]     0     0     1     0     0     0     0     0     0     0
#> [27,]     0     0     0     1     0     0     0     0     0     0
#> [28,]     0     0     0     0     1     0     0     0     0     0
#> [29,]     0     0     0     0     0     1     0     0     0     0
#> [30,]     0     0     0     0     0     0     0     0     0     0
#> [31,]     0     0     0     0     0     0     1     0     0     0
#> [32,]     0     0     0     0     0     0     0     0     0     0
#> [33,]     0     0     0     0     0     0     0     1     0     0
#> [34,]     0     0     0     0     0     0     0     0     1     0
#> [35,]     0     0     0     0     0     0     0     0     0     0
#> [36,]     0     0     0     0     0     0     0     0     0     1
#> 
#> $class
#> [1] "armamod"
#> 
#> $order
#> [1] 2 2 3 3
#> 
#> $n.par
#> [1] 23
#> 
#> $nu
#> [1] 3 2
#> 
# Number of columns of matrix H in affine mapping = number of free parameters
(n_par_deep = dim(tmpl_obj$H)[2])
#> [1] 23

fill_template(rnorm(n_par_deep), tmpl_obj)
#> ARMA model [2,2] with orders p = 3 and q = 3
#> AR polynomial a(z):
#>        z^0 [,1]  [,2]  z^1 [,1]       [,2] z^2 [,1]       [,2]   z^3 [,1]
#> [1,]  1.0000000     0 0.1586592  0.0000000 1.677256 -0.3179531 -0.3444367
#> [2,] -0.1285962     1 1.7350727 -0.2612972 1.044599 -0.2397774  0.0000000
#>            [,2]
#> [1,] -0.4538713
#> [2,]  0.0000000
#> MA polynomial b(z):
#>        z^0 [,1]  [,2]   z^1 [,1]      [,2]    z^2 [,1]       [,2]  z^3 [,1]
#> [1,]  1.0000000     0  0.9679175 0.1632936 0.433025297 0.07230028 -1.012777
#> [2,] -0.1285962     1 -0.5661814 0.1644222 0.007921939 0.79980420  0.000000
#>           [,2]
#> [1,] 0.3808726
#> [2,] 0.0000000
#> Left square root of noise covariance Sigma:
#>            u[1]       u[2]
#> u[1] -0.4935700  0.0000000
#> u[2]  0.6620541 -0.3514827
```

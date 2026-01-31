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
#>      z^0 [,1]  [,2]   z^1 [,1]       [,2]  z^2 [,1]      [,2]   z^3 [,1]
#> [1,]        1     0  0.4609162 -0.6868529 1.2240818 0.4007715 -0.5558411
#> [2,]        0     1 -1.2650612 -0.4456620 0.3598138 0.1106827  1.7869131
#>            [,2]
#> [1,]  0.4978505
#> [2,] -1.9666172
#> MA polynomial b(z):
#>      z^0 [,1]  [,2]   z^1 [,1]       [,2]
#> [1,]        1     0  0.7013559 -1.0678237
#> [2,]        0     1 -0.4727914 -0.2179749
#> Left square root of noise covariance Sigma:
#>          u[1]     u[2]
#> u[1] 1.564061 0.000000
#> u[2] 0.212037 1.703367
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
#>  [1]  0.4609162 -1.2650612 -0.6868529 -0.4456620  1.2240818  0.3598138
#>  [7]  0.4007715  0.1106827 -0.5558411  1.7869131  0.4978505 -1.9666172
#> [13]  0.7013559 -0.4727914 -1.0678237 -0.2179749  1.5640610  0.2120370
#> [19]  1.7033671


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
#>       z^0 [,1]  [,2]   z^1 [,1]      [,2]  z^2 [,1]      [,2]  z^3 [,1]
#> [1,]  1.000000     0 -0.7288912  0.000000 0.8377870 -1.138137 0.4264642
#> [2,] -1.026004     1 -0.6250393 -1.686693 0.1533731  1.253815 0.0000000
#>            [,2]
#> [1,] -0.2950715
#> [2,]  0.0000000
#> MA polynomial b(z):
#>       z^0 [,1]  [,2]  z^1 [,1]      [,2]    z^2 [,1]       [,2]  z^3 [,1]
#> [1,]  1.000000     0 0.8951257 0.8215811  0.55391765 -0.3059627 -0.694707
#> [2,] -1.026004     1 0.8781335 0.6886403 -0.06191171 -0.3804710  0.000000
#>            [,2]
#> [1,] -0.2079173
#> [2,]  0.0000000
#> Left square root of noise covariance Sigma:
#>           u[1]     u[2]
#> u[1] -1.265396 0.000000
#> u[2]  2.168956 1.207962
# Basic example
# Create a template
tmpl <- tmpl_stsp_ar(m = 2, p = 1)
# Generate a random model with this structure
th0 <- rnorm(tmpl$n.par, sd = 0.1)
model <- fill_template(th0, tmpl)
# Extract the "free" parameters from the model
result <- extract_theta(model, tmpl)
result
#> [1] -0.112310858 -0.040288484 -0.046665535  0.077996512 -0.008336907
#> [6]  0.025331851 -0.002854676
```

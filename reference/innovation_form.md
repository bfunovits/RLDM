# Innovation Form state space Model

**\[experimental\]** Convert a given state space model into innovation
form, i.e. the transformed model satisfies

- \\D=I\\

- the model is stable and minimum phase.

The procedure "flips" bad poles and zeroes by the helper functions
[`rationalmatrices::reflect_zeroes()`](https://bfunovits.github.io/rationalmatrices/reference/reflect_zeroes.html)
and
[`rationalmatrices::reflect_poles()`](https://bfunovits.github.io/rationalmatrices/reference/reflect_poles.html).
The transformed model is an equivalent description of the process in
terms of second order moments. This means that the spectral density is
not changed.

## Usage

``` r
innovation_form(model, echelon_form = TRUE, y = NULL)
```

## Arguments

- model:

  A state space model, i.e. an object of type
  [`stspmod()`](https://bfunovits.github.io/RLDM/reference/stspmod.md).

- echelon_form:

  boolean, if `TRUE` the innovation form model is in addition
  transformed to echelon form.

- y:

  `NULL` or a data sample. `as.matrix(y)` should return an
  (N,m)-dimensional numeric matrix. If not `NULL` the noise covariance
  matrix is estimated from this sample.

## Value

List with slots

- model:

  the state space model in innovation form.

- z:

  (complex) vector with the zeroes of the innovation form model.

- z_flipped:

  (boolean) vector which indicates which zeroes have been flipped.

- p:

  (complex) vector with the poles of the innovation form model.

- p_flipped:

  (boolean) vector which indicates which poles have been flipped.

## Examples

``` r
# in order to get reproducable results
set.seed(342)

model = r_model(tmpl_stsp_full(3, 3, 5))
print(model, digits = 4)
#> state space model [3,3] with s = 5 states
#>         s[1]    s[2]    s[3]    s[4]    s[5]    u[1]    u[2]    u[3]
#> s[1] -0.7701 -0.7534  0.4933 -1.2629  0.3608  0.0149  0.8278 -1.1295
#> s[2] -0.9100  2.2148 -0.4970 -2.2051 -1.0721 -0.9960  0.5980 -1.5275
#> s[3]  1.2751  0.9052 -0.7475 -0.5919 -1.6086 -1.1893  0.0486 -1.6221
#> s[4] -0.1740 -1.3902  1.4165  0.8076  0.5328 -0.9149 -0.2580 -0.4574
#> s[5]  0.5116 -0.0054 -0.6111  0.5597 -0.0738  0.0839 -0.3527 -0.6538
#> x[1]  0.1336  1.7917  0.5724  0.5907  0.1006  1.0000  0.0000  0.0000
#> x[2]  1.1534  1.2193  1.4340 -0.4257  0.6177  0.0000  1.0000  0.0000
#> x[3] -1.2113 -0.5287 -1.2432 -0.8186  0.8428  0.0000  0.0000  1.0000
#> Left square root of noise covariance Sigma:
#>            u[1]      u[2]     u[3]
#> u[1]  0.6968024 0.0000000 0.000000
#> u[2] -0.3246241 0.1801707 0.000000
#> u[3]  0.9493687 0.6805352 1.089781
# the model has two non-minimum phase zeroes and two non-stable poles.
z = zeroes(model, print_message = FALSE)
abs(z)
#> [1] 0.2430920 0.3822571 1.4472231 1.4472231 1.7947154
p = poles(model, print_message = FALSE)
abs(p)
#> [1] 0.3184163 0.3801766 1.1047541 1.1047541 7.2894009

# convert to innnovation form, by flipping the "bad" poles and zeroes.
out = innovation_form(model)
print(out$model, digits = 4)
#> state space model [3,3] with s = 5 states
#>         s[1]    s[2]    s[3]    s[4]    s[5]     u[1]     u[2]    u[3]
#> s[1]  0.0000  0.0000  0.0000  1.0000  0.0000  -2.1013   0.4187 -3.8366
#> s[2]  0.0000  0.0000  0.0000  0.0000  1.0000  -4.0003  -1.8191 -5.5912
#> s[3]  0.8764 -0.0799  1.0267  0.0449 -0.4833   6.2841   4.3642  4.3778
#> s[4] -1.0675 -0.1793 -2.2047 -2.5517  1.3806 -16.2235 -18.8212 -3.6213
#> s[5] -0.7212 -0.4139 -2.5557 -4.6404  2.3841 -17.7956 -23.0643 -1.8357
#> x[1]  1.0000  0.0000  0.0000  0.0000  0.0000   1.0000   0.0000  0.0000
#> x[2]  0.0000  1.0000  0.0000  0.0000  0.0000   0.0000   1.0000  0.0000
#> x[3]  0.0000  0.0000  1.0000  0.0000  0.0000   0.0000   0.0000  1.0000
#> Left square root of noise covariance Sigma:
#>             u[1]       u[2]     u[3]
#> u[1]  0.52953545  0.0000000 0.000000
#> u[2] -0.31786378  0.3060925 0.000000
#> u[3]  0.06273688 -1.0192496 1.099613
flip = function(x) {
  x[abs(x) < 1] = 1/x[abs(x) < 1]
  return(x)}
data.frame(poles.inno = out$p, flipped = out$p_flipped,
           poles.ori = p[match_vectors(out$p, p)],
           zeroes.inno = out$z, flipped = out$z_flipped,
           zeroes.ori = z[match_vectors(out$z, flip(z))])
#>              poles.inno flipped             poles.ori          zeroes.inno
#> 1  0.6456403-0.8964543i   FALSE  0.6456403-0.8964543i -1.7947154+0.000000i
#> 2  0.6456403+0.8964543i   FALSE  0.6456403+0.8964543i  0.8583822-1.165176i
#> 3 -7.2894009+0.0000000i   FALSE -7.2894009+0.0000000i  0.8583822+1.165176i
#> 4  3.1405429+0.0000000i    TRUE  0.3184163+0.0000000i -4.1136687+0.000000i
#> 5 -2.6303566+0.0000000i    TRUE -0.3801766+0.0000000i  2.6160404+0.000000i
#>   flipped.1           zeroes.ori
#> 1     FALSE -1.7947154+0.000000i
#> 2     FALSE  0.8583822-1.165176i
#> 3     FALSE  0.8583822+1.165176i
#> 4      TRUE -0.2430920+0.000000i
#> 5      TRUE  0.3822571+0.000000i

# check that the innovation form model describes the same process,
# by checking that the spectral density is not changed!
junk = spectrald(model, n.f = 128)
junk1 = spectrald(out$model, n.f = 128)
all.equal(junk, junk1)
#> [1] TRUE

# reset seed
set.seed(NULL)
```

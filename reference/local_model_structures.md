# Local Model Structures

Parametrization for "local" model classes, in particular, "Data Driven
Local Coordinates" as detailed in (McKelvey et al. 2004) and (Ribarits
et al. 2005) .

## Usage

``` r
tmpl_DDLC(
  model,
  balance = c("none", "lyapunov", "minimum phase"),
  sigma_L = c("chol", "symm", "identity")
)

tmpl_GRAM(
  model,
  balance = c("lyapunov", "minimum phase"),
  sigma_L = c("chol", "symm", "identity")
)
```

## Arguments

- model:

  [`stspmod()`](https://bfunovits.github.io/RLDM/reference/stspmod.md)
  object, which represents a state space model. Only the case \\m = n \>
  0\\ is implemented, i.e. the output process and the noise process must
  be of the same dimension.

- balance:

  (character string) For `balance = "lyapunov"` or
  `balance = "minimum phase"` the reference model is first balanced by
  the respective scheme.

- sigma_L:

  (character string) determines the form of the (left) square root of
  the noise covariance \\\Sigma\\. The choice `"chol"` gives a lower
  triangular matrix, `"symm"` gives a symmetric matrix and `"identity"`
  corresponds to the (fixed) identity matrix.

## Value

Model template, i.e. a list with slots

- `h`:

  \\((m+s)^2 + m^2)\\-dimensional vector,

- `H`:

  \\((m+s)^2 + m^2, k)\\-dimensional matrix,

- `class = "stspmod"`:

- `order = c(m,m,s)`:

  and

- `n.par`:

  number of free parameters \\=k\\.

See also
[`model structures()`](https://bfunovits.github.io/RLDM/reference/model_structures.md)
for more details.

## Details

The function `tmpl_DDLC` and `tmpl_GRAM` construct model templates which
describe models in a neighborhood of a given reference model. In a first
step the reference state space model is transformed to \\D=I\\ and
eventually (depending on the parameter `"balance"`) balanced.

state space models are described by a quadruple \\(A,B,C,D=I)\\ of
matrices which may be embedded into an \\(s^2+2ms)\\-dimensional
euclidean space. Note that the parameter matrices are not uniqely
determined from the ACF or the spectral density of the process, i.e.
there is an inherent non identifiablity problem. For minimal models the
"equivalence class" of models, which represent the same ACF is given by
the set of all models which may be obtained by a state transformation
\\(A,B,C,D) \rightarrow (TAT^{-1}, TB, CT^{-1}, D)\\.

The *DDLC* parametrization now considers models, \\(A,B,C,D=I)\\, which
are contained in the \\2ms\\-dimensional subspace, which is *orthogonal*
to the \\s^2\\-dimensional tangent space of the set of equivalent
models.

The routine `tmpl_GRAM` considers the \\2ms\\-dimensional subspace,
where models close to the reference models are "approximately" balanced.

Both schemes may fail for "non-generic" models. `tmpl_DDLC` issues a
warning message and `tmpl_GRAM` throws an error, in cases where the
\\2ms\\-dimensional subspace is not well defined.

Note that also the parametrization of the left square root `L=sigma_L`
of the noise covariance is "local", i.e. `th = 0` corresponds to the
(balanced) reference model.

## References

McKelvey T, Helmersson A, Ribarits T (2004). “Data driven local
coordinates for multivariable linear systems and their application to
system identification.” *Automatica*, **40**(9), 1629 - 1635.
[doi:10.1016/j.automatica.2004.04.015](https://doi.org/10.1016/j.automatica.2004.04.015)
.

Ribarits T, Deistler M, Hanzon B (2005). “An analysis of separable least
squares data driven local coordinates for maximum likelihood estimation
of linear systems.” *Automatica*, **41**(3), 531 - 544.
[doi:10.1016/j.automatica.2004.11.014](https://doi.org/10.1016/j.automatica.2004.11.014)
. .

## See also

For the computation of Grammians and for balancing of state space
models, see `rationalmatrices::balance()`.

## Examples

``` r
# create a random state space model with m outputs and s states
m = 3
s = 6
tmpl = tmpl_stsp_full(m, n = m, s, sigma_L = 'symm')
model = r_model(tmpl, bpoles = 1.1, bzeroes = 1.1, sd = 1/s)
model                              # note that sigma_L is symmetric
#> state space model [3,3] with s = 6 states
#>             s[1]        s[2]          s[3]         s[4]         s[5]
#> s[1] -0.06837539 -0.07765215 -0.1173139852  0.009878355 -0.005667642
#> s[2]  0.24411134 -0.10049338 -0.0696589081 -0.026992122 -0.366504332
#> s[3]  0.01849703  0.12483852  0.2544243832 -0.039176510 -0.054199056
#> s[4]  0.37406183 -0.21210549 -0.4072748906 -0.175718715  0.475811559
#> s[5] -0.03918026  0.06986259 -0.0997740674  0.062024192 -0.038264570
#> s[6]  0.34300122  0.22948153 -0.0009856264 -0.284862022  0.009643621
#> x[1]  0.13952692 -0.06613967  0.1054918253 -0.059415696 -0.092003476
#> x[2] -0.01087614  0.18230402 -0.3476755891 -0.007543803  0.231109079
#> x[3] -0.42554590  0.20387854  0.1684433703 -0.192112043  0.422032216
#>              s[6]        u[1]       u[2]        u[3]
#> s[1]  0.062751892 -0.11234615 0.17717855 -0.25093557
#> s[2] -0.056104427  0.11735589 0.10563535  0.40515858
#> s[3] -0.021586973 -0.05159080 0.39826474 -0.09013330
#> s[4]  0.062155844 -0.02746224 0.03065239  0.09599769
#> s[5] -0.107442839 -0.45898617 0.03320163 -0.15208060
#> s[6] -0.006567358  0.14293516 0.03991862 -0.02978739
#> x[1] -0.078594212  1.00000000 0.00000000  0.00000000
#> x[2] -0.175841601  0.00000000 1.00000000  0.00000000
#> x[3]  0.144933442  0.00000000 0.00000000  1.00000000
#> Left square root of noise covariance Sigma:
#>             u[1]       u[2]        u[3]
#> u[1]  0.11847317 -0.2704352 -0.06970479
#> u[2] -0.27043524  0.1201868 -0.20886560
#> u[3] -0.06970479 -0.2088656  0.29766034
model$sigma_L %*% t(model$sigma_L) # noise covariance Sigma
#>             [,1]        [,2]        [,3]
#> [1,]  0.09202987 -0.04998312  0.02747812
#> [2,] -0.04998312  0.13120492 -0.06842326
#> [3,]  0.02747812 -0.06842326  0.13708528

# tmpl_DDLC #############################################

# create a DDLC parametrization of a neighborhood of this model
tmpl = tmpl_DDLC(model, balance = 'lyapunov', sigma_L = 'chol')
# for th = 0, we get the original model (in balanced form)
model = fill_template(numeric(tmpl$n.par), tmpl)
model                                # note that sigma_L is lower triangular
#> state space model [3,3] with s = 6 states
#>             s[1]        s[2]        s[3]        s[4]         s[5]         s[6]
#> s[1] -0.04674028 -0.17261259  0.39066853  0.05101132  0.159228023 -0.003990457
#> s[2] -0.27599581  0.39541020  0.56954044 -0.19514673 -0.034168460 -0.002471917
#> s[3] -0.07491692 -0.12036024 -0.14717153  0.03323060  0.118664347  0.062964364
#> s[4]  0.07354901  0.39209590 -0.16077216 -0.48253761 -0.140622491  0.160024761
#> s[5] -0.05249468  0.06616410  0.20245703  0.14703738 -0.273924131 -0.318682819
#> s[6] -0.02222974 -0.01472977  0.06252047 -0.10773188  0.415848099  0.419968329
#> x[1] -0.10951170  0.07166656 -0.14201865 -0.13093253  0.058457617 -0.029815084
#> x[2]  0.36346815 -0.16194885  0.05683984 -0.09965656  0.001219085 -0.007883164
#> x[3]  0.25347210  0.32414936 -0.01683194  0.05861268  0.028981624 -0.003096200
#>               u[1]         u[2]        u[3]
#> s[1] -0.3115666795 -0.227320179  0.22361976
#> s[2] -0.0737555209  0.279041552  0.10537148
#> s[3]  0.1934560634 -0.001507147  0.27719141
#> s[4]  0.0009333888 -0.071846781  0.02748191
#> s[5]  0.0502772031 -0.053436728 -0.01911524
#> s[6]  0.0144344216 -0.011584255 -0.01110399
#> x[1]  1.0000000000  0.000000000  0.00000000
#> x[2]  0.0000000000  1.000000000  0.00000000
#> x[3]  0.0000000000  0.000000000  1.00000000
#> Left square root of noise covariance Sigma:
#>             u[1]       u[2]      u[3]
#> u[1]  0.30336425  0.0000000 0.0000000
#> u[2] -0.16476273  0.3225805 0.0000000
#> u[3]  0.09057799 -0.1658482 0.3183949
model$sigma_L %*% t(model$sigma_L)  # however Sigma is the same as above
#>             [,1]        [,2]        [,3]
#> [1,]  0.09202987 -0.04998312  0.02747812
#> [2,] -0.04998312  0.13120492 -0.06842326
#> [3,]  0.02747812 -0.06842326  0.13708528

#' apply a "small" state transformation T = (diag(s)+eps*X)
eps = sqrt(.Machine$double.eps)
sys = model$sys
d_sys = state_trafo(sys, diag(s) + matrix(rnorm(s^2, sd = eps), nrow = s, ncol = s))
d_pi = (as.vector(unclass(d_sys) - unclass(sys)))/eps
# The vector d_pi is (close to) an element of the tangent space
# of the set of models, which are generated by a state transformation
# of the reference model

# by construction d_pi is (close to) orthogonal to tmpl$H
max(abs(d_pi %*% tmpl$H[1:((m+s)^2), , drop = FALSE]))
#> [1] 5.040806e-08

# the tmpl_DDLC routine may fail in some exceptional cases
m = 1
s = 3
model = stspmod(sys = stsp(A = matrix(0, nrow = s, ncol = s),
                           B = matrix(rnorm(m*s), nrow = s, ncol = m),
                           C = matrix(rnorm(m*s), nrow = m, ncol = s),
                           D = diag(m)),
                sigma_L = diag(m))

# For this model "tmpl_DLLC" issues a warning.
junk = tmpl_DDLC(model, sigma_L = 'chol', balance = 'none')
#> Warning: The tangent space of the equivalence class does not have dimension s^2=9 (sv[1]=2.37294538950255, sv[9]=1.31089551872483e-17)

# tmpl_GRAM #############################################
model = fill_template(numeric(tmpl$n.par), tmpl)

tmpl = tmpl_GRAM(model, sigma_L = 'chol')
model = fill_template(numeric(tmpl$n.par), tmpl)

# check grammians
gr = grammians(model$sys, 'lyapunov')
P = gr$P
Q = gr$Q
# P=Q=diag() should hold!
print(round(cbind(P, P-Q), 6))
#>          [,1]    [,2]     [,3]    [,4]     [,5]     [,6] [,7] [,8] [,9] [,10]
#> [1,] 0.223614 0.00000 0.000000 0.00000 0.000000 0.000000    0    0    0     0
#> [2,] 0.000000 0.18092 0.000000 0.00000 0.000000 0.000000    0    0    0     0
#> [3,] 0.000000 0.00000 0.121045 0.00000 0.000000 0.000000    0    0    0     0
#> [4,] 0.000000 0.00000 0.000000 0.05018 0.000000 0.000000    0    0    0     0
#> [5,] 0.000000 0.00000 0.000000 0.00000 0.014839 0.000000    0    0    0     0
#> [6,] 0.000000 0.00000 0.000000 0.00000 0.000000 0.005145    0    0    0     0
#>      [,11] [,12]
#> [1,]     0     0
#> [2,]     0     0
#> [3,]     0     0
#> [4,]     0     0
#> [5,]     0     0
#> [6,]     0     0

# now consider a model close to the reference model
d_th = rnorm(tmpl$n.par, sd = eps)
d_model = fill_template(d_th, tmpl)
d_sys = d_model$sys
gr = grammians(d_sys, 'lyapunov')
d_P = gr$P - P
d_Q = gr$Q - Q

# the "disturbed" system should still be approximately balanced!
print(round(cbind(d_P, d_P - d_Q)/eps, 6) )
#>          [,1]     [,2]     [,3]     [,4]    [,5]     [,6] [,7] [,8] [,9] [,10]
#> [1,] 0.311185 0.000000 0.000000 0.000000 0.00000 0.000000    0    0    0     0
#> [2,] 0.000000 0.374576 0.000000 0.000000 0.00000 0.000000    0    0    0     0
#> [3,] 0.000000 0.000000 0.457341 0.000000 0.00000 0.000000    0    0    0     0
#> [4,] 0.000000 0.000000 0.000000 0.465702 0.00000 0.000000    0    0    0     0
#> [5,] 0.000000 0.000000 0.000000 0.000000 0.05869 0.000000    0    0    0     0
#> [6,] 0.000000 0.000000 0.000000 0.000000 0.00000 0.036142    0    0    0     0
#>      [,11] [,12]
#> [1,]     0     0
#> [2,]     0     0
#> [3,]     0     0
#> [4,]     0     0
#> [5,]     0     0
#> [6,]     0     0
```

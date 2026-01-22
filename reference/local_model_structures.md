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
#>             s[1]        s[2]         s[3]          s[4]        s[5]        s[6]
#> s[1] -0.22265095  0.50272019  0.083513362  0.2571544759 -0.15906105 -0.27229490
#> s[2] -0.11927751 -0.16415220  0.001178752  0.0960543961  0.20591719 -0.01219423
#> s[3] -0.11021934  0.05313687  0.091211239 -0.0025940641 -0.03484008  0.12386439
#> s[4] -0.07729205  0.15424939 -0.006312946 -0.2471635357  0.09182998 -0.05297297
#> s[5]  0.18620731  0.54056846 -0.074042125  0.0207133786 -0.11031516 -0.12790172
#> s[6] -0.04158061  0.06799590  0.240351977  0.0582696231  0.09531564 -0.15093084
#> x[1]  0.35927919  0.29838600 -0.086336419  0.1441695219  0.28723100  0.09486512
#> x[2]  0.16698048  0.08323002 -0.063497533 -0.1934404090 -0.08559753  0.03810958
#> x[3] -0.10010093 -0.13964260  0.209183518 -0.0005764036  0.16099128 -0.17881458
#>             u[1]        u[2]        u[3]
#> s[1]  0.16462478  0.31937121  0.25967546
#> s[2]  0.15654489  0.18129501  0.04779198
#> s[3] -0.25792220  0.02887677  0.30754427
#> s[4] -0.04657263  0.29725757  0.12495017
#> s[5] -0.15681295  0.19306445 -0.27792576
#> s[6] -0.11656068 -0.13340313  0.04655536
#> x[1]  1.00000000  0.00000000  0.00000000
#> x[2]  0.00000000  1.00000000  0.00000000
#> x[3]  0.00000000  0.00000000  1.00000000
#> Left square root of noise covariance Sigma:
#>             u[1]        u[2]        u[3]
#> u[1]  0.14463601  0.01497182 -0.28498768
#> u[2]  0.01497182 -0.17696139  0.30163687
#> u[3] -0.28498768  0.30163687  0.08626554
model$sigma_L %*% t(model$sigma_L) # noise covariance Sigma
#>             [,1]        [,2]        [,3]
#> [1,]  0.10236171 -0.08644676 -0.06128804
#> [2,] -0.08644676  0.12252429 -0.03162400
#> [3,] -0.06128804 -0.03162400  0.17964452

# tmpl_DDLC #############################################

# create a DDLC parametrization of a neighborhood of this model
tmpl = tmpl_DDLC(model, balance = 'lyapunov', sigma_L = 'chol')
# for th = 0, we get the original model (in balanced form)
model = fill_template(numeric(tmpl$n.par), tmpl)
model                                # note that sigma_L is lower triangular
#> state space model [3,3] with s = 6 states
#>             s[1]        s[2]         s[3]          s[4]         s[5]
#> s[1]  0.22749707 -0.09706120 -0.376418415 -0.2081900642  0.054668903
#> s[2]  0.13130247 -0.36833737  0.011801937  0.0340802642 -0.067264809
#> s[3]  0.30960298 -0.25986571 -0.555903533  0.5208816630  0.007463684
#> s[4] -0.06039423 -0.07845019 -0.104719211  0.1102402592 -0.142683665
#> s[5] -0.06439903 -0.17966764  0.054101790  0.0178486443 -0.054539103
#> s[6] -0.02066912  0.02215459 -0.005762771  0.1407234833 -0.189952759
#> x[1] -0.51782150 -0.08054127 -0.067977961 -0.0008951252  0.006844787
#> x[2] -0.03901319  0.19629373 -0.125363793 -0.0367534161 -0.056625548
#> x[3]  0.04949881 -0.26791384  0.017674281 -0.0588825594 -0.034572369
#>             s[6]         u[1]          u[2]         u[3]
#> s[1]  0.01002193 -0.168232498 -0.4847356383 -0.063645698
#> s[2]  0.07421788  0.315846871 -0.0963310647  0.109632249
#> s[3] -0.02324560 -0.056813276  0.0993283768  0.010965345
#> s[4]  0.05292509  0.060633747  0.0030335678 -0.190016927
#> s[5] -0.10497537 -0.025337571 -0.0001285062  0.006827089
#> s[6] -0.16295878  0.001324394 -0.0071282729  0.005706124
#> x[1]  0.00125503  1.000000000  0.0000000000  0.000000000
#> x[2] -0.00384912  0.000000000  1.0000000000  0.000000000
#> x[3] -0.01732466  0.000000000  0.0000000000  1.000000000
#> Left square root of noise covariance Sigma:
#>            u[1]       u[2]       u[3]
#> u[1]  0.3199402  0.0000000 0.00000000
#> u[2] -0.2701967  0.2225265 0.00000000
#> u[3] -0.1915610 -0.3747109 0.05040472
model$sigma_L %*% t(model$sigma_L)  # however Sigma is the same as above
#>             [,1]        [,2]        [,3]
#> [1,]  0.10236171 -0.08644676 -0.06128804
#> [2,] -0.08644676  0.12252429 -0.03162400
#> [3,] -0.06128804 -0.03162400  0.17964452

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
#> [1] 3.344665e-08

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
#> Warning: The tangent space of the equivalence class does not have dimension s^2=9 (sv[1]=2.11291112328716, sv[9]=4.22480766864052e-18)

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
#>          [,1]     [,2]    [,3]     [,4]     [,5]     [,6] [,7] [,8] [,9] [,10]
#> [1,] 0.299115 0.000000 0.00000 0.000000 0.000000 0.000000    0    0    0     0
#> [2,] 0.000000 0.146146 0.00000 0.000000 0.000000 0.000000    0    0    0     0
#> [3,] 0.000000 0.000000 0.09197 0.000000 0.000000 0.000000    0    0    0     0
#> [4,] 0.000000 0.000000 0.00000 0.043465 0.000000 0.000000    0    0    0     0
#> [5,] 0.000000 0.000000 0.00000 0.000000 0.006966 0.000000    0    0    0     0
#> [6,] 0.000000 0.000000 0.00000 0.000000 0.000000 0.001438    0    0    0     0
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
#>           [,1]     [,2]      [,3]    [,4]      [,5]      [,6] [,7] [,8] [,9]
#> [1,] -0.174217 0.000000  0.000000  0.0000  0.000000  0.000000    0    0    0
#> [2,]  0.000000 0.062264  0.000000  0.0000  0.000000  0.000000    0    0    0
#> [3,]  0.000000 0.000000 -0.372991  0.0000  0.000000  0.000000    0    0    0
#> [4,]  0.000000 0.000000  0.000000 -0.3011  0.000000  0.000000    0    0    0
#> [5,]  0.000000 0.000000  0.000000  0.0000 -0.037533  0.000000    0    0    0
#> [6,]  0.000000 0.000000  0.000000  0.0000  0.000000 -0.008466    0    0    0
#>      [,10] [,11] [,12]
#> [1,]     0     0     0
#> [2,]     0     0     0
#> [3,]     0     0     0
#> [4,]     0     0     0
#> [5,]     0     0     0
#> [6,]     0     0     0
```

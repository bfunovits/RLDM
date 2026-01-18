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
#>             s[1]        s[2]         s[3]        s[4]         s[5]        s[6]
#> s[1]  0.01321050  0.31123733 -0.009441997 -0.09718823  0.098663861  0.15103045
#> s[2] -0.09212396 -0.24171562  0.201285623  0.06670158  0.461532132 -0.01896857
#> s[3] -0.03836733  0.03447621 -0.078772912 -0.39696322 -0.168941640 -0.15171288
#> s[4] -0.10120253  0.16971301 -0.119165128  0.15511276 -0.099681676 -0.13453779
#> s[5] -0.27816388 -0.01823292  0.170094401 -0.15004522  0.072039933 -0.10276843
#> s[6]  0.19514689  0.19281588  0.006588536  0.09156063  0.127785700  0.13922648
#> x[1]  0.05765478  0.08948316 -0.015088562 -0.04971894 -0.101702893  0.10983098
#> x[2] -0.06704662  0.45932249 -0.096246648  0.13594879 -0.019748603  0.10850459
#> x[3]  0.09701793  0.09426726 -0.188310146  0.17387984  0.003516328 -0.12204231
#>               u[1]        u[2]        u[3]
#> s[1]  0.2392809195 -0.14529553 -0.17589796
#> s[2]  0.0271601047  0.17081427  0.03018147
#> s[3]  0.1858441523 -0.06684824 -0.34365827
#> s[4] -0.1170055161 -0.29034523  0.16466124
#> s[5]  0.2913418368  0.23268445  0.15136692
#> s[6] -0.0007903322  0.03729747 -0.04250798
#> x[1]  1.0000000000  0.00000000  0.00000000
#> x[2]  0.0000000000  1.00000000  0.00000000
#> x[3]  0.0000000000  0.00000000  1.00000000
#> Left square root of noise covariance Sigma:
#>             u[1]       u[2]        u[3]
#> u[1] -0.35636364 0.02519447  0.12560560
#> u[2]  0.02519447 0.08273967  0.02693086
#> u[3]  0.12560560 0.02693086 -0.31980116
model$sigma_L %*% t(model$sigma_L) # noise covariance Sigma
#>              [,1]         [,2]         [,3]
#> [1,]  0.143406569 -0.003511144 -0.084251576
#> [2,] -0.003511144  0.008205885 -0.003219702
#> [3,] -0.084251576 -0.003219702  0.118774819

# tmpl_DDLC #############################################

# create a DDLC parametrization of a neighborhood of this model
tmpl = tmpl_DDLC(model, balance = 'lyapunov', sigma_L = 'chol')
# for th = 0, we get the original model (in balanced form)
model = fill_template(numeric(tmpl$n.par), tmpl)
model                                # note that sigma_L is lower triangular
#> state space model [3,3] with s = 6 states
#>             s[1]        s[2]        s[3]         s[4]         s[5]         s[6]
#> s[1] -0.01648914 -0.48227471  0.31563376 -0.065616269 -0.074004828 -0.008931280
#> s[2] -0.26148205 -0.22750392  0.38241885  0.331408414  0.056781725  0.023382218
#> s[3] -0.05050684  0.03964702  0.28414112 -0.269525477  0.084887279 -0.046339198
#> s[4]  0.01846800 -0.09502178  0.53136447  0.243060472  0.034518652 -0.040168945
#> s[5] -0.06189414 -0.03044651 -0.05610569  0.157839698 -0.019773594  0.292895352
#> s[6]  0.02372034  0.04524106  0.01123395  0.083903181 -0.045606496 -0.204333789
#> x[1] -0.02524009  0.12191751 -0.02148435  0.113429905 -0.015047579 -0.010977896
#> x[2]  0.30329879  0.11155597  0.11506967 -0.001310153 -0.002854065  0.003859815
#> x[3]  0.20770309 -0.18088368 -0.11623649  0.050838070  0.009918857 -0.004657147
#>              u[1]         u[2]         u[3]
#> s[1] -0.162079637  0.059570234  0.287526328
#> s[2] -0.129839902  0.152593509 -0.171949970
#> s[3]  0.178612895  0.168350234  0.047665038
#> s[4]  0.046270200 -0.094024289 -0.009699282
#> s[5]  0.012618269 -0.002026864  0.007599202
#> s[6] -0.003635076  0.005888528  0.004167705
#> x[1]  1.000000000  0.000000000  0.000000000
#> x[2]  0.000000000  1.000000000  0.000000000
#> x[3]  0.000000000  0.000000000  1.000000000
#> Left square root of noise covariance Sigma:
#>              u[1]        u[2]      u[3]
#> u[1]  0.378690598  0.00000000 0.0000000
#> u[2] -0.009271802  0.09011059 0.0000000
#> u[3] -0.222481299 -0.05862247 0.2565936
model$sigma_L %*% t(model$sigma_L)  # however Sigma is the same as above
#>              [,1]         [,2]         [,3]
#> [1,]  0.143406569 -0.003511144 -0.084251576
#> [2,] -0.003511144  0.008205885 -0.003219702
#> [3,] -0.084251576 -0.003219702  0.118774819

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
#> [1] 2.567456e-08

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
#> Warning: The tangent space of the equivalence class does not have dimension s^2=9 (sv[1]=2.60272101527713, sv[9]=4.98973089849112e-17)

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
#>          [,1]     [,2]     [,3]     [,4]    [,5]     [,6] [,7] [,8] [,9] [,10]
#> [1,] 0.142761 0.000000 0.000000 0.000000 0.00000 0.000000    0    0    0     0
#> [2,] 0.000000 0.098768 0.000000 0.000000 0.00000 0.000000    0    0    0     0
#> [3,] 0.000000 0.000000 0.071289 0.000000 0.00000 0.000000    0    0    0     0
#> [4,] 0.000000 0.000000 0.000000 0.034166 0.00000 0.000000    0    0    0     0
#> [5,] 0.000000 0.000000 0.000000 0.000000 0.00199 0.000000    0    0    0     0
#> [6,] 0.000000 0.000000 0.000000 0.000000 0.00000 0.000628    0    0    0     0
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
#>         [,1]     [,2]     [,3]     [,4]     [,5]     [,6] [,7] [,8] [,9] [,10]
#> [1,] 0.62452 0.000000 0.000000 0.000000 0.000000 0.000000    0    0    0     0
#> [2,] 0.00000 0.432345 0.000000 0.000000 0.000000 0.000000    0    0    0     0
#> [3,] 0.00000 0.000000 0.189838 0.000000 0.000000 0.000000    0    0    0     0
#> [4,] 0.00000 0.000000 0.000000 0.054704 0.000000 0.000000    0    0    0     0
#> [5,] 0.00000 0.000000 0.000000 0.000000 0.013332 0.000000    0    0    0     0
#> [6,] 0.00000 0.000000 0.000000 0.000000 0.000000 0.001229    0    0    0     0
#>      [,11] [,12]
#> [1,]     0     0
#> [2,]     0     0
#> [3,]     0     0
#> [4,]     0     0
#> [5,]     0     0
#> [6,]     0     0
```

# Kullback–Leibler divergence

Compute the Kullback-Leibler divergence between a "true" state space
model and an estimated state space model. The function only works for
square systems, where the "true" model is stable and the estimate is
(strictly) miniphase.

## Usage

``` r
KL_divergence(model, model_hat)
```

## Arguments

- model:

  [`stspmod()`](https://bfunovits.github.io/RLDM/reference/stspmod.md)
  object, true model

- model_hat:

  [`stspmod()`](https://bfunovits.github.io/RLDM/reference/stspmod.md)
  object, estimated model

## Value

KL divergence.

## Details

The KL divergence is computed as follows. Suppose \\y_t = k(z) u_t\\,
with \\\mathbf{E} u_t u_t' = \Sigma\\ is the true model, and let \\y_t =
h(z) u_t\\, with \\\mathbf{E} u_t u_t'=\Omega\\ denote the estimate.
W.l.o.g. we assume that the models are in innovation form, i.e. \\k(0) =
h(0) = I\\. The procedure computes the covariance matrix, \\\Delta =
\mathbf{E} e_t e_t'\\ say, of the *one-step-ahead prediction errors*
\\e_t = h^{-1}(z) k(z) u_t\\ and then the KL divergence \$\$ \mathrm{KL}
= (1/2)(\mathrm{tr}(\Omega^{-1}\Delta) - m -
\ln\det(\Omega^{-1}\Delta))) \$\$ Note that this procedure breaks down
if the transfer function \\h^{-1}(z) k(z)\\ is not stable. Therefore the
true models has to be stable and the estimated model has to strictly
miniphase.

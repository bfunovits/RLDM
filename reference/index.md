# Package index

## Basic Classes - Models

Models for rational transfer functions including input covariance.

- [`RLDM`](https://bfunovits.github.io/RLDM/reference/RLDM-package.md)
  [`RLDM-package`](https://bfunovits.github.io/RLDM/reference/RLDM-package.md)
  : A Collection of Tools for VARMA and State Space Processes
- [`armamod()`](https://bfunovits.github.io/RLDM/reference/armamod.md) :
  Constructor for LMFD (ARMA) Models
- [`rmfdmod()`](https://bfunovits.github.io/RLDM/reference/rmfdmod.md)
  **\[experimental\]** : Constructor for RMFD Models
- [`stspmod()`](https://bfunovits.github.io/RLDM/reference/stspmod.md) :
  Creator for stspmod class

## Templates and Parametrization

Parameter templates and model structures.

- [`fill_template()`](https://bfunovits.github.io/RLDM/reference/fill_template.md)
  [`extract_theta()`](https://bfunovits.github.io/RLDM/reference/fill_template.md)
  : Connect Deep Parameters with a Model
- [`model2template()`](https://bfunovits.github.io/RLDM/reference/model_structures.md)
  [`tmpl_arma_pq()`](https://bfunovits.github.io/RLDM/reference/model_structures.md)
  [`tmpl_arma_echelon()`](https://bfunovits.github.io/RLDM/reference/model_structures.md)
  [`tmpl_rmfd_pq()`](https://bfunovits.github.io/RLDM/reference/model_structures.md)
  [`tmpl_rmfd_echelon()`](https://bfunovits.github.io/RLDM/reference/model_structures.md)
  [`tmpl_stsp_full()`](https://bfunovits.github.io/RLDM/reference/model_structures.md)
  [`tmpl_stsp_ar()`](https://bfunovits.github.io/RLDM/reference/model_structures.md)
  [`tmpl_stsp_echelon()`](https://bfunovits.github.io/RLDM/reference/model_structures.md)
  : Model Structures
- [`tmpl_DDLC()`](https://bfunovits.github.io/RLDM/reference/local_model_structures.md)
  [`tmpl_GRAM()`](https://bfunovits.github.io/RLDM/reference/local_model_structures.md)
  : Local Model Structures
- [`tmpl_llm()`](https://bfunovits.github.io/RLDM/reference/STSmodel.md)
  [`tmpl_lltm()`](https://bfunovits.github.io/RLDM/reference/STSmodel.md)
  [`tmpl_cycle()`](https://bfunovits.github.io/RLDM/reference/STSmodel.md)
  [`tmpl_season()`](https://bfunovits.github.io/RLDM/reference/STSmodel.md)
  [`cbind_templates()`](https://bfunovits.github.io/RLDM/reference/STSmodel.md)
  : Structural Time Series Models
- [`is.template()`](https://bfunovits.github.io/RLDM/reference/is.template.md)
  : Check templates
- [`tmpl_sigma_L()`](https://bfunovits.github.io/RLDM/reference/tmpl_sigma_L.md)
  : sigma_L Structure

## Test Objects and Model Conversion

- [`r_model()`](https://bfunovits.github.io/RLDM/reference/r_model.md) :
  Generate a Random Model
- [`as.stspmod()`](https://bfunovits.github.io/RLDM/reference/as.stspmod.md)
  : Coerce to State Space Model
- [`innovation_form()`](https://bfunovits.github.io/RLDM/reference/innovation_form.md)
  **\[experimental\]** : Innovation Form state space Model
- [`test_armamod()`](https://bfunovits.github.io/RLDM/reference/test_armamod.md)
  : Create Test ARMA model
- [`test_stspmod()`](https://bfunovits.github.io/RLDM/reference/test_stspmod.md)
  : Create Test state space Model

## Derived Properties

Compute properties of rational models.

- [`autocov()`](https://bfunovits.github.io/RLDM/reference/autocov.md) :
  Autocovariance, Autocorelation and Partial Autocorrelation Function
- [`fevardec()`](https://bfunovits.github.io/RLDM/reference/fevardec.md)
  : Forecast Error Variance Decomposition
- [`freqresp()`](https://bfunovits.github.io/RLDM/reference/freqresp.md)
  : Frequency Response Function
- [`impresp()`](https://bfunovits.github.io/RLDM/reference/impresp.md) :
  Impulse Response Function
- [`spectrald()`](https://bfunovits.github.io/RLDM/reference/spectrald.md)
  : Spectral Density
- [`dft_3D()`](https://bfunovits.github.io/RLDM/reference/dft_3D.md) :
  Discrete Time Fourier Transform

## Output and Visualization

Display and format model information.

- [`str(`*`<armamod>`*`)`](https://bfunovits.github.io/RLDM/reference/str.md)
  [`str(`*`<stspmod>`*`)`](https://bfunovits.github.io/RLDM/reference/str.md)
  [`str(`*`<impresp>`*`)`](https://bfunovits.github.io/RLDM/reference/str.md)
  [`str(`*`<autocov>`*`)`](https://bfunovits.github.io/RLDM/reference/str.md)
  [`str(`*`<fevardec>`*`)`](https://bfunovits.github.io/RLDM/reference/str.md)
  [`str(`*`<freqresp>`*`)`](https://bfunovits.github.io/RLDM/reference/str.md)
  [`str(`*`<spectrald>`*`)`](https://bfunovits.github.io/RLDM/reference/str.md)
  : Display the Structure of Objects
- [`print(`*`<armamod>`*`)`](https://bfunovits.github.io/RLDM/reference/print.md)
  [`print(`*`<rmfdmod>`*`)`](https://bfunovits.github.io/RLDM/reference/print.md)
  [`print(`*`<stspmod>`*`)`](https://bfunovits.github.io/RLDM/reference/print.md)
  [`print(`*`<impresp>`*`)`](https://bfunovits.github.io/RLDM/reference/print.md)
  [`print(`*`<autocov>`*`)`](https://bfunovits.github.io/RLDM/reference/print.md)
  [`print(`*`<fevardec>`*`)`](https://bfunovits.github.io/RLDM/reference/print.md)
  [`print(`*`<freqresp>`*`)`](https://bfunovits.github.io/RLDM/reference/print.md)
  [`print(`*`<spectrald>`*`)`](https://bfunovits.github.io/RLDM/reference/print.md)
  : Print Methods
- [`plot(`*`<impresp>`*`)`](https://bfunovits.github.io/RLDM/reference/plot.md)
  [`plot(`*`<autocov>`*`)`](https://bfunovits.github.io/RLDM/reference/plot.md)
  [`plot(`*`<freqresp>`*`)`](https://bfunovits.github.io/RLDM/reference/plot.md)
  [`plot(`*`<spectrald>`*`)`](https://bfunovits.github.io/RLDM/reference/plot.md)
  : Plot Methods
- [`plot(`*`<fevardec>`*`)`](https://bfunovits.github.io/RLDM/reference/plot.fevardec.md)
  : Plot Forecast Error Variance Decomposition
- [`plot_prediction()`](https://bfunovits.github.io/RLDM/reference/plot_prediction.md)
  : Plot Forecasts
- [`predict(`*`<armamod>`*`)`](https://bfunovits.github.io/RLDM/reference/predict.md)
  [`predict(`*`<stspmod>`*`)`](https://bfunovits.github.io/RLDM/reference/predict.md)
  [`evaluate_prediction()`](https://bfunovits.github.io/RLDM/reference/predict.md)
  : Model Predictions
- [`zeroes(`*`<armamod>`*`)`](https://bfunovits.github.io/RLDM/reference/poles_and_zeroes.md)
  [`poles(`*`<armamod>`*`)`](https://bfunovits.github.io/RLDM/reference/poles_and_zeroes.md)
  [`zeroes(`*`<rmfdmod>`*`)`](https://bfunovits.github.io/RLDM/reference/poles_and_zeroes.md)
  [`poles(`*`<rmfdmod>`*`)`](https://bfunovits.github.io/RLDM/reference/poles_and_zeroes.md)
  [`zeroes(`*`<stspmod>`*`)`](https://bfunovits.github.io/RLDM/reference/poles_and_zeroes.md)
  [`poles(`*`<stspmod>`*`)`](https://bfunovits.github.io/RLDM/reference/poles_and_zeroes.md)
  : Poles and Zeroes

## Time Series Operations

Solve difference equations and simulate.

- [`solve_de()`](https://bfunovits.github.io/RLDM/reference/solve_de.md)
  [`solve_inverse_de()`](https://bfunovits.github.io/RLDM/reference/solve_de.md)
  : Solve (linear) Difference Equations
- [`sim()`](https://bfunovits.github.io/RLDM/reference/sim.md) :
  Simulate from a State Space or VARMA Model
- [`solve_RMFD_R()`](https://bfunovits.github.io/RLDM/reference/solve_RMFD_R.md)
  : Solve RMFD system for given inputs
- [`solve_inverse_RMFD_R()`](https://bfunovits.github.io/RLDM/reference/solve_inverse_RMFD_R.md)
  : Obtain Inputs of RMFD System for Given Data
- [`solve_ARMA_R()`](https://bfunovits.github.io/RLDM/reference/solve_ARMA.md)
  : Solve ARMA system
- [`toepl_fwd()`](https://bfunovits.github.io/RLDM/reference/toepl_fwd.md)
  [`toepl_inv()`](https://bfunovits.github.io/RLDM/reference/toepl_fwd.md)
  : Toeplitz Calculations

## RcppArmadillo Implementations - Solving

- [`outputs_ARMA_cpp()`](https://bfunovits.github.io/RLDM/reference/outputs_ARMA_cpp.md)
  : Outputs of an ARMA systems
- [`outputs_STSP_cpp()`](https://bfunovits.github.io/RLDM/reference/outputs_STSP_cpp.md)
  : Outputs of a statespace system
- [`residuals_ARMA_cpp()`](https://bfunovits.github.io/RLDM/reference/residuals_ARMA_cpp.md)
  : Residuals of an ARMA system
- [`residuals_STSP_cpp()`](https://bfunovits.github.io/RLDM/reference/residuals_STSP_cpp.md)
  : Residuals of a statespace system
- [`cll_theta_ARMA_cpp()`](https://bfunovits.github.io/RLDM/reference/cll_theta_ARMA_cpp.md)
  : Compute the (concentrated) conditional log likelihood for ARMA
  models described by a model template.
- [`cll_theta_STSP_cpp()`](https://bfunovits.github.io/RLDM/reference/cll_theta_STSP_cpp.md)
  : Compute the (concentrated) conditional log likelihood for a
  statespace system described by a model template.
- [`fbsolve_STSP_cpp()`](https://bfunovits.github.io/RLDM/reference/fbsolve_STSP_cpp.md)
  : Forward-backward solution of statespace systems
- [`solve_rmfd_cpp()`](https://bfunovits.github.io/RLDM/reference/solve_rmfd_cpp.md)
  : Simulating Output from an RMFD Model (or obtain residuals)

## Estimation Methods - AR

- [`est_ar()`](https://bfunovits.github.io/RLDM/reference/est_ar.md)
  [`est_ar_yw()`](https://bfunovits.github.io/RLDM/reference/est_ar.md)
  [`est_ar_dlw()`](https://bfunovits.github.io/RLDM/reference/est_ar.md)
  [`est_ar_ols()`](https://bfunovits.github.io/RLDM/reference/est_ar.md)
  : Estimate Autoregressive Models

## Estimation Methods - ARMA

- [`est_arma_hrk()`](https://bfunovits.github.io/RLDM/reference/est_arma_hrk.md)
  **\[deprecated\]** : Hannan, Rissanen, Kavalieris estimation procedure
- [`est_arma_hrk3()`](https://bfunovits.github.io/RLDM/reference/est_arma_hrk3.md)
  : Different version of HRK Procedure

## Estimation Methods - State Space

- [`est_stsp_aoki()`](https://bfunovits.github.io/RLDM/reference/subspace-helpers.md)
  [`est_stsp_cca()`](https://bfunovits.github.io/RLDM/reference/subspace-helpers.md)
  [`est_stsp_cca_sample()`](https://bfunovits.github.io/RLDM/reference/subspace-helpers.md)
  : Subspace Helper Methods
- [`est_stsp_ss()`](https://bfunovits.github.io/RLDM/reference/subspace-methods.md)
  : Estimate State Space Models with Subspace Methods
- [`estorder_SVC()`](https://bfunovits.github.io/RLDM/reference/subspace-order-estimates.md)
  [`estorder_IVC()`](https://bfunovits.github.io/RLDM/reference/subspace-order-estimates.md)
  [`estorder_max()`](https://bfunovits.github.io/RLDM/reference/subspace-order-estimates.md)
  [`estorder_rkH()`](https://bfunovits.github.io/RLDM/reference/subspace-order-estimates.md)
  [`estorder_MOE()`](https://bfunovits.github.io/RLDM/reference/subspace-order-estimates.md)
  : Helper Functions for Order Estimation

## Estimation Methods - Likelihood Based

- [`est_ML()`](https://bfunovits.github.io/RLDM/reference/est_ML.md)
  **\[superseded\]** : Maximum Likelihood Estimation
- [`ll()`](https://bfunovits.github.io/RLDM/reference/ll.md) : Log
  Likelihood Methods
- [`ll_FUN()`](https://bfunovits.github.io/RLDM/reference/ll_FUN.md) :
  Log Likelihood Function Factory
- [`ll_theta()`](https://bfunovits.github.io/RLDM/reference/ll_theta.md)
  **\[superseded\]** : Log-likelihood Given Deep Parameters
- [`ll_kf()`](https://bfunovits.github.io/RLDM/reference/ll_kf.md)
  [`ll_kf_cpp()`](https://bfunovits.github.io/RLDM/reference/ll_kf.md)
  [`ll_kf2_cpp()`](https://bfunovits.github.io/RLDM/reference/ll_kf.md)
  : Gaussian log Likelihood of a State Space Model
- [`ll_kf_theta_cpp()`](https://bfunovits.github.io/RLDM/reference/ll_kf_theta_cpp.md)
  : Compute the log likelihood for a statespace system described by a
  model template.

## Estimation Methods - Recursive Least Squares

- [`arx_rls_core()`](https://bfunovits.github.io/RLDM/reference/arx_rls_core.md)
  : RLS function
- [`rls_exp_cpp()`](https://bfunovits.github.io/RLDM/reference/rls_exp_cpp.md)
  [`rls_window_cpp()`](https://bfunovits.github.io/RLDM/reference/rls_exp_cpp.md)
  : (Multivariate) Recursive Least Squares

## Model Comparison

Compare and evaluate estimated models.

- [`compare_estimates()`](https://bfunovits.github.io/RLDM/reference/compare_estimates.md)
  : Compare Estimated Models
- [`pm_test()`](https://bfunovits.github.io/RLDM/reference/pm_test.md) :
  Portmanteau Test for Serial Correlation
- [`KL_divergence()`](https://bfunovits.github.io/RLDM/reference/KL_divergence.md)
  : Kullback–Leibler divergence

## Kalman Filtering and Riccati

Kalman filter and matrix equation solvers.

- [`kf()`](https://bfunovits.github.io/RLDM/reference/kf.md)
  [`kf_cpp()`](https://bfunovits.github.io/RLDM/reference/kf.md)
  [`kf2_cpp()`](https://bfunovits.github.io/RLDM/reference/kf.md) :
  Kalman Filter
- [`riccati()`](https://bfunovits.github.io/RLDM/reference/riccati.md) :
  Solve a discrete time, algebraic Riccati equation

## Datasets

- [`BQdata`](https://bfunovits.github.io/RLDM/reference/BQdata.md)
  [`BQdata_ts`](https://bfunovits.github.io/RLDM/reference/BQdata.md)
  [`BQdata_xts`](https://bfunovits.github.io/RLDM/reference/BQdata.md) :
  Blanchard/Quah (1989) dataset
- [`RSdata`](https://bfunovits.github.io/RLDM/reference/RSdata.md)
  [`RSdata_ts`](https://bfunovits.github.io/RLDM/reference/RSdata.md)
  [`RSdata_xts`](https://bfunovits.github.io/RLDM/reference/RSdata.md) :
  Ramey/Shapiro dataset

test_that("est_ar works", {

  for (isim in (1:10)) {
    m = sample(1:5, 1)
    p = sample(1:4, 1)
    n.obs = sample(c(100, 200, 500), 1)
    mean_estimate = sample(c('zero', 'sample', 'intercept'), 1)
    p.max = 10
    tmpl = tmpl_arma_pq(m = m, n = m, p = p, q = 0)
    model = r_model(tmpl, bpoles = 1, sd = 0.25)
    # make sure that the diagonal entries of sigma_L are non negative
    model$sigma_L = model$sigma_L %*% diag(x = sign(diag(model$sigma_L)), nrow = m, ncol = m)

    ###############################################################
    # reconstruct the true AR model from the population ACF

    true_acf = autocov(model, lag.max = 12, type = 'covariance')
    ARest = est_ar(true_acf, p.max = p.max, method = 'yule-walker', penalty = 1e-6)
    expect_equivalent(model, ARest$model)

    ###############################################################
    # simulate a sample

    y = sim(model, n.obs = n.obs, start = list(s1 = NA))$y

    ###############################################################
    # estimate the AR(p) model with the true order p

    # OLS
    ARest = est_ar(y, ic = 'max', p.max = p, method = 'ols', mean_estimate = mean_estimate)
    # check the log Likelihood
    p.opt = ARest$p
    expect_equivalent(ll(ARest$model, scale(y, center = ARest$y.mean, scale = FALSE),
                         'conditional', skip = p.opt), ARest$ll)

    # Yule-Walker and Durbin-Levinson-Whittle are equivalent (up to numerical errors)
    ARest = est_ar(y, ic = 'max', p.max = p, method = 'yule-walker', mean_estimate = mean_estimate)
    junk = est_ar(y, ic = 'max', p.max = p, method = 'durbin-levinson-whittle', mean_estimate = mean_estimate)
    expect_equivalent(ARest$model, junk$model)

    # alternatively we may first estimate the sample autocovariance function
    # note that the 'type' of the ACF is irrelevant
    sample_acf = autocov(y, type = 'correlation', demean = !(mean_estimate == 'zero'))
    junk = est_ar(sample_acf, ic = 'max', p.max = p, method = 'yule-walker')
    expect_equivalent(ARest$model, junk$model)

    ###############################################################
    # estimate the order of the AR model with 'AIC', estimate a model with intercept
    ARest = est_ar(y, ic = 'AIC', p.max = p.max, method = 'ols', mean_estimate = "intercept")
    # print(ARest$model)

    # compare with the stats::ar function
    ARest2 = stats::ar(y, aic = TRUE, order.max = p.max, method = 'ols', intercept = TRUE)
    if (all(is.finite(ARest2$aic))) {
      # sometimes the stats:ar routine fails / returns Inf values?
      # the estimated coefficients are equal
      expect_equivalent(unclass(ARest$model$sys$a)[,,-1,drop = FALSE],
                        -aperm(ARest2$ar, c(2,3,1)), tolerance = 1e-6, scale = 1)
      # also the AIC values are up to scaling equivalent
      expect_equivalent( ARest$stats[,'ic'] - min(ARest$stats[,'ic']), ARest2$aic/n.obs)
    }

    # compare with est_arma_hrk
    ARest = est_ar(y, ic = 'max', p.max = p, method = 'ols', mean_estimate = mean_estimate)
    junk = est_arma_hrk(y, NULL, tmpl, trace = FALSE, mean_estimate = mean_estimate)
    expect_equivalent(ARest$model, junk$model)
    expect_equivalent(ARest$ll, junk$ll)
  }


})

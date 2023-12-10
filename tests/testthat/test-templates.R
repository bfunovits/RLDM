test_that("template tools work", {
  for (itest in (1:10)) {

    m = sample(1:4, 1)
    n = sample(1:4, 1)
    nu = sample(1:3, m, replace = TRUE)
    p = sample(0:4, 1)
    q = sample(0:4, 1)
    sigma_L = sample(c('chol','symm','identity'), 1)

    # full arma(p,q) model #########################################################
    tmpl = tmpl_arma_pq(m, n, p, q, sigma_L = sigma_L)
    th = round(rnorm(tmpl$n.par, sd = 0.1), 3)
    model = fill_template(th, tmpl)
    # model
    testthat::expect_equivalent(th, extract_theta(model, tmpl, on_error = 'stop'))

    k = impresp(model, 5)
    k0 = matrix(unclass(k$irf)[,,1], nrow = m, ncol = n)
    testthat::expect_equivalent(k0[1:min(m,n), 1:n, drop = FALSE], diag(x = 1, nrow = min(m,n), ncol = n))

    # echelon arma model #########################################################
    tmpl = tmpl_arma_echelon(nu, n = n, sigma_L = sigma_L)
    th = round(rnorm(tmpl$n.par, sd = 0.1), 3)
    model = fill_template(th, tmpl)
    # model
    testthat::expect_equivalent(th, extract_theta(model, tmpl, on_error = 'stop'))

    k = impresp(model, lag.max = 2*sum(nu) + 1)
    k0 = matrix(unclass(k$irf)[,,1], nrow = m, ncol = n)
    testthat::expect_equivalent(k0[1:min(m,n), 1:n, drop = FALSE], diag(x = 1, nrow = min(m,n), ncol = n))

    testthat::expect_equivalent(pseries2nu(k$irf), nu)

    # full state space model #########################################################
    tmpl = tmpl_stsp_full(m, n, s = sum(nu), sigma_L = sigma_L)
    th = round(rnorm(tmpl$n.par, sd = 0.1), 3)
    model = fill_template(th, tmpl)
    # model
    testthat::expect_equivalent(th, extract_theta(model, tmpl, on_error = 'stop'))

    k = impresp(model, lag.max = 2*sum(nu) + 1)
    k0 = matrix(unclass(k$irf)[,,1], nrow = m, ncol = n)
    testthat::expect_equivalent(k0[1:min(m,n), 1:n, drop = FALSE], diag(x = 1, nrow = min(m,n), ncol = n))

    # ar state space model #########################################################
    tmpl = tmpl_stsp_ar(m, p, sigma_L = sigma_L)
    th = round(rnorm(tmpl$n.par, sd = 0.1), 3)
    model = fill_template(th, tmpl)
    # model
    testthat::expect_equivalent(th, extract_theta(model, tmpl, on_error = 'stop'))

    k = impresp(model, lag.max = 2*sum(nu) + 1)
    k0 = matrix(unclass(k$irf)[,,1], nrow = m, ncol = m)
    testthat::expect_equivalent(k0, diag(m))

    # echelon state space model #########################################################
    tmpl = tmpl_stsp_echelon(nu, n = n, sigma_L = sigma_L)
    th = round(rnorm(tmpl$n.par, sd = 0.1), 3)
    model = fill_template(th, tmpl)
    # model
    testthat::expect_equivalent(th, extract_theta(model, tmpl, on_error = 'stop'))

    k = impresp(model, lag.max = 2*sum(nu) + 1)
    k0 = matrix(unclass(k$irf)[,,1], nrow = m, ncol = n)
    testthat::expect_equivalent(k0[1:min(m,n), 1:n, drop = FALSE], diag(x = 1, nrow = min(m,n), ncol = n))

    testthat::expect_equivalent(pseries2nu(k$irf), nu)

    # full rmfd(p,q) model #########################################################
    # swap m,n
    tmpl = tmpl_rmfd_pq(m = n, n = m, p, q, sigma_L = sigma_L)
    th = round(rnorm(tmpl$n.par, sd = 0.1), 3)
    model = fill_template(th, tmpl)
    # model
    testthat::expect_equivalent(th, extract_theta(model, tmpl, on_error = 'stop'))

    k = impresp(model, lag.max = 2*sum(nu) + 1)
    k0 = matrix(unclass(k$irf)[,,1], nrow = n, ncol = m)
    testthat::expect_equivalent(k0[1:min(m,n), 1:m, drop = FALSE], diag(x = 1, nrow = min(m,n), ncol = m))

    # echelon rmfd(p,q) model #########################################################
    # swap m,n
    tmpl = tmpl_rmfd_echelon(nu, m = n, sigma_L = sigma_L)
    th = round(rnorm(tmpl$n.par, sd = 0.1), 3)
    model = fill_template(th, tmpl)
    # model
    testthat::expect_equivalent(th, extract_theta(model, tmpl, on_error = 'stop'))

    k = impresp(model, lag.max = 2*sum(nu) + 1)
    k0 = matrix(unclass(k$irf)[,,1], nrow = n, ncol = m)
    testthat::expect_equivalent(k0[1:min(m,n), 1:m, drop = FALSE], diag(x = 1, nrow = min(m,n), ncol = m))

    testthat::expect_equivalent(pseries2nu(t(k$irf)), nu)
  }


})

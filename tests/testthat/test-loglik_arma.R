test_that("loglik_arma works", {

  for (itest in (1:100)) {
    m = sample(1:3, 1)
    p = sample(0:3, 1)
    q = sample(0:3, 1)
    T = 50

    k0 = matrix(NA_real_, nrow =m, ncol = m)
    k0[upper.tri(k0)] = 0
    sys = matrix(NA_real_, nrow = m, ncol = m*(p+q+2))
    sys[1:m, 1:m] = diag(m)

    # k[0] = lower triangular
    sys[1:m, (m*(p+1)+1):(m*(p+2))] = k0
    model1 = armamod(structure(sys, order = c(m,m,p,q), class = c('lmfd','ratm')),
                     matrix(NA_real_, nrow = m, ncol = m))
    tmpl1 = model2template(model1, sigma_L = 'chol')

    # k[0] = I
    sys[1:m, (m*(p+1)+1):(m*(p+2))] = diag(m)
    model2 = armamod(structure(sys, order = c(m,m,p,q), class = c('lmfd','ratm')),
                     matrix(NA_real_, nrow = m, ncol = m))
    tmpl2 = model2template(model2, sigma_L = 'chol')

    model2 = r_model(template = tmpl2,
                     bpoles = 1, bzeroes = 1, ntrials.max = 500, sd = 0.25)
    diag(model2$sigma_L) = 1 # scale the diagonal entries of sigma_L
    th2 = extract_theta(model2, tmpl2, on_error = 'stop')


    # generate a sample
    y = sim(model2, n.obs = T, n.burn_in = 100)$y

    # compute sample covariance of residuals and set Sigma_L correspondingly
    model2$sigma_L = t(chol(crossprod(solve_inverse_de(model2$sys, y)$u) / T))
    th2 = extract_theta(model2, tmpl2, on_error = 'stop')

    k0 = matrix(rnorm(m*m), nrow = m, ncol = m)
    k0[upper.tri(k0)] = 0
    model1 = model2
    model1$sys = lmfd(a = model1$sys$a, b = model1$sys$b %r% k0)
    model1$sigma_L = solve(k0, model1$sigma_L)
    # model1
    th1 = extract_theta(model1, tmpl1, on_error = 'stop')

    llfun1 = ll_FUN(tmpl1, y, 'conditional')
    cllfun1 = ll_FUN(tmpl1, y, 'concentrated')
    llfun2 = ll_FUN(tmpl2, y, 'conditional')
    cllfun2 = ll_FUN(tmpl2, y, 'concentrated')

    out = c(ll(model1, y, 'conditional'), ll(model1, y, 'concentrated'), llfun1(th1), cllfun1(th1),
            ll(model2, y, 'conditional'), ll(model2, y, 'concentrated'), llfun2(th2), cllfun2(th2))
    testthat::expect_equivalent(diff(range(out)), 0, scale = 1, tolerance = 1e-6)
  }
})

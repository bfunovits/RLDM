test_that("ll works", {
  for (isim in (1:10)) {
    n.obs = 250
    m = sample(1:4,1)
    p = sample(0:2,1)
    q = sample(0:2,1)
    arma_tmpl = tmpl_arma_pq(m, m, p, q)
    arma_model = r_model(arma_tmpl, bpoles = 1, bzeroes = 1,  sd = 0.25)
    # scale sigma_L
    diag(arma_model$sigma_L) = 1

    # data
    y = sim(arma_model, n.obs = n.obs)$y

    arma_model$sigma_L = t(chol(crossprod(solve_inverse_de(arma_model$sys, y)$u) / n.obs))
    arma_th = extract_theta(arma_model, arma_tmpl, on_error = 'stop')

    stsp_model = stspmod(sys = as.stsp(arma_model$sys), sigma_L = arma_model$sigma_L)
    s = dim(stsp_model$sys)[3]
    stsp_tmpl = tmpl_stsp_full(m, m, s)
    stsp_th = extract_theta(stsp_model, stsp_tmpl)

    arma_ll = ll_FUN(arma_tmpl, y, skip = 0, 'conditional')
    arma_cll = ll_FUN(arma_tmpl, y, skip = 0, 'concentrated')
    stsp_ll = ll_FUN(stsp_tmpl, y, skip = 0, 'conditional')
    stsp_cll = ll_FUN(stsp_tmpl, y, skip = 0, 'concentrated')

    ll0 = ll(arma_model, y, 'conditional')


    expect_equal(ll0, ll(arma_model, y, 'conditional'))
    expect_equal(ll0, ll(arma_model, y, 'concentrated'))
    expect_equal(ll0, ll_theta(arma_th, arma_tmpl, y, 'conditional'))
    expect_equal(ll0, ll_theta(arma_th, arma_tmpl, y, 'concentrated'))
    expect_equal(ll0, arma_ll(arma_th))
    expect_equal(ll0, arma_cll(arma_th))

    expect_equal(ll0, ll(stsp_model, y, 'conditional'))
    expect_equal(ll0, ll(stsp_model, y, 'concentrated'))
    expect_equal(ll0, ll_theta(stsp_th, stsp_tmpl, y, 'conditional'))
    expect_equal(ll0, ll_theta(stsp_th, stsp_tmpl, y, 'concentrated'))
    expect_equal(ll0, stsp_ll(stsp_th))
    expect_equal(ll0, stsp_cll(stsp_th))

  }
})

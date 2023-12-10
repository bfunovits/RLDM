# test-freqresp #########################################

n.obs = 20
n.f = 20
n.lags = 100
scenarios = expand.grid(m = 1:2, n = 0:2, p = 0:3, q = 0:3)

set.seed(12345)

test_that("freqresp methods work", {
  for (i in (1:nrow(scenarios))) {
    # for (i in (1:1)) {
    m = scenarios$m[i]
    n = scenarios$n[i]
    p = scenarios$p[i]
    q = scenarios$q[i]

    # cat('m =', m, 'n =', n, 'p =', p, 'q =', q, 'n.obs =', n.obs,'\n')

    model = test_armamod(dim = c(m,n), degrees = c(p,q), bpoles = 1.2)

    frr = freqresp(model, n.f = n.f)

    frr1 = freqresp(impresp(model, lag.max = n.lags), n.f = n.f)
    expect_equal(frr, frr1)

    frr1 = freqresp(as.stspmod(model), n.f = n.f)
    expect_equal(frr, frr1)
  }
})

set.seed(NULL)

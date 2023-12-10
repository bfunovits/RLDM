# test-solve_de_lmfd #########################################

n.obs = 20
n.obs1 = 10
scenarios = expand.grid(m = 1:2, n = 0:2, p = 0:2, q = -1:1)

test_that("solve_de.lmfd works", {
  for (i in (1:nrow(scenarios))) {
    # for (i in (1:1)) {
    m = scenarios$m[i]
    n = scenarios$n[i]
    p = scenarios$p[i]
    q = scenarios$q[i]

    # cat('m =', m, 'n =', n, 'p =', p, 'q =', q, 'n.obs =', n.obs,'\n')

    model = test_armamod(dim = c(m,n), degrees = c(p,q), bpoles = 1)
    sys = model$sys

    # check sim()
    data = sim(model, n.obs = n.obs, n.burnin = 10)
    expect_equal(dim(data$y), c(n.obs, m))
    expect_equal(dim(data$u), c(n.obs, n))

    u = data$u
    data1 = solve_de(sys, data$u[(n.obs1+1):n.obs, , drop = FALSE],
                     u0 = data$u[1:n.obs1, , drop = FALSE],
                     y0 = data$y[1:n.obs1, , drop = FALSE])
    expect_equal(data1$y, data$y[(n.obs1+1):n.obs, , drop = FALSE])

    # check initial values
    pp = sample((0:(p+1)), 1)
    y0 = matrix(rnorm((p+1)*m), nrow = p+1, ncol = m)
    y0[iseq(1,pp), ] = 0
    qq = sample((0:(q+2)), 1)
    u0 = matrix(rnorm((q+2)*n), nrow = q+2, ncol = n)
    u0[iseq(1,qq), ] = 0

    data = solve_de(sys, u = u, y0 = y0, u0 = u0)
    data1 = solve_de(sys, u = u, y0 = y0[iseq(pp+1, p+1), , drop = FALSE], u0 = u0[iseq(qq+1, q+2), , drop = FALSE])
    expect_equal(data, data1)


    # check inverse system
    if ((n > 0) && (q>=0)) {
      # square system
      sys = test_lmfd(dim = c(n,n), degrees = c(p,q), bpoles = 1, bzeroes = 1)

      pp = sample((0:p), 1)
      y0 = matrix(rnorm(pp*n), nrow = pp, ncol = n)
      qq = sample((0:q), 1)
      u0 = matrix(rnorm(qq*n), nrow = qq, ncol = n)

      data = solve_de(sys, u = u, y0 = y0, u0 = u0)
      data1 = solve_inverse_de(sys, y = data$y, y0 = y0, u0 = u0)
      expect_equal(data, data1)
    }
  }
})

# test-solve_de_stsp #########################################

n.obs = 20
n.obs1 = 10
scenarios = expand.grid(m = 1:2, n = 0:2, s = 0:2)

test_that("solve_de.stsp works", {
  for (i in (1:nrow(scenarios))) {
    # for (i in (1:1)) {
    m = scenarios$m[i]
    n = scenarios$n[i]
    s = scenarios$s[i]

    # cat('m =', m, 'n =', n, 's =', s, 'n.obs =', n.obs,'\n')

    model = test_stspmod(dim = c(m,n), s = s, bpoles = 1)
    sys = model$sys

    # check sim()
    data = sim(model, n.obs = n.obs, n.burnin = 10)
    expect_equal(dim(data$y), c(n.obs, m))
    expect_equal(dim(data$u), c(n.obs, n))
    expect_equal(dim(data$a), c(n.obs+1, s))

    u = data$u
    data1 = solve_de(sys, data$u[(n.obs1+1):n.obs, , drop = FALSE],
                     a1 = data$a[n.obs1+1, ])
    expect_equal(data1$y, data$y[(n.obs1+1):n.obs, , drop = FALSE])

    # check inverse system
    if (n > 0) {
      # square system
      sys = test_stsp(dim = c(n,n), s = s, bpoles = 1, bzeroes = 1)

      a1 = rnorm(s)

      data = solve_de(sys, u = u, a1 = a1)
      data1 = solve_inverse_de(sys, y = data$y, a1 = a1)
      expect_equal(data, data1)
    }
  }
})

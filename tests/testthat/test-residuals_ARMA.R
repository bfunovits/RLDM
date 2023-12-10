tmpl_arma_full = function(m,p,q) {
  if (m)
  if ( (m<=0) || (min(m,p) < 0) || ((p+q)<=0) ) stop('illegal inputs')
  n.par = (m^2)*(p+q+2)
  h = array(0, dim =c(m,m,p+q+2))
  h[,,1] = diag(m)
  h[,,p+2] = diag(m)
  tmpl = list(h = c(h, diag(m)),
              H = rbind(diag(n.par), matrix(0, nrow = m^2, ncol = n.par)),
              class = "armamod",
              order = c(m,m,p,q), n.par = n.par)
}

res_gr = function(th, tmpl, y) {
  order = tmpl$order
  m = order[1]
  p = order[3]
  q = order[4]
  model = fill_template(th, tmpl)
  a = unclass(model$sys$a)
  b = unclass(model$sys$b)
  ib0 = solve(matrix(b[,,1], nrow = m, ncol = m))
  B1 = -ib0 %*% matrix(b[,,rev(iseq(2,q+1))], nrow = m, ncol = m*q)
  A = ib0 %*% matrix(a, nrow = m, ncol = m*(p+1))
  y = t(y)
  u = matrix(0, nrow = m, ncol = ncol(y))
  dU = matrix(0, nrow = m*ncol(y), ncol = (m^2)*(p+q+2))

  residuals_ARMA_cpp(ib0, B1, A, 1, y, u, dU)
  return(list(u = t(u), dU = dU))
}

eps = sqrt(.Machine$double.eps)


test_that("resARMA works", {
  for ( i in (1:10)) {
    n.obs = sample(c(50, 100), 1)
    m = sample(1:2, 1)
    p = 0
    q = 0
    while ((p+q) == 0) {
      p = sample(0:2, 1)
      q = sample(0:2, 1)
    }
    tmpl = tmpl_arma_full(m, p , q)

    model = r_model(tmpl, bzeroes = 1, bpoles = 1, sd = 0.5,  ntrials.max = 1000)
    th = extract_theta(model, tmpl, on_error = 'stop')

    data = sim(model, n.obs = n.obs)
    out = res_gr(th, tmpl, data$y)
    u = out$u
    dU = out$dU

    # check residuals
    expect_equivalent(u, solve_inverse_de(model$sys, data$y)$u)

    # check derivative with respect to the parameters "th"
    dth = rnorm(length(th), sd = eps)

    # finite difference
    model1 = fill_template(th+dth, tmpl)
    u1 = solve_inverse_de(model1$sys, data$y)$u
    du = u1 - u

    # analytic derivative
    du0 = dU %*% dth
    du0 = t(matrix(du0, nrow = m, ncol = n.obs))

    expect_equivalent(du, du0, scale = mean(abs(u)), tolerance = 10*eps*sqrt(10*eps))
    # print(mean(abs(de -de0)))
  }
})

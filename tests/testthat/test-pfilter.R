test_that("pfilter runs without errors", {
  set.seed(123)
  m <- 1
  s <- 2
  n <- s
  tmpl <- tmpl_stsp_full(m, n, s, sigma_L = "chol")
  model <- r_model(tmpl, bpoles = 1, sd = 0.5)
  data <- sim(model, n.obs = 20)
  y <- data$y

  # Test SIR filter
  expect_no_error({
    pf_sir <- pfilter(model, y, N_particles = 100, method = "sir")
  })
  expect_true(inherits(pf_sir, "pfilter"))
  expect_named(pf_sir, c("filtered_states", "predicted_states", "particles", "weight_trajectories", "weights",
                         "log_likelihood", "log_likelihood_contributions",
                         "effective_sample_size", "N_particles", "resampling_method"))
  expect_equal(pf_sir$N_particles, 100)
  expect_equal(length(pf_sir$weights), 100)
  expect_equal(dim(pf_sir$filtered_states), c(21, s))
  expect_equal(dim(pf_sir$particles), c(s, 100, 21))

  # Test APF filter
  expect_no_error({
    pf_apf <- pfilter(model, y, N_particles = 100, method = "apf")
  })
  expect_true(inherits(pf_apf, "pfilter"))
  expect_equal(pf_apf$N_particles, 100)
})

test_that("ll_pfilter works", {
  set.seed(456)
  m <- 1
  s <- 2
  n <- s
  tmpl <- tmpl_stsp_full(m, n, s, sigma_L = "chol")
  model <- r_model(tmpl, bpoles = 1, sd = 0.5)
  data <- sim(model, n.obs = 30)
  y <- data$y

  # Test SIR likelihood approximation
  expect_no_error({
    ll_sir <- ll_pfilter(model, y, N_particles = 500, filter_type = "sir", N_runs = 2)
  })
  expect_type(ll_sir, "double")
  expect_length(ll_sir, 1)
  expect_false(is.na(ll_sir))

  # Test APF likelihood approximation
  expect_no_error({
    ll_apf <- ll_pfilter(model, y, N_particles = 500, filter_type = "apf", N_runs = 2)
  })
  expect_type(ll_apf, "double")
  expect_length(ll_apf, 1)
  expect_false(is.na(ll_apf))
})

test_that("plot.pfilter works without errors", {
  set.seed(789)
  m <- 1
  s <- 2
  n <- s
  tmpl <- tmpl_stsp_full(m, n, s, sigma_L = "chol")
  model <- r_model(tmpl, bpoles = 1, sd = 0.5)
  data <- sim(model, n.obs = 15)
  y <- data$y

  pf <- pfilter(model, y, N_particles = 50, method = "sir")

  # Test each plot type
  expect_no_error(plot(pf, type = "states"))
  expect_no_error(plot(pf, type = "ess"))
  expect_no_error(plot(pf, type = "weights"))
  expect_no_error(plot(pf, type = "likelihood"))
})

test_that("different resampling methods work", {
  set.seed(101)
  m <- 1
  s <- 2
  n <- s
  tmpl <- tmpl_stsp_full(m, n, s, sigma_L = "chol")
  model <- r_model(tmpl, bpoles = 1, sd = 0.5)
  data <- sim(model, n.obs = 10)
  y <- data$y

  for (resamp in c("systematic", "multinomial", "stratified")) {
    expect_no_error({
      pfilter(model, y, N_particles = 100, method = "sir", resampling = resamp)
    })
    expect_no_error({
      pfilter(model, y, N_particles = 100, method = "apf", resampling = resamp)
    })
  }
})

test_that("pfilter handles edge cases", {
  set.seed(202)
  m <- 1
  s <- 0  # zero-dimensional state space
  n <- m
  tmpl <- tmpl_stsp_full(m, n, s, sigma_L = "chol")
  model <- r_model(tmpl, bpoles = 1, sd = 0.5)
  data <- sim(model, n.obs = 5)
  y <- data$y

  # Should work with zero states
  expect_no_error({
    pf <- pfilter(model, y, N_particles = 50, method = "sir")
  })
  expect_equal(dim(pf$filtered_states), c(6, 0))
  expect_equal(dim(pf$particles), c(0, 50, 6))
})

test_that("pfilter returns consistent filtered states", {
  set.seed(303)
  m <- 1
  s <- 2
  n <- s
  tmpl <- tmpl_stsp_full(m, n, s, sigma_L = "chol")
  model <- r_model(tmpl, bpoles = 1, sd = 0.5)
  data <- sim(model, n.obs = 10)
  y <- data$y

  pf <- pfilter(model, y, N_particles = 200, method = "sir")

  # Filtered states should be weighted average of particles
  for (t in 1:11) {
    particles_t <- pf$particles[,,t]
    weights_t <- pf$weight_trajectories[t,]
    weighted_avg <- particles_t %*% weights_t
    expect_equal(as.numeric(pf$filtered_states[t,]), as.numeric(weighted_avg),
                 tolerance = 1e-10)
  }
})
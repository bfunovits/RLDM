test_that("R and Rcpp versions of solve_de and solve_inverse_de get same result (square case,  positive degrees)", {
  dim_in = 2;
  dim_out = 2;
  deg_c = 1;
  deg_d = 1;
  n_obs = 10;

  set.seed(1000)
  (robj = test_rmfd(dim = c(dim_out, dim_in), degrees = c(deg_c, deg_d)))
  (data_in = matrix(rnorm(dim_in*n_obs), nrow = dim_in, ncol = n_obs))
  (data_out = matrix(rnorm(dim_out*n_obs), nrow = dim_out, ncol = n_obs))
  robj
  robj$c
  data_in
  data_out

  #debugonce(solve_RMFD_R)
  (y1 = solve_RMFD_R(robj$c, robj$d, data_input = data_in, t0 = 1))
  #debugonce(solve_de)
  (y2 = solve_de(robj, u = t(data_in))$y %>% t)

  expect_equivalent(y1, y2)

  #debugonce(solve_inverse_RMFD_R)
  (u1 = solve_inverse_RMFD_R(robj$c, robj$d, data_output = data_out, t0 = 1))
  #debugonce(solve_inverse_de)
  (u2 = solve_inverse_de(robj, y = t(data_out))$u %>% t())

  expect_equivalent(u1, u2)
} )

test_that("R and Rcpp versions of solve_de and solve_inverse_de get same result (non-square case, positive degrees)", {
dim_in = 2;
dim_out = 3;
deg_c = 1;
deg_d = 1;
n_obs = 10;

set.seed(1000)
(robj = test_rmfd(dim = c(dim_out, dim_in), degrees = c(deg_c, deg_d)))
(data_in = matrix(rnorm(dim_in*n_obs), nrow = dim_in, ncol = n_obs))
(data_out = matrix(rnorm(dim_out*n_obs), nrow = dim_out, ncol = n_obs))
robj
robj$c
data_in
data_out

#debugonce(solve_RMFD_R)
(y1 = solve_RMFD_R(robj$c, robj$d, data_input = data_in, t0 = 1))
#debugonce(solve_de)
(y2 = solve_de(robj, u = t(data_in))$y %>% t)

expect_equivalent(y1, y2)

#debugonce(solve_inverse_RMFD_R)
(u1 = solve_inverse_RMFD_R(robj$c, robj$d, data_output = data_out, t0 = 1))
#debugonce(solve_inverse_de)
(u2 = solve_inverse_de(robj, y = t(data_out))$u %>% t())

expect_equivalent(u1, u2)
})


test_that("R and Rcpp versions of solve_de and solve_inverse_de get same result (non-square case, deg_c == 0)", {
dim_in = 2;
dim_out = 3;
deg_c = 0;
deg_d = 1;
n_obs = 10;

set.seed(1000)
(robj = test_rmfd(dim = c(dim_out, dim_in), degrees = c(deg_c, deg_d)))
(data_in = matrix(rnorm(dim_in*n_obs), nrow = dim_in, ncol = n_obs))
(data_out = matrix(rnorm(dim_out*n_obs), nrow = dim_out, ncol = n_obs))
robj
robj$c
data_in
data_out

#debugonce(solve_RMFD_R)
(y1 = solve_RMFD_R(robj$c, robj$d, data_input = data_in, t0 = 1))
#debugonce(solve_de)
(y2 = solve_de(robj, u = t(data_in))$y %>% t)

expect_equivalent(y1, y2)

#debugonce(solve_inverse_RMFD_R)
(u1 = solve_inverse_RMFD_R(robj$c, robj$d, data_output = data_out, t0 = 1))
#debugonce(solve_inverse_de)
(u2 = solve_inverse_de(robj, y = t(data_out))$u %>% t())

expect_equivalent(u1, u2)
})

test_that("R and Rcpp versions of solve_de and solve_inverse_de get same result (non-square case, deg_d == 0)", {

dim_in = 2;
dim_out = 3;
deg_c = 1;
deg_d = 0;
n_obs = 10;

set.seed(1000)
(robj = test_rmfd(dim = c(dim_out, dim_in), degrees = c(deg_c, deg_d)))
(data_in = matrix(rnorm(dim_in*n_obs), nrow = dim_in, ncol = n_obs))
(data_out = matrix(rnorm(dim_out*n_obs), nrow = dim_out, ncol = n_obs))
robj
robj$c
data_in
data_out

#debugonce(solve_RMFD_R)
(y1 = solve_RMFD_R(robj$c, robj$d, data_input = data_in, t0 = 1))
#debugonce(solve_de)
(y2 = solve_de(robj, u = t(data_in))$y %>% t)

expect_equivalent(y1, y2)

#debugonce(solve_inverse_RMFD_R)
(u1 = solve_inverse_RMFD_R(robj$c, robj$d, data_output = data_out, t0 = 1))
#debugonce(solve_inverse_de)
(u2 = solve_inverse_de(robj, y = t(data_out))$u %>% t())

expect_equivalent(u1, u2)
})

devtools::load_all()
set.seed(123)
tmpl <- tmpl_stsp_full(1,1,2,sigma_L='chol')
model <- r_model(tmpl, bpoles=1, sd=0.5)
data <- sim(model, n.obs=50, a1=NA)
y <- data$y
kf_result <- kf(model, y)
cat('KF log-likelihood:', kf_result$ll, '\n')

particle_counts <- c(100, 500, 1000, 5000, 10000)
results <- data.frame(N = particle_counts, rmse = NA, ll_diff = NA, ess_min = NA)
for (N in particle_counts) {
  cat('\nTesting N =', N, '\n')
  pf_result <- pfilter(model, y, N_particles = N, method = 'sir')
  rmse <- sqrt(mean((kf_result$a[1:50,] - pf_result$filtered_states[1:50,])^2))
  ll_diff <- pf_result$log_likelihood - kf_result$ll
  ess_min <- min(pf_result$effective_sample_size, na.rm=TRUE)
  results[results$N == N, ] <- c(N, rmse, ll_diff, ess_min)
  cat('RMSE:', rmse, 'LL diff:', ll_diff, 'ESS min:', ess_min, '\n')
}
print(results)
# Plot convergence
plot(results$N, results$rmse, type='b', log='x', xlab='Number of particles', ylab='RMSE', main='PF convergence to KF')
plot(results$N, results$ll_diff, type='b', log='x', xlab='Number of particles', ylab='LL difference')
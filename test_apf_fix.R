#!/usr/bin/env Rscript
devtools::load_all()
set.seed(123)
m <- 1
s <- 2
n <- s
tmpl <- tmpl_stsp_full(m, n, s, sigma_L = "chol")
model <- r_model(tmpl, bpoles = 1, sd = 0.5)
data <- sim(model, n.obs=20, a1=NA)
y <- data$y
cat('Model generated\n')

# Kalman filter baseline
kf_result <- kf(model, y)
cat('KF log-likelihood:', kf_result$ll, '\n')

# Test SIR
cat('\n--- SIR particle filter ---\n')
pf_sir <- pfilter(model, y, N_particles = 1000, method = 'sir')
cat('SIR log-likelihood:', pf_sir$log_likelihood, '\n')
cat('ESS min:', min(pf_sir$effective_sample_size, na.rm=TRUE), '\n')
cat('ESS max:', max(pf_sir$effective_sample_size, na.rm=TRUE), '\n')
cat('Final weight summary:\n')
print(summary(pf_sir$weights))

# Test APF
cat('\n--- APF particle filter ---\n')
pf_apf <- pfilter(model, y, N_particles = 1000, method = 'apf')
cat('APF log-likelihood:', pf_apf$log_likelihood, '\n')
cat('ESS min:', min(pf_apf$effective_sample_size, na.rm=TRUE), '\n')
cat('ESS max:', max(pf_apf$effective_sample_size, na.rm=TRUE), '\n')
cat('Final weight summary:\n')
print(summary(pf_apf$weights))

# Compare RMSE
kf_states <- kf_result$a[1:20,]
sir_states <- pf_sir$filtered_states[1:20,]
apf_states <- pf_apf$filtered_states[1:20,]
rmse_sir <- sqrt(mean((kf_states - sir_states)^2))
rmse_apf <- sqrt(mean((kf_states - apf_states)^2))
cat('\nRMSE vs KF: SIR =', rmse_sir, 'APF =', rmse_apf, '\n')

# Plot weight distributions
pdf(file=NULL)
par(mfrow=c(1,2))
hist(pf_sir$weights, breaks=30, main='SIR weights', xlab='Weight')
abline(v=1/1000, col='red', lty=2)
hist(pf_apf$weights, breaks=30, main='APF weights', xlab='Weight')
abline(v=1/1000, col='red', lty=2)
dev.off()

# Test likelihood approximation with multiple runs
cat('\n--- Likelihood approximation (averaged) ---\n')
ll_sir <- ll_pfilter(model, y, N_particles = 2000, N_runs = 5, filter_type = 'sir')
ll_apf <- ll_pfilter(model, y, N_particles = 2000, N_runs = 5, filter_type = 'apf')
cat('SIR averaged LL:', ll_sir, '\n')
cat('APF averaged LL:', ll_apf, '\n')
cat('Difference from KF:', ll_sir - kf_result$ll, '(SIR)', ll_apf - kf_result$ll, '(APF)\n')

cat('\nTest completed.\n')
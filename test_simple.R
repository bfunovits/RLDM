#!/usr/bin/env Rscript
devtools::load_all()
set.seed(123)
tmpl <- tmpl_stsp_full(1,1,2,sigma_L='chol')
model <- r_model(tmpl, bpoles=1, sd=0.5)
data <- sim(model, n.obs=50, a1=NA)
y <- data$y
cat('y dim:', dim(y), '\n')
cat('y class:', class(y), '\n')
# Kalman filter
kf_result <- kf(model, y)
cat('KF log-likelihood:', kf_result$ll, '\n')
# Particle filter
pf_result <- pfilter(model, y, N_particles=500, method='sir')
cat('PF log-likelihood:', pf_result$log_likelihood, '\n')
cat('Success!\n')
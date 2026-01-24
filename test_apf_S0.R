devtools::load_all()
set.seed(123)
m <- 1
s <- 2
n <- s
# Create model with diagonal sigma_L (no cross-correlation)
tmpl <- tmpl_stsp_full(m, n, s, sigma_L = "diag")
model <- r_model(tmpl, bpoles = 1, sd = 0.5)
data <- sim(model, n.obs=20, a1=NA)
y <- data$y
cat('Model with diagonal sigma_L (S=0)\n')
cat('S matrix (cross-covariance):\n')
sigma = tcrossprod(model$sigma_L)
S = model$sys$B %*% sigma %*% t(model$sys$D)
print(S)
# Kalman filter baseline
kf_result <- kf(model, y)
cat('KF log-likelihood:', kf_result$ll, '\n')
# Test SIR
pf_sir <- pfilter(model, y, N_particles = 1000, method = 'sir')
cat('SIR log-likelihood:', pf_sir$log_likelihood, '\n')
# Test APF
pf_apf <- pfilter(model, y, N_particles = 1000, method = 'apf')
cat('APF log-likelihood:', pf_apf$log_likelihood, '\n')
# Compare
cat('Difference from KF:', pf_sir$log_likelihood - kf_result$ll, '(SIR)', pf_apf$log_likelihood - kf_result$ll, '(APF)\n')
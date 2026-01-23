devtools::load_all()
set.seed(123)
m <- 1
s <- 2
n <- s
tmpl <- tmpl_stsp_full(m, n, s, sigma_L = "chol")
model <- r_model(tmpl, bpoles = 1, sd = 0.5)
data <- sim(model, n.obs = 5, a1 = NA)
y <- data$y

# Extract matrices
A <- model$sys$A
C <- model$sys$C
sigma <- tcrossprod(model$sigma_L)
R <- model$sys$D %*% sigma %*% t(model$sys$D)
S <- model$sys$B %*% sigma %*% t(model$sys$D)
Q <- model$sys$B %*% sigma %*% t(model$sys$B)

cat("A:\n"); print(A)
cat("C:\n"); print(C)
cat("Q:\n"); print(Q)
cat("R:\n"); print(R)
cat("S:\n"); print(S)

# Compute predictive covariance for optimal proposal
S_cov <- C %*% Q %*% t(C) + R + C %*% S + t(S) %*% t(C)
cat("S_cov (CQC' + R + CS + S'C'):\n"); print(S_cov)

# Compute Kalman gain
K <- (Q %*% t(C) + S) %*% solve(S_cov)
cat("Kalman gain K:\n"); print(K)

# Compute proposal covariance Sigma
Sigma <- Q - K %*% (C %*% Q + t(S))
cat("Proposal covariance Sigma:\n"); print(Sigma)

# Regularization epsilon
eps <- 1e-10 * norm(S_cov, "F")
cat("Regularization eps:", eps, "\n")
S_cov_reg <- S_cov + eps * diag(m)
cat("Regularized S_cov:\n"); print(S_cov_reg)

# Compare with Kalman filter innovation covariance from kf()
kf_result <- kf(model, y)
cat("\nKalman filter innovation covariances (first few):\n")
for (t in 1:min(3, nrow(y))) {
  cat("t=", t, "F_t (innovation covariance):\n")
  print(kf_result$F[,,t])
}

# Compute steady-state solution via DARE? Not needed.

# Run optimal filter with 5 particles and print internal matrices from C++?
# We'll compute predictive likelihood for each particle manually
pf <- pfilter(model, y, N_particles = 5, method = "optimal")
cat("\nParticle filter log-likelihood:", pf$log_likelihood, "\n")
cat("Kalman filter log-likelihood:", ll_kf(model, y), "\n")

# Let's compute predictive likelihood for first particle at t=1
particles_t0 <- pf$particles[,,1]  # t=0 initial particles
x_pred <- A %*% particles_t0[,1]   # prediction for particle 1
y_pred <- C %*% x_pred
residual <- y[1,] - t(y_pred)
scaled_residual <- solve(t(chol(S_cov_reg)), residual)
log_pdf <- -0.5 * (m * log(2*pi) + 2*sum(log(diag(chol(S_cov_reg)))) + sum(scaled_residual^2))
cat("\nManual predictive log pdf for particle 1, t=1:", log_pdf, "\n")
cat("Residual:", residual, "\n")
cat("S_cov_reg determinant:", det(S_cov_reg), "\n")
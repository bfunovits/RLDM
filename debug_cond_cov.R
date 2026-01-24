devtools::load_all()
set.seed(123)
m <- 1
s <- 2
n <- s
tmpl <- tmpl_stsp_full(m, n, s, sigma_L = "chol")
model <- r_model(tmpl, bpoles = 1, sd = 0.5)
sigma <- tcrossprod(model$sigma_L)
Q <- model$sys$B %*% sigma %*% t(model$sys$B)
R <- model$sys$D %*% sigma %*% t(model$sys$D)
S <- model$sys$B %*% sigma %*% t(model$sys$D)
C <- model$sys$C
cat('Q:\n'); print(Q)
cat('R:\n'); print(R)
cat('S:\n'); print(S)
cat('C:\n'); print(C)
# Regularize
eps_q <- 1e-10 * norm(Q, 'F')
if (eps_q == 0) eps_q <- 1e-10
Q_reg <- Q + eps_q * diag(s)
eps_r <- 1e-10 * norm(R, 'F')
if (eps_r == 0) eps_r <- 1e-10
R_reg <- R + eps_r * diag(m)
Q_inv <- solve(Q_reg)
cond_cov <- R - t(S) %*% Q_inv %*% S
cond_cov2 <- R_reg - t(S) %*% Q_inv %*% S
cat('cond_cov:\n'); print(cond_cov)
cat('cond_cov2 (regularized R):\n'); print(cond_cov2)
cat('Eigenvalues of cond_cov:\n'); print(eigen(cond_cov)$values)
cat('Eigenvalues of cond_cov2:\n'); print(eigen(cond_cov2)$values)
# Compute S_cov with cross terms
S_cov <- C %*% Q %*% t(C) + R + C %*% S + t(S) %*% t(C)
cat('S_cov:\n'); print(S_cov)
cat('Eigenvalues of S_cov:\n'); print(eigen(S_cov)$values)
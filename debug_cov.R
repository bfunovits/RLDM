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
S_cov_wrong <- C %*% Q %*% t(C) + R
S_cov_correct <- C %*% Q %*% t(C) + R + C %*% S + t(S) %*% t(C)
cat('S_cov (wrong):\n'); print(S_cov_wrong)
cat('S_cov (correct):\n'); print(S_cov_correct)
cat('Difference:\n'); print(S_cov_correct - S_cov_wrong)
# Check if difference is positive definite
eig <- eigen(S_cov_correct - S_cov_wrong)
cat('Eigenvalues of diff:\n'); print(eig$values)
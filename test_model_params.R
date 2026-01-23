#!/usr/bin/env Rscript
devtools::load_all()

set.seed(123)
m <- 1
s <- 2
n <- s
tmpl <- tmpl_stsp_full(m, n, s, sigma_L = "chol")
model <- r_model(tmpl, bpoles = 1, sd = 0.5)

cat("Model parameters:\n")
cat("A:\n"); print(model$sys$A)
cat("B:\n"); print(model$sys$B)
cat("C:\n"); print(model$sys$C)
cat("D:\n"); print(model$sys$D)
cat("sigma_L:\n"); print(model$sigma_L)

sigma <- tcrossprod(model$sigma_L)
R <- model$sys$D %*% sigma %*% t(model$sys$D)
S <- model$sys$B %*% sigma %*% t(model$sys$D)
Q <- model$sys$B %*% sigma %*% t(model$sys$B)

cat("\nDerived matrices:\n")
cat("Q (state noise cov):\n"); print(Q)
cat("R (obs noise cov):\n"); print(R)
cat("S (cross cov):\n"); print(S)

# Check if S is zero
cat("\nIs S zero?", all(abs(S) < 1e-10), "\n")

# Compute predictive covariance
pred_cov <- C %*% Q %*% t(C) + R + C %*% S + t(S) %*% t(C)
cat("\nPredictive covariance CQC' + R + CS + S'C':\n"); print(pred_cov)

# Compute Kalman gain for optimal proposal
K_opt <- (Q %*% t(C) + S) %*% solve(pred_cov)
cat("\nOptimal Kalman gain (Q C' + S) (C Q C' + R + C S + S' C')^{-1}:\n"); print(K_opt)

# Check eigenvalues
cat("\nEigenvalues of Q:", eigen(Q)$values, "\n")
cat("Eigenvalues of R:", eigen(R)$values, "\n")
cat("Eigenvalues of predictive covariance:", eigen(pred_cov)$values, "\n")
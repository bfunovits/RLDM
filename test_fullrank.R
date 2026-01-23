devtools::load_all()
set.seed(123)
m <- 1
s <- 2
n <- s  # inputs = states
tmpl <- tmpl_stsp_full(m, n, s, sigma_L = "chol")
model <- r_model(tmpl, bpoles = 1, sd = 0.5)
print(model$sys$B)
sigma <- tcrossprod(model$sigma_L)
Q <- model$sys$B %*% sigma %*% t(model$sys$B)
cat('Q eigenvalues:', eigen(Q, symmetric=TRUE, only.values=TRUE)$values, '\n')
data <- sim(model, n.obs=50, a1=NA)
y <- data$y
kf_result <- kf(model, y)
cat('KF log-likelihood:', kf_result$ll, '\n')
pf_result <- pfilter(model, y, N_particles=5000, method='sir')
cat('PF log-likelihood:', pf_result$log_likelihood, '\n')
rmse <- sqrt(mean((kf_result$a[1:50,] - pf_result$filtered_states[1:50,])^2))
cat('RMSE:', rmse, '\n')
cat('ESS min:', min(pf_result$effective_sample_size, na.rm=TRUE), '\n')
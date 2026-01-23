#!/usr/bin/env Rscript
devtools::load_all()
set.seed(123)
m <- 1
s <- 2
n <- s
tmpl <- tmpl_stsp_full(m, n, s, sigma_L = "chol")
model <- r_model(tmpl, bpoles = 1, sd = 0.5)
data <- sim(model, n.obs=10, a1=NA)  # small sample for debugging
y <- data$y
cat('Model parameters:\n')
print(model$sys$A)
print(model$sys$B)
print(model$sys$C)
print(model$sys$D)
sigma <- tcrossprod(model$sigma_L)
Q <- model$sys$B %*% sigma %*% t(model$sys$B)
R <- model$sys$D %*% sigma %*% t(model$sys$D)
S <- model$sys$B %*% sigma %*% t(model$sys$D)
cat('\nQ eigenvalues:', eigen(Q, symmetric=TRUE, only.values=TRUE)$values, '\n')
cat('R:', R, '\n')

# Run particle filter
N_particles <- 100
pf_result <- pfilter(model, y, N_particles = N_particles, method = 'sir')
cat('\nParticle filter results:\n')
cat('Log-likelihood:', pf_result$log_likelihood, '\n')
cat('ESS min:', min(pf_result$effective_sample_size, na.rm=TRUE), '\n')
cat('ESS max:', max(pf_result$effective_sample_size, na.rm=TRUE), '\n')

# Examine weights at final time
weights <- pf_result$weights
cat('\nFinal weights summary:\n')
print(summary(weights))
cat('Weight sum:', sum(weights), '\n')
cat('Number of zero weights:', sum(weights == 0), '\n')
cat('Effective sample size (1/sum(w^2)):', 1/sum(weights^2), '\n')

# Plot weight distribution
hist(weights, breaks=30, main='Final particle weights', xlab='Weight')
abline(v=1/N_particles, col='red', lty=2)

# Let's manually compute weights for the first time step to verify
# Get particles at time 0 (initial)
particles_t0 <- pf_result$particles[,,1]  # s x N_particles
cat('\nInitial particles shape:', dim(particles_t0), '\n')

# We can't easily recompute weights because we need the internal random seeds
# But we can examine the particle trajectories
cat('\nParticle trajectories shape:', dim(pf_result$particles), '\n')
cat('Filtered states shape:', dim(pf_result$filtered_states), '\n')

# Compute mean and variance of particles over time
particle_means <- apply(pf_result$particles, c(1,3), mean)
cat('\nParticle means vs filtered states (first few time steps):\n')
print(head(cbind(particle_means[1,], pf_result$filtered_states[,1])))

# Check if filtered states equal weighted average of particles
weighted_avg <- t(pf_result$particles[,,11]) %*% weights
cat('\nWeighted avg of particles at final time:', weighted_avg, '\n')
cat('Filtered state at final time:', pf_result$filtered_states[11,], '\n')

# Examine log-likelihood contributions
cat('\nLog-likelihood contributions:\n')
print(pf_result$log_likelihood_contributions)

# Test with more particles to see if degeneracy persists
cat('\n\n--- Testing with 5000 particles ---\n')
pf_result2 <- pfilter(model, y, N_particles = 5000, method = 'sir')
weights2 <- pf_result2$weights
cat('ESS min:', min(pf_result2$effective_sample_size, na.rm=TRUE), '\n')
cat('Weight summary:\n')
print(summary(weights2))
cat('Number of zero weights:', sum(weights2 == 0), '\n')
hist(weights2, breaks=30, main='Final particle weights (N=5000)', xlab='Weight')
abline(v=1/5000, col='red', lty=2)
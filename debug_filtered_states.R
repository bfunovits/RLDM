#!/usr/bin/env Rscript
devtools::load_all()

cat("Debugging filtered state consistency issue...\n")

set.seed(123)
m <- 1
s <- 2
n <- s
tmpl <- tmpl_stsp_full(m, n, s, sigma_L = "chol")
model <- r_model(tmpl, bpoles = 1, sd = 0.5)
data <- sim(model, n.obs = 5, a1 = NA)  # Very small sample for debugging
y <- data$y

cat("Model created with s =", s, "states\n")
cat("Data: n_obs =", nrow(y), "\n\n")

# Run particle filter with small N for debugging
N_particles <- 10
cat("Running SIR filter with N =", N_particles, "particles...\n")
pf <- pfilter(model, y, N_particles = N_particles, method = "sir")

cat("\n=== Checking filtered state consistency ===\n")

# Check each time point
for (t in 1:6) {  # t=1 initial, t=6 final (5 observations + 1 initial)
  particles_t <- pf$particles[,,t]

  # Use weight trajectories
  weights_t <- pf$weight_trajectories[t,]
  cat(sprintf("\nt=%d: using weight_trajectories\n", t-1))

  # Compute weighted average
  weighted_avg <- as.numeric(particles_t %*% weights_t)
  filtered_state <- as.numeric(pf$filtered_states[t,])

  cat(sprintf("  Filtered state from output: [%s]\n",
              paste(formatC(filtered_state, digits=4), collapse=", ")))
  cat(sprintf("  Weighted avg of particles: [%s]\n",
              paste(formatC(weighted_avg, digits=4), collapse=", ")))

  diff <- abs(filtered_state - weighted_avg)
  cat(sprintf("  Difference: [%s]\n",
              paste(formatC(diff, digits=4), collapse=", ")))

  if (any(diff > 1e-10)) {
    cat("  ** MISMATCH DETECTED **\n")

    # Examine particles and weights
    cat("  Particles (columns):\n")
    print(particles_t)
    cat("  Weights:\n")
    print(weights_t)
    cat("  Weight sum:", sum(weights_t), "\n")
  }
}

cat("\n=== Examining C++ output structure ===\n")
cat("Names in pf result:", names(pf), "\n")
cat("Filtered states dim:", dim(pf$filtered_states), "\n")
cat("Particles dim:", dim(pf$particles), "\n")
cat("Weights length:", length(pf$weights), "\n")
cat("ESS length:", length(pf$effective_sample_size), "\n")

cat("\n=== Hypothesis: weights are only final weights ===\n")
cat("The C++ code may only return final weights, not weights at each time.\n")
cat("If true, filtered_states[t,] should use weights at time t, but we only have final weights.\n")

# Let's manually compute what filtered states should be using only final weights
cat("\n=== Recomputing filtered states using final weights ===\n")
for (t in 1:6) {
  particles_t <- pf$particles[,,t]
  recomputed <- as.numeric(particles_t %*% pf$weights)
  original <- as.numeric(pf$filtered_states[t,])
  diff <- abs(recomputed - original)

  cat(sprintf("t=%d: Original [%s], Recomp [%s], Diff [%s]\n",
              t-1,
              paste(formatC(original, digits=4), collapse=", "),
              paste(formatC(recomputed, digits=4), collapse=", "),
              paste(formatC(diff, digits=4), collapse=", ")))
}

cat("\n=== Checking if filtered_states are particle means (not weighted) ===\n")
for (t in 1:6) {
  particles_t <- pf$particles[,,t]
  particle_mean <- as.numeric(rowMeans(particles_t))
  filtered_state <- as.numeric(pf$filtered_states[t,])
  diff <- abs(particle_mean - filtered_state)

  cat(sprintf("t=%d: Filtered [%s], Mean [%s], Diff [%s]\n",
              t-1,
              paste(formatC(filtered_state, digits=4), collapse=", "),
              paste(formatC(particle_mean, digits=4), collapse=", "),
              paste(formatC(diff, digits=4), collapse=", ")))
}

cat("\n=== Conclusion ===\n")
cat("Based on this debugging, the issue appears to be:\n")
cat("1. C++ returns only final weights (not weights at each time)\n")
cat("2. Filtered_states[t,] should be weighted average using weights at time t\n")
cat("3. But we're comparing with final weights, causing mismatch\n")
cat("\nSolution: Either:\n")
cat("1. Store weights at each time in C++ and return them\n")
cat("2. Or document that filtered_states use final weights (incorrect)\n")
cat("3. Or fix C++ to compute filtered_states correctly\n")
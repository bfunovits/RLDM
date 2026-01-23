#!/usr/bin/env Rscript

library(RLDM)

# Check if pf from RLDM is different from stats::pf
cat("Checking function origins...\n")
cat("pf from RLDM:", environmentName(environment(pf)), "\n")
cat("pf from stats:", environmentName(environment(stats::pf)), "\n")

# Check if our generic exists
cat("\nChecking generic...\n")
cat("isGeneric('pf'):", isGeneric("pf"), "\n")
cat("isS3generic('pf'):", isS3generic("pf"), "\n")
cat("isS3method('pf'):", isS3method("pf"), "\n")

# Check methods
cat("\nMethods for pf:\n")
print(methods("pf"))

# Create a simple model
cat("\nCreating test model...\n")
tmpl <- tmpl_stsp_full(1, 1, 2, sigma_L = "chol")
model <- r_model(tmpl, bpoles = 1, sd = 0.5)
data <- sim(model, n.obs = 20)
y <- data$y

# Try calling with namespace
cat("\nTrying RLDM::pf...\n")
result1 <- try(RLDM::pf(model, y, N_particles = 200, method = "sir"))
if (inherits(result1, "try-error")) {
  cat("Error with RLDM::pf:", result1, "\n")
} else {
  cat("Success! Log-likelihood:", result1$log_likelihood, "\n")
}

# Try calling without namespace (should use RLDM's pf after loading)
cat("\nTrying pf (unqualified)...\n")
result2 <- try(pf(model, y, N_particles = 200, method = "sir"))
if (inherits(result2, "try-error")) {
  cat("Error with pf:", result2, "\n")
} else {
  cat("Success! Log-likelihood:", result2$log_likelihood, "\n")
}

# Test ll_pf
cat("\nTesting ll_pf...\n")
ll <- try(ll_pf(model, y, N_particles = 500, N_runs = 2))
if (inherits(ll, "try-error")) {
  cat("Error with ll_pf:", ll, "\n")
} else {
  cat("Success! Log-likelihood:", ll, "\n")
}

cat("\nDone.\n")
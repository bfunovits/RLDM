#!/usr/bin/env Rscript

# Test with explicit namespace
library(RLDM)

cat("Testing with namespace...\n")

# Check what's exported
cat("Exported names containing 'pfilter':\n")
exports <- getNamespaceExports(getNamespace("RLDM"))
print(exports[grep("pfilter", exports)])

# Check if RLDM::pfilter works
cat("\nTrying RLDM::pfilter...\n")
tmpl <- tmpl_stsp_full(1, 1, 1, sigma_L = "chol")
model <- r_model(tmpl, bpoles = 1, sd = 0.5)
data <- sim(model, n.obs = 10)
y <- data$y

result <- try(RLDM::pfilter(model, y, N_particles = 100))
if (inherits(result, "try-error")) {
  cat("Error:", result[1], "\n")
} else {
  cat("Success!\n")
}

# Check if it's in search path
cat("\nSearch path check:\n")
for (env in search()) {
  if (exists("pfilter", where=env, inherits=FALSE)) {
    cat("Found pfilter in:", env, "\n")
  }
}
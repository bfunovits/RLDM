#!/usr/bin/env Rscript

cat("Debugging export issues...\n")

# Clean load
devtools::load_all()

cat("\n1. Checking namespace...\n")
ns <- getNamespace("RLDM")
cat("Objects in namespace containing 'pfilter':\n")
print(ls(ns)[grep("pfilter", ls(ns))])

cat("\n2. Checking exports...\n")
exports <- getNamespaceExports(ns)
cat("Exports containing 'pfilter':\n")
print(exports[grep("pfilter", exports)])

cat("\n3. Checking if generic is S3 generic...\n")
cat("isGeneric('pfilter'):", isGeneric("pfilter"), "\n")
cat("isS3Generic('pfilter'):", tryCatch(isS3generic("pfilter"), error=function(e) "error"), "\n")

cat("\n4. Checking generic definition...\n")
gen <- get("pfilter", envir = ns)
cat("Class:", class(gen), "\n")
cat("Is function:", is.function(gen), "\n")
cat("Body starts with:", as.character(body(gen))[1], "\n")

cat("\n5. Trying to call from namespace...\n")
tmpl <- tmpl_stsp_full(1, 1, 1, sigma_L = "chol")
model <- r_model(tmpl, bpoles = 1, sd = 0.5)
data <- sim(model, n.obs = 5)
y <- data$y

cat("   Attempt 1: RLDM::pfilter...")
result1 <- try(RLDM::pfilter(model, y, N_particles = 50))
if (inherits(result1, "try-error")) {
  cat(" ERROR:", result1[1], "\n")
} else {
  cat(" SUCCESS\n")
}

cat("   Attempt 2: pfilter (unqualified)...")
result2 <- try(pfilter(model, y, N_particles = 50))
if (inherits(result2, "try-error")) {
  cat(" ERROR:", result2[1], "\n")
} else {
  cat(" SUCCESS\n")
}

cat("\n6. Checking search path...\n")
cat("Search path:\n")
print(search())
cat("\nWhere is pfilter found?\n")
for (env in search()) {
  if (exists("pfilter", where=env, inherits=FALSE)) {
    cat("  Found in:", env, "\n")
  }
}

cat("\n7. Checking if package is attached...\n")
cat("Is RLDM in search path?", any(grepl("package:RLDM", search())), "\n")
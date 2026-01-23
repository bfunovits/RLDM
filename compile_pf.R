#!/usr/bin/env Rscript

cat("Compiling particle filter extensions for RLDM package...\n")

# Load required packages
if (!requireNamespace("Rcpp", quietly = TRUE)) {
  stop("Rcpp package is required")
}
if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("devtools package is required")
}

# Set working directory to package root
pkg_root <- getwd()
cat("Package root:", pkg_root, "\n")

# 1. Compile Rcpp attributes
cat("Running Rcpp::compileAttributes()...\n")
Rcpp::compileAttributes(pkgdir = pkg_root, verbose = TRUE)

# 2. Load the package to test compilation
cat("Loading package with devtools::load_all()...\n")
devtools::load_all(pkg_root, quiet = FALSE)

# 3. Test that functions are available
cat("Testing function availability...\n")
if (exists(".Call")) {
  cat("Checking .Call wrappers...\n")
  # Check if Rcpp exports were generated
  rcpp_exports <- file.path(pkg_root, "R", "RcppExports.R")
  if (file.exists(rcpp_exports)) {
    cat("RcppExports.R exists\n")
  } else {
    cat("WARNING: RcppExports.R not found\n")
  }

  # Try to find our functions
  env <- asNamespace("RLDM")
  funcs <- ls(env)
  pf_funcs <- grep("pf|ll_pf", funcs, value = TRUE)
  cat("Found functions:", paste(pf_funcs, collapse = ", "), "\n")
}

cat("Done!\n")
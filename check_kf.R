devtools::load_all(quiet=TRUE)
set.seed(123)
m <- 1; s <- 2; n <- s
tmpl <- tmpl_stsp_full(m, n, s, sigma_L = "chol")
model <- r_model(tmpl, bpoles = 1, sd = 0.5)
data <- sim(model, n.obs = 3, a1 = NA)
y <- data$y

kf_result <- kf(model, y)
cat("Names in kf_result:", names(kf_result), "\n")
cat("Class:", class(kf_result), "\n")
cat("Structure:\n")
str(kf_result, max.level=2)

# Try to access F
if ("F" %in% names(kf_result)) {
  cat("F dim:", dim(kf_result$F), "\n")
  cat("F[,,1]:", kf_result$F[,,1], "\n")
  cat("F[1,1,1]:", kf_result$F[1,1,1], "\n")
}

# Check ll_kf
cat("\nll_kf:", ll_kf(model, y), "\n")
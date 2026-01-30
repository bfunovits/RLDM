#!/usr/bin/env Rscript
# Script to identify functions missing examples

library(stringr)

# Get all R files (excluding RcppExports.R)
r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
r_files <- r_files[!grepl("RcppExports\\.R$", r_files)]

# Parse NAMESPACE to get exported functions
namespace_lines <- readLines("NAMESPACE")
exported_functions <- character(0)

for (line in namespace_lines) {
  if (grepl("^export\\(", line)) {
    # Extract function name from export(...)
    func_name <- gsub("^export\\(([^)]+)\\).*", "\\1", line)
    # Remove quotes if present
    func_name <- gsub('"', '', func_name)
    func_name <- gsub("'", '', func_name)
    # Handle special case for %>%
    if (func_name == "%>%") next  # Skip magrittr pipe
    exported_functions <- c(exported_functions, func_name)
  }
}

cat("Total exported functions (excluding %>%):", length(exported_functions), "\n")

# Function to check if a function has examples
check_function_examples <- function(file_path, func_name) {
  lines <- readLines(file_path)

  # Find the function definition
  func_pattern <- paste0("^", func_name, "\\s*<-\\s*function|^", func_name, "\\s*=")
  func_line <- which(grepl(func_pattern, lines))

  if (length(func_line) == 0) {
    # Try S3 method pattern
    func_pattern <- paste0("\\.", func_name, "\\s*<-\\s*function|\\.", func_name, "\\s*=")
    func_line <- which(grepl(func_pattern, lines))
  }

  if (length(func_line) == 0) {
    return(list(has_examples = FALSE, reason = "Function not found in file"))
  }

  # Find the start of the roxygen block (go backwards from function)
  start_line <- func_line[1]
  while (start_line > 1 && grepl("^#'", lines[start_line - 1])) {
    start_line <- start_line - 1
  }

  # Check if this is actually a roxygen block (might be just comments)
  if (start_line == func_line[1] || !any(grepl("^#'", lines[start_line:func_line[1]]))) {
    return(list(has_examples = FALSE, reason = "No roxygen block found"))
  }

  # Extract roxygen block
  roxygen_block <- lines[start_line:(func_line[1] - 1)]

  # Check for @examples tag
  has_examples_tag <- any(grepl("^#'\\s*@examples", roxygen_block))

  # Check for @example tag (singular)
  has_example_tag <- any(grepl("^#'\\s*@example", roxygen_block))

  # Check for actual example code (lines starting with #' #)
  has_example_code <- any(grepl("^#'\\s*#", roxygen_block))

  has_examples <- has_examples_tag || has_example_tag || has_example_code

  if (!has_examples) {
    return(list(has_examples = FALSE, reason = "No @examples/@example tag or example code"))
  }

  # Check example quality
  example_quality <- "good"
  issues <- character(0)

  if (has_examples_tag) {
    # Find the @examples line
    examples_line <- which(grepl("^#'\\s*@examples", roxygen_block))[1]
    # Check next few lines for example code
    example_start <- examples_line + 1
    example_end <- min(examples_line + 20, length(roxygen_block))
    example_lines <- roxygen_block[example_start:example_end]

    # Check for set.seed in stochastic functions
    is_stochastic <- func_name %in% c("sim", "r_model", "pfilter", "ll_pfilter", "est_ML",
                                      "est_ar", "est_arma_hrk", "est_arma_hrk3", "est_stsp_cca",
                                      "est_stsp_cca_sample", "est_stsp_aoki", "est_stsp_ss")

    if (is_stochastic) {
      example_text <- paste(example_lines, collapse = " ")
      if (!grepl("set\\.seed", example_text)) {
        issues <- c(issues, "Missing set.seed() for stochastic function")
        example_quality <- "needs_improvement"
      }
    }

    # Check for actual R code (not just comments)
    has_r_code <- any(grepl("^#'\\s*[^#]", example_lines))
    if (!has_r_code) {
      issues <- c(issues, "No executable R code in examples")
      example_quality <- "needs_improvement"
    }
  }

  return(list(
    has_examples = TRUE,
    quality = example_quality,
    issues = if (length(issues) > 0) issues else "none"
  ))
}

# Check each exported function
results <- data.frame(
  function_name = character(),
  file = character(),
  has_examples = logical(),
  quality = character(),
  issues = character(),
  stringsAsFactors = FALSE
)

for (func in exported_functions) {
  cat("Checking:", func, "\n")

  found <- FALSE
  for (r_file in r_files) {
    # Quick check if function might be in this file
    file_content <- readLines(r_file, warn = FALSE)
    if (any(grepl(paste0("\\b", func, "\\b"), file_content))) {
      check_result <- check_function_examples(r_file, func)
      if (!is.null(check_result$reason) && check_result$reason == "Function not found in file") {
        next
      }

      results <- rbind(results, data.frame(
        function_name = func,
        file = basename(r_file),
        has_examples = check_result$has_examples,
        quality = if (check_result$has_examples) check_result$quality else "none",
        issues = if (check_result$has_examples) paste(check_result$issues, collapse = "; ") else check_result$reason,
        stringsAsFactors = FALSE
      ))
      found <- TRUE
      break
    }
  }

  if (!found) {
    results <- rbind(results, data.frame(
      function_name = func,
      file = "unknown",
      has_examples = FALSE,
      quality = "none",
      issues = "Function not found in any R file",
      stringsAsFactors = FALSE
    ))
  }
}

# Summary statistics
total_functions <- nrow(results)
functions_with_examples <- sum(results$has_examples)
functions_needing_improvement <- sum(results$quality == "needs_improvement")
coverage_percentage <- round(functions_with_examples / total_functions * 100, 1)

cat("\n=== SUMMARY ===\n")
cat("Total exported functions:", total_functions, "\n")
cat("Functions with examples:", functions_with_examples, "\n")
cat("Functions needing example improvement:", functions_needing_improvement, "\n")
cat("Example coverage:", coverage_percentage, "%\n")

# Functions missing examples
missing_examples <- results[!results$has_examples, ]
cat("\n=== FUNCTIONS MISSING EXAMPLES (", nrow(missing_examples), ") ===\n", sep = "")
if (nrow(missing_examples) > 0) {
  for (i in 1:nrow(missing_examples)) {
    cat(i, ". ", missing_examples$function_name[i], " (", missing_examples$file[i], "): ",
        missing_examples$issues[i], "\n", sep = "")
  }
} else {
  cat("None!\n")
}

# Functions needing improvement
needs_improvement <- results[results$quality == "needs_improvement", ]
cat("\n=== FUNCTIONS NEEDING EXAMPLE IMPROVEMENT (", nrow(needs_improvement), ") ===\n", sep = "")
if (nrow(needs_improvement) > 0) {
  for (i in 1:nrow(needs_improvement)) {
    cat(i, ". ", needs_improvement$function_name[i], " (", needs_improvement$file[i], "): ",
        needs_improvement$issues[i], "\n", sep = "")
  }
} else {
  cat("None!\n")
}

# Save detailed results
write.csv(results, "example_check_results.csv", row.names = FALSE)
cat("\nDetailed results saved to: example_check_results.csv\n")

# Print files that need modification
files_to_modify <- unique(results$file[!results$has_examples | results$quality == "needs_improvement"])
files_to_modify <- files_to_modify[files_to_modify != "unknown"]
cat("\n=== FILES NEEDING MODIFICATION ===\n")
cat(paste(files_to_modify, collapse = "\n"), "\n")
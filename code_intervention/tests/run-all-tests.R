## Run all intervention-side test files in this directory.
##
## Must be run with this tests/ directory as the working directory, or
## sourced from there. Each test file manages its own setwd("../") to
## reach code_intervention/ for sourcing intervention code.
##
##   cd tests && Rscript run-all-tests.R
##   source("run-all-tests.R")   # from R with tests/ as wd

library(testthat)


options(intervention_every_gene_cores = parallel::detectCores())
## For laptop
if (system("hostname", intern = TRUE) == "Triturus")
  intervention_every_gene_cores <- 4

tests_files <- sort(dir(pattern = glob2rx("*-TESTS.R")))

if (length(tests_files) == 0)
    stop("No *-TESTS.R files found. Is the working directory tests/?")

for (f in tests_files) {
    cat("\n#######################################\n")
    cat("#### Running", f, "\n")
    cat("#######################################\n\n")
    results <- test_file(f)
    n_failures <- sum(as.data.frame(results)$failed)
    if (n_failures > 0) {
        stop(n_failures, " failure(s) in ", f, ". Aborting.")
    }
}

cat("\n\n#######################################\n")
cat("#### ALL INTERVENTION-SIDE TESTS PASSED\n")
cat("#######################################\n")


## Interpreting warnings:
## Check that there are no failures. That is the key. Some warnings are
## intentional (they ARE the point of certain tests):
##   - "weighted_fgraph contains unreachable destinations":
##       kill-gene-equivalences-TESTS.R,
##       intervene-fitness-landscapes-TESTS.R,
##       miscellanea-and-former-issues-TESTS.R.
##   - "No accessible genotypes": same.
##   - "Model has 0 rows": in miscellanea-and-former-issues-TESTS.R under
##       "Issue 1 is solved", intentional.
##   - "Issue 1 is solved": my own warning() call used as a label/marker
##       for a test context.

## Run all intervention-side test files in this directory.
##
## Must be run with this tests/ directory as the working directory, or
## sourced from there. Each test file manages its own setwd("../") to
## reach code_intervention/ for sourcing intervention code.
##

##   cd tests &&  R --vanilla -f run-all-tests.R &> run-all-tests.Rout
##   from R with tests/ as wd



## This is probably not needed, since in each individual file
## but I've tripped on this so often, that it is worth adding it here
## too.
## options(intervention_every_gene_cores = parallel::detectCores())
## markovchain, Armadillo/BLAS solve, and genots_from_trm, etc,
## can be multithreaded and we are using mclapply with detectCores()
## The Sys.setenv call will affect programs called from R via system()
## or system2(). But this won't affect the BLAS loaded when R starts
## and to limit those threads we use RhpcBLASctl.
Sys.setenv(OMP_NUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)

version
date()

## Make sure all packages we need exist. Otherwise, we will abort
## immediately

source("../dependencies.R", echo = TRUE)


options(intervention_every_gene_cores = parallel::detectCores())
## For my laptop (Triturus), use only 4 cores. This must set the option
## (not a variable): the code reads it with
## getOption("intervention_every_gene_cores").
## Sys.info() is base R and, unlike system("hostname"), does not need
## the hostname program to exist (it might not, e.g., in a minimal
## Docker image).
if (Sys.info()[["nodename"]] == "Triturus")
  options(intervention_every_gene_cores = 4)

tests_files <- sort(dir(pattern = glob2rx("*-TESTS.R")))

if (length(tests_files) == 0)
    stop("No *-TESTS.R files found. Is the working directory tests/?")

for (f in tests_files) {
    cat("\n#######################################\n")
    cat("#### Running", f, "\n")
    cat("#######################################\n\n")
    results <- test_file(f)
    ## Why we check more than the "failed" column.
    ##
    ## test_file() returns a table with one row per test_that block.
    ## Its "failed" column counts only failed expectations (e.g., an
    ## expect_true() that got FALSE). There are several ways in which
    ## a test file can be broken and still have failed = 0. If we only
    ## checked "failed", those broken files would pass silently. So we
    ## check each of these cases:
    ##
    ## 1. A package loaded with library() is not installed. testthat
    ##    does not report this as a failure or as an error. Instead, it
    ##    records no tests at all, so the table has zero rows. Zero
    ##    rows means nothing was tested, so we abort. (dependencies.R,
    ##    sourced above, should already catch the packages we know
    ##    about; this check catches any we missed.)
    ##
    ## 2. An error happens outside any test_that block; for example,
    ##    suppose `source("intervention.R")` at the top of a test file
    ##    fails. testthat records a single row with error = TRUE, failed =
    ##    0 and passed = 0. But none of the tests in that file ran.
    ##
    ## 3. An error happens inside a test_that block. That block's row
    ##    gets error = TRUE, but failed stays 0, because an error is
    ##    not a failed expectation. And the rest of that block does not
    ##    run.
    ##
    ##    For cases 2 and 3, we abort if any row has error = TRUE.
    ##
    ## 4. Skipped tests give skipped = TRUE. A test_that block with no
    ##    expectations also counts as skipped. Because no tests in our
    ##    test suite are skipped on purpose, if a test is skipped that
    ##    means something is broken. Thus, we abort as soon as a test
    ##    is skipped.
    ##
    ## 5. As a last safety net, we abort if no expectation passed in
    ##    the file.

    res_df <- as.data.frame(results)
    if (nrow(res_df) == 0)
        stop("No test results recorded in ", f, ". Aborting.")
    n_failures <- sum(res_df$failed)
    n_errors <- sum(res_df$error)
    n_skipped <- sum(res_df$skipped)
    n_passed <- sum(res_df$passed)
    if (n_failures > 0 || n_errors > 0 || n_skipped > 0 || n_passed == 0) {
        stop("In ", f, ": ", n_failures, " failure(s), ", n_errors,
             " error(s), ", n_skipped, " skip(s), ", n_passed,
             " passed expectation(s). Aborting.")
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
##   - "No acessible genotypes": same.
##   - "Model has 0 rows": in miscellanea-and-former-issues-TESTS.R under
##       "Issue 1 is solved", intentional.
##   - "Issue 1 is solved": my own warning() call used as a label/marker
##       for a test context.

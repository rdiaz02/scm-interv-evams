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

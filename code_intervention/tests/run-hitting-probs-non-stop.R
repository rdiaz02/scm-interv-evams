## Run, non-stop, only the tests of the hitting probabilities and of
## HyperHMM. These are the three test files:
##
##   - hitting-probs-direct-and-drop-unreachable-TESTS.R
##   - kill-HyperHMM-TESTS.R
##   - trm-TESTS.R
##
## Why: all the numerical problems found so far by the full non-stop
## loop (run-tests-non-stop.R) were in this code. The full loop spends
## most of its time elsewhere, so this one tries many more random seeds
## of just this code per hour (about 1.3 minutes per loop).
##
## This just runs run-tests-non-stop.R, in a new R process (the same R
## as the one running this file), giving it the names of the three test
## files. Everything else (output in non-stop-logs/, the FAILED-loop
## files, summary.log, how to replay a failure) is as explained in
## run-tests-non-stop.R.
##
## Must be run with this tests/ directory as the working directory:
##
##   cd tests && R --vanilla -f run-hitting-probs-non-stop.R \
##        &> run-hitting-probs-non-stop.Rout

tests_to_run <- c("hitting-probs-direct-and-drop-unreachable-TESTS.R",
                  "kill-HyperHMM-TESTS.R",
                  "trm-TESTS.R")

R_bin <- file.path(R.home("bin"), "R")
system2(R_bin, c("--vanilla", "-f", "run-tests-non-stop.R",
                 "--args", tests_to_run))

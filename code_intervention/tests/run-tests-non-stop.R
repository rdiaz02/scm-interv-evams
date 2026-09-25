## Run run-all-tests.R over and over, until killed (Ctrl-C, or kill).
##
## Why: many tests use random data, with a new seed each run (see
## set_and_print_seed in ../intervention.R). A test can pass most of the
## time and fail only for some seeds. Running the whole suite for, say,
## a day, is a way to catch those.
##
## Must be run with this tests/ directory as the working directory:
##
##   cd tests && R --vanilla -f run-tests-non-stop.R &> run-tests-non-stop.Rout
##
## Each loop runs run-all-tests.R in a new R process (the same R as the
## one running this file), exactly as in the header of run-all-tests.R,
## so that no state (seeds, sourced code, options) carries over from one
## loop to the next.
##
## Output, in directory non-stop-logs/:
##
##   - current.Rout: output of the loop now running (overwritten at the
##     start of each loop).
##
##   - FAILED-loop-<n>-<date-time>.Rout: output of each loop that failed,
##     kept. To replay a failure, look for the last "Seed used was"
##     before the failure (see set_and_print_seed in ../intervention.R).
##
##   - summary.log: one line per finished loop (loop number, start time,
##     minutes, PASSED/FAILED), plus running totals. It is written as we
##     go, so the count of loops done is there even after this is
##     killed.
##
## A failed loop does not stop this: we keep going, to see how often
## failures happen. Set stop_on_failure to TRUE to stop at the first one.

stop_on_failure <- FALSE

if (!file.exists("run-all-tests.R"))
    stop("run-all-tests.R not found. Is the working directory tests/?")

log_dir <- "non-stop-logs"
dir.create(log_dir, showWarnings = FALSE)
current_out <- file.path(log_dir, "current.Rout")
summary_log <- file.path(log_dir, "summary.log")

R_bin <- file.path(R.home("bin"), "R")

log_line <- function(...) {
    line <- paste0(...)
    cat(line, "\n")
    cat(line, "\n", file = summary_log, append = TRUE)
}

log_line("#### Started ", format(Sys.time()), " with ", R.version.string,
         " on ", Sys.info()[["nodename"]])

n_loops <- 0
n_failed <- 0

repeat {
    n_loops <- n_loops + 1
    start <- Sys.time()
    status <- system2(R_bin, c("--vanilla", "-f", "run-all-tests.R"),
                      stdout = current_out, stderr = current_out)
    minutes <- round(as.numeric(difftime(Sys.time(), start,
                                         units = "mins")), 1)
    ## A non-zero exit status is how run-all-tests.R signals a failure
    ## (it calls stop()). As a second check, the final message must be
    ## there too: this catches, e.g., the R process being killed.
    passed <- (status == 0) &&
        any(grepl("ALL INTERVENTION-SIDE TESTS PASSED",
                  readLines(current_out, warn = FALSE), fixed = TRUE))
    if (!passed) {
        n_failed <- n_failed + 1
        failed_out <- file.path(log_dir,
                                paste0("FAILED-loop-", n_loops, "-",
                                       format(start, "%Y-%m-%d_%H-%M-%S"),
                                       ".Rout"))
        file.copy(current_out, failed_out)
    }
    log_line("Loop ", n_loops, " started ", format(start), ", ", minutes,
             " min: ", if (passed) "PASSED" else "FAILED",
             if (!passed) paste0(" (exit status ", status, "; see ",
                                 failed_out, ")"),
             ". Totals: ", n_loops, " loops, ", n_loops - n_failed,
             " passed, ", n_failed, " failed.")
    if (!passed && stop_on_failure)
        stop("Stopping at the first failure (stop_on_failure is TRUE).")
}

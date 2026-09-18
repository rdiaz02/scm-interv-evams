## This file shows, via simple tests, a few features and equivalences
pwd <- getwd()
setwd("../")

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

source("generate_all_fitness_landscape.R")
setwd(pwd)

## Set a new random seed and print it. Used in the tests instead of
## set.seed(NULL), so that a failing run can be replayed: look for the
## last "Seed used was" printed before the failure, and replace the
## call to set_and_print_seed() there by set.seed(<that number>).
set_and_print_seed <- function() {
  ## First, set.seed(NULL): it re-seeds from the clock and the process
  ## ID. Without it, if an earlier line did, say, set.seed(1), the
  ## "random" seed drawn below would be the same number in every run.
  set.seed(NULL)
  ## sample.int gives a whole number between 1 and 1e9, a valid seed
  ## (the largest integer R allows is about 2.1e9).
  seed <- sample.int(1e9, 1)
  set.seed(seed)
  cat("\n Seed used was ", seed, "\n")
  return(invisible(seed))
}

### Generate fitness lanscapes using DAGs

#### About random numbers and replaying a failing run with landscapes
##
## generate_n_f_landscape_requir() draws all its random numbers inside
## mclapply(). Here we always call it with a first argument (n, the
## number of landscapes) of 1. With a single item, mclapply() does not
## fork: it just runs lapply() in the current R process. So all
## random numbers are drawn in the current process, and the seed set
## (and printed) by set_and_print_seed() is enough to replay a run.
##
## But if n were 2 or more, and more than one core were used,
## mclapply() would fork one child process per landscape. With R's
## default random number generator (Mersenne-Twister), each child
## throws away the seed it inherited and re-seeds itself from the
## clock and its process ID. So the landscapes would change from run
## to run, whatever seed we set here.
##
## Using RNGkind("L'Ecuyer-CMRG") would make those forked runs
## reproducible, but with care, because of a trap: forking does not advance
## the random state of the parent process. So, if using
## RNGkind("L'Ecuyer-CMRG"), two calls to mclapply() in a row, without
## drawing any random number in the parent between them, give each child
## the same random numbers as in the first call; e.g., two calls to
## generate_n_f_landscape_requir() would return the same landscapes.
## The non-advancing seed I think is not widely documented or known
## but has been mentioned here
## https://irudnyts.github.io//setting-a-seed-in-r-when-using-parallel-simulation/
## (also in R-bloggers: https://www.r-bloggers.com/2018/07/%F0%9F%8C%B1-setting-a-seed-in-r-when-using-parallel-simulation/)
## and is easy to check
## Expect differences between the mc.preschedule TRUE/FALSE
## set.seed(123, "L'Ecuyer")
## x22b <- unlist(mclapply(1:20,
##                         function(x) sum(runif(5)),
##                         mc.preschedule = FALSE,
##                         mc.cores = 3))
## x23b <- unlist(mclapply(1:20,
##                         function(x) sum(runif(5)),
##                         mc.preschedule = FALSE,
##                         mc.cores = 3))
## x24b <- unlist(mclapply(1:20,
##                         function(x) sum(runif(5)),
##                         mc.preschedule = TRUE,
##                         mc.cores = 3))
## x25b <- unlist(mclapply(1:20,
##                         function(x) sum(runif(5)),
##                         mc.preschedule = TRUE,
##                         mc.cores = 3))


#### The transition rate matrices based on the DAG and the fitness landscape are identical
##   Also, the scaled trm = trm * a * c  (now, a * c is = 1)

test_that("The transition rate matrices based on the DAG and the fitness landscape are identical", {
  set_and_print_seed()
  ## run a few times if you want
  for (i in 1:10) {
    cat("############### Doing i = ", i, "\n")
    rr_CBN <- generate_n_f_landscape_requir(1, 7, "CBN")
    ## Check the transformed trm = original trm * a * c
    ## where a = 0.006

    ## rr_CBN[[1]]$other$CBN_model
    ## rr_CBN[[1]]$c
    ## rr_CBN[[1]][["trm_scaled"]][1:5, 1:5]
    ## rr_CBN[[1]][["other"]][["CBN_trans_rate_mat"]][1:5, 1:5]
    ## rr_CBN[[1]][["trm_scaled"]][1:8,1:8]/rr_CBN[[1]][["other"]][["CBN_trans_rate_mat"]][1:8,1:8]

    a_c <- as.vector(rr_CBN[[1]][["trm_scaled"]]/
                       rr_CBN[[1]][["other"]][["CBN_trans_rate_mat"]])
    a_c <- as.vector(a_c[!is.na(a_c)])
    ## 0.006 is the value used to multiply lambda to get fitness
    ## c * 0.006 should be one
    expect_true(all.equal(a_c, rep(rr_CBN[[1]]$c * 0.006, length(a_c)), check.attributes = FALSE))

    expect_true(check_frequencies_of_genotypes_in_trm(rr_CBN[[1]][["trm_scaled"]], LETTERS[1:7])$test_OK)
    expect_true(check_frequencies_of_genotypes_in_trm(rr_CBN[[1]][["other"]][["CBN_trans_rate_mat"]], LETTERS[1:7])$test_OK)
    expect_equal(
      check_frequencies_of_genotypes_in_trm(rr_CBN[[1]][["trm_scaled"]], LETTERS[1:7])$genots_freq_from_check,
      check_frequencies_of_genotypes_in_trm(rr_CBN[[1]][["other"]][["CBN_trans_rate_mat"]], LETTERS[1:7])$genots_freq_from_check)

    ## Yes, expect this to take a bunch of iterations, as we require
    ## all three relationships to be present
    rr_HESBCN <- generate_n_f_landscape_requir(1, 7, "HESBCN",
                                               hesbcn_relations = c("AND", "OR", "XOR"))
    ## Check the transformed trm = original trm * a * c
    ## where a = 0.006

    ## rr_HESBCN[[1]]$other$CBN_model
    ## rr_HESBCN[[1]]$c
    ## rr_HESBCN[[1]][["trm_scaled"]][1:5, 1:5]
    ## rr_HESBCN[[1]][["other"]][["CBN_trans_rate_mat"]][1:5, 1:5]
    ## rr_HESBCN[[1]][["trm_scaled"]][1:8,1:8]/rr_HESBCN[[1]][["other"]][["CBN_trans_rate_mat"]][1:8,1:8]

    a_ch <- as.vector(rr_HESBCN[[1]][["trm_scaled"]]/
                        rr_HESBCN[[1]][["other"]][["HESBCN_trans_rate_mat"]])
    a_ch <- as.vector(a_ch[!is.na(a_ch)])
    ## 0.006 is the value used to multiply lambda to get fitness
    ## c * 0.006 should be one
    expect_true(all.equal(a_ch, rep(rr_HESBCN[[1]]$c * 0.006, length(a_ch)), check.attributes = FALSE))

    expect_true(check_frequencies_of_genotypes_in_trm(rr_HESBCN[[1]][["trm_scaled"]], LETTERS[1:7])$test_OK)
    expect_true(check_frequencies_of_genotypes_in_trm(rr_HESBCN[[1]][["other"]][["HESBCN_trans_rate_mat"]], LETTERS[1:7])$test_OK)

    expect_equal(
      check_frequencies_of_genotypes_in_trm(rr_HESBCN[[1]][["trm_scaled"]], LETTERS[1:7])$genots_freq_from_check,
      check_frequencies_of_genotypes_in_trm(rr_HESBCN[[1]][["other"]][["HESBCN_trans_rate_mat"]], LETTERS[1:7])$genots_freq_from_check)

    rr_HESBCN_2 <- generate_n_f_landscape_requir(1, 7, "HESBCN",
                                                 hesbcn_relations = c("AND", "OR"))
    ## Check the transformed trm = original trm * a * c
    ## where a = 0.006

    ## rr_HESBCN_2[[1]]$other$CBN_model
    ## rr_HESBCN_2[[1]]$c
    ## rr_HESBCN_2[[1]][["trm_scaled"]][1:5, 1:5]
    ## rr_HESBCN_2[[1]][["other"]][["CBN_trans_rate_mat"]][1:5, 1:5]
    ## rr_HESBCN_2[[1]][["trm_scaled"]][1:8,1:8]/rr_HESBCN_2[[1]][["other"]][["CBN_trans_rate_mat"]][1:8,1:8]

    a_ch2 <- as.vector(rr_HESBCN_2[[1]][["trm_scaled"]]/
                         rr_HESBCN_2[[1]][["other"]][["HESBCN_trans_rate_mat"]])
    a_ch2 <- as.vector(a_ch2[!is.na(a_ch2)])
    ## 0.006 is the value used to multiply lambda to get fitness
    ## c * 0.006 should be one
    expect_true(all.equal(a_ch, rep(rr_HESBCN_2[[1]]$c * 0.006, length(a_ch)), check.attributes = FALSE))

    expect_true(check_frequencies_of_genotypes_in_trm(rr_HESBCN_2[[1]][["trm_scaled"]], LETTERS[1:7])$test_OK)
    expect_true(check_frequencies_of_genotypes_in_trm(rr_HESBCN_2[[1]][["other"]][["HESBCN_trans_rate_mat"]], LETTERS[1:7])$test_OK)
    expect_equal(
      check_frequencies_of_genotypes_in_trm(rr_HESBCN_2[[1]][["trm_scaled"]], LETTERS[1:7])$genots_freq_from_check,
      check_frequencies_of_genotypes_in_trm(rr_HESBCN_2[[1]][["other"]][["HESBCN_trans_rate_mat"]], LETTERS[1:7])$genots_freq_from_check)
  }
})

set.seed(NULL)

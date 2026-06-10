## Testing: this code has been tested in evamtools.  I add just one manual
## test and a bunch of random tests of the relative fitness difference by
## comparing with the transition matrix.

library(testthat)

pwd <- getwd()
setwd("../")
source("rfitness_to_trm.R")
setwd(pwd)

set.seed(NULL)

## show all columns
options(sparse.colnames = TRUE)


test_that("no_evam_genots_2_fgraph_and_trans_mat_rf works OK", {

  ## relative diff. vector w.r.t x, with min set to 0
  rel_f_diff <- function(y, x) {
    tmp <- (y - x)/x
    pmax(tmp, 0)
  }
  ## Yes, this seed is set on purpose, as the test below
  ## is for this specific fitness landscape.
  ## Other fitness landscapes could have very different sets
  ## of accessible genotypes
  ## Note test below does not fix seed.
  set.seed(1)
  u <- rfitness(4)
  ul <- no_evam_rfitness_to_letter(u)
  rul <- no_evam_genots_2_fgraph_and_trans_mat_rf(ul)
  expect_equal(rul$relative_fitness_differences[1, 2:5],
               rel_f_diff(rul$accessible_genotypes, 1)[1:4])
  expect_equal(rul$relative_fitness_differences[2, 6:7],
               rel_f_diff(rul$accessible_genotypes,
                          rul$accessible_genotypes["A"])[5:6])
  expect_equal(rul$relative_fitness_differences[3, 8:9],
               rel_f_diff(rul$accessible_genotypes,
                          rul$accessible_genotypes["B"])[7:8])
  expect_equal(rul$relative_fitness_differences[4, 10],
               unname(rel_f_diff(rul$accessible_genotypes,
                                 rul$accessible_genotypes["D"])[9]))
  expect_equal(rul$relative_fitness_differences[5, 11],
               unname(rel_f_diff(rul$accessible_genotypes,
                                 rul$accessible_genotypes["A, C"])[10]))
  expect_equal(rul$relative_fitness_differences[6, 12],
               unname(rel_f_diff(rul$accessible_genotypes,
                                 rul$accessible_genotypes["B, C"])[11]))
  expect_equal(rul$relative_fitness_differences[7, 12],
               unname(rel_f_diff(rul$accessible_genotypes,
                                 rul$accessible_genotypes["B, D"])[11]))
  expect_equal(rul$relative_fitness_differences[8, 12],
               unname(rel_f_diff(rul$accessible_genotypes,
                                 rul$accessible_genotypes["C, D"])[11]))

  m2 <- matrix(0L, nrow = 8, ncol = 12)
  colnames(m2) <- c("WT", names(rul$accessible_genotypes))
  rownames(m2) <- colnames(m2)[-c(4, 7, 11, 12)]
  m2[1, 2:5] <- 1L
  m2[2, 6:7] <- 1L
  m2[3, 8:9] <- 1L
  m2[4, 10]  <- 1L
  m2[5, 11]  <- 1L
  m2[6, 12]  <- 1L
  m2[7, 12]  <- 1L
  m2[8, 12]  <- 1L
  expect_equal(as.matrix(rul$fitness_graph), m2)
  expect_true(all(as.matrix(rul$fitness_graph - m2) == 0L))


  tm_from_rel_fit_diff <- function(x) {
    ## this sweep puts 0s in the sparse matrix
    if (!is.null(dim(x)))
      sweep(x, 1, rowSums(x), "/")
    else ## so a vector
      x/sum(x)
  }

  set.seed(NULL)
  ## Tested with much larger i and nn up to 9.
  for (i in 1:20) {
    nn <- sample(3:6, size = 1)
    f1 <- rfitness(nn)
    f1l <- no_evam_rfitness_to_letter(f1)
    rf1l <- no_evam_genots_2_fgraph_and_trans_mat_rf(f1l)
    if (length(OncoSimulR:::wrap_accessibleGenotypes(f1, 0)) > 1)
      expect_equal(as.matrix(tm_from_rel_fit_diff(rf1l$relative_fitness_differences)),
                   as.matrix(rf1l$transition_matrix))
  }


})


test_that("no_evam_genots_2_fgraph_and_trans_mat_rf works with WT fitness different from 1", {
  rel_f_diff <- function(y, x) {
    tmp <- (y - x)/x
    pmax(tmp, 0)
  }

  set.seed(1)
  u <- rfitness(4, model = "RMF", wt_is_1 = "no", sd=0.5)
  fitness_birth_column <- ifelse("Fitness" %in% colnames(u),
                                 "Fitness", "Birth") # Compatible with prevoius and new versions of OncoSimulR
  ul <- no_evam_rfitness_to_letter(u)
  rul <- no_evam_genots_2_fgraph_and_trans_mat_rf(ul)
  expect_equal(rul$relative_fitness_differences[1, 2:5],
               rel_f_diff(rul$accessible_genotypes, u[1, fitness_birth_column])[1:4])
  expect_equal(rul$relative_fitness_differences[2,6],
               unname(rel_f_diff(rul$accessible_genotypes,
                                 rul$accessible_genotypes["A"])[5]))
  expect_equal(rul$relative_fitness_differences[3, 7],
               unname(rel_f_diff(rul$accessible_genotypes,
                                 rul$accessible_genotypes["B"])[6]))
  expect_equal(rul$relative_fitness_differences[4, 8:9],
               rel_f_diff(rul$accessible_genotypes,
                          rul$accessible_genotypes["D"])[7:8])
  expect_equal(rul$relative_fitness_differences[5, 10],
               unname(rel_f_diff(rul$accessible_genotypes,
                                 rul$accessible_genotypes["A, C"])[9]))
  expect_equal(rul$relative_fitness_differences[6, 11],
               unname(rel_f_diff(rul$accessible_genotypes,
                                 rul$accessible_genotypes["B, D"])[10]))
  expect_equal(rul$relative_fitness_differences[7, 11],
               unname(rel_f_diff(rul$accessible_genotypes,
                                 rul$accessible_genotypes["C, D"])[10]))

  m2 <- matrix(0L, nrow = 7, ncol = 11)
  colnames(m2) <- c("WT", names(rul$accessible_genotypes))
  rownames(m2) <- colnames(m2)[-c(4, 7, 10, 11)]
  m2[1, 2:5] <- 1L
  m2[2, 6] <- 1L
  m2[3, 7] <- 1L
  m2[4, 8:9]  <- 1L
  m2[5, 10]  <- 1L
  m2[6, 11]  <- 1L
  m2[7, 11]  <- 1L
  expect_equal(as.matrix(rul$fitness_graph), m2)
  expect_true(all(as.matrix(rul$fitness_graph - m2) == 0L))


  tm_from_rel_fit_diff <- function(x) {
    ## this sweep puts 0s in the sparse matrix
    if (!is.null(dim(x)))
      sweep(x, 1, rowSums(x), "/")
    else ## so a vector
      x/sum(x)
  }

  set.seed(NULL)
  ## Tested with much larger i(i.e 2000) and nn up to 9.
  for (i in 1:20) {
    nn <- sample(3:6, size = 1)
    f1 <- rfitness(nn, wt_is_1 = "no")
    f1l <- no_evam_rfitness_to_letter(f1)
    rf1l <- no_evam_genots_2_fgraph_and_trans_mat_rf(f1l)
    if (length(OncoSimulR:::wrap_accessibleGenotypes(f1, 0)) > 1)
      expect_equal(as.matrix(tm_from_rel_fit_diff(rf1l$relative_fitness_differences)),
                   as.matrix(rf1l$transition_matrix))
  }
})

test_that("Stops if any fitness is 0 or less than 0",{
  set.seed(1)
  fit <- rfitness(4, truncate_at_0 = FALSE)
  fit_letter <- no_evam_rfitness_to_letter(fit)
  expect_error(no_evam_genots_2_fgraph_and_trans_mat_rf(fit_letter),
               "All fitness values must be above 0.")

  fit_letter[fit_letter < 0] <- 0 # If no negative, only 0
  expect_error(no_evam_genots_2_fgraph_and_trans_mat_rf(fit_letter),
               "All fitness values must be above 0.")
  set.seed(1)
  fit_2 <- rfitness(4)
  fit_letter_2 <- no_evam_rfitness_to_letter(fit_2)
  expect_error(no_evam_genots_2_fgraph_and_trans_mat_rf(fit_letter_2), NA)

})


set.seed(NULL)

## options(intervention_every_gene_cores = parallel::detectCores())
library(testthat)
pwd <- getwd()
setwd("../")
source("intervention.R")
setwd(pwd)

set.seed(NULL)


test_that("genots_to_bin transforms genotypes into binary correctly",{
    str_genots <- c("WT", "A", "A, B", "B", "A, C", "A, B, C")
    result <- matrix(
        c(0,0,0, # WT
          1,0,0, # A
          1,1,0, # A,B
          0,1,0, # B
          1,0,1, #A, C
          1,1,1 # A, B, C
          ),
        ncol = 3,
        byrow = TRUE
    )
    colnames(result) <- c("A", "B", "C")
    expect_equal(genots_to_bin(str_genots, 3), result)
    res2 <- matrix(c(0, 1), nrow = 2)
    colnames(res2) <- c("A")
    expect_equal(genots_to_bin(c("WT", "A"), 1), res2)

})

test_that("get_square_matrix returns a square matrix", {
  matrix_2x3 <- get_square_matrix(matrix(runif(6), nrow = 2,
                       dimnames = list(c("A", "B"), c("A", "B", "C"))))
  matrix_2x3_sparse <- get_square_matrix(Matrix(runif(6), nrow = 2,
                       dimnames = list(c("A", "B"), c("A", "B", "C")), sparse = TRUE))
  expect_equal(nrow(matrix_2x3), ncol(matrix_2x3))
  expect_equal(nrow(matrix_2x3), 3)
  expect_equal(nrow(matrix_2x3_sparse), ncol(matrix_2x3_sparse))
  expect_equal(nrow(matrix_2x3_sparse), 3)
})

test_that("get_square_matrix returns expected ouput", {
  mtx <- Matrix(c(1,2,3,4,5,6), nrow = 2, byrow = TRUE,
                dimnames = list(c("A", "B"), c("A", "B", "C")))
  expect_equal(get_square_matrix(mtx, 0), Matrix(c(1,2,3,4,5,6,0,0,0), nrow = 3,
                                                 byrow = TRUE, sparse = TRUE,
                                                 dimnames = list(c("A", "B", "C"), c("A", "B", "C"))))
  mtx2 <- Matrix(c(1,2,3,4,5,6), nrow = 2, byrow = TRUE,
                 dimnames = list(c("A", "D"), c("A", "B", "C")))
  expect_equal(get_square_matrix(mtx2, 0), Matrix(c(1,2,3,0,0,0,0,0,0,0,0,0,4,5,6,0), nrow = 4,
                                                  byrow = TRUE, sparse = TRUE,
                                                  dimnames = list(c("A", "B", "C", "D"),
                                                                  c("A", "B", "C", "D"))))
})

set.seed(NULL)

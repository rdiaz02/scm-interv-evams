## Copyright 2026 Ramon Diaz-Uriarte

## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU Affero General Public License (AGPLv3.0) as
## published by the Free Software Foundation, either version 3 of the
## License, or (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU Affero General Public License for more details.

## You should have received a copy of the GNU Affero General Public License
## along with this program. If not, see <http://www.gnu.org/licenses/>.


### What this is

## The question that this file answers (2026-09-21).
##
## When we intervene on gene g, the code modifies the model first (the
## DAG of restrictions, or the theta matrix for MHN) and then builds the
## transition matrix from the modified model. We never take the
## transition matrix of the model without intervention and cut pieces
## out of it.
##
## So: would we get the same transition matrix if we did it the other
## way round? That is, if we took the transition matrix of the model
## without intervention, kept only the genotypes that do not have gene
## g, and then divided each row by its sum (renormalized, so that each
## row adds up to 1 again)?
##
## The answer is yes, and that is what this file checks, for OT,
## OncoBN, CBN, H-ESBCN and MHN.

#### Why it is so?

## In the weighted fitness graph (function cpm2tm in evamtools), the
## weight of the arrow from genotype x to genotype x + h depends only on
## the gene that is gained, h: the conditional probability for OT and
## OncoBN, the lambda for CBN and H-ESBCN. It does not depend on the
## gene we kill, g. (For MHN the rate depends on the genes already
## present in x, but x never has g, so again it does not depend on g.)
##
## So intervening on g does two things, and only two: it removes the
## genotypes that have g (and the genes that can only appear together
## with g, whose genotypes also have g), and it leaves the weights of
## all the surviving arrows untouched. The transition matrix is then
## obtained by dividing each row by its sum.
##
## Going the other way round we divide exactly the same surviving
## weights by exactly the same sum. So we get the same matrix. This
## file checks that this is indeed so, and that it is so for the five
## models.

#### Why we care, and what this is not.

## For OT and OncoBN there is no transition rate matrix: these models are
## untimed. We still obtain hitting probabilities for them, using the "as
## if" transition matrix between genotypes built by cpm2tm from the DAG
## (this is the heuristic procedure of Diaz-Colunga and Diaz-Uriarte, 2021,
## based on the same one use in Diaz-Uriarte and Vasallo, 2019). So it is
## worth knowing what that matrix is, after an intervention, in terms of
## the matrix before the intervention.
##

##### How is this related to the rescaling of genotype freqs. and "conditioning ain't intervention"?

## This is not the rescaling of genotype frequencies that the paper says is
## wrong (the OT example figure: setting to 0 the frequencies of the
## genotypes that have the killed gene, and rescaling the rest so they add
## up to 1; that is conditioning, not intervening). Here we renormalize the
## rows of a transition matrix, which is a different thing: for each
## genotype, the probability of the step that is now blocked is shared out
## among the steps that are still possible. And, for OT and OncoBN, the
## predicted genotype frequencies do not come from this matrix anyway: they
## come from the conditional probabilities of the DAG (function
## OT_model_2_predict_genots in evamtools). Only the hitting probabilities
## use this matrix.

#### What about HyperHMM? Didn't you do a totally different thing there?

## HyperHMM is the one case where the code does remove rows and columns
## from a transition matrix (function kill_gene_HyperHMM, and then
## kill_gene_HyperHMM_drop_unreachable, in file
## kill-gene-and-output-from-cpm.R), and it does not renormalize the
## rows: the probability of the blocked moves is put in the diagonal of
## the genotypes that survive.
##
## Why the difference? Because of what is primary in each model. In
## HyperHMM there are no parameters: a step is an attempt at mutation, the
## entries of $R$ are probabilities per step, killing sets one of them to 0
## and leaves the others alone, and the remainder stays on the diagonal
## because the row must sum to 1. In OT and OncoBN there are parameters
## (the pi of OT, for example), so we set to zero the affected one, and
## then compute the rest much like we do for CBN.

## Yes, modularity, in both OT/OncoBN and HyperHMM is not as tenable as for
## CBN/H-ESBCN. But, for predictions under interventions, the most sensible
## procedures, without invoking extra assumptions, are to do what we do for
## OT/OncoBN on the one hand, and HyperHMM on the other (i.e., actually do
## different things).


### Thread control

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

library(testthat)
pwd <- getwd()
setwd("../")
source("intervention.R")
setwd(pwd)


### Utility functions

## Transition matrix, names of the genotypes to keep -> that matrix with
## only those genotypes, and each row divided by its sum. Rows that add
## up to 0 (genotypes from which nothing else can be reached) are left
## as they are, all zeros.
restrict_and_renormalize <- function(trans_mat, keep) {
  m <- as.matrix(trans_mat)[keep, keep, drop = FALSE]
  row_sums <- rowSums(m)
  m[row_sums > 0, ] <- m[row_sums > 0, ] / row_sums[row_sums > 0]
  return(m)
}

## Names of the genotypes that do not have the gene. The genotype names
## are the genes separated by ", " (e.g., "A, B, D"), and "WT" for the
## wildtype, so we split each name and look for an exact match.
genots_without_gene <- function(genot_names, gene) {
  has_gene <- vapply(strsplit(genot_names, ", ", fixed = TRUE),
                     function(z) gene %in% z,
                     logical(1))
  return(genot_names[!has_gene])
}

## Model (a data frame for all models except MHN; the log-theta matrix
## for MHN), method, gene to kill -> nothing; it runs the checks.
##
## The two matrices that we compare are:
##   - the one we get the usual way: kill the gene in the model, and
##     build the transition matrix from the killed model;
##   - the one we get the other way round: take the transition matrix of
##     the model without intervention, keep only the genotypes without
##     the gene, and divide each row by its sum.
check_kill_equals_restrict_renormalize <- function(model, method, gene,
                                                   label = "") {
  trans_mat_name <- paste0(method, "_trans_mat")
  mat_before <- as.matrix(
      suppressWarnings(get_full_output(model))[[trans_mat_name]])
  killed_model <- suppressWarnings(kill_gene(model, gene))
  mat_after <- as.matrix(
      suppressWarnings(get_full_output(killed_model))[[trans_mat_name]])
  keep <- rownames(mat_after)

  what <- paste0(method, " ", label, " kill ", gene)
  ## The genotypes that are left must be exactly the genotypes of the
  ## original matrix that do not have the killed gene
  expect_true(setequal(keep, genots_without_gene(rownames(mat_before), gene)),
              label = paste0(what, ": genotypes that are left"))
  ## And the two matrices must be the same
  expect_true(all.equal(restrict_and_renormalize(mat_before, keep),
                        mat_after),
              label = paste0(what, ": matrices"))
  return(invisible(NULL))
}

## Gene names of a model given as a data frame with From and To columns
genes_of_model <- function(model) {
  return(sort(setdiff(unique(c(model[, "From"], model[, "To"])), "Root")))
}


### Hand-made models
##  One per model type, so that the checks below do not depend on a
##  random model, and so that killing any single gene still leaves
##  something to compare.

##  OT: Root -> A, Root -> B, A -> C, B -> D
ot_model <- data.frame(From = c("Root", "Root", "A", "B"),
                       To   = c("A",    "B",    "C", "D"),
                       OT_edgeWeight = c(0.6, 0.4, 0.7, 0.3))

##  OncoBN, DBN version (the relation is OR): Root -> A, Root -> B,
##  A -> C, B -> C, C -> D
oncobn_dbn_model <- data.frame(
    From = c("Root", "Root", "A",  "B",  "C"),
    To   = c("A",    "B",    "C",  "C",  "D"),
    theta = c(0.6, 0.5, 0.4, 0.4, 0.7),
    Relation = c("Single", "Single", "OR", "OR", "Single"))

##  OncoBN, CBN version (the relation is AND)
oncobn_cbn_model <- data.frame(
    From = c("Root", "Root", "A",   "B",   "C"),
    To   = c("A",    "B",    "C",   "C",   "D"),
    theta = c(0.6, 0.5, 0.4, 0.4, 0.7),
    Relation = c("Single", "Single", "AND", "AND", "Single"))

##  CBN: Root -> A, Root -> B, A -> C, B -> C, C -> D
cbn_model <- data.frame(From = c("Root", "Root", "A", "B", "C"),
                        To   = c("A",    "B",    "C", "C", "D"),
                        rerun_lambda = c(1, 2, 3, 3, 4))

##  H-ESBCN, with a mixture of AND, OR and XOR
hesbcn_model <- data.frame(
    From = c("Root", "Root", "Root", "A",   "B",   "C",   "D"),
    To   = c("A",    "B",    "C",    "D",   "D",   "E",   "E"),
    Lambdas = c(1, 2, 3, 4, 4, 5, 5),
    Relation = c("Single", "Single", "Single", "OR", "OR", "XOR", "XOR"))

##  MHN: the log-theta matrix of a 4-gene model. The diagonal holds the
##  baseline log-rates; theta[i, j] is the effect of gene j on gene i.
mhn_theta <- matrix(
    c(-2.0,  0.5, -0.3,  0.0,
       0.3, -1.5,  0.4, -0.2,
      -0.1,  0.2, -1.0,  0.3,
       0.0, -0.3,  0.1, -3.0),
    nrow = 4, ncol = 4, byrow = TRUE,
    dimnames = list(c("A", "B", "C", "D"), c("A", "B", "C", "D")))


### The checks with the hand-made models

test_that("OT: killing in the DAG is restricting and renormalizing", {
  local_edition(3)
  cat("\n OT: hand-made model\n")
  for (gene in genes_of_model(ot_model)) {
    check_kill_equals_restrict_renormalize(ot_model, "OT", gene,
                                           label = "hand-made")
  }
})

test_that("OncoBN, DBN: killing in the DAG is restricting and renormalizing", {
  local_edition(3)
  cat("\n OncoBN, DBN version: hand-made model\n")
  for (gene in genes_of_model(oncobn_dbn_model)) {
    check_kill_equals_restrict_renormalize(oncobn_dbn_model, "OncoBN", gene,
                                           label = "hand-made, DBN")
  }
})

test_that("OncoBN, CBN: killing in the DAG is restricting and renormalizing", {
  local_edition(3)
  cat("\n OncoBN, CBN version: hand-made model\n")
  for (gene in genes_of_model(oncobn_cbn_model)) {
    check_kill_equals_restrict_renormalize(oncobn_cbn_model, "OncoBN", gene,
                                           label = "hand-made, CBN")
  }
})

test_that("CBN: killing in the DAG is restricting and renormalizing", {
  local_edition(3)
  cat("\n CBN: hand-made model\n")
  for (gene in genes_of_model(cbn_model)) {
    check_kill_equals_restrict_renormalize(cbn_model, "CBN", gene,
                                           label = "hand-made")
  }
})

test_that("H-ESBCN: killing in the DAG is restricting and renormalizing", {
  local_edition(3)
  cat("\n H-ESBCN: hand-made model\n")
  for (gene in genes_of_model(hesbcn_model)) {
    check_kill_equals_restrict_renormalize(hesbcn_model, "HESBCN", gene,
                                           label = "hand-made")
  }
})

test_that("MHN: killing in theta is restricting and renormalizing", {
  local_edition(3)
  cat("\n MHN: hand-made model\n")
  for (gene in colnames(mhn_theta)) {
    check_kill_equals_restrict_renormalize(mhn_theta, "MHN", gene,
                                           label = "hand-made")
  }
})


### The same checks with random models
##  The hand-made models above are small and have a fixed shape. Here we
##  use random models, to see the same thing with many different shapes.
##  Killing a gene can leave a model with no rows at all (e.g., if every
##  other gene depends on the killed one); there is then no matrix to
##  compare, and we skip that gene.

test_that("Random models: killing is restricting and renormalizing", {
  local_edition(3)
  set_and_print_seed()
  total_iters <- 5
  for (i in 1:total_iters) {
    cat("\n #################### Doing iteration ", i, "\n\n")
    random_models <- list(
        OT = evamtools::random_evam(5, model = "OT",
                                    ot_oncobn_epos = 0)$OT_model,
        OncoBN_DBN = evamtools::random_evam(5, model = "OncoBN",
                                            ot_oncobn_epos = 0,
                                            oncobn_model = "DBN")$OncoBN_model,
        OncoBN_CBN = evamtools::random_evam(5, model = "OncoBN",
                                            ot_oncobn_epos = 0,
                                            oncobn_model = "CBN")$OncoBN_model,
        CBN = evamtools::random_evam(5, model = "CBN")$CBN_model,
        HESBCN = evamtools::random_evam(5, model = "HESBCN")$HESBCN_model)
    methods <- c(OT = "OT", OncoBN_DBN = "OncoBN", OncoBN_CBN = "OncoBN",
                 CBN = "CBN", HESBCN = "HESBCN")
    for (nm in names(random_models)) {
      model <- random_models[[nm]]
      for (gene in genes_of_model(model)) {
        ## Nothing left to compare if the killed model has no rows
        if (nrow(suppressWarnings(kill_gene(model, gene))) == 0) next
        check_kill_equals_restrict_renormalize(model, methods[[nm]], gene,
                                               label = paste0("random ", nm))
      }
    }
    ## MHN: a random log-theta matrix of 5 genes
    random_theta <- matrix(runif(25, -3, 3), ncol = 5)
    colnames(random_theta) <- rownames(random_theta) <- LETTERS[1:5]
    for (gene in colnames(random_theta)) {
      check_kill_equals_restrict_renormalize(random_theta, "MHN", gene,
                                             label = "random MHN")
    }
  }
})

set.seed(NULL)

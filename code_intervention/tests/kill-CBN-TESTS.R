## Copyright 2022 Ramon Diaz-Uriarte

## This program is free software: you can redistribute it and/or modify it under
## the terms of the GNU Affero General Public License (AGPLv3.0) as published by
## the Free Software Foundation, either version 3 of the License, or (at your
## option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU Affero General Public License for more details.

## You should have received a copy of the GNU Affero General Public License along
## with this program.  If not, see <http://www.gnu.org/licenses/>.


### What this is

## First, a simple and quick test. Then, thorough tests of CBN
## interventions using six larger hand-crafted DAG structures
## (DAG_2 through DAG_7; same structures as in kill-HESBCN-DAG-2-7-TESTS.R).
##
## CBN has only conjunctive (AND) dependencies, so each DAG gives exactly
## one model — no relation combinations to enumerate.
## Total: 6 CBN models, one per DAG.
##
## For each model, four procedures are compared:
##   1. Standard kill_gene (modifies DAG edges)
##   2. kill_gene_by_params_to_0 (zeroes parameters)
##   3. intervene_cpm_trm_rm_every_gene (removes genotypes from pre-computed TRM;
##      does NOT call get_full_output for the intervention step)
##   4. intervene_fitness_landscape_every_gene (fitness landscape derived from model)
##
##   Having shown all 4 are identical, we then check
##   (Full ground truth sections)  that manual
##   structure removal is identical to the first one.
##
## Plus: procedure 1 with a row-permuted model must equal procedure 1 with
##       the original model order.
##
## Procedures 3 and 4 share no key code with procedures 1 and 2
## (which go through get_full_output / get_genotype_freqs_cpm).
## This makes these tests a genuine independent check.
##
## Gene names in each DAG are sequential from A (sorted alphabetically from
## the original names used when designing the DAG topology). A is intentionally
## NOT the first structural gene or a root gene in most DAGs, so gene-name
## ordering is not trivially sequential.
##


## options(intervention_every_gene_cores = parallel::detectCores())
library(testthat)
pwd <- getwd()
setwd("../")
source("intervention.R")
setwd(pwd)

set.seed(NULL)

## Extract a single named intervention from the full intervention output as a
## list(genot_freqs, hitting_probs_from_WT), for direct comparison with
## preds_from_model / preds_from_hesbcn_model.
get_interv <- function(res, name) {
    return(res[[name]])
}

test_that("Simple CBN checks of structure (standard) killing and equiv of procedures", {
  local_edition(3)
  m1 <- data.frame(From = c("Root", "Root", "B", "B", "C", "C"),
                   To   = c("A",    "B",    "C", "D",  "E", "F"),
                   rerun_lambda = 1:6)

  m2 <- data.frame(From = c("Root", "Root", "A", "B", "B", "C", "C"),
                   To   = c("A",    "B",    "C", "C", "D",  "E", "F"),
                   rerun_lambda = c(1, 2, 3, 3, 4, 5, 6))

  m3 <- data.frame(From = c("Root", "Root", "A", "B", "B", "C", "C"),
                   To   = c("A",    "B",    "E", "C", "D",  "E", "F"),
                   rerun_lambda = c(1, 2, 3, 4, 5, 3, 6) )

  kills_m1 <- lapply(LETTERS[1:6], function(v) kill_gene(m1, v))
  kills_m2 <- lapply(LETTERS[1:6], function(v) kill_gene(m2, v))
  kills_m3 <- lapply(LETTERS[1:6], function(v) kill_gene(m3, v))

  names(kills_m1)  <- names(kills_m2) <-
    names(kills_m3) <- LETTERS[1:6]

  m1_A <- m1[-c(1), ]
  m1_B <- m1[1, ]
  m1_C <- m1[-c(3, 5, 6), ]
  m1_D <- m1[-c(4), ]
  m1_E <- m1[-c(5), ]
  m1_F <- m1[-c(6), ]

  expect_equal(kills_m1[["A"]], m1_A)
  expect_equal(kills_m1[["B"]], m1_B)
  expect_equal(kills_m1[["C"]], m1_C)
  expect_equal(kills_m1[["D"]], m1_D)
  expect_equal(kills_m1[["E"]], m1_E)
  expect_equal(kills_m1[["F"]], m1_F)

  m2_A <- m2[-c(1, 3, 4, 6, 7), ]
  m2_B <- m2[c(1), ]
  m2_C <- m2[-c(3, 4, 6, 7), ]
  m2_D <- m2[-c(5), ]
  m2_E <- m2[-c(6), ]
  m2_F <- m2[-c(7), ]

  expect_equal(kills_m2[["A"]], m2_A)
  expect_equal(kills_m2[["B"]], m2_B)
  expect_equal(kills_m2[["C"]], m2_C)
  expect_equal(kills_m2[["D"]], m2_D)
  expect_equal(kills_m2[["E"]], m2_E)
  expect_equal(kills_m2[["F"]], m2_F)

  m3_A <- m3[-c(1, 3, 6), ]
  m3_B <- m3[1, ]
  m3_C <- m3[-c(3, 4, 6, 7), ]
  m3_D <- m3[-c(5), ]
  m3_E <- m3[-c(3, 6), ]
  m3_F <- m3[-c(7), ]

  expect_equal(kills_m3[["A"]], m3_A)
  expect_equal(kills_m3[["B"]], m3_B)
  expect_equal(kills_m3[["C"]], m3_C)
  expect_equal(kills_m3[["D"]], m3_D)
  expect_equal(kills_m3[["E"]], m3_E)
  expect_equal(kills_m3[["F"]], m3_F)

  preds_from_model <- function(x) {
    if (nrow(x) == 0) {
      return(list(genot_freqs = c(WT = 1), hitting_probs_from_WT = c(WT = 1.0)))
    }
    tmp <- suppressWarnings(get_full_output(x))
    f <- tmp$CBN_predicted_genotype_freqs
    hp <- tmp$CBN_hitting_probs_from_WT
    list(genot_freqs = f[f > 0], hitting_probs_from_WT = filter_hp_keep_wt(hp))
  }

  m1_int <- intervene_cpm_every_gene(list(CBN_model = m1), "CBN")
  m1_int_p0 <- intervene_cpm_every_gene(list(CBN_model = m1), "CBN",
                                        kill_gene_funct = kill_gene_by_params_to_0)
  m1_int_trm <- intervene_cpm_trm_rm_every_gene(c(list(CBN_model = m1),
                                                  get_full_output(m1)), "CBN")


  m2_int <- intervene_cpm_every_gene(list(CBN_model = m2), "CBN")
  m2_int_p0 <- intervene_cpm_every_gene(list(CBN_model = m2), "CBN",
                                        kill_gene_funct = kill_gene_by_params_to_0)
  m2_int_trm <- intervene_cpm_trm_rm_every_gene(c(list(CBN_model = m2),
                                                  get_full_output(m2)), "CBN")


  m3_int <- intervene_cpm_every_gene(list(CBN_model = m3), "CBN")
  m3_int_p0 <- intervene_cpm_every_gene(list(CBN_model = m3), "CBN",
                                        kill_gene_funct = kill_gene_by_params_to_0)
  m3_int_trm <- intervene_cpm_trm_rm_every_gene(c(list(CBN_model = m3),
                                                  get_full_output(m3)), "CBN")


  expect_equal(m1_int, m1_int_p0)
  expect_equal(m2_int, m2_int_p0)
  expect_equal(m3_int, m3_int_p0)
  expect_equal(m1_int, m1_int_trm)
  expect_equal(m2_int, m2_int_trm)
  expect_equal(m3_int, m3_int_trm)

  gene_names <- sort(setdiff(union(m1$From, m1$To), "Root"))
  for (g in gene_names) {
    message("Doing gene ", g)
    expect_equal(get_interv(m1_int, paste0("I:", g)),
                 preds_from_model(get(paste0("m1_", g))))
  }

  gene_names <- sort(setdiff(union(m2$From, m2$To), "Root"))
  for (g in gene_names) {
    message("Doing gene ", g)
    expect_equal(get_interv(m2_int, paste0("I:", g)),
                 preds_from_model(get(paste0("m2_", g))))
  }

  gene_names <- sort(setdiff(union(m3$From, m3$To), "Root"))
  for (g in gene_names) {
    message("Doing gene ", g)
    expect_equal(get_interv(m3_int, paste0("I:", g)),
                 preds_from_model(get(paste0("m3_", g))))
  }

  ## Equivalence to fitness landscape killing
  ## First, utility to create the landscape from the bare model
  ## Based on generate_n_f_landscape_requir.
  cbn_to_landscape_obj <- function(model, n_genes) {
    a_cpm_2_si <- 0.006
    genot_fitness <- ev2_cpm_to_fitness_genots(model, a = a_cpm_2_si)
    ## Newer OncoSimulR uses "Birth"; older uses "Fitness"
    fitness_col <- ifelse("Birth" %in% colnames(genot_fitness), "Birth", "Fitness")
    rfo <- cbind(genots_to_bin(genot_fitness$Genotype, n_genes),
                 Fitness = genot_fitness[[fitness_col]])
    rfo <- rfo[rfo[, "Fitness"] > 0, ]
    class(rfo) <- c("matrix", "array")
    trm_and_c <- suppressMessages(get_scaled_trm_adaptive(rfo, c = 1/a_cpm_2_si))
    list(fitness_landscape = rfo,
         c = trm_and_c$c,
         trm_scaled = trm_and_c$trm_scaled)
  }

  m1_fl <- cbn_to_landscape_obj(m1, 6)
  m2_fl <- cbn_to_landscape_obj(m2, 6)
  m3_fl <- cbn_to_landscape_obj(m3, 6)
  m1_fl_i <- intervene_fitness_landscape_every_gene(m1_fl)
  m2_fl_i <- intervene_fitness_landscape_every_gene(m2_fl)
  m3_fl_i <- intervene_fitness_landscape_every_gene(m3_fl)

  expect_equal(m1_int, m1_fl_i)
  expect_equal(m2_int, m2_fl_i)
  expect_equal(m3_int, m3_fl_i)

})


### Comprehensive tests


## Same DAG structures as in kill-HESBCN-TESTS.R.
## CBN models use rerun_lambda; no Relation column.
##
## Original gene names mapped to sequential letters (alphabetical sort -> A,B,...):
##   DAG_2: A->A, B->B, C->C, D->D, F->E, G->F, M->G, S->H, W->I, Y->J
##   DAG_3: A->A, B->B, G->C, K->D, M->E, N->F, T->G, W->H, X->I
##   DAG_4: G->A, H->B, S->C, T->D, V->E, W->F, X->G, Y->H, Z->I
##   DAG_5: C->A, D->B, F->C, G->D, I->E, K->F, N->G, T->H, Z->I
##   DAG_6: B->A, C->B, D->C, E->D, F->E, G->F, H->G
##   DAG_7: B->A, C->B, D->C, E->D, F->E, G->F

## DAG_2: 10 genes (A-J)
## Multi-parent nodes: H (5 parents: F,D,G,J,B),
##                     I (5 parents: F,D,G,J,B),
##                     A (4 parents: D,G,J,B)
## Root genes: F, D, G, J, B
DAG_2 <- data.frame(
    From = c("Root", "Root", "Root", "Root", "Root",
             "F", "D", "G", "J", "B",
             "F", "D", "G", "J", "B",
             "D", "G", "J", "B",
             "H", "H"),
    To = c("F", "D", "G", "J", "B",
           rep("H", 5),
           rep("I", 5),
           rep("A", 4),
           "C", "E"),
    rerun_lambda = c(2, 1, .5, .35, .8,
                     2.7, 2.7, 2.7, 2.7, 2.7,
                     1.4, 1.4, 1.4, 1.4, 1.4,
                     0.6, 0.6, 0.6, 0.6,
                     0.9, 2.6))

## DAG_3: 9 genes (A-I)
## Multi-parent node: C (4 parents: H,B,F,D)
## Root genes: I, H
DAG_3 <- data.frame(
    From =         c("Root", "Root", "I", "I", "H", "E", "B", "G", "G", "F", "A", "D"),
    To =           c("I",    "H",    "E", "B", "C", "G", "C", "F", "A", "C", "D", "C"),
    rerun_lambda = c(.35,     .8,    .7,  .6,  .5,  .9,  .5,  1,   2,   .5,  .4,  .5)
)

## DAG_4: 9 genes (A-I)
## Multi-parent nodes: E (3 parents: C,A,B),
##                     G (3 parents: F,E,D),
##                     H (2 parents: I,G)
## Root genes: C, A, B, I
DAG_4 <- data.frame(
    From = c("Root", "Root", "Root", "Root",
             "C", "C", "A", "B", "B", "I",
             "F", "E", "D", "G"),
    To   = c("C", "A", "B", "I",
             "F", "E", "E", "E", "D", "H",
             "G", "G", "G", "H"),
    rerun_lambda = c(.8, 2, .35, .4,
                     .5, 1, 1, 1, .6, .9,
                     .7, .7, .7, .9)
)

## DAG_5: 9 genes (A-I)
## Multi-parent nodes: C (4 parents: A,F,G,I),
##                     D (3 parents: A,F,G),
##                     B (2 parents: E,D)
## Root genes: A, F, G, I
DAG_5 <- data.frame(
    From = c("Root", "Root", "Root", "Root",
             "A", "F", "G", "I",
             "A", "F", "G",
             "C", "D", "E", "D"),
    To = c("A", "F", "G", "I",
           "C", "C", "C", "C",
           "D", "D", "D",
           "E", "H", "B", "B"),
    rerun_lambda = c(2, .7, .8, .10, .4, .4, .4, .4, .5, .5, .5, .6, .9, 2.3, 2.3)
)

## DAG_6: 7 genes (A-G)
## Multi-parent nodes: A (2 parents: C,D),
##                     B (2 parents: F,G),
##                     E (2 parents: A,B)
## Root genes: C, D, F, G
DAG_6 <- data.frame(
    From = c(rep("Root", 4),     "C", "D", "F", "G", "A", "B"),
    To   = c("C", "D", "F", "G", "A", "A", "B", "B", "E", "E"),
    rerun_lambda = c(1, 2, 2.3, .4, .5, .5, .6, .6, .7, .7)
)

## DAG_7: 6 genes (A-F)
## Multi-parent nodes: A (3 parents: C,D,F),
##                     B (2 parents: D,F),
##                     E (2 parents: A,B)
## Root genes: C, D, F
DAG_7 <- data.frame(
    From = c(rep("Root", 3),   "C", "D", "F", "D", "F", "A", "B"),
    To   = c("C", "D", "F",    "A", "A", "A", "B", "B", "E", "E"),
    rerun_lambda = c(1, 2, 2.3, .4, .4, .4, .5, .5, .6, .6)
)

## You can get an idea of what they look like by doing
## evamtools:::DAG_plot_graphAM(DAG_2, "DAG_2")
## evamtools:::DAG_plot_graphAM(DAG_3, "DAG_3")
## evamtools:::DAG_plot_graphAM(DAG_4, "DAG_4")
## evamtools:::DAG_plot_graphAM(DAG_5, "DAG_5")
## evamtools:::DAG_plot_graphAM(DAG_6, "DAG_6")
## evamtools:::DAG_plot_graphAM(DAG_7, "DAG_7")


library(testthat)

set.seed(NULL)


### Utility functions

## CBN model data frame, number of genes ->
##   list(fitness_landscape, c, trm_scaled)
## suitable as input to intervene_fitness_landscape_every_gene.
## The model must have genes named A, B, C, ... (sequential from 1).
## a_cpm_2_si = 0.006 as for CBN/HESBCN in generate_f_landscape.
## c = 1/0.006 as used in generate_n_f_landscape_requir.
cbn_to_landscape_obj <- function(model, n_genes) {
    a_cpm_2_si <- 0.006
    genot_fitness <- ev2_cpm_to_fitness_genots(model, a = a_cpm_2_si)
    ## Newer OncoSimulR uses "Birth"; older uses "Fitness"
    fitness_col <- ifelse("Birth" %in% colnames(genot_fitness), "Birth", "Fitness")
    rfo <- cbind(genots_to_bin(genot_fitness$Genotype, n_genes),
                 Fitness = genot_fitness[[fitness_col]])
    rfo <- rfo[rfo[, "Fitness"] > 0, ]
    class(rfo) <- c("matrix", "array")
    trm_and_c <- suppressMessages(get_scaled_trm_adaptive(rfo, c = 1/a_cpm_2_si))
    list(fitness_landscape = rfo,
         c = trm_and_c$c,
         trm_scaled = trm_and_c$trm_scaled)
}

## CBN model data frame (with rerun_lambda column) ->
##   run four procedures + row-permuted model and assert all equal.
## Genes must be named sequentially A, B, C, ... (no gaps).
check_cbn_four_procedures <- function(model, label = "") {
    m <- model
    n <- length(setdiff(unique(c(m$From, m$To)), "Root"))
    cat("\n CBN", if (nchar(label) > 0) label else "", ": n =", n, "genes\n")

    ## Procedure 1: standard kill_gene via DAG edge removal
    i1 <- intervene_cpm_every_gene(list(CBN_model = m), "CBN")

    ## Procedure 2: kill by setting parameters to 0
    i2 <- suppressWarnings(
        intervene_cpm_every_gene(list(CBN_model = m), "CBN",
                                 kill_gene_funct = kill_gene_by_params_to_0))

    ## Procedure 3: remove genotypes from pre-computed TRM.
    ## Requires a list with CBN_model and CBN_trans_rate_mat.
    ## get_full_output provides the trans_rate_mat; we add CBN_model
    ## explicitly (in case CBN_model_2_output does not include it).
    full_out <- suppressWarnings(get_full_output(m))
    full_out$CBN_model <- m
    i3 <- intervene_cpm_trm_rm_every_gene(full_out, "CBN")

    ## Procedure 4: fitness landscape derived from model.
    landscape_obj <- cbn_to_landscape_obj(m, n)
    i4 <- intervene_fitness_landscape_every_gene(landscape_obj)

    ## Procedure 1 with row-permuted model (order must not matter)
    m_perm <- m[sample(1:nrow(m)), ]
    i1_perm <- intervene_cpm_every_gene(list(CBN_model = m_perm), "CBN")

    lbl <- if (nchar(label) > 0) paste0(" [", label, "]") else ""
    expect_equal(i1, i2, label = paste0("i1 vs i2", lbl))
    expect_equal(i1, i3, label = paste0("i1 vs i3", lbl))
    expect_equal(i1, i4, label = paste0("i1 vs i4", lbl))
    expect_equal(i1, i1_perm, label = paste0("i1 vs i1_perm", lbl))
}


### Test blocks: one per DAG
##
## DAG structures (multi-parent nodes, all AND in CBN):
##   DAG_2: H (5 parents: F,D,G,J,B), I (5 parents: F,D,G,J,B), A (4 parents: D,G,J,B)
##   DAG_3: C (4 parents: H,B,F,D)
##   DAG_4: E (3 parents: C,A,B), G (3 parents: F,E,D), H (2 parents: I,G)
##   DAG_5: C (4 parents: A,F,G,I), D (3 parents: A,F,G), B (2 parents: E,D)
##   DAG_6: A (2 parents: C,D), B (2 parents: F,G), E (2 parents: A,B)
##   DAG_7: A (3 parents: C,D,F), B (2 parents: D,F), E (2 parents: A,B)

test_that("CBN DAG_2", {
    local_edition(3)
    check_cbn_four_procedures(DAG_2, label = "DAG_2")
})

test_that("CBN DAG_3", {
    local_edition(3)
    check_cbn_four_procedures(DAG_3, label = "DAG_3")
})

test_that("CBN DAG_4", {
    local_edition(3)
    check_cbn_four_procedures(DAG_4, label = "DAG_4")
})

test_that("CBN DAG_5", {
    local_edition(3)
    check_cbn_four_procedures(DAG_5, label = "DAG_5")
})

test_that("CBN DAG_6", {
    local_edition(3)
    check_cbn_four_procedures(DAG_6, label = "DAG_6")
})

test_that("CBN DAG_7", {
    local_edition(3)
    check_cbn_four_procedures(DAG_7, label = "DAG_7")
})


### Cascade tests with ground-truth comparison
##
## CBN is always AND, so killing any parent of a multi-parent node
## always cascades to the child (and its descendants).
##
## Two cascade models:
##
## (a) 2-parent AND cascade:
##     Root->A, Root->B, A->C, B->C, C->D  (rows 1-5)
##     Kill A: C needs both A and B; C unreachable -> D cascades.
##             Surviving rows: 2 (Root->B).
##     Kill B: symmetric. Surviving rows: 1 (Root->A).
##     Kill C: D cascades. Surviving rows: 1,2.
##     Kill D: leaf, no cascade. Surviving rows: 1,2,3,4.
##
## (b) Chain cascade:
##     Root->A, A->B, B->C  (rows 1-3)
##     Kill A: B and C both unreachable (chain). Surviving: 0 rows -> WT = 1.
##     Kill B: C cascades. Surviving rows: 1 (Root->A).
##     Kill C: leaf. Surviving rows: 1,2.
##
## Surviving models specified by row index, not by kill_gene, to avoid
## circularity.

test_that("CBN cascade: 2-parent AND — killing one parent removes child and its descendants", {
    local_edition(3)

    m <- data.frame(
        From        = c("Root", "Root", "A", "B", "C"),
        To          = c("A",    "B",    "C", "C", "D"),
        rerun_lambda = c(2,     0.8,    1,   1,   0.5),
        stringsAsFactors = FALSE
    )

    preds_from_model <- function(x) {
        if (nrow(x) == 0) {
            return(list(genot_freqs = c(WT = 1), hitting_probs_from_WT = c(WT = 1.0)))
        }
        tmp <- suppressWarnings(get_full_output(x))
        f <- tmp$CBN_predicted_genotype_freqs
        hp <- tmp$CBN_hitting_probs_from_WT
        list(genot_freqs = f[f > 0], hitting_probs_from_WT = filter_hp_keep_wt(hp))
    }

    res <- intervene_cpm_every_gene(list(CBN_model = m), "CBN")

    ## Kill A: surviving rows 2 only
    expect_equal(get_interv(res, "I:A"), preds_from_model(m[2, ]))

    ## Kill B: surviving rows 1 only
    expect_equal(get_interv(res, "I:B"), preds_from_model(m[1, ]))

    ## Kill C: surviving rows 1,2
    expect_equal(get_interv(res, "I:C"), preds_from_model(m[c(1, 2), ]))

    ## Kill D: surviving rows 1,2,3,4
    expect_equal(get_interv(res, "I:D"), preds_from_model(m[c(1, 2, 3, 4), ]))

    ## Also verify all four procedures agree on this model
    check_cbn_four_procedures(m, label = "cascade 2-parent AND")
})


test_that("CBN cascade: chain — killing a gene removes all descendants in the chain", {
    local_edition(3)

    m <- data.frame(
        From        = c("Root", "A",  "B"),
        To          = c("A",    "B",  "C"),
        rerun_lambda = c(2,     0.8,  0.5),
        stringsAsFactors = FALSE
    )

    preds_from_model <- function(x) {
        if (nrow(x) == 0) {
            return(list(genot_freqs = c(WT = 1), hitting_probs_from_WT = c(WT = 1.0)))
        }
        tmp <- suppressWarnings(get_full_output(x))
        f <- tmp$CBN_predicted_genotype_freqs
        hp <- tmp$CBN_hitting_probs_from_WT
        list(genot_freqs = f[f > 0], hitting_probs_from_WT = filter_hp_keep_wt(hp))
    }

    res <- intervene_cpm_every_gene(list(CBN_model = m), "CBN")

    ## Kill A: B and C both cascade; no rows survive -> WT = 1
    expect_equal(get_interv(res, "I:A"), preds_from_model(m[integer(0), ]))

    ## Kill B: C cascades; surviving rows 1 (Root->A)
    expect_equal(get_interv(res, "I:B"), preds_from_model(m[1, ]))

    ## Kill C: leaf; surviving rows 1,2
    expect_equal(get_interv(res, "I:C"), preds_from_model(m[c(1, 2), ]))

    ## Also verify all four procedures agree on this model
    check_cbn_four_procedures(m, label = "cascade chain")
})


### No-over-cascade test
##
## Model: Root->A, Root->B, A->C, B->C, B->D  (rows 1-5)
## C needs both A and B (AND). D needs only B.
##
## Kill A: C cascades (needs A+B, no A). D does NOT cascade (B still present).
##         Surviving: rows 2,5 (Root->B, B->D).
## Kill B: C cascades (needs B). D also cascades (needs B). Surviving: row 1.
## Kill C: D unaffected. Surviving: rows 1,2,5.
## Kill D: leaf. Surviving: rows 1,2,3,4.
##
## The key check is kill A: verifies the cascade stops at C and does not
## propagate to D, which has a surviving parent (B).

test_that("CBN no-over-cascade: cascade stops at node that still has all parents alive", {
    local_edition(3)

    m <- data.frame(
        From        = c("Root", "Root", "A", "B", "B"),
        To          = c("A",    "B",    "C", "C", "D"),
        rerun_lambda = c(2,     0.8,    1,   1,   0.5),
        stringsAsFactors = FALSE
    )

    preds_from_model <- function(x) {
        if (nrow(x) == 0) {
            return(list(genot_freqs = c(WT = 1), hitting_probs_from_WT = c(WT = 1.0)))
        }
        tmp <- suppressWarnings(get_full_output(x))
        f <- tmp$CBN_predicted_genotype_freqs
        hp <- tmp$CBN_hitting_probs_from_WT
        list(genot_freqs = f[f > 0], hitting_probs_from_WT = filter_hp_keep_wt(hp))
    }

    res <- intervene_cpm_every_gene(list(CBN_model = m), "CBN")

    ## Kill A: C cascades (AND needs A+B), D survives (only needs B).
    expect_equal(get_interv(res, "I:A"), preds_from_model(m[c(2, 5), ]))

    ## Kill B: C cascades (needs B), D also cascades (needs B).
    expect_equal(get_interv(res, "I:B"), preds_from_model(m[1, ]))

    ## Kill C: D unaffected.
    expect_equal(get_interv(res, "I:C"), preds_from_model(m[c(1, 2, 5), ]))

    ## Kill D: leaf.
    expect_equal(get_interv(res, "I:D"), preds_from_model(m[c(1, 2, 3, 4), ]))

    ## Also verify all four procedures agree on this model
    check_cbn_four_procedures(m, label = "no-over-cascade")
})


### Full ground-truth test: DAG_6, every gene killed
##
## DAG_6 rows:
##   1: Root -> C
##   2: Root -> D
##   3: Root -> F
##   4: Root -> G
##   5: C -> A  (AND, needs C)
##   6: D -> A  (AND, needs D)
##   7: F -> B  (AND, needs F)
##   8: G -> B  (AND, needs G)
##   9: A -> E  (AND, needs A)
##  10: B -> E  (AND, needs B)
##
## Multi-parent nodes: A (parents C,D), B (parents F,G), E (parents A,B).
##
## Surviving rows after each kill (worked by hand):
##   Kill C: A unreachable (needs C+D) -> cascade rows 5,6; E unreachable -> cascade rows 9,10. B ok.
##           Surviving: 2,3,4,7,8
##   Kill D: A unreachable -> cascade rows 5,6; E unreachable -> cascade rows 9,10. B ok.
##           Surviving: 1,3,4,7,8
##   Kill F: B unreachable (needs F+G) -> cascade rows 7,8; E unreachable -> cascade rows 9,10. A ok.
##           Surviving: 1,2,4,5,6
##   Kill G: B unreachable -> cascade rows 7,8; E unreachable -> cascade rows 9,10. A ok.
##           Surviving: 1,2,3,5,6
##   Kill A: E unreachable (needs A+B) -> cascade rows 9,10. B ok.
##           Surviving: 1,2,3,4,7,8
##   Kill B: E unreachable -> cascade rows 9,10. A ok.
##           Surviving: 1,2,3,4,5,6
##   Kill E: leaf, no cascade.
##           Surviving: 1,2,3,4,5,6,7,8

test_that("CBN DAG_6: ground-truth comparison for all gene kills", {
    local_edition(3)

    preds_from_model <- function(x) {
        if (nrow(x) == 0) {
            return(list(genot_freqs = c(WT = 1), hitting_probs_from_WT = c(WT = 1.0)))
        }
        tmp <- suppressWarnings(get_full_output(x))
        f <- tmp$CBN_predicted_genotype_freqs
        hp <- tmp$CBN_hitting_probs_from_WT
        list(genot_freqs = f[f > 0], hitting_probs_from_WT = filter_hp_keep_wt(hp))
    }

    res <- intervene_cpm_every_gene(list(CBN_model = DAG_6), "CBN")

    expect_equal(get_interv(res, "I:C"), preds_from_model(DAG_6[c(2, 3, 4, 7, 8), ]))
    expect_equal(get_interv(res, "I:D"), preds_from_model(DAG_6[c(1, 3, 4, 7, 8), ]))
    expect_equal(get_interv(res, "I:F"), preds_from_model(DAG_6[c(1, 2, 4, 5, 6), ]))
    expect_equal(get_interv(res, "I:G"), preds_from_model(DAG_6[c(1, 2, 3, 5, 6), ]))
    expect_equal(get_interv(res, "I:A"), preds_from_model(DAG_6[c(1, 2, 3, 4, 7, 8), ]))
    expect_equal(get_interv(res, "I:B"), preds_from_model(DAG_6[c(1, 2, 3, 4, 5, 6), ]))
    expect_equal(get_interv(res, "I:E"), preds_from_model(DAG_6[c(1, 2, 3, 4, 5, 6, 7, 8), ]))
})


### Full ground-truth test: DAG_7, every gene killed
##
## DAG_7 rows:
##   1: Root -> C
##   2: Root -> D
##   3: Root -> F
##   4: C -> A  (AND, needs C+D+F)
##   5: D -> A  (AND)
##   6: F -> A  (AND)
##   7: D -> B  (AND, needs D+F)
##   8: F -> B  (AND)
##   9: A -> E  (AND, needs A+B)
##  10: B -> E  (AND)
##
## Multi-parent nodes: A (parents C,D,F), B (parents D,F), E (parents A,B).
## Key structural difference from DAG_6: D is a shared parent of both A and B,
## so killing D cascades TWO separate multi-parent nodes simultaneously.
##
## Surviving rows after each kill (worked by hand):
##   Kill C: A unreachable (needs C+D+F) -> cascade rows 4,5,6;
##           B still reachable (D and F survive);
##           E unreachable (needs A+B, A gone) -> cascade rows 9,10.
##           Surviving: 2,3,7,8
##   Kill D: A unreachable -> cascade rows 4,5,6;
##           B unreachable (needs D+F) -> cascade rows 7,8;
##           E unreachable -> cascade rows 9,10.
##           Surviving: 1,3
##   Kill F: A unreachable -> cascade rows 4,5,6;
##           B unreachable -> cascade rows 7,8;
##           E unreachable -> cascade rows 9,10.
##           Surviving: 1,2
##   Kill A: E unreachable (needs A+B) -> cascade rows 9,10. B ok.
##           Surviving: 1,2,3,7,8
##   Kill B: E unreachable -> cascade rows 9,10. A ok.
##           Surviving: 1,2,3,4,5,6
##   Kill E: leaf, no cascade.
##           Surviving: 1,2,3,4,5,6,7,8

test_that("CBN DAG_7: ground-truth comparison for all gene kills", {
    ## Recall all four procedures are compared above on this very DAG
    ## so running them here again adds nothing.
    local_edition(3)

    preds_from_model <- function(x) {
        if (nrow(x) == 0) {
            return(list(genot_freqs = c(WT = 1), hitting_probs_from_WT = c(WT = 1.0)))
        }
        tmp <- suppressWarnings(get_full_output(x))
        f <- tmp$CBN_predicted_genotype_freqs
        hp <- tmp$CBN_hitting_probs_from_WT
        list(genot_freqs = f[f > 0], hitting_probs_from_WT = filter_hp_keep_wt(hp))
    }

    res <- intervene_cpm_every_gene(list(CBN_model = DAG_7), "CBN")

    expect_equal(get_interv(res, "I:C"), preds_from_model(DAG_7[c(2, 3, 7, 8), ]))
    expect_equal(get_interv(res, "I:D"), preds_from_model(DAG_7[c(1, 3), ]))
    expect_equal(get_interv(res, "I:F"), preds_from_model(DAG_7[c(1, 2), ]))
    expect_equal(get_interv(res, "I:A"), preds_from_model(DAG_7[c(1, 2, 3, 7, 8), ]))
    expect_equal(get_interv(res, "I:B"), preds_from_model(DAG_7[c(1, 2, 3, 4, 5, 6), ]))
    expect_equal(get_interv(res, "I:E"), preds_from_model(DAG_7[c(1, 2, 3, 4, 5, 6, 7, 8), ]))
})


### Full ground-truth test: DAG_2, every gene killed
##
## DAG_2 rows:
##   1: Root->F   2: Root->D   3: Root->G   4: Root->J   5: Root->B
##   6: F->H   7: D->H   8: G->H   9: J->H   10: B->H   (H needs all 5: F,D,G,J,B)
##  11: F->I  12: D->I  13: G->I  14: J->I  15: B->I    (I needs all 5: F,D,G,J,B)
##  16: D->A  17: G->A  18: J->A  19: B->A              (A needs 4:    D,G,J,B)
##  20: H->C  21: H->E                                   (C,E: single parent H)
##
## Kill F: H unreachable (needs all 5), I unreachable; A unaffected (F not its parent);
##         C,E cascade. Surviving: 2,3,4,5,16,17,18,19
## Kill D: H, I, A all unreachable; C,E cascade. Surviving: 1,3,4,5
## Kill G: H, I, A all unreachable; C,E cascade. Surviving: 1,2,4,5
## Kill J: H, I, A all unreachable; C,E cascade. Surviving: 1,2,3,5
## Kill B: H, I, A all unreachable; C,E cascade. Surviving: 1,2,3,4
## Kill H: C,E cascade (only parent H); I and A unaffected. Surviving: 1-5,11-19
## Kill I: leaf (no children). Surviving: 1-10,16-21
## Kill A: leaf. Surviving: 1-15,20,21
## Kill C: leaf. Surviving: 1-19,21
## Kill E: leaf. Surviving: 1-20

test_that("CBN DAG_2: ground-truth comparison for all gene kills", {
    local_edition(3)

    preds_from_model <- function(x) {
        if (nrow(x) == 0) {
            return(list(genot_freqs = c(WT = 1), hitting_probs_from_WT = c(WT = 1.0)))
        }
        tmp <- suppressWarnings(get_full_output(x))
        f <- tmp$CBN_predicted_genotype_freqs
        hp <- tmp$CBN_hitting_probs_from_WT
        list(genot_freqs = f[f > 0], hitting_probs_from_WT = filter_hp_keep_wt(hp))
    }

    res <- intervene_cpm_every_gene(list(CBN_model = DAG_2), "CBN")

    expect_equal(get_interv(res, "I:F"), preds_from_model(DAG_2[c(2, 3, 4, 5, 16, 17, 18, 19), ]))
    expect_equal(get_interv(res, "I:D"), preds_from_model(DAG_2[c(1, 3, 4, 5), ]))
    expect_equal(get_interv(res, "I:G"), preds_from_model(DAG_2[c(1, 2, 4, 5), ]))
    expect_equal(get_interv(res, "I:J"), preds_from_model(DAG_2[c(1, 2, 3, 5), ]))
    expect_equal(get_interv(res, "I:B"), preds_from_model(DAG_2[c(1, 2, 3, 4), ]))
    expect_equal(get_interv(res, "I:H"), preds_from_model(DAG_2[c(1:5, 11:19), ]))
    expect_equal(get_interv(res, "I:I"), preds_from_model(DAG_2[c(1:10, 16:21), ]))
    expect_equal(get_interv(res, "I:A"), preds_from_model(DAG_2[c(1:15, 20, 21), ]))
    expect_equal(get_interv(res, "I:C"), preds_from_model(DAG_2[c(1:19, 21), ]))
    expect_equal(get_interv(res, "I:E"), preds_from_model(DAG_2[c(1:20), ]))
})


### Full ground-truth test: DAG_3, every gene killed
##
## DAG_3 rows:
##   1: Root->I   2: Root->H
##   3: I->E   4: I->B
##   5: H->C   (AND parent of C)
##   6: E->G
##   7: B->C   (AND parent of C)
##   8: G->F   9: G->A
##  10: F->C   (AND parent of C)
##  11: A->D
##  12: D->C   (AND parent of C)
## C needs all 4 parents: H, B, F, D.
##
## Kill I: E,B cascade; G cascades (only parent E); F,A cascade (only parent G);
##         D cascades (only parent A); C unreachable (H survives but B,F,D gone).
##         Surviving: 2
## Kill H: C unreachable (H gone, needs H+B+F+D). Rest of chain unaffected.
##         Surviving: 1,3,4,6,8,9,11
## Kill E: G cascades; F,A cascade; D cascades; C unreachable (B survives, F,D gone).
##         Surviving: 1,2,4
## Kill B: C unreachable (B gone). Rest unaffected.
##         Surviving: 1,2,3,6,8,9,11
## Kill G: F,A cascade; D cascades; C unreachable (H,B survive, F,D gone).
##         Surviving: 1,2,3,4
## Kill F: C unreachable (F gone). A,D chain unaffected.
##         Surviving: 1,2,3,4,6,9,11
## Kill A: D cascades (only parent A); C unreachable (D gone).
##         Surviving: 1,2,3,4,6,8
## Kill D: C unreachable (D gone). A unaffected.
##         Surviving: 1,2,3,4,6,8,9
## Kill C: leaf. Surviving: 1,2,3,4,6,8,9,11

test_that("CBN DAG_3: ground-truth comparison for all gene kills", {
    local_edition(3)

    preds_from_model <- function(x) {
        if (nrow(x) == 0) {
            return(list(genot_freqs = c(WT = 1), hitting_probs_from_WT = c(WT = 1.0)))
        }
        tmp <- suppressWarnings(get_full_output(x))
        f <- tmp$CBN_predicted_genotype_freqs
        hp <- tmp$CBN_hitting_probs_from_WT
        list(genot_freqs = f[f > 0], hitting_probs_from_WT = filter_hp_keep_wt(hp))
    }

    res <- intervene_cpm_every_gene(list(CBN_model = DAG_3), "CBN")

    expect_equal(get_interv(res, "I:I"), preds_from_model(DAG_3[2, ]))
    expect_equal(get_interv(res, "I:H"), preds_from_model(DAG_3[c(1, 3, 4, 6, 8, 9, 11), ]))
    expect_equal(get_interv(res, "I:E"), preds_from_model(DAG_3[c(1, 2, 4), ]))
    expect_equal(get_interv(res, "I:B"), preds_from_model(DAG_3[c(1, 2, 3, 6, 8, 9, 11), ]))
    expect_equal(get_interv(res, "I:G"), preds_from_model(DAG_3[c(1, 2, 3, 4), ]))
    expect_equal(get_interv(res, "I:F"), preds_from_model(DAG_3[c(1, 2, 3, 4, 6, 9, 11), ]))
    expect_equal(get_interv(res, "I:A"), preds_from_model(DAG_3[c(1, 2, 3, 4, 6, 8), ]))
    expect_equal(get_interv(res, "I:D"), preds_from_model(DAG_3[c(1, 2, 3, 4, 6, 8, 9), ]))
    expect_equal(get_interv(res, "I:C"), preds_from_model(DAG_3[c(1, 2, 3, 4, 6, 8, 9, 11), ]))
})


### Full ground-truth test: DAG_4, every gene killed
##
## DAG_4 rows:
##   1: Root->C   2: Root->A   3: Root->B   4: Root->I
##   5: C->F
##   6: C->E   7: A->E   8: B->E    (E needs all 3: C,A,B)
##   9: B->D
##  10: I->H                          (H needs both: I,G)
##  11: F->G  12: E->G  13: D->G    (G needs all 3: F,E,D)
##  14: G->H                         (H needs both: I,G)
##
## Kill C: F cascades (only parent C); E unreachable (C gone, needs C+A+B);
##         G unreachable (F,E gone, needs F+E+D); H unreachable (G gone).
##         D unaffected (parent B). Surviving: 2,3,4,9
## Kill A: E unreachable; G unreachable; H unreachable. F,D unaffected.
##         Surviving: 1,3,4,5,9
## Kill B: E unreachable; D cascades (only parent B); G unreachable; H unreachable.
##         F unaffected. Surviving: 1,2,4,5
## Kill I: H unreachable (needs I+G; I gone). Surviving: 1,2,3,5,6,7,8,9,11,12,13
## Kill F: G unreachable (needs F+E+D; F gone); H unreachable.
##         Surviving: 1,2,3,4,6,7,8,9
## Kill E: G unreachable; H unreachable. F,D unaffected.
##         Surviving: 1,2,3,4,5,9
## Kill D: G unreachable; H unreachable. F,E unaffected.
##         Surviving: 1,2,3,4,5,6,7,8
## Kill G: H unreachable. Surviving: 1,2,3,4,5,6,7,8,9
## Kill H: leaf. Surviving: 1,2,3,4,5,6,7,8,9,11,12,13

test_that("CBN DAG_4: ground-truth comparison for all gene kills", {
    local_edition(3)

    preds_from_model <- function(x) {
        if (nrow(x) == 0) {
            return(list(genot_freqs = c(WT = 1), hitting_probs_from_WT = c(WT = 1.0)))
        }
        tmp <- suppressWarnings(get_full_output(x))
        f <- tmp$CBN_predicted_genotype_freqs
        hp <- tmp$CBN_hitting_probs_from_WT
        list(genot_freqs = f[f > 0], hitting_probs_from_WT = filter_hp_keep_wt(hp))
    }

    res <- intervene_cpm_every_gene(list(CBN_model = DAG_4), "CBN")

    expect_equal(get_interv(res, "I:C"), preds_from_model(DAG_4[c(2, 3, 4, 9), ]))
    expect_equal(get_interv(res, "I:A"), preds_from_model(DAG_4[c(1, 3, 4, 5, 9), ]))
    expect_equal(get_interv(res, "I:B"), preds_from_model(DAG_4[c(1, 2, 4, 5), ]))
    expect_equal(get_interv(res, "I:I"), preds_from_model(DAG_4[c(1, 2, 3, 5, 6, 7, 8, 9, 11, 12, 13), ]))
    expect_equal(get_interv(res, "I:F"), preds_from_model(DAG_4[c(1, 2, 3, 4, 6, 7, 8, 9), ]))
    expect_equal(get_interv(res, "I:E"), preds_from_model(DAG_4[c(1, 2, 3, 4, 5, 9), ]))
    expect_equal(get_interv(res, "I:D"), preds_from_model(DAG_4[c(1, 2, 3, 4, 5, 6, 7, 8), ]))
    expect_equal(get_interv(res, "I:G"), preds_from_model(DAG_4[c(1, 2, 3, 4, 5, 6, 7, 8, 9), ]))
    expect_equal(get_interv(res, "I:H"), preds_from_model(DAG_4[c(1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13), ]))
})


### Full ground-truth test: DAG_5, every gene killed
##
## DAG_5 rows:
##   1: Root->A   2: Root->F   3: Root->G   4: Root->I
##   5: A->C   6: F->C   7: G->C   8: I->C    (C needs all 4: A,F,G,I)
##   9: A->D  10: F->D  11: G->D              (D needs all 3: A,F,G)
##  12: C->E
##  13: D->H
##  14: E->B  15: D->B                        (B needs both: E,D)
##
## Kill A: C unreachable (needs A,F,G,I); D unreachable (needs A,F,G);
##         E,H,B all cascade. Surviving: 2,3,4
## Kill F: C unreachable; D unreachable; E,H,B cascade. Surviving: 1,3,4
## Kill G: C unreachable; D unreachable; E,H,B cascade. Surviving: 1,2,4
## Kill I: C unreachable (I not a parent of D, so D survives);
##         E cascades; B unreachable (E gone, D alive but needs both).
##         D,H survive. Surviving: 1,2,3,9,10,11,13
## Kill C: E cascades; B unreachable (E gone). D,H unaffected.
##         Surviving: 1,2,3,4,9,10,11,13
## Kill D: H cascades (only parent D); B unreachable (D gone).
##         C,E unaffected. Surviving: 1,2,3,4,5,6,7,8,12
## Kill E: B unreachable (E gone, needs E+D). D,H unaffected.
##         Surviving: 1,2,3,4,5,6,7,8,9,10,11,13
## Kill H: leaf. Surviving: 1,2,3,4,5,6,7,8,9,10,11,12,14,15
## Kill B: leaf. Surviving: 1,2,3,4,5,6,7,8,9,10,11,12,13

test_that("CBN DAG_5: ground-truth comparison for all gene kills", {
    local_edition(3)

    preds_from_model <- function(x) {
        if (nrow(x) == 0) {
            return(list(genot_freqs = c(WT = 1), hitting_probs_from_WT = c(WT = 1.0)))
        }
        tmp <- suppressWarnings(get_full_output(x))
        f <- tmp$CBN_predicted_genotype_freqs
        hp <- tmp$CBN_hitting_probs_from_WT
        list(genot_freqs = f[f > 0], hitting_probs_from_WT = filter_hp_keep_wt(hp))
    }

    res <- intervene_cpm_every_gene(list(CBN_model = DAG_5), "CBN")

    expect_equal(get_interv(res, "I:A"), preds_from_model(DAG_5[c(2, 3, 4), ]))
    expect_equal(get_interv(res, "I:F"), preds_from_model(DAG_5[c(1, 3, 4), ]))
    expect_equal(get_interv(res, "I:G"), preds_from_model(DAG_5[c(1, 2, 4), ]))
    expect_equal(get_interv(res, "I:I"), preds_from_model(DAG_5[c(1, 2, 3, 9, 10, 11, 13), ]))
    expect_equal(get_interv(res, "I:C"), preds_from_model(DAG_5[c(1, 2, 3, 4, 9, 10, 11, 13), ]))
    expect_equal(get_interv(res, "I:D"), preds_from_model(DAG_5[c(1, 2, 3, 4, 5, 6, 7, 8, 12), ]))
    expect_equal(get_interv(res, "I:E"), preds_from_model(DAG_5[c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13), ]))
    expect_equal(get_interv(res, "I:H"), preds_from_model(DAG_5[c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15), ]))
    expect_equal(get_interv(res, "I:B"), preds_from_model(DAG_5[c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13), ]))
})

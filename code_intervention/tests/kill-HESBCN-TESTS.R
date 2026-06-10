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

## First, some basic tests of HESBCN-killing. Next, thorough tests of
## HESBCN interventions using six hand-crafted DAG structures
## focused on what we really care about: probs. of predicted genotypes
## (DAG_2 through DAG_7 from new-dags-for-trm-test.R).
##
## For each DAG, we exhaustively enumerate all AND/OR/XOR combinations for
## multi-parent nodes, building 27+3+27+27+27+27 = 138 HESBCN variants.
##
## For each variant, four procedures are compared:
##   1. Standard kill_gene (modifies DAG edges)
##   2. kill_gene_by_params_to_0 (zeroes parameters)
##   3. intervene_cpm_trm_rm_every_gene (removes genotypes from pre-computed TRM;
##      does NOT call get_full_output for the intervention step)
##   4. intervene_fitness_landscape_every_gene (fitness landscape derived from model)
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
## These tests are slow. Run on a machine with spare time.

## options(intervention_every_gene_cores = parallel::detectCores())

library(testthat)
pwd <- getwd()
setwd("../")
source("intervention.R")
setwd(pwd)

set.seed(NULL)

## Extract a single named intervention from the full intervention output as a
## list(genot_freqs, hitting_probs_from_WT), for direct comparison with
## preds_from_hesbcn_model / preds_from_model.
get_interv <- function(res, name) {
    return(res[[name]])
}


test_that("HESBCN, mixture of relations", {
    local_edition(3)

    m1x <- data.frame(
        From = c("Root", "Root", "Root", "Root", "A", "B", "B", "D", "B", "F"),
        To =   c("A",    "B",    "D",    "F",    "C", "C", "E", "E", "G", "G"),
        Lambdas = c(1, 2, 3, 4, 5, 5, 6, 6, 7, 7),
        Relation = c(rep("Single", 4), "XOR", "XOR", "AND", "AND", "OR", "OR"))

    ## Permutations on similar theme
    m2x <- data.frame(
        From = c("Root", "Root", "Root", "A", "B", "D", "C", "E"),
        To =   c("A",    "B",    "D",    "C", "C", "E", "E", "F"),
        Lambdas = c(1, 2, 3, 4, 4, 5, 5, 6),
        Relation = c(rep("Single", 3), "AND", "AND", "OR", "OR", "Single"))

    m3x <- data.frame(
        From = c("Root", "Root", "Root", "A", "B", "D", "C", "E"),
        To =   c("A",    "B",    "D",    "C", "C", "E", "E", "F"),
        Lambdas = c(1, 2, 3, 4, 4, 5, 5, 6),
        Relation = c(rep("Single", 3), "AND", "AND", "XOR", "XOR", "Single"))

    m4x <- data.frame(
        From = c("Root", "Root", "Root", "A", "B", "D", "C", "E"),
        To =   c("A",    "B",    "D",    "C", "C", "E", "E", "F"),
        Lambdas = c(1, 2, 3, 4, 4, 5, 5, 6),
        Relation = c(rep("Single", 3), "OR", "OR", "AND", "AND", "Single"))


    kills_m1x <- lapply(LETTERS[1:7], function(v) kill_gene(m1x, v))
    names(kills_m1x) <- LETTERS[1:7]

    kills_m2x <- lapply(LETTERS[1:6], function(v) kill_gene(m2x, v))
    kills_m3x <- lapply(LETTERS[1:6], function(v) kill_gene(m3x, v))
    kills_m4x <- lapply(LETTERS[1:6], function(v) kill_gene(m4x, v))
    names(kills_m2x) <- names(kills_m3x) <- names(kills_m4x) <- LETTERS[1:6]


    expect_equal(kills_m1x[["A"]], m1x[-c(1, 5), ])
    expect_equal(kills_m1x[["B"]], m1x[-c(2, 6, 7, 8, 9), ])
    expect_equal(kills_m1x[["C"]], m1x[-c(5, 6), ])
    expect_equal(kills_m1x[["D"]], m1x[-c(3, 7, 8), ])
    expect_equal(kills_m1x[["E"]], m1x[-c(7, 8), ])
    expect_equal(kills_m1x[["F"]], m1x[-c(4, 10), ])
    expect_equal(kills_m1x[["G"]], m1x[-c(9, 10), ])

    expect_equal(kills_m2x[["A"]], m2x[-c(1, 4, 5, 7), ])
    expect_equal(kills_m2x[["B"]], m2x[-c(2, 4, 5, 7), ])
    expect_equal(kills_m2x[["C"]], m2x[-c(4, 5, 7), ])
    expect_equal(kills_m2x[["D"]], m2x[-c(3, 6), ])
    expect_equal(kills_m2x[["E"]], m2x[-c(6:8), ])
    expect_equal(kills_m2x[["F"]], m2x[-c(8), ])

    expect_equal(kills_m3x[["A"]], m3x[-c(1, 4, 5, 7), ])
    expect_equal(kills_m3x[["B"]], m3x[-c(2, 4, 5, 7), ])
    expect_equal(kills_m3x[["C"]], m3x[-c(4, 5, 7), ])
    expect_equal(kills_m3x[["D"]], m3x[-c(3, 6), ])
    expect_equal(kills_m3x[["E"]], m3x[-c(6:8), ])
    expect_equal(kills_m3x[["F"]], m3x[-c(8), ])

    expect_equal(kills_m4x[["A"]], m4x[-c(1, 4), ])
    expect_equal(kills_m4x[["B"]], m4x[-c(2, 5), ])
    expect_equal(kills_m4x[["C"]], m4x[-c(4:8), ])
    expect_equal(kills_m4x[["D"]], m4x[-c(3, 6:8), ])
    expect_equal(kills_m4x[["E"]], m4x[-c(6:8), ])
    expect_equal(kills_m4x[["F"]], m4x[-c(8), ])

    ## From a m4d before
    m4dx <- data.frame(From = c("Root", "Root", "A", "B", "B", "C", "D"),
                       To   = c("A",    "B",    "C", "C", "D",  "E", "F"),
                       Lambdas = c(1, 2, 3, 3, 4, 5, 6),
                       Relation = c("Single", "Single", "XOR", "XOR", "Single",
                                    "Single", "Single"))

    kills_m4dx <- lapply(LETTERS[1:6], function(v) kill_gene(m4dx, v))
    names(kills_m4dx) <- LETTERS[1:6]

    expect_equal(kills_m4dx[["A"]], m4dx[-c(1, 3), ])
    expect_equal(kills_m4dx[["B"]], m4dx[-c(2, 4, 5, 7), ])
    expect_equal(kills_m4dx[["C"]], m4dx[-c(3, 4, 6), ])
    expect_equal(kills_m4dx[["D"]], m4dx[-c(5, 7), ])
    expect_equal(kills_m4dx[["E"]], m4dx[-c(6), ])
    expect_equal(kills_m4dx[["F"]], m4dx[-c(7), ])
})



### Tests with DAG_2 to DAG_7


## The DAGs are the same as used for CBN and OncoBN

## You can see what they look like by doing
## evamtools:::DAG_plot_graphAM(DAG_2, "DAG_2")
## evamtools:::DAG_plot_graphAM(DAG_3, "DAG_3")
## evamtools:::DAG_plot_graphAM(DAG_4, "DAG_4")
## evamtools:::DAG_plot_graphAM(DAG_5, "DAG_5")
## evamtools:::DAG_plot_graphAM(DAG_6, "DAG_6")
## evamtools:::DAG_plot_graphAM(DAG_7, "DAG_7")



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


#### Utility functions

## HESBCN model data frame, number of genes ->
##   list(fitness_landscape, c, trm_scaled)
## suitable as input to intervene_fitness_landscape_every_gene.
## The model must have genes named A, B, C, ... (sequential from 1).
## a_cpm_2_si = 0.006 as for CBN/HESBCN in generate_f_landscape.
## c = 1/0.006 as used in generate_n_f_landscape_requir.
hesbcn_to_landscape_obj <- function(model, n_genes) {
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

## Independent reimplementation of HESBCN cascade rules for ground-truth testing.
## Returns indices of rows in model that survive after killing gene `killed`.
##
## Rules (applied iteratively until no more genes go dead):
##   AND node: unreachable if ANY non-Root parent is dead.
##   OR/XOR/Single node: unreachable only if ALL non-Root parents are dead.
##
## Root genes (those with no non-Root parents) never cascade; they can only
## be removed by being the explicitly killed gene.
surviving_rows_hesbcn <- function(model, killed) {
    genes <- setdiff(unique(c(model$From, model$To)), "Root")
    dead <- killed
    changed <- TRUE
    while (changed) {
        changed <- FALSE
        for (g in genes) {
            if (g %in% dead) next
            pars <- model$From[model$To == g & model$From != "Root"]
            if (length(pars) == 0) next  ## root gene: never cascades
            alive_pars <- setdiff(pars, dead)
            rel <- model$Relation[model$To == g & model$From != "Root"][1]
            goes_dead <- if (rel == "AND") length(alive_pars) < length(pars)
                         else length(alive_pars) == 0
            if (goes_dead) {
                dead <- c(dead, g)
                changed <- TRUE
            }
        }
    }
    which(!(model$From %in% dead) & !(model$To %in% dead))
}

## Predicted genotype frequencies and hitting probs from a surviving HESBCN sub-model.
preds_from_hesbcn_model <- function(x) {
    if (nrow(x) == 0) {
        return(list(genot_freqs = c(WT = 1), hitting_probs_from_WT = c(WT = 1.0)))
    }
    tmp <- suppressWarnings(get_full_output(x))
    f <- tmp$HESBCN_predicted_genotype_freqs
    hp <- tmp$HESBCN_hitting_probs_from_WT
    list(genot_freqs = f[f > 0], hitting_probs_from_WT = filter_hp_keep_wt(hp))
}

## HESBCN model data frame (with Lambdas and Relation columns) ->
##   run four procedures + row-permuted model and assert all equal,
##   then compare procedure 1 against independent cascade ground truth
##   for every gene kill.
## Genes must be named sequentially A, B, C, ... (no gaps).
check_hesbcn_four_procedures <- function(model, label = "") {
    m <- model
    n <- length(setdiff(unique(c(m$From, m$To)), "Root"))
    genes <- setdiff(unique(c(m$From, m$To)), "Root")

    ## Procedure 1: standard kill_gene via DAG edge removal
    i1 <- intervene_cpm_every_gene(list(HESBCN_model = m), "HESBCN")

    ## Procedure 2: kill by setting parameters to 0
    i2 <- suppressWarnings(
        intervene_cpm_every_gene(list(HESBCN_model = m), "HESBCN",
                                 kill_gene_funct = kill_gene_by_params_to_0))

    ## Procedure 3: remove genotypes from pre-computed TRM.
    ## Requires a list with HESBCN_model and HESBCN_trans_rate_mat.
    ## get_full_output provides the trans_rate_mat; we add HESBCN_model
    ## explicitly (in case HESBCN_model_2_output does not include it).
    full_out <- suppressWarnings(get_full_output(m))
    full_out$HESBCN_model <- m
    i3 <- intervene_cpm_trm_rm_every_gene(full_out, "HESBCN")

    ## Procedure 4: fitness landscape derived from model.
    landscape_obj <- hesbcn_to_landscape_obj(m, n)
    i4 <- intervene_fitness_landscape_every_gene(landscape_obj)

    ## Procedure 1 with row-permuted model (order must not matter)
    m_perm <- m[sample(1:nrow(m)), ]
    i1_perm <- intervene_cpm_every_gene(list(HESBCN_model = m_perm), "HESBCN")

    lbl <- if (nchar(label) > 0) paste0(" [", label, "]") else ""
    expect_equal(i1, i2, label = paste0("i1 vs i2", lbl))
    expect_equal(i1, i3, label = paste0("i1 vs i3", lbl))
    expect_equal(i1, i4, label = paste0("i1 vs i4", lbl))
    expect_equal(i1, i1_perm, label = paste0("i1 vs i1_perm", lbl))

    ## Ground-truth: for each gene kill, compare procedure 1 against the
    ## surviving-row prediction from the independent cascade reimplementation.
    for (gene in genes) {
        sr <- surviving_rows_hesbcn(m, gene)
        expected <- preds_from_hesbcn_model(m[sr, ])
        expect_equal(get_interv(i1, paste0("I:", gene)), expected,
                     label = paste0("ground truth kill ", gene, lbl))
    }
}


### Test blocks: one per DAG

## DAG_2: multi-parent nodes H (5 parents: F,D,G,J,B),
##                            I (5 parents: F,D,G,J,B),
##                            A (4 parents: D,G,J,B)
## Rows:
##   1-5:   Root -> F, D, G, J, B  (Single)
##   6-10:  F,D,G,J,B -> H         (rel_H)
##   11-15: F,D,G,J,B -> I         (rel_I)
##   16-19: D,G,J,B -> A           (rel_A)
##   20-21: H -> C, H -> E         (Single)
## 3^3 = 27 combinations.

test_that("HESBCN DAG_2: 27 AND/OR/XOR combinations for H, I, A", {
    local_edition(3)
    relations <- c("AND", "OR", "XOR")
    combos <- expand.grid(rel_H = relations, rel_I = relations, rel_A = relations,
                          stringsAsFactors = FALSE)

    for (k in 1:nrow(combos)) {
        rel_H <- combos$rel_H[k]
        rel_I <- combos$rel_I[k]
        rel_A <- combos$rel_A[k]
        cat("\n DAG_2 combo", k, "of", nrow(combos),
            ": rel_H =", rel_H, ", rel_I =", rel_I, ", rel_A =", rel_A, "\n")

        model <- DAG_2
        colnames(model)[colnames(model) == "rerun_lambda"] <- "Lambdas"
        model$Relation <- c(rep("Single", 5),   ## Root -> F, D, G, J, B
                            rep(rel_H, 5),       ## F,D,G,J,B -> H
                            rep(rel_I, 5),       ## F,D,G,J,B -> I
                            rep(rel_A, 4),       ## D,G,J,B -> A
                            rep("Single", 2))    ## H -> C, H -> E

        check_hesbcn_four_procedures(model,
            label = paste0("DAG_2 k=", k, " rel_H=", rel_H,
                           " rel_I=", rel_I, " rel_A=", rel_A))
    }
})


## DAG_3: multi-parent node C only (4 parents: H,B,F,D)
## Rows:
##   1: Root -> I   (Single)
##   2: Root -> H   (Single)
##   3: I -> E      (Single)
##   4: I -> B      (Single)
##   5: H -> C      (rel_C)
##   6: E -> G      (Single)
##   7: B -> C      (rel_C)
##   8: G -> F      (Single)
##   9: G -> A      (Single)
##  10: F -> C      (rel_C)
##  11: A -> D      (Single)
##  12: D -> C      (rel_C)
## 3 combinations (single multi-parent node).

test_that("HESBCN DAG_3: 3 AND/OR/XOR combinations for C", {
    local_edition(3)
    for (rel_C in c("AND", "OR", "XOR")) {
        cat("\n DAG_3: rel_C =", rel_C, "\n")

        model <- DAG_3
        colnames(model)[colnames(model) == "rerun_lambda"] <- "Lambdas"
        model$Relation <- c("Single",  ##  1: Root -> I
                            "Single",  ##  2: Root -> H
                            "Single",  ##  3: I -> E
                            "Single",  ##  4: I -> B
                            rel_C,     ##  5: H -> C
                            "Single",  ##  6: E -> G
                            rel_C,     ##  7: B -> C
                            "Single",  ##  8: G -> F
                            "Single",  ##  9: G -> A
                            rel_C,     ## 10: F -> C
                            "Single",  ## 11: A -> D
                            rel_C)     ## 12: D -> C

        check_hesbcn_four_procedures(model,
            label = paste0("DAG_3 rel_C=", rel_C))
    }
})


## DAG_4: multi-parent nodes E (3 parents: C,A,B),
##                             G (3 parents: F,E,D),
##                             H (2 parents: I,G)
## Rows:
##   1:  Root -> C  (Single)
##   2:  Root -> A  (Single)
##   3:  Root -> B  (Single)
##   4:  Root -> I  (Single)
##   5:  C -> F     (Single)
##   6:  C -> E     (rel_E)
##   7:  A -> E     (rel_E)
##   8:  B -> E     (rel_E)
##   9:  B -> D     (Single)
##  10:  I -> H     (rel_H)
##  11:  F -> G     (rel_G)
##  12:  E -> G     (rel_G)
##  13:  D -> G     (rel_G)
##  14:  G -> H     (rel_H)
## 3^3 = 27 combinations.

test_that("HESBCN DAG_4: 27 AND/OR/XOR combinations for E, G, H", {
    local_edition(3)
    relations <- c("AND", "OR", "XOR")
    combos <- expand.grid(rel_E = relations, rel_G = relations, rel_H = relations,
                          stringsAsFactors = FALSE)

    for (k in 1:nrow(combos)) {
        rel_E_val <- combos$rel_E[k]
        rel_G_val <- combos$rel_G[k]
        rel_H_val <- combos$rel_H[k]
        cat("\n DAG_4 combo", k, "of", nrow(combos),
            ": rel_E =", rel_E_val, ", rel_G =", rel_G_val,
            ", rel_H =", rel_H_val, "\n")

        model <- DAG_4
        colnames(model)[colnames(model) == "rerun_lambda"] <- "Lambdas"
        model$Relation <- c("Single",        ##  1: Root -> C
                            "Single",        ##  2: Root -> A
                            "Single",        ##  3: Root -> B
                            "Single",        ##  4: Root -> I
                            "Single",        ##  5: C -> F
                            rel_E_val, rel_E_val, rel_E_val,  ##  6-8: C,A,B -> E
                            "Single",        ##  9: B -> D
                            rel_H_val,       ## 10: I -> H
                            rel_G_val, rel_G_val, rel_G_val,  ## 11-13: F,E,D -> G
                            rel_H_val)       ## 14: G -> H

        check_hesbcn_four_procedures(model,
            label = paste0("DAG_4 k=", k, " rel_E=", rel_E_val,
                           " rel_G=", rel_G_val, " rel_H=", rel_H_val))
    }
})


## DAG_5: multi-parent nodes C (4 parents: A,F,G,I),
##                             D (3 parents: A,F,G),
##                             B (2 parents: E,D)
## Rows:
##   1:  Root -> A  (Single)
##   2:  Root -> F  (Single)
##   3:  Root -> G  (Single)
##   4:  Root -> I  (Single)
##   5:  A -> C     (rel_C)
##   6:  F -> C     (rel_C)
##   7:  G -> C     (rel_C)
##   8:  I -> C     (rel_C)
##   9:  A -> D     (rel_D)
##  10:  F -> D     (rel_D)
##  11:  G -> D     (rel_D)
##  12:  C -> E     (Single)
##  13:  D -> H     (Single)
##  14:  E -> B     (rel_B)
##  15:  D -> B     (rel_B)
## 3^3 = 27 combinations.

test_that("HESBCN DAG_5: 27 AND/OR/XOR combinations for C, D, B", {
    local_edition(3)
    relations <- c("AND", "OR", "XOR")
    combos <- expand.grid(rel_C = relations, rel_D = relations, rel_B = relations,
                          stringsAsFactors = FALSE)

    for (k in 1:nrow(combos)) {
        rel_C_val <- combos$rel_C[k]
        rel_D_val <- combos$rel_D[k]
        rel_B_val <- combos$rel_B[k]
        cat("\n DAG_5 combo", k, "of", nrow(combos),
            ": rel_C =", rel_C_val, ", rel_D =", rel_D_val,
            ", rel_B =", rel_B_val, "\n")

        model <- DAG_5
        colnames(model)[colnames(model) == "rerun_lambda"] <- "Lambdas"
        model$Relation <- c("Single", "Single", "Single", "Single",  ##  1-4: Root -> A,F,G,I
                            rel_C_val, rel_C_val, rel_C_val, rel_C_val,  ##  5-8: A,F,G,I -> C
                            rel_D_val, rel_D_val, rel_D_val,             ##  9-11: A,F,G -> D
                            "Single",   ## 12: C -> E
                            "Single",   ## 13: D -> H
                            rel_B_val, rel_B_val)  ## 14-15: E,D -> B

        check_hesbcn_four_procedures(model,
            label = paste0("DAG_5 k=", k, " rel_C=", rel_C_val,
                           " rel_D=", rel_D_val, " rel_B=", rel_B_val))
    }
})


## DAG_6: multi-parent nodes A (2 parents: C,D),
##                             B (2 parents: F,G),
##                             E (2 parents: A,B)
## Rows:
##   1:  Root -> C  (Single)
##   2:  Root -> D  (Single)
##   3:  Root -> F  (Single)
##   4:  Root -> G  (Single)
##   5:  C -> A     (rel_A)
##   6:  D -> A     (rel_A)
##   7:  F -> B     (rel_B)
##   8:  G -> B     (rel_B)
##   9:  A -> E     (rel_E)
##  10:  B -> E     (rel_E)
## 3^3 = 27 combinations.

test_that("HESBCN DAG_6: 27 AND/OR/XOR combinations for A, B, E", {
    local_edition(3)
    relations <- c("AND", "OR", "XOR")
    combos <- expand.grid(rel_A = relations, rel_B = relations, rel_E = relations,
                          stringsAsFactors = FALSE)

    for (k in 1:nrow(combos)) {
        rel_A_val <- combos$rel_A[k]
        rel_B_val <- combos$rel_B[k]
        rel_E_val <- combos$rel_E[k]
        cat("\n DAG_6 combo", k, "of", nrow(combos),
            ": rel_A =", rel_A_val, ", rel_B =", rel_B_val,
            ", rel_E =", rel_E_val, "\n")

        model <- DAG_6
        colnames(model)[colnames(model) == "rerun_lambda"] <- "Lambdas"
        model$Relation <- c("Single", "Single", "Single", "Single",  ## 1-4: Root -> C,D,F,G
                            rel_A_val, rel_A_val,   ##  5-6: C,D -> A
                            rel_B_val, rel_B_val,   ##  7-8: F,G -> B
                            rel_E_val, rel_E_val)   ##  9-10: A,B -> E

        check_hesbcn_four_procedures(model,
            label = paste0("DAG_6 k=", k, " rel_A=", rel_A_val,
                           " rel_B=", rel_B_val, " rel_E=", rel_E_val))
    }
})


## DAG_7: multi-parent nodes A (3 parents: C,D,F),
##                             B (2 parents: D,F),
##                             E (2 parents: A,B)
## Rows:
##   1:  Root -> C  (Single)
##   2:  Root -> D  (Single)
##   3:  Root -> F  (Single)
##   4:  C -> A     (rel_A)
##   5:  D -> A     (rel_A)
##   6:  F -> A     (rel_A)
##   7:  D -> B     (rel_B)
##   8:  F -> B     (rel_B)
##   9:  A -> E     (rel_E)
##  10:  B -> E     (rel_E)
## 3^3 = 27 combinations.

test_that("HESBCN DAG_7: 27 AND/OR/XOR combinations for A, B, E", {
    local_edition(3)
    relations <- c("AND", "OR", "XOR")
    combos <- expand.grid(rel_A = relations, rel_B = relations, rel_E = relations,
                          stringsAsFactors = FALSE)

    for (k in 1:nrow(combos)) {
        rel_A_val <- combos$rel_A[k]
        rel_B_val <- combos$rel_B[k]
        rel_E_val <- combos$rel_E[k]
        cat("\n DAG_7 combo", k, "of", nrow(combos),
            ": rel_A =", rel_A_val, ", rel_B =", rel_B_val,
            ", rel_E =", rel_E_val, "\n")

        model <- DAG_7
        colnames(model)[colnames(model) == "rerun_lambda"] <- "Lambdas"
        model$Relation <- c("Single", "Single", "Single",  ## 1-3: Root -> C,D,F
                            rel_A_val, rel_A_val, rel_A_val,  ## 4-6: C,D,F -> A
                            rel_B_val, rel_B_val,             ## 7-8: D,F -> B
                            rel_E_val, rel_E_val)             ## 9-10: A,B -> E

        check_hesbcn_four_procedures(model,
            label = paste0("DAG_7 k=", k, " rel_A=", rel_A_val,
                           " rel_B=", rel_B_val, " rel_E=", rel_E_val))
    }
})


### Cascade tests with ground-truth comparison
##
## Model: Root->A (Single), Root->B (Single), A->C (rel), B->C (rel), C->D (Single)
## Rows: 1=Root->A, 2=Root->B, 3=A->C, 4=B->C, 5=C->D
##
## AND case:
##   Kill A: C requires both A and B; without A, C unreachable, D cascades.
##           Surviving rows: 2 only (Root->B).
##   Kill B: symmetric. Surviving rows: 1 only (Root->A).
##   Kill C: D cascades (only parent is C). Surviving rows: 1,2.
##   Kill D: leaf, no cascade. Surviving rows: 1,2,3,4.
##
## OR/XOR case:
##   Kill A: C still reachable via B alone; no cascade to C or D.
##           Surviving rows: 2,4,5 (Root->B, B->C, C->D).
##   Kill B: symmetric. Surviving rows: 1,3,5.
##   Kill C: D cascades (only parent is C). Surviving rows: 1,2.
##   Kill D: leaf. Surviving rows: 1,2,3,4.
##
## Surviving models specified by row index, not by kill_gene, to avoid
## circularity.

test_that("HESBCN cascade AND: killing one AND-parent removes child and its descendants", {
    local_edition(3)

    m <- data.frame(
        From     = c("Root", "Root", "A", "B", "C"),
        To       = c("A",    "B",    "C", "C", "D"),
        Lambdas  = c(2,      0.8,    1,   1,   0.5),
        Relation = c("Single", "Single", "AND", "AND", "Single"),
        stringsAsFactors = FALSE
    )

    ## predicted genotype freqs from a hand-specified surviving model
    preds_from_model <- function(x) {
        tmp <- suppressWarnings(get_full_output(x))
        method <- if (!is.null(tmp$HESBCN_predicted_genotype_freqs)) "HESBCN" else "CBN"
        f <- tmp[[paste0(method, "_predicted_genotype_freqs")]]
        hp <- tmp[[paste0(method, "_hitting_probs_from_WT")]]
        list(genot_freqs = f[f > 0], hitting_probs_from_WT = filter_hp_keep_wt(hp))
    }

    res <- intervene_cpm_every_gene(list(HESBCN_model = m), "HESBCN")

    ## Kill A: surviving rows 2 only
    expect_equal(get_interv(res, "I:A"), preds_from_model(m[2, ]))

    ## Kill B: surviving rows 1 only
    expect_equal(get_interv(res, "I:B"), preds_from_model(m[1, ]))

    ## Kill C: surviving rows 1,2
    expect_equal(get_interv(res, "I:C"), preds_from_model(m[c(1, 2), ]))

    ## Kill D: surviving rows 1,2,3,4
    expect_equal(get_interv(res, "I:D"), preds_from_model(m[c(1, 2, 3, 4), ]))
})


test_that("HESBCN cascade OR/XOR: killing one OR/XOR-parent does NOT cascade to child", {
    local_edition(3)

    for (rel in c("OR", "XOR")) {
        m <- data.frame(
            From     = c("Root", "Root", "A", "B", "C"),
            To       = c("A",    "B",    "C", "C", "D"),
            Lambdas  = c(2,      0.8,    1,   1,   0.5),
            Relation = c("Single", "Single", rel, rel, "Single"),
            stringsAsFactors = FALSE
        )

        preds_from_model <- function(x) {
            tmp <- suppressWarnings(get_full_output(x))
            method <- if (!is.null(tmp$HESBCN_predicted_genotype_freqs)) "HESBCN" else "CBN"
            f <- tmp[[paste0(method, "_predicted_genotype_freqs")]]
            hp <- tmp[[paste0(method, "_hitting_probs_from_WT")]]
            list(genot_freqs = f[f > 0], hitting_probs_from_WT = filter_hp_keep_wt(hp))
        }

        res <- intervene_cpm_every_gene(list(HESBCN_model = m), "HESBCN")

        ## Kill A: C still reachable via B; surviving rows 2,4,5
        expect_equal(get_interv(res, "I:A"), preds_from_model(m[c(2, 4, 5), ]))

        ## Kill B: symmetric; surviving rows 1,3,5
        expect_equal(get_interv(res, "I:B"), preds_from_model(m[c(1, 3, 5), ]))

        ## Kill C: D cascades; surviving rows 1,2
        expect_equal(get_interv(res, "I:C"), preds_from_model(m[c(1, 2), ]))

        ## Kill D: leaf; surviving rows 1,2,3,4
        expect_equal(get_interv(res, "I:D"), preds_from_model(m[c(1, 2, 3, 4), ]))
    }
})

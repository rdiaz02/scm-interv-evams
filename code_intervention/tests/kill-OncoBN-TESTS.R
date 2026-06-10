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

## First, some minimal tests that killing gives identical
## structures to hand removal. Then, thorough tests of OncoBN
## interventions using six hand-crafted DAG structures
## (DAG_2 through DAG_7; same topologies as in kill-CBN-DAG-2-7-TESTS.R).
## where we compare the key output for us: probabilities of genotypes

## OncoBN models have two variants per DAG:
##   - all multi-parent relations AND  (equivalent to CBN mode in OncoBN)
##   - all multi-parent relations OR   (equivalent to DBN mode in OncoBN)
## Total: 6 x 2 = 12 OncoBN models.
##
## For each model, two procedures are compared:
##   1. Standard kill_gene (modifies DAG edges / zeroes affected rows)
##   2. kill_gene_by_params_to_0 (sets theta to 0 explicitly)
##
## OncoBN has no TRM-based procedure and no fitness landscape procedure,
## so only these two are available.
##
## Plus: procedure 1 with a row-permuted model must equal the original.
##
## Additional ground-truth tests compare the intervention output against
## manually identified surviving model rows (computed by hand).  The key
## difference between AND and OR cascades is:
##   AND: killing any one parent of a multi-parent AND-node makes that node
##        unreachable; all of its incoming edges cascade away, then downstream.
##   OR:  killing one parent of a multi-parent OR-node does NOT make it
##        unreachable (other parents still live); only that parent's direct
##        edges are removed; downstream is unaffected unless the node loses
##        ALL its parents.
##
## Gene names are sequential A, B, C, ... (same mapping as in
## kill-CBN-DAG-TESTS.R).
##
## theta values are derived from the CBN lambda values of the same DAGs
## via round(lambda/4 + 0.10, 2), mapping the [0.1, 2.7] lambda range
## to approximately [0.13, 0.78].  Multi-parent edges into the same child
## node always share the same theta value.

## options(intervention_every_gene_cores = parallel::detectCores())
library(testthat)

pwd <- getwd()
setwd("../")
source("intervention.R")
setwd(pwd)

set.seed(NULL)

## Extract a single named intervention from the full intervention output as a
## list(genot_freqs, hitting_probs_from_WT), for direct comparison with
## preds_from_oncobn.
get_interv <- function(res, name) {
    return(res[[name]])
}


test_that("OncoBN, conjuntive", {
    local_edition(3)
    m1 <- data.frame(From = c("Root", "Root", "B", "B", "C", "C"),
                     To   = c("A",    "B",    "C", "D",  "E", "F"),
                     theta = 0.5,
                     Relation = c("Single"))

    m2 <- data.frame(From = c("Root", "Root", "A", "B", "B", "C", "C"),
                     To   = c("A",    "B",    "C", "C", "D",  "E", "F"),
                     theta = 0.5,
                     Relation = c("Single", "Single", "AND", "AND", rep("Single", 3)))

    m3 <- data.frame(From = c("Root", "Root", "A", "B", "B", "C", "C"),
                     To   = c("A",    "B",    "E", "C", "D",  "E", "F"),
                     theta = 0.5,
                     Relation = c("Single", "Single", "AND", "Single", "Single",
                                  "AND", "Single"))

    kills_m1 <- lapply(LETTERS[1:6], function(v) kill_gene(m1, v))
    kills_m2 <- lapply(LETTERS[1:6], function(v) kill_gene(m2, v))
    kills_m3 <- lapply(LETTERS[1:6], function(v) kill_gene(m3, v))

    names(kills_m1) <- names(kills_m2) <- names(kills_m3) <- LETTERS[1:6]

    expect_equal(kills_m1[["A"]], m1[-c(1), ])
    expect_equal(kills_m1[["B"]], m1[1, ])
    expect_equal(kills_m1[["C"]], m1[-c(3, 5, 6), ])
    expect_equal(kills_m1[["D"]], m1[-c(4), ])
    expect_equal(kills_m1[["E"]], m1[-c(5), ])
    expect_equal(kills_m1[["F"]], m1[-c(6), ])

    expect_equal(kills_m2[["A"]], m2[-c(1, 3, 4, 6, 7), ])
    expect_equal(kills_m2[["B"]], m2[c(1), ])
    expect_equal(kills_m2[["C"]], m2[-c(3, 4, 6, 7), ])
    expect_equal(kills_m2[["D"]], m2[-c(5), ])
    expect_equal(kills_m2[["E"]], m2[-c(6), ])
    expect_equal(kills_m2[["F"]], m2[-c(7), ])

    expect_equal(kills_m3[["A"]], m3[-c(1, 3, 6), ])
    expect_equal(kills_m3[["B"]], m3[1, ])
    expect_equal(kills_m3[["C"]], m3[-c(3, 4, 6, 7), ])
    expect_equal(kills_m3[["D"]], m3[-c(5), ])
    expect_equal(kills_m3[["E"]], m3[-c(3, 6), ])
    expect_equal(kills_m3[["F"]], m3[-c(7), ])
})



test_that("OncoBN, disjunctive", {
    local_edition(3)
    m2d <- data.frame(From = c("Root", "Root", "A", "B", "B", "C", "C"),
                      To   = c("A",    "B",    "C", "C", "D",  "E", "F"),
                      theta = 0.5,
                      Relation = c("Single", "Single", "OR", "OR",
                                   rep("Single", 3)))

    m3d <- data.frame(From = c("Root", "Root", "A", "B", "B", "C", "C"),
                      To   = c("A",    "B",    "E", "C", "D",  "E", "F"),
                      theta = 0.5,
                      Relation = c("Single", "Single", "OR", "Single", "Single",
                                   "OR", "Single"))

    m4d <- data.frame(From = c("Root", "Root", "A", "B", "B", "C", "D"),
                      To   = c("A",    "B",    "C", "C", "D",  "E", "F"),
                      theta = 0.5,
                      Relation = c("Single", "Single", "OR", "OR", "Single",
                                   "Single", "Single"))

    ## Next three: permute and change types of joint dep
    m5d <- data.frame(From = c("Root", "Root", "A", "B", "Root", "C", "D", "E"),
                      To   = c("A",    "B",     "C", "C", "D",    "E", "E", "F"),
                      theta = 0.5,
                      Relation = c("Single", "Single", "AND", "AND", "Single",
                                   "OR", "OR", "Single"))

    m6d <- data.frame(From = c("Root", "Root", "A", "B", "Root", "C", "D", "E"),
                      To   = c("A",    "B",     "C", "C", "D",    "E", "E", "F"),
                      theta = 0.5,
                      Relation = c("Single", "Single", "OR", "OR", "Single",
                                   "AND", "AND", "Single"))

    m7d <- data.frame(From = c("Root", "Root", "A", "B", "Root", "C", "D", "E"),
                      To   = c("A",    "B",     "C", "C", "D",    "E", "E", "F"),
                      theta = 0.5,
                      Relation = c("Single", "Single", "OR", "OR", "Single",
                                   "OR", "OR", "Single"))

    kills_m2d <- lapply(LETTERS[1:6], function(v) kill_gene(m2d, v))
    kills_m3d <- lapply(LETTERS[1:6], function(v) kill_gene(m3d, v))
    kills_m4d <- lapply(LETTERS[1:6], function(v) kill_gene(m4d, v))
    kills_m5d <- lapply(LETTERS[1:6], function(v) kill_gene(m5d, v))
    kills_m6d <- lapply(LETTERS[1:6], function(v) kill_gene(m6d, v))
    kills_m7d <- lapply(LETTERS[1:6], function(v) kill_gene(m7d, v))

    names(kills_m2d) <- names(kills_m3d) <- names(kills_m4d) <-
        names(kills_m5d) <- names(kills_m6d) <- names(kills_m7d) <- LETTERS[1:6]

    expect_equal(kills_m2d[["A"]], m2d[-c(1, 3), ])
    expect_equal(kills_m2d[["B"]], m2d[-c(2, 4, 5), ])
    expect_equal(kills_m2d[["C"]], m2d[-c(3, 4, 6, 7), ])
    expect_equal(kills_m2d[["D"]], m2d[-c(5), ])
    expect_equal(kills_m2d[["E"]], m2d[-c(6), ])
    expect_equal(kills_m2d[["F"]], m2d[-c(7), ])

    expect_equal(kills_m3d[["A"]], m3d[-c(1, 3), ])
    expect_equal(kills_m3d[["B"]], m3d[-c(2, 4, 5, 6, 7), ])
    expect_equal(kills_m3d[["C"]], m3d[-c(4, 6, 7), ])
    expect_equal(kills_m3d[["D"]], m3d[-c(5), ])
    expect_equal(kills_m3d[["E"]], m3d[-c(3, 6), ])
    expect_equal(kills_m3d[["F"]], m3d[-c(7), ])

    expect_equal(kills_m4d[["A"]], m4d[-c(1, 3), ])
    expect_equal(kills_m4d[["B"]], m4d[-c(2, 4, 5, 7), ])
    expect_equal(kills_m4d[["C"]], m4d[-c(3, 4, 6), ])
    expect_equal(kills_m4d[["D"]], m4d[-c(5, 7), ])
    expect_equal(kills_m4d[["E"]], m4d[-c(6), ])
    expect_equal(kills_m4d[["F"]], m4d[-c(7), ])

    expect_equal(kills_m5d[["A"]], m5d[-c(1, 3, 4, 6), ])
    expect_equal(kills_m5d[["B"]], m5d[-c(2, 4, 3, 6), ])
    expect_equal(kills_m5d[["C"]], m5d[-c(3, 4, 6), ])
    expect_equal(kills_m5d[["D"]], m5d[-c(5, 7), ])
    expect_equal(kills_m5d[["E"]], m5d[-c(6, 7, 8), ])
    expect_equal(kills_m5d[["F"]], m5d[-c(8), ])

    expect_equal(kills_m6d[["A"]], m6d[-c(1, 3), ])
    expect_equal(kills_m6d[["B"]], m6d[-c(2, 4), ])
    expect_equal(kills_m6d[["C"]], m6d[-c(3, 4, 6, 7, 8), ])
    expect_equal(kills_m6d[["D"]], m6d[-c(5, 6, 7, 8), ])
    expect_equal(kills_m6d[["E"]], m6d[-c(6, 7, 8), ])
    expect_equal(kills_m6d[["F"]], m6d[-c(8), ])

    expect_equal(kills_m7d[["A"]], m7d[-c(1, 3), ])
    expect_equal(kills_m7d[["B"]], m7d[-c(2, 4), ])
    expect_equal(kills_m7d[["C"]], m7d[-c(3, 4, 6), ])
    expect_equal(kills_m7d[["D"]], m7d[-c(5, 7), ])
    expect_equal(kills_m7d[["E"]], m7d[-c(6:8), ])
    expect_equal(kills_m7d[["F"]], m7d[-c(8), ])
})


## DAG topologies identical to kill-CBN-DAG-TESTS.R.
## theta column: round(rerun_lambda / 4 + 0.10, 2).
## Relation column is NOT included here; it is added inside each test_that
## block by the loop over AND / OR.

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
    theta = c(.60, .35, .23, .19, .30,
              .78, .78, .78, .78, .78,
              .45, .45, .45, .45, .45,
              .25, .25, .25, .25,
              .33, .75),
    stringsAsFactors = FALSE
)

## DAG_3: 9 genes (A-I)
## Multi-parent node: C (4 parents: H,B,F,D)
## Root genes: I, H
DAG_3 <- data.frame(
    From  = c("Root", "Root", "I",  "I",  "H",  "E",  "B",  "G",  "G",  "F",  "A",  "D"),
    To    = c("I",    "H",    "E",  "B",  "C",  "G",  "C",  "F",  "A",  "C",  "D",  "C"),
    theta = c(.19,    .30,    .28,  .25,  .23,  .33,  .23,  .35,  .60,  .23,  .20,  .23),
    stringsAsFactors = FALSE
)

## DAG_4: 9 genes (A-I)
## Multi-parent nodes: E (3 parents: C,A,B),
##                     G (3 parents: F,E,D),
##                     H (2 parents: I,G)
## Root genes: C, A, B, I
DAG_4 <- data.frame(
    From  = c("Root", "Root", "Root", "Root",
              "C",  "C",  "A",  "B",  "B",  "I",
              "F",  "E",  "D",  "G"),
    To    = c("C",  "A",  "B",  "I",
              "F",  "E",  "E",  "E",  "D",  "H",
              "G",  "G",  "G",  "H"),
    theta = c(.30,  .60,  .19,  .20,
              .23,  .35,  .35,  .35,  .25,  .33,
              .28,  .28,  .28,  .33),
    stringsAsFactors = FALSE
)

## DAG_5: 9 genes (A-I)
## Multi-parent nodes: C (4 parents: A,F,G,I),
##                     D (3 parents: A,F,G),
##                     B (2 parents: E,D)
## Root genes: A, F, G, I
DAG_5 <- data.frame(
    From  = c("Root", "Root", "Root", "Root",
              "A",  "F",  "G",  "I",
              "A",  "F",  "G",
              "C",  "D",  "E",  "D"),
    To    = c("A",  "F",  "G",  "I",
              "C",  "C",  "C",  "C",
              "D",  "D",  "D",
              "E",  "H",  "B",  "B"),
    theta = c(.60,  .28,  .30,  .13,
              .20,  .20,  .20,  .20,
              .23,  .23,  .23,
              .25,  .33,  .68,  .68),
    stringsAsFactors = FALSE
)

## DAG_6: 7 genes (A-G)
## Multi-parent nodes: A (2 parents: C,D),
##                     B (2 parents: F,G),
##                     E (2 parents: A,B)
## Root genes: C, D, F, G
DAG_6 <- data.frame(
    From  = c(rep("Root", 4),    "C",  "D",  "F",  "G",  "A",  "B"),
    To    = c("C",  "D",  "F",  "G",  "A",  "A",  "B",  "B",  "E",  "E"),
    theta = c(.35,  .60,  .68,  .20,  .23,  .23,  .25,  .25,  .28,  .28),
    stringsAsFactors = FALSE
)

## DAG_7: 6 genes (A-F)
## Multi-parent nodes: A (3 parents: C,D,F),
##                     B (2 parents: D,F),
##                     E (2 parents: A,B)
## Root genes: C, D, F
DAG_7 <- data.frame(
    From  = c(rep("Root", 3),    "C",  "D",  "F",  "D",  "F",  "A",  "B"),
    To    = c("C",  "D",  "F",   "A",  "A",  "A",  "B",  "B",  "E",  "E"),
    theta = c(.35,  .60,  .68,   .20,  .20,  .20,  .23,  .23,  .25,  .25),
    stringsAsFactors = FALSE
)

## You can see what they look like by doing
## evamtools:::DAG_plot_graphAM(DAG_2, "DAG_2")
## evamtools:::DAG_plot_graphAM(DAG_3, "DAG_3")
## evamtools:::DAG_plot_graphAM(DAG_4, "DAG_4")
## evamtools:::DAG_plot_graphAM(DAG_5, "DAG_5")
## evamtools:::DAG_plot_graphAM(DAG_6, "DAG_6")
## evamtools:::DAG_plot_graphAM(DAG_7, "DAG_7")





#### Utility functions

## OncoBN model -> predicted genotype frequencies (positive entries only).
## Used for ground-truth row-removal comparisons.
preds_from_oncobn <- function(x) {
    if (nrow(x) == 0) {
        return(list(genot_freqs = c(WT = 1), hitting_probs_from_WT = c(WT = 1.0)))
    }
    tmp <- suppressWarnings(get_full_output(x))
    f <- tmp$OncoBN_predicted_genotype_freqs
    hp <- tmp$OncoBN_hitting_probs_from_WT
    list(genot_freqs = f[f > 0], hitting_probs_from_WT = filter_hp_keep_wt(hp))
}

## OncoBN model data frame (with theta and Relation columns) ->
##   run two procedures + row-permuted model and assert all equal.
## Genes must be named sequentially A, B, C, ... (no gaps).
check_oncobn_two_procedures <- function(model, label = "") {
    cat("\n OncoBN", if (nchar(label) > 0) label else "",
        ": n =", length(setdiff(unique(c(model$From, model$To)), "Root")),
        "genes\n")

    ## Procedure 1: standard kill_gene
    i1 <- intervene_cpm_every_gene(list(OncoBN_model = model), "OncoBN")

    ## Procedure 2: kill by setting parameters to 0
    i2 <- suppressWarnings(
        intervene_cpm_every_gene(list(OncoBN_model = model), "OncoBN",
                                 kill_gene_funct = kill_gene_by_params_to_0))

    ## Procedure 1 with row-permuted model (order must not matter)
    model_perm <- model[sample(1:nrow(model)), ]
    i1_perm <- intervene_cpm_every_gene(list(OncoBN_model = model_perm), "OncoBN")

    lbl <- if (nchar(label) > 0) paste0(" [", label, "]") else ""
    expect_equal(i1, i2,      label = paste0("i1 vs i2",      lbl))
    expect_equal(i1, i1_perm, label = paste0("i1 vs i1_perm", lbl))
}


### Test blocks: one per DAG, looping over AND / OR
##
## DAG structures (multi-parent nodes):
##   DAG_2: H (5 parents: F,D,G,J,B), I (5 parents: F,D,G,J,B), A (4 parents: D,G,J,B)
##   DAG_3: C (4 parents: H,B,F,D)
##   DAG_4: E (3 parents: C,A,B), G (3 parents: F,E,D), H (2 parents: I,G)
##   DAG_5: C (4 parents: A,F,G,I), D (3 parents: A,F,G), B (2 parents: E,D)
##   DAG_6: A (2 parents: C,D), B (2 parents: F,G), E (2 parents: A,B)
##   DAG_7: A (3 parents: C,D,F), B (2 parents: D,F), E (2 parents: A,B)

## DAG_2 rows:
##   1-5:   Root -> F,D,G,J,B   Single
##   6-10:  F,D,G,J,B -> H      rel
##   11-15: F,D,G,J,B -> I      rel
##   16-19: D,G,J,B -> A        rel
##   20-21: H -> C, H -> E      Single

test_that("OncoBN DAG_2: AND and OR variants", {
    local_edition(3)
    for (rel in c("AND", "OR")) {
        cat("\n DAG_2:", rel, "\n")
        model <- DAG_2
        model$Relation <- c(rep("Single", 5),
                            rep(rel, 5),    ## F,D,G,J,B -> H
                            rep(rel, 5),    ## F,D,G,J,B -> I
                            rep(rel, 4),    ## D,G,J,B -> A
                            rep("Single", 2))
        check_oncobn_two_procedures(model, label = paste("DAG_2", rel))
    }
})


## DAG_3 rows:
##  1: Root -> I    Single
##  2: Root -> H    Single
##  3: I -> E       Single
##  4: I -> B       Single
##  5: H -> C       rel   (C parents: H,B,F,D)
##  6: E -> G       Single
##  7: B -> C       rel
##  8: G -> F       Single
##  9: G -> A       Single
## 10: F -> C       rel
## 11: A -> D       Single
## 12: D -> C       rel

test_that("OncoBN DAG_3: AND and OR variants", {
    local_edition(3)
    for (rel in c("AND", "OR")) {
        cat("\n DAG_3:", rel, "\n")
        model <- DAG_3
        model$Relation <- c("Single",  ##  1: Root -> I
                            "Single",  ##  2: Root -> H
                            "Single",  ##  3: I -> E
                            "Single",  ##  4: I -> B
                            rel,       ##  5: H -> C
                            "Single",  ##  6: E -> G
                            rel,       ##  7: B -> C
                            "Single",  ##  8: G -> F
                            "Single",  ##  9: G -> A
                            rel,       ## 10: F -> C
                            "Single",  ## 11: A -> D
                            rel)       ## 12: D -> C
        check_oncobn_two_procedures(model, label = paste("DAG_3", rel))
    }
})


## DAG_4 rows:
##  1: Root -> C    Single
##  2: Root -> A    Single
##  3: Root -> B    Single
##  4: Root -> I    Single
##  5: C -> F       Single
##  6: C -> E       rel   (E parents: C,A,B)
##  7: A -> E       rel
##  8: B -> E       rel
##  9: B -> D       Single
## 10: I -> H       rel   (H parents: I,G)
## 11: F -> G       rel   (G parents: F,E,D)
## 12: E -> G       rel
## 13: D -> G       rel
## 14: G -> H       rel

test_that("OncoBN DAG_4: AND and OR variants", {
    local_edition(3)
    for (rel in c("AND", "OR")) {
        cat("\n DAG_4:", rel, "\n")
        model <- DAG_4
        model$Relation <- c("Single",        ##  1: Root -> C
                            "Single",        ##  2: Root -> A
                            "Single",        ##  3: Root -> B
                            "Single",        ##  4: Root -> I
                            "Single",        ##  5: C -> F
                            rel, rel, rel,   ##  6-8: C,A,B -> E
                            "Single",        ##  9: B -> D
                            rel,             ## 10: I -> H
                            rel, rel, rel,   ## 11-13: F,E,D -> G
                            rel)             ## 14: G -> H
        check_oncobn_two_procedures(model, label = paste("DAG_4", rel))
    }
})


## DAG_5 rows:
##  1-4: Root -> A,F,G,I   Single
##  5-8: A,F,G,I -> C      rel   (C parents: A,F,G,I)
##  9-11: A,F,G -> D       rel   (D parents: A,F,G)
## 12: C -> E              Single
## 13: D -> H              Single
## 14: E -> B              rel   (B parents: E,D)
## 15: D -> B              rel

test_that("OncoBN DAG_5: AND and OR variants", {
    local_edition(3)
    for (rel in c("AND", "OR")) {
        cat("\n DAG_5:", rel, "\n")
        model <- DAG_5
        model$Relation <- c(rep("Single", 4),   ##  1-4: Root -> A,F,G,I
                            rep(rel, 4),         ##  5-8: A,F,G,I -> C
                            rep(rel, 3),         ##  9-11: A,F,G -> D
                            "Single",            ## 12: C -> E
                            "Single",            ## 13: D -> H
                            rel, rel)            ## 14-15: E,D -> B
        check_oncobn_two_procedures(model, label = paste("DAG_5", rel))
    }
})


## DAG_6 rows:
##  1-4: Root -> C,D,F,G   Single
##  5-6: C,D -> A          rel   (A parents: C,D)
##  7-8: F,G -> B          rel   (B parents: F,G)
##  9-10: A,B -> E         rel   (E parents: A,B)

test_that("OncoBN DAG_6: AND and OR variants", {
    local_edition(3)
    for (rel in c("AND", "OR")) {
        cat("\n DAG_6:", rel, "\n")
        model <- DAG_6
        model$Relation <- c(rep("Single", 4),   ## 1-4: Root -> C,D,F,G
                            rel, rel,            ## 5-6: C,D -> A
                            rel, rel,            ## 7-8: F,G -> B
                            rel, rel)            ## 9-10: A,B -> E
        check_oncobn_two_procedures(model, label = paste("DAG_6", rel))
    }
})


## DAG_7 rows:
##  1-3: Root -> C,D,F     Single
##  4-6: C,D,F -> A        rel   (A parents: C,D,F)
##  7-8: D,F -> B          rel   (B parents: D,F)
##  9-10: A,B -> E         rel   (E parents: A,B)

test_that("OncoBN DAG_7: AND and OR variants", {
    local_edition(3)
    for (rel in c("AND", "OR")) {
        cat("\n DAG_7:", rel, "\n")
        model <- DAG_7
        model$Relation <- c(rep("Single", 3),   ## 1-3: Root -> C,D,F
                            rel, rel, rel,       ## 4-6: C,D,F -> A
                            rel, rel,            ## 7-8: D,F -> B
                            rel, rel)            ## 9-10: A,B -> E
        check_oncobn_two_procedures(model, label = paste("DAG_7", rel))
    }
})


### Cascade tests with ground-truth row-removal comparison
##
## Model: Root->A(1), Root->B(2), A->C(3), B->C(4), C->D(5)
## C has two parents A and B.
##
## AND mode — kill one parent cascades through C to D:
##   Kill A: C unreachable (needs A+B); D cascades.
##           Remove rows 1(Root->A), 3(A->C), 4(B->C cascade), 5(C->D cascade).
##           Surviving: row 2.
##   Kill B: symmetric. Surviving: row 1.
##   Kill C: D cascades. Remove rows 3,4,5. Surviving: rows 1,2.
##   Kill D: leaf. Remove row 5. Surviving: rows 1,2,3,4.
##
## OR mode — kill one parent does NOT cascade through C (C still reachable):
##   Kill A: removes rows 1(Root->A) and 3(A->C-OR). C alive via B. D ok.
##           Surviving: rows 2,4,5.
##   Kill B: removes rows 2,4. C alive via A. D ok.
##           Surviving: rows 1,3,5.
##   Kill C: C is the only parent of D (Single). D cascades. Remove rows 3,4,5.
##           Surviving: rows 1,2.
##   Kill D: leaf. Remove row 5. Surviving: rows 1,2,3,4.

test_that("OncoBN cascade: 2-parent AND vs OR", {
    local_edition(3)

    m_base <- data.frame(
        From  = c("Root", "Root", "A",  "B",  "C"),
        To    = c("A",    "B",    "C",  "C",  "D"),
        theta = c(.60,    .40,    .50,  .50,  .40),
        stringsAsFactors = FALSE
    )

    for (rel in c("AND", "OR")) {
        m <- m_base
        m$Relation <- c("Single", "Single", rel, rel, "Single")

        res <- intervene_cpm_every_gene(list(OncoBN_model = m), "OncoBN")

        if (rel == "AND") {
            ## Kill A: C unreachable, D cascades. Surviving: row 2.
            expect_equal(get_interv(res, "I:A"), preds_from_oncobn(m[2, ]),
                         label = "AND kill A")
            ## Kill B: symmetric. Surviving: row 1.
            expect_equal(get_interv(res, "I:B"), preds_from_oncobn(m[1, ]),
                         label = "AND kill B")
            ## Kill C: D cascades. Surviving: rows 1,2.
            expect_equal(get_interv(res, "I:C"), preds_from_oncobn(m[c(1, 2), ]),
                         label = "AND kill C")
            ## Kill D: leaf. Surviving: rows 1,2,3,4.
            expect_equal(get_interv(res, "I:D"), preds_from_oncobn(m[c(1, 2, 3, 4), ]),
                         label = "AND kill D")
        } else {
            ## Kill A: C alive via B. Surviving: rows 2,4,5.
            expect_equal(get_interv(res, "I:A"), preds_from_oncobn(m[c(2, 4, 5), ]),
                         label = "OR kill A")
            ## Kill B: C alive via A. Surviving: rows 1,3,5.
            expect_equal(get_interv(res, "I:B"), preds_from_oncobn(m[c(1, 3, 5), ]),
                         label = "OR kill B")
            ## Kill C: D cascades (C is D's sole parent). Surviving: rows 1,2.
            expect_equal(get_interv(res, "I:C"), preds_from_oncobn(m[c(1, 2), ]),
                         label = "OR kill C")
            ## Kill D: leaf. Surviving: rows 1,2,3,4.
            expect_equal(get_interv(res, "I:D"), preds_from_oncobn(m[c(1, 2, 3, 4), ]),
                         label = "OR kill D")
        }

        ## Also verify both procedures agree on this model
        check_oncobn_two_procedures(m, label = paste("cascade 2-parent", rel))
    }
})


### Chain cascade test
##
## Model: Root->A(1), A->B(2), B->C(3)
## Every node has exactly one parent, so AND vs OR makes no difference.
##
##   Kill A: B and C both unreachable (chain). No rows survive -> WT = 1.
##   Kill B: C cascades. Surviving: row 1 (Root->A).
##   Kill C: leaf. Surviving: rows 1,2.

test_that("OncoBN cascade: chain (single-parent; AND and OR identical)", {
    local_edition(3)

    m_base <- data.frame(
        From  = c("Root", "A",   "B"),
        To    = c("A",    "B",   "C"),
        theta = c(.60,    .40,   .40),
        stringsAsFactors = FALSE
    )

    for (rel in c("AND", "OR")) {
        m <- m_base
        ## All edges are single-parent; Relation = "Single" in both cases.
        m$Relation <- rep("Single", 3)

        res <- intervene_cpm_every_gene(list(OncoBN_model = m), "OncoBN")

        ## Kill A: full chain cascade -> WT only.
        expect_equal(get_interv(res, "I:A"), preds_from_oncobn(m[integer(0), ]),
                     label = paste(rel, "chain kill A"))
        ## Kill B: C cascades. Surviving: row 1.
        expect_equal(get_interv(res, "I:B"), preds_from_oncobn(m[1, ]),
                     label = paste(rel, "chain kill B"))
        ## Kill C: leaf. Surviving: rows 1,2.
        expect_equal(get_interv(res, "I:C"), preds_from_oncobn(m[c(1, 2), ]),
                     label = paste(rel, "chain kill C"))

        check_oncobn_two_procedures(m, label = paste("chain", rel))
    }
})


### No-over-cascade test
##
## Model: Root->A(1), Root->B(2), A->C(3), B->C(4), B->D(5)
## C has two parents A and B.  D has one parent B.
##
## AND mode:
##   Kill A: C unreachable (needs A+B). D ok (only needs B).
##           Remove rows 1,3,4(cascade A->C and B->C since AND).
##           Surviving: rows 2,5.
##   Kill B: C unreachable (needs B). D unreachable (only needs B).
##           Remove rows 2,3,4,5. Surviving: row 1.
##   Kill C: D unaffected. Remove rows 3,4. Surviving: rows 1,2,5.
##   Kill D: leaf. Remove row 5. Surviving: rows 1,2,3,4.
##
## OR mode:
##   Kill A: C alive via B. D unaffected. Remove rows 1,3.
##           Surviving: rows 2,4,5.
##   Kill B: C alive via A. D unreachable (B was sole parent). Remove rows 2,4,5.
##           Surviving: rows 1,3.
##   Kill C: D unaffected. Remove rows 3,4. Surviving: rows 1,2,5.
##   Kill D: leaf. Remove row 5. Surviving: rows 1,2,3,4.
##
## The key check: kill A with AND removes C (cascade stops at C because D's
## parent B is still alive).  With OR, kill A leaves C alive entirely.

test_that("OncoBN no-over-cascade: AND vs OR", {
    local_edition(3)

    m_base <- data.frame(
        From  = c("Root", "Root", "A",  "B",  "B"),
        To    = c("A",    "B",    "C",  "C",  "D"),
        theta = c(.60,    .40,    .50,  .50,  .40),
        stringsAsFactors = FALSE
    )

    for (rel in c("AND", "OR")) {
        m <- m_base
        m$Relation <- c("Single", "Single", rel, rel, "Single")

        res <- intervene_cpm_every_gene(list(OncoBN_model = m), "OncoBN")

        if (rel == "AND") {
            ## Kill A: C cascades (AND needs A+B), D survives (only needs B).
            expect_equal(get_interv(res, "I:A"), preds_from_oncobn(m[c(2, 5), ]),
                         label = "AND no-over kill A")
            ## Kill B: C and D both cascade.
            expect_equal(get_interv(res, "I:B"), preds_from_oncobn(m[1, ]),
                         label = "AND no-over kill B")
            ## Kill C: D unaffected.
            expect_equal(get_interv(res, "I:C"), preds_from_oncobn(m[c(1, 2, 5), ]),
                         label = "AND no-over kill C")
            ## Kill D: leaf.
            expect_equal(get_interv(res, "I:D"), preds_from_oncobn(m[c(1, 2, 3, 4), ]),
                         label = "AND no-over kill D")
        } else {
            ## Kill A: C alive via B, D alive via B.
            expect_equal(get_interv(res, "I:A"), preds_from_oncobn(m[c(2, 4, 5), ]),
                         label = "OR no-over kill A")
            ## Kill B: C alive via A; D unreachable (sole parent B gone).
            expect_equal(get_interv(res, "I:B"), preds_from_oncobn(m[c(1, 3), ]),
                         label = "OR no-over kill B")
            ## Kill C: D unaffected.
            expect_equal(get_interv(res, "I:C"), preds_from_oncobn(m[c(1, 2, 5), ]),
                         label = "OR no-over kill C")
            ## Kill D: leaf.
            expect_equal(get_interv(res, "I:D"), preds_from_oncobn(m[c(1, 2, 3, 4), ]),
                         label = "OR no-over kill D")
        }

        check_oncobn_two_procedures(m, label = paste("no-over-cascade", rel))
    }
})


### Full ground-truth test: DAG_6, every gene, AND and OR
##
## DAG_6 rows:
##   1: Root -> C  (Single)
##   2: Root -> D  (Single)
##   3: Root -> F  (Single)
##   4: Root -> G  (Single)
##   5: C -> A     (rel)
##   6: D -> A     (rel)
##   7: F -> B     (rel)
##   8: G -> B     (rel)
##   9: A -> E     (rel)
##  10: B -> E     (rel)
##
## Multi-parent nodes: A (parents C,D), B (parents F,G), E (parents A,B).
##
## AND mode — surviving rows after each kill (worked by hand):
##   Kill C: A unreachable (needs C+D) -> cascade rows 5,6; E unreachable -> cascade rows 9,10.
##           Surviving: 2,3,4,7,8
##   Kill D: A unreachable -> cascade rows 5,6; E unreachable -> cascade rows 9,10.
##           Surviving: 1,3,4,7,8
##   Kill F: B unreachable (needs F+G) -> cascade rows 7,8; E unreachable -> cascade rows 9,10.
##           Surviving: 1,2,4,5,6
##   Kill G: B unreachable -> cascade rows 7,8; E unreachable -> cascade rows 9,10.
##           Surviving: 1,2,3,5,6
##   Kill A: E unreachable (needs A+B) -> cascade rows 9,10.
##           Surviving: 1,2,3,4,7,8
##   Kill B: E unreachable -> cascade rows 9,10.
##           Surviving: 1,2,3,4,5,6
##   Kill E: leaf. Surviving: 1,2,3,4,5,6,7,8
##
## OR mode — surviving rows after each kill (worked by hand):
##   Kill C: A alive via D (row 6 survives). E alive. Remove rows 1,5.
##           Surviving: 2,3,4,6,7,8,9,10
##   Kill D: A alive via C (row 5 survives). E alive. Remove rows 2,6.
##           Surviving: 1,3,4,5,7,8,9,10
##   Kill F: B alive via G (row 8 survives). E alive. Remove rows 3,7.
##           Surviving: 1,2,4,5,6,8,9,10
##   Kill G: B alive via F (row 7 survives). E alive. Remove rows 4,8.
##           Surviving: 1,2,3,5,6,7,9,10
##   Kill A: E alive via B (row 10 survives). Remove rows 5,6,9.
##           Surviving: 1,2,3,4,7,8,10
##   Kill B: E alive via A (row 9 survives). Remove rows 7,8,10.
##           Surviving: 1,2,3,4,5,6,9
##   Kill E: leaf. Remove rows 9,10. Surviving: 1,2,3,4,5,6,7,8

### Full ground-truth test: DAG_7, every gene killed
##
## DAG_7 rows:
##   1: Root -> C   2: Root -> D   3: Root -> F
##   4: C -> A  (rel)   5: D -> A  (rel)   6: F -> A  (rel)   (A parents: C,D,F)
##   7: D -> B  (rel)   8: F -> B  (rel)                       (B parents: D,F)
##   9: A -> E  (rel)  10: B -> E  (rel)                       (E parents: A,B)
##
## AND mode surviving rows (same as CBN DAG_7):
##   Kill C: A unreachable (needs C+D+F) -> cascade 4,5,6; E -> 9,10. B ok.
##           Surviving: 2,3,7,8
##   Kill D: A unreachable -> cascade 4,5,6; B unreachable -> cascade 7,8; E -> 9,10.
##           Surviving: 1,3
##   Kill F: A,B both unreachable; E -> 9,10. Surviving: 1,2
##   Kill A: E unreachable -> cascade 9,10. B ok. Surviving: 1,2,3,7,8
##   Kill B: E unreachable -> cascade 9,10. A ok. Surviving: 1,2,3,4,5,6
##   Kill E: leaf. Surviving: 1,2,3,4,5,6,7,8
##
## OR mode surviving rows:
##   Kill C: remove 1(Root->C),4(C->A). A alive via D,F. B,E alive.
##           Surviving: 2,3,5,6,7,8,9,10
##   Kill D: remove 2,5(D->A),7(D->B). A alive via C,F. B alive via F. E alive.
##           Surviving: 1,3,4,6,8,9,10
##   Kill F: remove 3,6(F->A),8(F->B). A alive via C,D. B alive via D. E alive.
##           Surviving: 1,2,4,5,7,9,10
##   Kill A: remove 4,5,6(->A),9(A->E). E alive via B (OR). Row 10 stays.
##           Surviving: 1,2,3,7,8,10
##   Kill B: remove 7,8(->B),10(B->E). E alive via A (OR). Row 9 stays.
##           Surviving: 1,2,3,4,5,6,9
##   Kill E: leaf; same as AND. Surviving: 1,2,3,4,5,6,7,8

test_that("OncoBN DAG_7: ground-truth comparison for all gene kills, AND and OR", {
    local_edition(3)

    for (rel in c("AND", "OR")) {
        cat("\n DAG_7 ground-truth:", rel, "\n")
        m <- DAG_7
        m$Relation <- c(rep("Single", 3),   ## 1-3: Root -> C,D,F
                        rel, rel, rel,       ## 4-6: C,D,F -> A
                        rel, rel,            ## 7-8: D,F -> B
                        rel, rel)            ## 9-10: A,B -> E

        res <- intervene_cpm_every_gene(list(OncoBN_model = m), "OncoBN")

        if (rel == "AND") {
            expect_equal(get_interv(res, "I:C"), preds_from_oncobn(m[c(2, 3, 7, 8), ]),
                         label = "AND kill C")
            expect_equal(get_interv(res, "I:D"), preds_from_oncobn(m[c(1, 3), ]),
                         label = "AND kill D")
            expect_equal(get_interv(res, "I:F"), preds_from_oncobn(m[c(1, 2), ]),
                         label = "AND kill F")
            expect_equal(get_interv(res, "I:A"), preds_from_oncobn(m[c(1, 2, 3, 7, 8), ]),
                         label = "AND kill A")
            expect_equal(get_interv(res, "I:B"), preds_from_oncobn(m[c(1, 2, 3, 4, 5, 6), ]),
                         label = "AND kill B")
            expect_equal(get_interv(res, "I:E"), preds_from_oncobn(m[c(1, 2, 3, 4, 5, 6, 7, 8), ]),
                         label = "AND kill E")
        } else {
            expect_equal(get_interv(res, "I:C"), preds_from_oncobn(m[c(2, 3, 5, 6, 7, 8, 9, 10), ]),
                         label = "OR kill C")
            expect_equal(get_interv(res, "I:D"), preds_from_oncobn(m[c(1, 3, 4, 6, 8, 9, 10), ]),
                         label = "OR kill D")
            expect_equal(get_interv(res, "I:F"), preds_from_oncobn(m[c(1, 2, 4, 5, 7, 9, 10), ]),
                         label = "OR kill F")
            expect_equal(get_interv(res, "I:A"), preds_from_oncobn(m[c(1, 2, 3, 7, 8, 10), ]),
                         label = "OR kill A")
            expect_equal(get_interv(res, "I:B"), preds_from_oncobn(m[c(1, 2, 3, 4, 5, 6, 9), ]),
                         label = "OR kill B")
            expect_equal(get_interv(res, "I:E"), preds_from_oncobn(m[c(1, 2, 3, 4, 5, 6, 7, 8), ]),
                         label = "OR kill E")
        }
    }
})


### Full ground-truth test: DAG_2, every gene killed
##
## DAG_2 rows:
##   1: Root->F   2: Root->D   3: Root->G   4: Root->J   5: Root->B
##   6: F->H   7: D->H   8: G->H   9: J->H  10: B->H   (H parents: F,D,G,J,B)
##  11: F->I  12: D->I  13: G->I  14: J->I  15: B->I   (I parents: F,D,G,J,B)
##  16: D->A  17: G->A  18: J->A  19: B->A             (A parents: D,G,J,B)
##  20: H->C  21: H->E                                  (C,E: sole parent H)
##
## AND mode (same as CBN DAG_2):
##   Kill F: H,I unreachable (need F); A unaffected (F not a parent). C,E cascade.
##           Surviving: 2,3,4,5,16,17,18,19
##   Kill D: H,I,A all unreachable. C,E cascade. Surviving: 1,3,4,5
##   Kill G: H,I,A unreachable. Surviving: 1,2,4,5
##   Kill J: H,I,A unreachable. Surviving: 1,2,3,5
##   Kill B: H,I,A unreachable. Surviving: 1,2,3,4
##   Kill H: C,E cascade (sole parent H). I,A unaffected. Surviving: 1-5,11-19
##   Kill I: leaf. Surviving: 1-10,16-21
##   Kill A: leaf. Surviving: 1-15,20,21
##   Kill C: leaf. Surviving: 1-19,21
##   Kill E: leaf. Surviving: 1-20
##
## OR mode (H,I,A are OR nodes; C,E have sole parent H so AND=OR for kill H,I,A,C,E):
##   Kill F: remove 1(Root->F),6(F->H),11(F->I). H alive via D,G,J,B. I alive. A unaffected.
##           C,E via H alive. Surviving: 2,3,4,5,7,8,9,10,12,13,14,15,16,17,18,19,20,21
##   Kill D: remove 2,7,12,16. H,I,A alive. C,E alive.
##           Surviving: 1,3,4,5,6,8,9,10,11,13,14,15,17,18,19,20,21
##   Kill G: remove 3,8,13,17. H,I,A alive.
##           Surviving: 1,2,4,5,6,7,9,10,11,12,14,15,16,18,19,20,21
##   Kill J: remove 4,9,14,18. H,I,A alive.
##           Surviving: 1,2,3,5,6,7,8,10,11,12,13,15,16,17,19,20,21
##   Kill B: remove 5,10,15,19. H,I,A alive.
##           Surviving: 1,2,3,4,6,7,8,9,11,12,13,14,16,17,18,20,21
##   Kill H: same as AND (C,E have sole parent H). Surviving: 1,2,3,4,5,11-19
##   Kill I,A,C,E: leaf behaviour same as AND.

test_that("OncoBN DAG_2: ground-truth comparison for all gene kills, AND and OR", {
    local_edition(3)

    for (rel in c("AND", "OR")) {
        cat("\n DAG_2 ground-truth:", rel, "\n")
        m <- DAG_2
        m$Relation <- c(rep("Single", 5),   ## 1-5: Root -> F,D,G,J,B
                        rep(rel, 5),         ## 6-10: F,D,G,J,B -> H
                        rep(rel, 5),         ## 11-15: F,D,G,J,B -> I
                        rep(rel, 4),         ## 16-19: D,G,J,B -> A
                        rep("Single", 2))    ## 20-21: H->C, H->E

        res <- intervene_cpm_every_gene(list(OncoBN_model = m), "OncoBN")

        if (rel == "AND") {
            expect_equal(get_interv(res, "I:F"), preds_from_oncobn(m[c(2, 3, 4, 5, 16, 17, 18, 19), ]),
                         label = "AND kill F")
            expect_equal(get_interv(res, "I:D"), preds_from_oncobn(m[c(1, 3, 4, 5), ]),
                         label = "AND kill D")
            expect_equal(get_interv(res, "I:G"), preds_from_oncobn(m[c(1, 2, 4, 5), ]),
                         label = "AND kill G")
            expect_equal(get_interv(res, "I:J"), preds_from_oncobn(m[c(1, 2, 3, 5), ]),
                         label = "AND kill J")
            expect_equal(get_interv(res, "I:B"), preds_from_oncobn(m[c(1, 2, 3, 4), ]),
                         label = "AND kill B")
            expect_equal(get_interv(res, "I:H"), preds_from_oncobn(m[c(1:5, 11:19), ]),
                         label = "AND kill H")
            expect_equal(get_interv(res, "I:I"), preds_from_oncobn(m[c(1:10, 16:21), ]),
                         label = "AND kill I")
            expect_equal(get_interv(res, "I:A"), preds_from_oncobn(m[c(1:15, 20, 21), ]),
                         label = "AND kill A")
            expect_equal(get_interv(res, "I:C"), preds_from_oncobn(m[c(1:19, 21), ]),
                         label = "AND kill C")
            expect_equal(get_interv(res, "I:E"), preds_from_oncobn(m[c(1:20), ]),
                         label = "AND kill E")
        } else {
            expect_equal(get_interv(res, "I:F"),
                         preds_from_oncobn(m[c(2,3,4,5,7,8,9,10,12,13,14,15,16,17,18,19,20,21), ]),
                         label = "OR kill F")
            expect_equal(get_interv(res, "I:D"),
                         preds_from_oncobn(m[c(1,3,4,5,6,8,9,10,11,13,14,15,17,18,19,20,21), ]),
                         label = "OR kill D")
            expect_equal(get_interv(res, "I:G"),
                         preds_from_oncobn(m[c(1,2,4,5,6,7,9,10,11,12,14,15,16,18,19,20,21), ]),
                         label = "OR kill G")
            expect_equal(get_interv(res, "I:J"),
                         preds_from_oncobn(m[c(1,2,3,5,6,7,8,10,11,12,13,15,16,17,19,20,21), ]),
                         label = "OR kill J")
            expect_equal(get_interv(res, "I:B"),
                         preds_from_oncobn(m[c(1,2,3,4,6,7,8,9,11,12,13,14,16,17,18,20,21), ]),
                         label = "OR kill B")
            ## H,I,A,C,E: same as AND (C,E sole parent H; H killed; I,A,C,E leaves)
            expect_equal(get_interv(res, "I:H"), preds_from_oncobn(m[c(1:5, 11:19), ]),
                         label = "OR kill H")
            expect_equal(get_interv(res, "I:I"), preds_from_oncobn(m[c(1:10, 16:21), ]),
                         label = "OR kill I")
            expect_equal(get_interv(res, "I:A"), preds_from_oncobn(m[c(1:15, 20, 21), ]),
                         label = "OR kill A")
            expect_equal(get_interv(res, "I:C"), preds_from_oncobn(m[c(1:19, 21), ]),
                         label = "OR kill C")
            expect_equal(get_interv(res, "I:E"), preds_from_oncobn(m[c(1:20), ]),
                         label = "OR kill E")
        }
    }
})


### Full ground-truth test: DAG_3, every gene killed
##
## DAG_3 rows:
##   1: Root->I   2: Root->H
##   3: I->E   4: I->B
##   5: H->C   (rel; C parent H)
##   6: E->G
##   7: B->C   (rel; C parent B)
##   8: G->F   9: G->A
##  10: F->C   (rel; C parent F)
##  11: A->D
##  12: D->C   (rel; C parent D)
##  C has 4 parents: H,B,F,D.
##
## AND mode (same as CBN DAG_3):
##   Kill I: E,B cascade; G cascades; F,A cascade; D cascades; C unreachable. Surviving: 2
##   Kill H: C unreachable. Surviving: 1,3,4,6,8,9,11
##   Kill E: G,F,A,D cascade; C unreachable. Surviving: 1,2,4
##   Kill B: C unreachable. Surviving: 1,2,3,6,8,9,11
##   Kill G: F,A,D cascade; C unreachable. Surviving: 1,2,3,4
##   Kill F: C unreachable. Surviving: 1,2,3,4,6,9,11
##   Kill A: D cascades; C unreachable. Surviving: 1,2,3,4,6,8
##   Kill D: C unreachable. Surviving: 1,2,3,4,6,8,9
##   Kill C: leaf. Surviving: 1,2,3,4,6,8,9,11
##
## OR mode:
##   Kill I: E,B,G,F,A,D all cascade (all single-parent chains). C alive via H only!
##           Surviving: 2,5
##   Kill H: remove 2,5(H->C). C alive via B,F,D (OR). All other rows intact.
##           Surviving: 1,3,4,6,7,8,9,10,11,12
##   Kill E: G,F,A,D cascade (single chains). C alive via H,B (B unaffected).
##           Surviving: 1,2,4,5,7
##   Kill B: remove 4(I->B),7(B->C). C alive via H,F,D.
##           Surviving: 1,2,3,5,6,8,9,10,11,12
##   Kill G: F,A,D cascade. C alive via H,B (row 5,7 survive).
##           Surviving: 1,2,3,4,5,7
##   Kill F: remove 8(G->F),10(F->C). C alive via H,B,D.
##           Surviving: 1,2,3,4,5,6,7,9,11,12
##   Kill A: D cascades (sole parent A); remove 12(D->C). C alive via H,B,F (rows 5,7,10).
##           Surviving: 1,2,3,4,5,6,7,8,10
##   Kill D: remove 11(A->D),12(D->C). C alive via H,B,F.
##           Surviving: 1,2,3,4,5,6,7,8,9,10
##   Kill C: leaf; same as AND. Surviving: 1,2,3,4,6,8,9,11

test_that("OncoBN DAG_3: ground-truth comparison for all gene kills, AND and OR", {
    local_edition(3)

    for (rel in c("AND", "OR")) {
        cat("\n DAG_3 ground-truth:", rel, "\n")
        m <- DAG_3
        m$Relation <- c("Single",  ##  1: Root -> I
                        "Single",  ##  2: Root -> H
                        "Single",  ##  3: I -> E
                        "Single",  ##  4: I -> B
                        rel,       ##  5: H -> C
                        "Single",  ##  6: E -> G
                        rel,       ##  7: B -> C
                        "Single",  ##  8: G -> F
                        "Single",  ##  9: G -> A
                        rel,       ## 10: F -> C
                        "Single",  ## 11: A -> D
                        rel)       ## 12: D -> C

        res <- intervene_cpm_every_gene(list(OncoBN_model = m), "OncoBN")

        if (rel == "AND") {
            expect_equal(get_interv(res, "I:I"), preds_from_oncobn(m[2, ]),
                         label = "AND kill I")
            expect_equal(get_interv(res, "I:H"), preds_from_oncobn(m[c(1, 3, 4, 6, 8, 9, 11), ]),
                         label = "AND kill H")
            expect_equal(get_interv(res, "I:E"), preds_from_oncobn(m[c(1, 2, 4), ]),
                         label = "AND kill E")
            expect_equal(get_interv(res, "I:B"), preds_from_oncobn(m[c(1, 2, 3, 6, 8, 9, 11), ]),
                         label = "AND kill B")
            expect_equal(get_interv(res, "I:G"), preds_from_oncobn(m[c(1, 2, 3, 4), ]),
                         label = "AND kill G")
            expect_equal(get_interv(res, "I:F"), preds_from_oncobn(m[c(1, 2, 3, 4, 6, 9, 11), ]),
                         label = "AND kill F")
            expect_equal(get_interv(res, "I:A"), preds_from_oncobn(m[c(1, 2, 3, 4, 6, 8), ]),
                         label = "AND kill A")
            expect_equal(get_interv(res, "I:D"), preds_from_oncobn(m[c(1, 2, 3, 4, 6, 8, 9), ]),
                         label = "AND kill D")
            expect_equal(get_interv(res, "I:C"), preds_from_oncobn(m[c(1, 2, 3, 4, 6, 8, 9, 11), ]),
                         label = "AND kill C")
        } else {
            expect_equal(get_interv(res, "I:I"), preds_from_oncobn(m[c(2, 5), ]),
                         label = "OR kill I")
            expect_equal(get_interv(res, "I:H"), preds_from_oncobn(m[c(1, 3, 4, 6, 7, 8, 9, 10, 11, 12), ]),
                         label = "OR kill H")
            expect_equal(get_interv(res, "I:E"), preds_from_oncobn(m[c(1, 2, 4, 5, 7), ]),
                         label = "OR kill E")
            expect_equal(get_interv(res, "I:B"), preds_from_oncobn(m[c(1, 2, 3, 5, 6, 8, 9, 10, 11, 12), ]),
                         label = "OR kill B")
            expect_equal(get_interv(res, "I:G"), preds_from_oncobn(m[c(1, 2, 3, 4, 5, 7), ]),
                         label = "OR kill G")
            expect_equal(get_interv(res, "I:F"), preds_from_oncobn(m[c(1, 2, 3, 4, 5, 6, 7, 9, 11, 12), ]),
                         label = "OR kill F")
            expect_equal(get_interv(res, "I:A"), preds_from_oncobn(m[c(1, 2, 3, 4, 5, 6, 7, 8, 10), ]),
                         label = "OR kill A")
            expect_equal(get_interv(res, "I:D"), preds_from_oncobn(m[c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), ]),
                         label = "OR kill D")
            expect_equal(get_interv(res, "I:C"), preds_from_oncobn(m[c(1, 2, 3, 4, 6, 8, 9, 11), ]),
                         label = "OR kill C")
        }
    }
})


### Full ground-truth test: DAG_4, every gene killed
##
## DAG_4 rows:
##   1: Root->C   2: Root->A   3: Root->B   4: Root->I
##   5: C->F  (Single)
##   6: C->E   7: A->E   8: B->E   (E parents: C,A,B)
##   9: B->D  (Single)
##  10: I->H   (rel; H parent I)
##  11: F->G  12: E->G  13: D->G   (G parents: F,E,D)
##  14: G->H   (rel; H parent G)
##
## AND mode (same as CBN DAG_4):
##   Kill C: F cascades; E unreachable; G,H unreachable. D ok. Surviving: 2,3,4,9
##   Kill A: E unreachable; G,H cascade. F,D ok. Surviving: 1,3,4,5,9
##   Kill B: E unreachable; D cascades; G,H cascade. F ok. Surviving: 1,2,4,5
##   Kill I: H unreachable (needs I+G). Surviving: 1,2,3,5,6,7,8,9,11,12,13
##   Kill F: G unreachable; H unreachable. Surviving: 1,2,3,4,6,7,8,9
##   Kill E: G unreachable; H unreachable. Surviving: 1,2,3,4,5,9
##   Kill D: G unreachable; H unreachable. Surviving: 1,2,3,4,5,6,7,8
##   Kill G: H unreachable. Surviving: 1,2,3,4,5,6,7,8,9
##   Kill H: leaf. Surviving: 1,2,3,4,5,6,7,8,9,11,12,13
##
## OR mode:
##   Kill C: remove 1,5(C->F),6(C->E). F gone (sole parent C); remove 11(F->G).
##           E alive via A,B. G alive via E,D. H alive.
##           Surviving: 2,3,4,7,8,9,10,12,13,14
##   Kill A: remove 2,7(A->E). E alive via C,B. G,H alive.
##           Surviving: 1,3,4,5,6,8,9,10,11,12,13,14
##   Kill B: remove 3,8(B->E),9(B->D). E alive via C,A. D gone; remove 13(D->G).
##           G alive via F,E. H alive.
##           Surviving: 1,2,4,5,6,7,10,11,12,14
##   Kill I: remove 4,10(I->H). H alive via G (OR). Surviving: 1,2,3,5,6,7,8,9,11,12,13,14
##   Kill F: remove 5(C->F),11(F->G). G alive via E,D. H alive.
##           Surviving: 1,2,3,4,6,7,8,9,10,12,13,14
##   Kill E: remove 6,7,8(->E),12(E->G). G alive via F,D. H alive.
##           Surviving: 1,2,3,4,5,9,10,11,13,14
##   Kill D: remove 9(B->D),13(D->G). G alive via F,E. H alive.
##           Surviving: 1,2,3,4,5,6,7,8,10,11,12,14
##   Kill G: remove 11,12,13(->G),14(G->H). H alive via I (OR).
##           Surviving: 1,2,3,4,5,6,7,8,9,10
##   Kill H: leaf; same as AND. Surviving: 1,2,3,4,5,6,7,8,9,11,12,13

test_that("OncoBN DAG_4: ground-truth comparison for all gene kills, AND and OR", {
    local_edition(3)

    for (rel in c("AND", "OR")) {
        cat("\n DAG_4 ground-truth:", rel, "\n")
        m <- DAG_4
        m$Relation <- c("Single",        ##  1: Root -> C
                        "Single",        ##  2: Root -> A
                        "Single",        ##  3: Root -> B
                        "Single",        ##  4: Root -> I
                        "Single",        ##  5: C -> F
                        rel, rel, rel,   ##  6-8: C,A,B -> E
                        "Single",        ##  9: B -> D
                        rel,             ## 10: I -> H
                        rel, rel, rel,   ## 11-13: F,E,D -> G
                        rel)             ## 14: G -> H

        res <- intervene_cpm_every_gene(list(OncoBN_model = m), "OncoBN")

        if (rel == "AND") {
            expect_equal(get_interv(res, "I:C"), preds_from_oncobn(m[c(2, 3, 4, 9), ]),
                         label = "AND kill C")
            expect_equal(get_interv(res, "I:A"), preds_from_oncobn(m[c(1, 3, 4, 5, 9), ]),
                         label = "AND kill A")
            expect_equal(get_interv(res, "I:B"), preds_from_oncobn(m[c(1, 2, 4, 5), ]),
                         label = "AND kill B")
            expect_equal(get_interv(res, "I:I"), preds_from_oncobn(m[c(1, 2, 3, 5, 6, 7, 8, 9, 11, 12, 13), ]),
                         label = "AND kill I")
            expect_equal(get_interv(res, "I:F"), preds_from_oncobn(m[c(1, 2, 3, 4, 6, 7, 8, 9), ]),
                         label = "AND kill F")
            expect_equal(get_interv(res, "I:E"), preds_from_oncobn(m[c(1, 2, 3, 4, 5, 9), ]),
                         label = "AND kill E")
            expect_equal(get_interv(res, "I:D"), preds_from_oncobn(m[c(1, 2, 3, 4, 5, 6, 7, 8), ]),
                         label = "AND kill D")
            expect_equal(get_interv(res, "I:G"), preds_from_oncobn(m[c(1, 2, 3, 4, 5, 6, 7, 8, 9), ]),
                         label = "AND kill G")
            expect_equal(get_interv(res, "I:H"), preds_from_oncobn(m[c(1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13), ]),
                         label = "AND kill H")
        } else {
            expect_equal(get_interv(res, "I:C"), preds_from_oncobn(m[c(2, 3, 4, 7, 8, 9, 10, 12, 13, 14), ]),
                         label = "OR kill C")
            expect_equal(get_interv(res, "I:A"), preds_from_oncobn(m[c(1, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14), ]),
                         label = "OR kill A")
            expect_equal(get_interv(res, "I:B"), preds_from_oncobn(m[c(1, 2, 4, 5, 6, 7, 10, 11, 12, 14), ]),
                         label = "OR kill B")
            expect_equal(get_interv(res, "I:I"), preds_from_oncobn(m[c(1, 2, 3, 5, 6, 7, 8, 9, 11, 12, 13, 14), ]),
                         label = "OR kill I")
            expect_equal(get_interv(res, "I:F"), preds_from_oncobn(m[c(1, 2, 3, 4, 6, 7, 8, 9, 10, 12, 13, 14), ]),
                         label = "OR kill F")
            expect_equal(get_interv(res, "I:E"), preds_from_oncobn(m[c(1, 2, 3, 4, 5, 9, 10, 11, 13, 14), ]),
                         label = "OR kill E")
            expect_equal(get_interv(res, "I:D"), preds_from_oncobn(m[c(1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 14), ]),
                         label = "OR kill D")
            expect_equal(get_interv(res, "I:G"), preds_from_oncobn(m[c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), ]),
                         label = "OR kill G")
            expect_equal(get_interv(res, "I:H"), preds_from_oncobn(m[c(1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13), ]),
                         label = "OR kill H")
        }
    }
})


### Full ground-truth test: DAG_5, every gene killed
##
## DAG_5 rows:
##   1: Root->A   2: Root->F   3: Root->G   4: Root->I
##   5: A->C   6: F->C   7: G->C   8: I->C   (C parents: A,F,G,I)
##   9: A->D  10: F->D  11: G->D             (D parents: A,F,G)
##  12: C->E  (Single)
##  13: D->H  (Single)
##  14: E->B  15: D->B                        (B parents: E,D)
##
## AND mode (same as CBN DAG_5):
##   Kill A: C,D unreachable; E,H,B cascade. Surviving: 2,3,4
##   Kill F: C,D unreachable; E,H,B cascade. Surviving: 1,3,4
##   Kill G: C,D unreachable; E,H,B cascade. Surviving: 1,2,4
##   Kill I: C unreachable; D,H ok; E cascades; B unreachable. Surviving: 1,2,3,9,10,11,13
##   Kill C: E cascades; B unreachable. D,H ok. Surviving: 1,2,3,4,9,10,11,13
##   Kill D: H cascades; B unreachable. C,E ok. Surviving: 1,2,3,4,5,6,7,8,12
##   Kill E: B unreachable. D,H ok. Surviving: 1,2,3,4,5,6,7,8,9,10,11,13
##   Kill H: leaf. Surviving: 1,2,3,4,5,6,7,8,9,10,11,12,14,15
##   Kill B: leaf. Surviving: 1,2,3,4,5,6,7,8,9,10,11,12,13
##
## OR mode:
##   Kill A: remove 1,5(A->C),9(A->D). C alive via F,G,I. D alive via F,G.
##           E alive. H alive. B alive. Surviving: 2,3,4,6,7,8,10,11,12,13,14,15
##   Kill F: remove 2,6,10. C alive. D alive. Surviving: 1,3,4,5,7,8,9,11,12,13,14,15
##   Kill G: remove 3,7,11. C alive. D alive. Surviving: 1,2,4,5,6,8,9,10,12,13,14,15
##   Kill I: remove 4,8(I->C). C alive via A,F,G. D unaffected.
##           Surviving: 1,2,3,5,6,7,9,10,11,12,13,14,15
##   Kill C: remove 5,6,7,8,12(C->E). E gone; remove 14(E->B). B alive via D (OR).
##           Surviving: 1,2,3,4,9,10,11,13,15
##   Kill D: remove 9,10,11,13(D->H),15(D->B). H gone (sole parent). B alive via E (OR).
##           Surviving: 1,2,3,4,5,6,7,8,12,14
##   Kill E: remove 12(C->E),14(E->B). B alive via D (OR).
##           Surviving: 1,2,3,4,5,6,7,8,9,10,11,13,15
##   Kill H: leaf; same as AND. Surviving: 1,2,3,4,5,6,7,8,9,10,11,12,14,15
##   Kill B: leaf; same as AND. Surviving: 1,2,3,4,5,6,7,8,9,10,11,12,13

test_that("OncoBN DAG_5: ground-truth comparison for all gene kills, AND and OR", {
    local_edition(3)

    for (rel in c("AND", "OR")) {
        cat("\n DAG_5 ground-truth:", rel, "\n")
        m <- DAG_5
        m$Relation <- c(rep("Single", 4),   ##  1-4: Root -> A,F,G,I
                        rep(rel, 4),         ##  5-8: A,F,G,I -> C
                        rep(rel, 3),         ##  9-11: A,F,G -> D
                        "Single",            ## 12: C -> E
                        "Single",            ## 13: D -> H
                        rel, rel)            ## 14-15: E,D -> B

        res <- intervene_cpm_every_gene(list(OncoBN_model = m), "OncoBN")

        if (rel == "AND") {
            expect_equal(get_interv(res, "I:A"), preds_from_oncobn(m[c(2, 3, 4), ]),
                         label = "AND kill A")
            expect_equal(get_interv(res, "I:F"), preds_from_oncobn(m[c(1, 3, 4), ]),
                         label = "AND kill F")
            expect_equal(get_interv(res, "I:G"), preds_from_oncobn(m[c(1, 2, 4), ]),
                         label = "AND kill G")
            expect_equal(get_interv(res, "I:I"), preds_from_oncobn(m[c(1, 2, 3, 9, 10, 11, 13), ]),
                         label = "AND kill I")
            expect_equal(get_interv(res, "I:C"), preds_from_oncobn(m[c(1, 2, 3, 4, 9, 10, 11, 13), ]),
                         label = "AND kill C")
            expect_equal(get_interv(res, "I:D"), preds_from_oncobn(m[c(1, 2, 3, 4, 5, 6, 7, 8, 12), ]),
                         label = "AND kill D")
            expect_equal(get_interv(res, "I:E"), preds_from_oncobn(m[c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13), ]),
                         label = "AND kill E")
            expect_equal(get_interv(res, "I:H"), preds_from_oncobn(m[c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15), ]),
                         label = "AND kill H")
            expect_equal(get_interv(res, "I:B"), preds_from_oncobn(m[c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13), ]),
                         label = "AND kill B")
        } else {
            expect_equal(get_interv(res, "I:A"),
                         preds_from_oncobn(m[c(2, 3, 4, 6, 7, 8, 10, 11, 12, 13, 14, 15), ]),
                         label = "OR kill A")
            expect_equal(get_interv(res, "I:F"),
                         preds_from_oncobn(m[c(1, 3, 4, 5, 7, 8, 9, 11, 12, 13, 14, 15), ]),
                         label = "OR kill F")
            expect_equal(get_interv(res, "I:G"),
                         preds_from_oncobn(m[c(1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14, 15), ]),
                         label = "OR kill G")
            expect_equal(get_interv(res, "I:I"),
                         preds_from_oncobn(m[c(1, 2, 3, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15), ]),
                         label = "OR kill I")
            expect_equal(get_interv(res, "I:C"), preds_from_oncobn(m[c(1, 2, 3, 4, 9, 10, 11, 13, 15), ]),
                         label = "OR kill C")
            expect_equal(get_interv(res, "I:D"), preds_from_oncobn(m[c(1, 2, 3, 4, 5, 6, 7, 8, 12, 14), ]),
                         label = "OR kill D")
            expect_equal(get_interv(res, "I:E"),
                         preds_from_oncobn(m[c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15), ]),
                         label = "OR kill E")
            expect_equal(get_interv(res, "I:H"),
                         preds_from_oncobn(m[c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15), ]),
                         label = "OR kill H")
            expect_equal(get_interv(res, "I:B"),
                         preds_from_oncobn(m[c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13), ]),
                         label = "OR kill B")
        }
    }
})


### Full ground-truth test: DAG_6, every gene, AND and OR
##
test_that("OncoBN DAG_6: ground-truth comparison for all gene kills, AND and OR", {
    local_edition(3)

    for (rel in c("AND", "OR")) {
        cat("\n DAG_6 ground-truth:", rel, "\n")
        m <- DAG_6
        m$Relation <- c(rep("Single", 4), rel, rel, rel, rel, rel, rel)

        res <- intervene_cpm_every_gene(list(OncoBN_model = m), "OncoBN")

        if (rel == "AND") {
            expect_equal(get_interv(res, "I:C"), preds_from_oncobn(m[c(2, 3, 4, 7, 8), ]),
                         label = "AND kill C")
            expect_equal(get_interv(res, "I:D"), preds_from_oncobn(m[c(1, 3, 4, 7, 8), ]),
                         label = "AND kill D")
            expect_equal(get_interv(res, "I:F"), preds_from_oncobn(m[c(1, 2, 4, 5, 6), ]),
                         label = "AND kill F")
            expect_equal(get_interv(res, "I:G"), preds_from_oncobn(m[c(1, 2, 3, 5, 6), ]),
                         label = "AND kill G")
            expect_equal(get_interv(res, "I:A"), preds_from_oncobn(m[c(1, 2, 3, 4, 7, 8), ]),
                         label = "AND kill A")
            expect_equal(get_interv(res, "I:B"), preds_from_oncobn(m[c(1, 2, 3, 4, 5, 6), ]),
                         label = "AND kill B")
            expect_equal(get_interv(res, "I:E"), preds_from_oncobn(m[c(1, 2, 3, 4, 5, 6, 7, 8), ]),
                         label = "AND kill E")
        } else {
            expect_equal(get_interv(res, "I:C"), preds_from_oncobn(m[c(2, 3, 4, 6, 7, 8, 9, 10), ]),
                         label = "OR kill C")
            expect_equal(get_interv(res, "I:D"), preds_from_oncobn(m[c(1, 3, 4, 5, 7, 8, 9, 10), ]),
                         label = "OR kill D")
            expect_equal(get_interv(res, "I:F"), preds_from_oncobn(m[c(1, 2, 4, 5, 6, 8, 9, 10), ]),
                         label = "OR kill F")
            expect_equal(get_interv(res, "I:G"), preds_from_oncobn(m[c(1, 2, 3, 5, 6, 7, 9, 10), ]),
                         label = "OR kill G")
            expect_equal(get_interv(res, "I:A"), preds_from_oncobn(m[c(1, 2, 3, 4, 7, 8, 10), ]),
                         label = "OR kill A")
            expect_equal(get_interv(res, "I:B"), preds_from_oncobn(m[c(1, 2, 3, 4, 5, 6, 9), ]),
                         label = "OR kill B")
            expect_equal(get_interv(res, "I:E"), preds_from_oncobn(m[c(1, 2, 3, 4, 5, 6, 7, 8), ]),
                         label = "OR kill E")
        }
    }
})


set.seed(NULL)

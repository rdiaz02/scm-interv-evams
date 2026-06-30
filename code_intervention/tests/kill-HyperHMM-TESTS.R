## Copyright 2025 Ramon Diaz-Uriarte

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

## This file tests the killing procedure of HyperHMM
## (similar to what kill-CBN/OT/OncoBN/HESBCN/MHN-TESTS.R do for
## other methods). Note that for HyperHMM, in contrast to
## other methods, there are no multiple equivalent ways of
## killing; see file kill-gene-equivalences-TESTS.R
## for CBN/OT/OncoBN/HESBCN/MHN equivalences.


## About testing predicted genotype frequencies:
## That is tested in evamtools. In evamtools:
## For OT, evamtools tests OT_model_2_output the same function used here
##         in get_full_output
## For OncoBN, evamtools tests OncoBN_model_2_output, the same function
##         used here in get_full_output
## For CBN, HESBCN, and MHN evamtools tests probs_from_trm, same function used here
##      (called in get_full_output), for any of CBN, HESBCN, MHN
##      And creating the transition rate matrix is also tested in evamtools.


library(testthat)

pwd <- getwd()
setwd("../")
source("intervention.R")
setwd(pwd)

set.seed(NULL)

test_that("HyperHMM, error if trying to use wrong kill_gene funct", {
    rmhn <- random_evam(model = "MHN", ngenes = 3)
    sample_mhn <- sample_evam(rmhn, N = 1000, obs_noise = 0.05)
    dd <- sample_mhn$MHN_sampled_genotype_counts_as_data
    h1 <- evam_like_HyperHMM(dd)

    expect_error(intervene_cpm_every_gene(h1,
                                          method = "HyperHMM",
                                          kill_gene_funct = kill_gene_by_params_to_0),
                 "HyperHMM can only use")
})


test_that("HyperHMM interventions and predicted genot calculations", {
  reps <- 10
    for (i in 1:reps) {
        cat("\n ##### Doing rep = ", i, "\n")
        cat("\n ##### Doing 3 genes\n")
        rmhn <- random_evam(model = "MHN", ngenes = 3)
        sample_mhn <- sample_evam(rmhn, N = round(runif(1, 500, 3000)),
                                  obs_noise = 0.05)
        dd <- sample_mhn$MHN_sampled_genotype_counts_as_data
        h1 <- evam_like_HyperHMM(dd)

        cat("used_prob.set:\n ")
        print(h1$HyperHMM_used_prob.set)
        o1 <- intervene_cpm_every_gene(h1, method = "HyperHMM")
        ## get_full_output for HyperHMM needs the transition matrix
        ## (with its method_output attribute), not the full output list.
        h1_full <- get_full_output(h1$HyperHMM_trans_mat)

        expect_equal(h1$HyperHMM_predicted_genotype_freqs,
                     o1[["no_intervention"]]$genot_freqs)
        local({hp <- h1_full$HyperHMM_hitting_probs_from_WT
               expect_equal(filter_hp_keep_wt(hp), o1[["no_intervention"]]$hitting_probs_from_WT)})

        ## I assume kill_gene is correct, as tested in
        ## already.
        ## Nevertheless, we check here again

        h1_tm <- as.matrix(h1$HyperHMM_trans_mat)

        ## And we use a different logic from the one in the code
        ## Here, matrix powers, not iterated multiplications

        p_0 <- c(1.0, rep(0, 7))
        names(p_0) <- allGenotypesLetter(3)

        ## Baseline, no interventions. This tests probs_from_HyperHMM
        ## And we use a different logic from the one in the code
        ## Here, matrix powers, not iterated multiplications

        P_2 <- h1_tm %^% 2
        P_3 <- h1_tm %^% 3

        p_1 <-  p_0 %*% h1_tm
        p_2 <-  p_0 %*% P_2
        p_3 <-  p_0 %*% P_3

        p_0_w <- p_0 * h1$HyperHMM_used_prob.set[1]
        p_1_w <- p_1 * h1$HyperHMM_used_prob.set[2]
        p_2_w <- p_2 * h1$HyperHMM_used_prob.set[3]
        p_3_w <- p_3 * h1$HyperHMM_used_prob.set[4]
        pp <- p_0_w + p_1_w + p_2_w + p_3_w

        expect_equal(pp[1, ], ## we want a named vector
                     h1$HyperHMM_predicted_genotype_freqs)

        expect_equal(pp[1, ], ## we want a named vector
                     o1[["no_intervention"]]$genot_freqs)

        ## Intervened
        o1_A <- h1_tm
        o1_B <- h1_tm
        o1_C <- h1_tm

        ## Change values manually
        o1_A["A", ] <- 0
        o1_A["A, B", ] <- 0
        o1_A["A, C", ] <- 0
        o1_A["WT", "WT"] <- h1_tm["WT", "A"]
        o1_A["WT", "A"] <- 0
        o1_A["B", "B"] <- h1_tm["B", "A, B"]
        o1_A["B", "A, B"] <- 0
        o1_A["C", "C"] <- h1_tm["C", "A, C"]
        o1_A["C", "A, C"] <- 0
        o1_A["B, C", "B, C"] <- 1
        o1_A["B, C", "A, B, C"] <- 0

        o1_B["B", ] <- 0
        o1_B["A, B", ] <- 0
        o1_B["B, C", ] <- 0
        o1_B["WT", "WT"] <- h1_tm["WT", "B"]
        o1_B["WT", "B"] <- 0
        o1_B["A", "A"] <- h1_tm["A", "A, B"]
        o1_B["A", "A, B"] <- 0
        o1_B["C", "C"] <- h1_tm["C", "B, C"]
        o1_B["C", "B, C"] <- 0
        o1_B["A, C", "A, C"] <- 1
        o1_B["A, C", "A, B, C"] <- 0

        o1_C["C", ] <- 0
        o1_C["B, C", ] <- 0
        o1_C["A, C", ] <- 0
        o1_C["WT", "WT"] <- h1_tm["WT", "C"]
        o1_C["WT", "C"] <- 0
        o1_C["A", "A"] <- h1_tm["A", "A, C"]
        o1_C["A", "A, C"] <- 0
        o1_C["B", "B"] <- h1_tm["B", "B, C"]
        o1_C["B", "B, C"] <- 0
        o1_C["A, B", "A, B"] <- 1
        o1_C["A, B", "A, B, C"] <- 0


        P_A_2 <- o1_A %^% 2
        P_A_3 <- o1_A %^% 3
        p_A_1 <- p_0 %*% o1_A
        p_A_2 <- p_0 %*% P_A_2
        p_A_3 <- p_0 %*% P_A_3

        p_A_0_w <- p_0 * h1$HyperHMM_used_prob.set[1]
        p_A_1_w <- p_A_1 * h1$HyperHMM_used_prob.set[2]
        p_A_2_w <- p_A_2 * h1$HyperHMM_used_prob.set[3]
        p_A_3_w <- p_A_3 * h1$HyperHMM_used_prob.set[4]
        pp_A <- p_A_0_w + p_A_1_w + p_A_2_w + p_A_3_w
        pp_A <- pp_A[1, ]
        (pp_A <- pp_A[pp_A > 0])

        P_B_2 <- o1_B %^% 2
        P_B_3 <- o1_B %^% 3
        p_B_1 <- p_0 %*% o1_B
        p_B_2 <- p_0 %*% P_B_2
        p_B_3 <- p_0 %*% P_B_3

        p_B_0_w <- p_0 * h1$HyperHMM_used_prob.set[1]
        p_B_1_w <- p_B_1 * h1$HyperHMM_used_prob.set[2]
        p_B_2_w <- p_B_2 * h1$HyperHMM_used_prob.set[3]
        p_B_3_w <- p_B_3 * h1$HyperHMM_used_prob.set[4]
        pp_B <- p_B_0_w + p_B_1_w + p_B_2_w + p_B_3_w
        pp_B <- pp_B[1, ]
        (pp_B <- pp_B[pp_B > 0])

        P_C_2 <- o1_C %^% 2
        P_C_3 <- o1_C %^% 3
        p_C_1 <- p_0 %*% o1_C
        p_C_2 <- p_0 %*% P_C_2
        p_C_3 <- p_0 %*% P_C_3

        p_C_0_w <- p_0 * h1$HyperHMM_used_prob.set[1]
        p_C_1_w <- p_C_1 * h1$HyperHMM_used_prob.set[2]
        p_C_2_w <- p_C_2 * h1$HyperHMM_used_prob.set[3]
        p_C_3_w <- p_C_3 * h1$HyperHMM_used_prob.set[4]
        pp_C <- p_C_0_w + p_C_1_w + p_C_2_w + p_C_3_w
        pp_C <- pp_C[1, ]
        (pp_C <- pp_C[pp_C > 0])

        expect_equal(pp_A, o1[["I:A"]]$genot_freqs)
        expect_equal(pp_B, o1[["I:B"]]$genot_freqs)
        expect_equal(pp_C, o1[["I:C"]]$genot_freqs)

        ## Explicit hitting-prob checks for the intervened cases (the
        ## genot_freqs checks above only verify WT mass implicitly via the
        ## o1_X["WT", "WT"] <- h1_tm["WT", X] redirection). The hand-built
        ## o1_X matrices equal kill_gene's output (verified in the
        ## test_that("HyperHMM", ...) block below), so feeding them straight
        ## to hitting_probs_from_WT reproduces the production path
        ## (kill_gene -> get_full_output -> hitting_probs_from_WT). WT is
        ## retained in the HP vector by filter_hp_keep_wt (its first-passage
        ## hitting prob is 0 by convention).
        expect_equal(filter_hp_keep_wt(hitting_probs_from_WT(o1_A)),
                     o1[["I:A"]]$hitting_probs_from_WT)
        expect_equal(filter_hp_keep_wt(hitting_probs_from_WT(o1_B)),
                     o1[["I:B"]]$hitting_probs_from_WT)
        expect_equal(filter_hp_keep_wt(hitting_probs_from_WT(o1_C)),
                     o1[["I:C"]]$hitting_probs_from_WT)
        ## WT must be present (kept by filter_hp_keep_wt) after intervention
        expect_true("WT" %in% names(o1[["I:A"]]$hitting_probs_from_WT))
        expect_true("WT" %in% names(o1[["I:B"]]$hitting_probs_from_WT))
        expect_true("WT" %in% names(o1[["I:C"]]$hitting_probs_from_WT))

        cat("\n ##### Doing 5 genes\n")
        ## And now, a 5 gene example
        rmhn <- random_evam(model = "MHN", ngenes = 5)
        sample_mhn <- sample_evam(rmhn, N = round(runif(1, 500, 3000)),
                                  obs_noise = 0.02)
        dd <- sample_mhn$MHN_sampled_genotype_counts_as_data
        h1 <- evam_like_HyperHMM(dd)
        cat("used_prob.set: \n")
        print(h1$HyperHMM_used_prob.set)
        o1 <- intervene_cpm_every_gene(h1, method = "HyperHMM")
        h1_full <- get_full_output(h1$HyperHMM_trans_mat)

        expect_equal(h1$HyperHMM_predicted_genotype_freqs,
                     o1[["no_intervention"]]$genot_freqs)
        local({hp <- h1_full$HyperHMM_hitting_probs_from_WT
          expect_equal(filter_hp_keep_wt(hp),
                       o1[["no_intervention"]]$hitting_probs_from_WT)})

        h1_tm <- as.matrix(h1$HyperHMM_trans_mat)

        p_0 <- c(1.0, rep(0, 31))
        names(p_0) <- allGenotypesLetter(5)

        ## Baseline, no interventions. This tests probs_from_HyperHMM
        ## And we use a different logic from the one in the code
        ## Here, matrix powers, not iterated multiplications

        P_2 <- h1_tm %^% 2
        P_3 <- h1_tm %^% 3
        P_4 <- h1_tm %^% 4
        P_5 <- h1_tm %^% 5

        p_1 <-  p_0 %*% h1_tm
        p_2 <-  p_0 %*% P_2
        p_3 <-  p_0 %*% P_3
        p_4 <-  p_0 %*% P_4
        p_5 <-  p_0 %*% P_5

        p_0_w <- p_0 * h1$HyperHMM_used_prob.set[1]
        p_1_w <- p_1 * h1$HyperHMM_used_prob.set[2]
        p_2_w <- p_2 * h1$HyperHMM_used_prob.set[3]
        p_3_w <- p_3 * h1$HyperHMM_used_prob.set[4]
        p_4_w <- p_4 * h1$HyperHMM_used_prob.set[5]
        p_5_w <- p_5 * h1$HyperHMM_used_prob.set[6]

        pp <- p_0_w + p_1_w + p_2_w + p_3_w + p_4_w + p_5_w

        expect_equal(pp[1, ], ## we want a named vector
                     h1$HyperHMM_predicted_genotype_freqs)

        expect_equal(pp[1, ], ## we want a named vector
                     o1[["no_intervention"]]$genot_freqs)

        ## Just on C for this test.
        o1_C <- h1_tm

        o1_C["C", ] <- 0
        o1_C["B, C", ] <- 0
        o1_C["A, C", ] <- 0
        o1_C["C, D", ] <- 0
        o1_C["C, E", ] <- 0
        o1_C["A, B, C", ] <- 0
        o1_C["A, C, D", ] <- 0
        o1_C["A, C, E", ] <- 0
        o1_C["B, C, D", ] <- 0
        o1_C["B, C, E", ] <- 0
        o1_C["C, D, E", ] <- 0
        o1_C["A, B, C, D", ] <- 0
    o1_C["A, B, C, E", ] <- 0
        o1_C["A, C, D, E", ] <- 0
        o1_C["B, C, D, E", ] <- 0

        o1_C["WT", "WT"] <- h1_tm["WT", "C"]
        o1_C["WT", "C"] <- 0
        o1_C["A", "A"] <- h1_tm["A", "A, C"]
        o1_C["A", "A, C"] <- 0
        o1_C["B", "B"] <- h1_tm["B", "B, C"]
        o1_C["B", "B, C"] <- 0
        o1_C["D", "D"] <- h1_tm["D", "C, D"]
        o1_C["D", "C, D"] <- 0
        o1_C["E", "E"] <- h1_tm["E", "C, E"]
        o1_C["E", "C, E"] <- 0

        o1_C["A, B", "A, B"] <- h1_tm["A, B", "A, B, C"]
        o1_C["A, B", "A, B, C"] <- 0

        o1_C["A, D", "A, D"] <- h1_tm["A, D", "A, C, D"]
        o1_C["A, D", "A, C, D"] <- 0

        o1_C["A, E", "A, E"] <- h1_tm["A, E", "A, C, E"]
        o1_C["A, E", "A, C, E"] <- 0

        o1_C["B, D", "B, D"] <- h1_tm["B, D", "B, C, D"]
        o1_C["B, D", "B, C, D"] <- 0

        o1_C["B, E", "B, E"] <- h1_tm["B, E", "B, C, E"]
        o1_C["B, E", "B, C, E"] <- 0

        o1_C["D, E", "D, E"] <- h1_tm["D, E", "C, D, E"]
        o1_C["D, E", "C, D, E"] <- 0

        o1_C["A, B, D", "A, B, D"] <- h1_tm["A, B, D", "A, B, C, D"]
        o1_C["A, B, D", "A, B, C, D"] <- 0

        o1_C["A, B, E", "A, B, E"] <- h1_tm["A, B, E", "A, B, C, E"]
        o1_C["A, B, E", "A, B, C, E"] <- 0

        o1_C["A, D, E", "A, D, E"] <- h1_tm["A, D, E", "A, C, D, E"]
        o1_C["A, D, E", "A, C, D, E"] <- 0

        o1_C["B, D, E", "B, D, E"] <- h1_tm["B, D, E", "B, C, D, E"]
        o1_C["B, D, E", "B, C, D, E"] <- 0

        o1_C["A, B, D, E", "A, B, D, E"] <- 1
        o1_C["A, B, D, E", "A, B, C, D, E"] <- 0

        P_C_2 <- o1_C %^% 2
        P_C_3 <- o1_C %^% 3
        P_C_4 <- o1_C %^% 4
        P_C_5 <- o1_C %^% 5

        p_C_1 <- p_0 %*% o1_C
        p_C_2 <- p_0 %*% P_C_2
        p_C_3 <- p_0 %*% P_C_3
        p_C_4 <- p_0 %*% P_C_4
        p_C_5 <- p_0 %*% P_C_5

        p_C_0_w <- p_0 * h1$HyperHMM_used_prob.set[1]
        p_C_1_w <- p_C_1 * h1$HyperHMM_used_prob.set[2]
        p_C_2_w <- p_C_2 * h1$HyperHMM_used_prob.set[3]
        p_C_3_w <- p_C_3 * h1$HyperHMM_used_prob.set[4]
        p_C_4_w <- p_C_4 * h1$HyperHMM_used_prob.set[5]
        p_C_5_w <- p_C_5 * h1$HyperHMM_used_prob.set[6]

        pp_C <- p_C_0_w + p_C_1_w + p_C_2_w + p_C_3_w + p_C_4_w + p_C_5_w
        pp_C <- pp_C[1, ]
        (pp_C <- pp_C[pp_C > 0])

        expect_equal(pp_C, o1[["I:C"]]$genot_freqs)

        ## Explicit hitting-prob check for I:C (5-gene), as above
        expect_equal(filter_hp_keep_wt(hitting_probs_from_WT(o1_C)),
                     o1[["I:C"]]$hitting_probs_from_WT)
        expect_true("WT" %in% names(o1[["I:C"]]$hitting_probs_from_WT))

        ## Kill E
        o1_E <- h1_tm

        o1_E["E", ] <- 0
        o1_E["A, E", ] <- 0
        o1_E["B, E", ] <- 0
        o1_E["C, E", ] <- 0
        o1_E["D, E", ] <- 0
        o1_E["A, B, E", ] <- 0
        o1_E["A, C, E", ] <- 0
        o1_E["A, D, E", ] <- 0
        o1_E["B, C, E", ] <- 0
        o1_E["B, D, E", ] <- 0
        o1_E["C, D, E", ] <- 0
        o1_E["A, B, C, E", ] <- 0
        o1_E["A, B, D, E", ] <- 0
        o1_E["A, C, D, E", ] <- 0
        o1_E["B, C, D, E", ] <- 0

        o1_E["WT", "WT"] <- h1_tm["WT", "E"]
        o1_E["WT", "E"] <- 0
        o1_E["A", "A"] <- h1_tm["A", "A, E"]
        o1_E["A", "A, E"] <- 0
        o1_E["B", "B"] <- h1_tm["B", "B, E"]
        o1_E["B", "B, E"] <- 0
        o1_E["C", "C"] <- h1_tm["C", "C, E"]
        o1_E["C", "C, E"] <- 0
        o1_E["D", "D"] <- h1_tm["D", "D, E"]
        o1_E["D", "D, E"] <- 0

        o1_E["A, B", "A, B"] <- h1_tm["A, B", "A, B, E"]
        o1_E["A, B", "A, B, E"] <- 0
        o1_E["A, C", "A, C"] <- h1_tm["A, C", "A, C, E"]
        o1_E["A, C", "A, C, E"] <- 0
        o1_E["A, D", "A, D"] <- h1_tm["A, D", "A, D, E"]
        o1_E["A, D", "A, D, E"] <- 0
        o1_E["B, C", "B, C"] <- h1_tm["B, C", "B, C, E"]
        o1_E["B, C", "B, C, E"] <- 0
        o1_E["B, D", "B, D"] <- h1_tm["B, D", "B, D, E"]
        o1_E["B, D", "B, D, E"] <- 0
        o1_E["C, D", "C, D"] <- h1_tm["C, D", "C, D, E"]
        o1_E["C, D", "C, D, E"] <- 0

        o1_E["A, B, C", "A, B, C"] <- h1_tm["A, B, C", "A, B, C, E"]
        o1_E["A, B, C", "A, B, C, E"] <- 0
        o1_E["A, B, D", "A, B, D"] <- h1_tm["A, B, D", "A, B, D, E"]
        o1_E["A, B, D", "A, B, D, E"] <- 0
        o1_E["A, C, D", "A, C, D"] <- h1_tm["A, C, D", "A, C, D, E"]
        o1_E["A, C, D", "A, C, D, E"] <- 0
        o1_E["B, C, D", "B, C, D"] <- h1_tm["B, C, D", "B, C, D, E"]
        o1_E["B, C, D", "B, C, D, E"] <- 0

        o1_E["A, B, C, D", "A, B, C, D"] <- 1
        o1_E["A, B, C, D", "A, B, C, D, E"] <- 0

        P_E_2 <- o1_E %^% 2
        P_E_3 <- o1_E %^% 3
        P_E_4 <- o1_E %^% 4
        P_E_5 <- o1_E %^% 5

        p_E_1 <- p_0 %*% o1_E
        p_E_2 <- p_0 %*% P_E_2
        p_E_3 <- p_0 %*% P_E_3
        p_E_4 <- p_0 %*% P_E_4
        p_E_5 <- p_0 %*% P_E_5

        p_E_0_w <- p_0 * h1$HyperHMM_used_prob.set[1]
        p_E_1_w <- p_E_1 * h1$HyperHMM_used_prob.set[2]
        p_E_2_w <- p_E_2 * h1$HyperHMM_used_prob.set[3]
        p_E_3_w <- p_E_3 * h1$HyperHMM_used_prob.set[4]
        p_E_4_w <- p_E_4 * h1$HyperHMM_used_prob.set[5]
        p_E_5_w <- p_E_5 * h1$HyperHMM_used_prob.set[6]

        pp_E <- p_E_0_w + p_E_1_w + p_E_2_w + p_E_3_w + p_E_4_w + p_E_5_w
        pp_E <- pp_E[1, ]
        (pp_E <- pp_E[pp_E > 0])

        expect_equal(pp_E, o1[["I:E"]]$genot_freqs)

        ## Explicit hitting-prob check for I:E (5-gene), as above
        expect_equal(filter_hp_keep_wt(hitting_probs_from_WT(o1_E)),
                     o1[["I:E"]]$hitting_probs_from_WT)
        expect_true("WT" %in% names(o1[["I:E"]]$hitting_probs_from_WT))
    }
})


test_that("HyperHMM", {
    ## Note: additional checks of kill_gene for HyperHMM
    ## in intervention-TESTS.R
    ## How are the tests in intervention-TESTS.R different?
    ## They are more comprehensive. Here we only check
    ## what happens to the transition matrix.
    ## In intervention-TESTS.R we check the predicted
    ## genotypes too.

    hm1 <- Matrix(0, nrow = 8, ncol = 8, sparse = TRUE)
    colnames(hm1) <- rownames(hm1) <- allGenotypesLetter(3)
    hm1[1, 2:4] <- c(0.1, 0.2, 0.7)
    hm1[2, 5:6] <- c(0.6, 0.4)
    hm1[3, c(5, 7)] <- c(0.2, 0.8)
    hm1[4, 6:7] <- c(0.1, 0.9)
    hm1[5:7, 8] <- 1
    attr(hm1, "method_output") <- "HyperHMM_trans_mat"
    attr(hm1, "num_prob.set") <- c(0.1, 0.2, 0.3, 0.4)
    attr(hm1, "num_features") <- 3

    hm1_A <- kill_gene(hm1, "A")
    hm1_B <- kill_gene(hm1, "B")
    hm1_C <- kill_gene(hm1, "C")

    e_A <- hm1
    e_A[1, 1] <- 0.1
    e_A[1, 2] <- 0
    e_A["A", ] <- 0
    e_A["B", "B"] <- 0.2
    e_A["B", "A, B"] <- 0
    e_A["C", "C"] <- 0.1
    e_A["C", "A, C"] <- 0
    e_A["A, B", ] <- 0
    e_A["A, C", ] <- 0
    e_A["B, C", "A, B, C"] <- 0
    e_A["B, C", "B, C"] <- 1.0

    ## Why not use expect_equal? Because
    ## expect_equal sucks with all the 2e, 3e, etc, arguments.
    ## I am sick of tracking an eternally moving target.
    expect_true(isTRUE(all.equal(e_A, hm1_A, check.attributes = FALSE)))

    e_B <- hm1
    e_B["WT", "WT"] <- 0.2
    e_B["WT", "B"] <- 0
    e_B["B", ] <- 0
    e_B["A", "A"] <- 0.6
    e_B["A", "A, B"] <- 0
    e_B["C", "C"] <- 0.9
    e_B["C", "B, C"] <- 0
    e_B["A, B", ] <- 0
    e_B["B, C", ] <- 0
    e_B["A, C", "A, B, C"] <- 0
    e_B["A, C", "A, C"] <- 1.0

    expect_true(isTRUE(all.equal(e_B, hm1_B, check.attributes = FALSE)))

    e_C <- hm1
    e_C["WT", "WT"] <- 0.7
    e_C["WT", "C"] <- 0
    e_C["C", ] <- 0
    e_C["A", "A"] <- 0.4
    e_C["B", "B"] <- 0.8
    e_C["A", "A, C"] <- 0
    e_C["B", "B, C"] <- 0
    e_C["A, C", ] <- 0
    e_C["B, C", ] <- 0
    e_C["A, B", "A, B, C"] <- 0
    e_C["A, B", "A, B"] <- 1.0

    expect_true(isTRUE(all.equal(e_C, hm1_C, check.attributes = FALSE)))

    expect_warning(get_genotype_freqs_cpm(hm1_B, t = 2),
                   "With methods")
})


set.seed(NULL)

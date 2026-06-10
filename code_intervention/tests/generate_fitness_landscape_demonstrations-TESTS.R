## This file shows, via simple tests, a few features and equivalences
pwd <- getwd()
setwd("../")
source("generate_all_fitness_landscape.R")
setwd(pwd)

### Generate fitness lanscapes using DAGs

#### The transition rate matrices based on the DAG and the fitness landscape are identical
##   Also, the scaled trm = trm * a * c  (now, a * c is = 1)

test_that("The transition rate matrices based on the DAG and the fitness landscape are identical", {
    set.seed(NULL)
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

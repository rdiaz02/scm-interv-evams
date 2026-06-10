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


## We put here together a collection of tests that
##  show that the kill-gene code is equivalent to:

##  - OT and OncoBN: our default kill gene (modify DAG) same as setting the
##         conditional probabilities = 0
##
##
##  - MHN: our default kill gene (remove gene) same as setting
##         \Theta_gg = 0, \Theta_xg = 0
##  - MHN: our default kill gene same are removing the genotypes from the
##         transition rate matrix. Note that function
##         intervene_cpm_trm_rm_every_gene
##         does not go through "get_full_output" so this procedure
##         shares no key code with the other two
##         (which go through get_full_output)
##         - intervene_cpm_every_gene calls get_genotype_freqs_cpm
##         which itself calls get_full_output
##         - intervene_cpm_trm_rm_every_gene does not call
##         get_genotype_freqs_cpm for the intervened genes
##         (nor makes any other call to get_full_output)
##         but rather directly calls evamtools:::probs_from_trm
##
##
##  - CBN and H-ESBCN: our default kill gene (remove gene, and any
##         descendant that cannot exist without removed gene, from the
##         DAG) is the same as setting \lambda_g = 0.
##  - CBN and H-ESBCN: our default kill gene is the same as removing
##         the genotypes from the transition rate matrix.
##         See MHN explanation for why the different paths share no
##         relevant code.
##  - CBN and H-ESBCN: our default kill gene is the same as intervening on
##         the fitness landscape: intervene_fitness_landscape_every_gene

## "same as": predicted genotype frequencies are the same

## However, *much more comprehensive* tests, including equivalence
## between all the possible ways of killing for each method, are
## available from the kill-CBN/OT/OncoBN/HESBCN/MHN/HyperHMM-TESTS.R
## files. Note that for HyperHMM there are no equivalences.

## File intervene-fitness-landscapes-TESTS.R contains
## tests of the intervention on fitness landscapes directly
## (its focus is not testing intervene_cpm_every_gene).


## About testing predicted genotype frequencies:
## That is tested in evamtools. In evamtools:
##
## For OT, evamtools tests OT_model_2_output the same function used here
##         in get_full_output
##
## For OncoBN, evamtools tests OncoBN_model_2_output, the same function
##         used here in get_full_output
##
## For CBN, HESBCN, and MHN evamtools tests probs_from_trm, same function
##      used here (called in get_full_output), for any of CBN, HESBCN, MHN
##      And creating the transition rate matrix is also tested in evamtools.



## options(intervention_every_gene_cores = parallel::detectCores())
library(testthat)
pwd <- getwd()
setwd("../")
source("intervention.R")
source("generate_all_fitness_landscape.R")
setwd(pwd)

set.seed(NULL)


### Utility functions

## Two intervention outputs -> stop unless identical
stop_unless_intervention_identical <- function(x, y) {
    stopifnot(all(unlist(lapply(1:length(x),
                                function(i)
                                    all.equal(x[[i]], y[[i]])))))
}



### An example of killing a gene from which all depend
##  This uses the redefined adjm_rm_no_access
##  This is a test case included in kill-gene.R
##  but here we also use killing via setting parameters to 0

test_that("kill gene from which all depend", {
    local_edition(3)
    ## CBN
    m1 <- data.frame(From = c("Root", "A", "A", "A"),
                     To   = c("A",    "B", "C", "D"),
                     rerun_lambda = 1:4)
    ## OT
    m2 <- data.frame(From = c("Root", "B", "B", "B"),
                     To   = c("B",    "A", "C", "D"),
                     OT_edgeWeight = rep(0.3, 4))
    ## OncoBN
    m3 <- data.frame(From = c("Root", "Z", "Z", "Z"),
                     To   = c("Z",    "B", "C", "D"),
                     theta = rep(0.2, 4),
                     Relation = "Single")
    ## H-ESBCN
    m4 <- data.frame(From = c("Root", "C", "C", "C"),
                     To   = c("C",    "B", "A", "D"),
                     Lambdas = 1:4,
                     Relation = "Single")
    expect_true(nrow(suppressWarnings(kill_gene(m1, "A"))) == 0)
    expect_true(nrow(suppressWarnings(kill_gene(m2, "B"))) == 0)
    expect_true(nrow(suppressWarnings(kill_gene(m3, "Z"))) == 0)
    expect_true(nrow(suppressWarnings(kill_gene(m4, "C"))) == 0)

    expect_equal(kill_gene(m1, "B"), m1[-2, ])

    just_wt <- c(WT = 1)
    just_wt_full <- list(genot_freqs = just_wt, hitting_probs_from_WT = c(WT = 1.0))
    ## Structural kill reduces model dimensions, params-to-0 keeps them.
    ## Filter > 0 before comparing when dimensions may differ.
    filter_preds <- function(x) list(
        genot_freqs = x$genot_freqs[x$genot_freqs > 0],
        hitting_probs_from_WT = x$hitting_probs_from_WT[x$hitting_probs_from_WT > 0])

    expect_identical(
        get_genotype_freqs_cpm(suppressWarnings(kill_gene(m1, "A"))),
        just_wt_full)

    expect_identical(
        get_genotype_freqs_cpm(suppressWarnings(kill_gene(m2, "B"))),
        just_wt_full)

    expect_identical(
        get_genotype_freqs_cpm(suppressWarnings(kill_gene(m3, "Z"))),
        just_wt_full)

    expect_identical(
      get_genotype_freqs_cpm(suppressWarnings(kill_gene(m4, "C"))),
      just_wt_full)

    ## When setting parameters to 0. For pedagogical purposes, do
    ## not suppres warnings nor messages

    kill_gene_by_params_to_0(m1, "A", verbose = TRUE)
    kill_gene_by_params_to_0(m2, "B", verbose = TRUE)
    kill_gene_by_params_to_0(m3, "Z", verbose = TRUE)
    kill_gene_by_params_to_0(m4, "C", verbose = TRUE)

    expect_warning(get_genotype_freqs_cpm(kill_gene_by_params_to_0(m1, "A", verbose = TRUE)),
                   "weighted_fgraph contains unreachable destinations")
    expect_warning(get_genotype_freqs_cpm(kill_gene_by_params_to_0(m2, "B", verbose = TRUE)),
                   "weighted_fgraph contains unreachable destinations")
    get_genotype_freqs_cpm(kill_gene_by_params_to_0(m3, "Z", verbose = TRUE))
    expect_warning(get_genotype_freqs_cpm(kill_gene_by_params_to_0(m4, "C", verbose = TRUE)),
                   "weighted_fgraph contains unreachable destinations")

    ## Same predictions
    pm1 <- suppressWarnings(get_genotype_freqs_cpm(kill_gene_by_params_to_0(m1, "A")))
    pm2 <- suppressWarnings(get_genotype_freqs_cpm(kill_gene_by_params_to_0(m2, "B")))
    pm3 <- suppressWarnings(get_genotype_freqs_cpm(kill_gene_by_params_to_0(m3, "Z")))
    pm4 <- suppressWarnings(get_genotype_freqs_cpm(kill_gene_by_params_to_0(m4, "C")))

    expect_equal(filter_preds(pm1), just_wt_full)
    expect_equal(filter_preds(pm2), just_wt_full)
    expect_equal(filter_preds(pm3), just_wt_full)
    expect_equal(filter_preds(pm4), just_wt_full)

    ## For CBN and HESBCN, two additional methods: landscape and trans.
    ## rate matrix
    ## First, utility to create the landscape from the bare model
    ## Based on generate_n_f_landscape_requir.
    ## Remember for fitness landscape, gene names must be sequential letters
    ## starting from A (though A need not be "the first" in the DAG).
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

    m1_tr_i <- intervene_cpm_trm_rm_every_gene(c(list(CBN_model = m1),
                                                 get_full_output(m1)),
                                               "CBN")
    m1_fl <- cbn_to_landscape_obj(m1, 4)
    m1_fl_i <- intervene_fitness_landscape_every_gene(m1_fl)
    m1_tr_i[["I:A"]]$genot_freqs
    expect_equal(m1_tr_i[["I:A"]]$genot_freqs, just_wt)
    expect_equal(m1_tr_i[["I:A"]]$hitting_probs_from_WT, just_wt_full$hitting_probs_from_WT)
    ## Full structure comparison (genot_freqs AND hitting_probs must match)
    expect_equal(m1_tr_i, m1_fl_i)

    m4_tr_i <- intervene_cpm_trm_rm_every_gene(c(list(HESBCN_model = m4),
                                                 get_full_output(m4)),
                                               "HESBCN")
    m4_fl <- cbn_to_landscape_obj(m4, 4)
    m4_fl_i <- intervene_fitness_landscape_every_gene(m4_fl)
    m4_tr_i[["I:C"]]$genot_freqs
    expect_equal(m4_tr_i[["I:C"]]$genot_freqs, just_wt)
    expect_equal(m4_tr_i[["I:C"]]$hitting_probs_from_WT, just_wt_full$hitting_probs_from_WT)
    ## Full structure comparison (genot_freqs AND hitting_probs must match)
    expect_equal(m4_tr_i, m4_fl_i)
})




### "Standard" DAG intervention identical to setting parameters to 0

## Run and check a few times
local({
    set.seed(NULL)
    total_iters <- 10
    cat("Standard DAG intervention identical to setting parameters to 0")
    for (i in 1:total_iters) {
        cat("\n #################### Doing iteration ", i, "\n\n")
#### Models and fitness landscapes

        ## CPM models simulated from scratch

        rcbn <- evamtools::random_evam(7, model = "CBN")
        rhes <- evamtools::random_evam(7, model = "HESBCN")
        rot <- evamtools::random_evam(7, model = "OT", ot_oncobn_epos = 0)
        rdbn_d <- evamtools::random_evam(7, model = "OncoBN", ot_oncobn_epos = 0,
                                         oncobn_model = "DBN")
        rdbn_c <- evamtools::random_evam(7, model = "OncoBN", ot_oncobn_epos = 0,
                                         oncobn_model = "CBN")
        rmhn <- evamtools::random_evam(7, model = "MHN")


        ## CPM models simulated ensuring restrictions fulfilled
        ## The model itself is in x[[1]][[11]][[1]]
        ## MHN cannot be generated this way (no model <-> fitness landscape)
        ## OT and OncoBN subsumed in CBN and/or HESBCN
        rcbn_f <- suppressMessages(generate_n_f_landscape_requir(1, 7, "CBN"))
        rhes_f <- suppressMessages(
            generate_n_f_landscape_requir(1, 7, "HESBCN",
                                          hesbcn_relations = c("AND", "OR", "XOR")))



#### "Standard" kill-gene procedure
        i_rcbn <- intervene_cpm_every_gene(rcbn, "CBN", verbose = TRUE)
        i_rhes <- intervene_cpm_every_gene(rhes, "HESBCN", verbose = TRUE)
        i_rcbn_f <- intervene_cpm_every_gene(rcbn_f[[1]]$other, "CBN", verbose = TRUE)
        i_rhes_f <- intervene_cpm_every_gene(rhes_f[[1]]$other, "HESBCN", verbose = TRUE)
        i_rot <- intervene_cpm_every_gene(rot, "OT", verbose = TRUE)
        i_rdbn_d <- intervene_cpm_every_gene(rdbn_d, "OncoBN", verbose = TRUE)
        i_rdbn_c <- intervene_cpm_every_gene(rdbn_c, "OncoBN", verbose = TRUE)
        i_rmhn <- intervene_cpm_every_gene(rmhn, "MHN", verbose = TRUE)

#### Kill by setting parameters to 0

        ## Yes, expect, as should be the case, warning for unreachable destinations
        i_0_rcbn <-
            intervene_cpm_every_gene(rcbn, "CBN",
                                     kill_gene_funct = kill_gene_by_params_to_0,
                                     verbose = TRUE)
        i_0_rhes <-
            intervene_cpm_every_gene(rhes, "HESBCN",
                                     kill_gene_funct = kill_gene_by_params_to_0,
                                     verbose = TRUE)

        i_0_rcbn_f <-
            intervene_cpm_every_gene(rcbn_f[[1]]$other, "CBN",
                                     kill_gene_funct = kill_gene_by_params_to_0,
                                     verbose = TRUE)
        i_0_rhes_f <-
            intervene_cpm_every_gene(rhes_f[[1]]$other, "HESBCN",
                                     kill_gene_funct = kill_gene_by_params_to_0,
                                     verbose = TRUE)
        i_0_rot <-
            intervene_cpm_every_gene(rot, "OT",
                                     kill_gene_funct = kill_gene_by_params_to_0,
                                     verbose = TRUE)
        i_0_rdbn_d <-
            intervene_cpm_every_gene(rdbn_d, "OncoBN",
                                     kill_gene_funct = kill_gene_by_params_to_0,
                                     verbose = TRUE)
        i_0_rdbn_c <-
            intervene_cpm_every_gene(rdbn_c, "OncoBN",
                                     kill_gene_funct = kill_gene_by_params_to_0,
                                     verbose = TRUE)
        i_0_rmhn <-
            intervene_cpm_every_gene(rmhn, "MHN",
                                     kill_gene_funct = kill_gene_by_params_to_0,
                                     verbose = TRUE)

        ## Yes, it detects failures. For example
        ## stop_unless_intervention_identical(i_0_rcbn, i_0_rot)
        ## stop_unless_intervention_identical(i_0_rdbn_c, i_0_rdbn_d)

        stop_unless_intervention_identical(i_0_rcbn, i_rcbn)
        stop_unless_intervention_identical(i_0_rhes, i_rhes)
        stop_unless_intervention_identical(i_0_rcbn_f, i_rcbn_f)
        stop_unless_intervention_identical(i_0_rhes_f, i_rhes_f)

        stop_unless_intervention_identical(i_0_rot, i_rot)
        stop_unless_intervention_identical(i_0_rdbn_d, i_rdbn_d)
        stop_unless_intervention_identical(i_0_rdbn_c, i_rdbn_c)
        stop_unless_intervention_identical(i_0_rmhn, i_rmhn)
    }
})


### CBN and H-ESBCN: intervene by modifying the fitness landscape identical to DAG intervention

local({
    set.seed(NULL)
    ## Run a few times
    total_iters <- 10
    cat("CBN and H-ESBCN: intervene by modifying the fitness landscape identical to DAG intervention")
    for (i in 1:total_iters) {
        cat("\n #################### Doing iteration ", i, "\n\n")

        ## Generate fitness landscapes
        rcbn_f_2 <- suppressMessages(generate_n_f_landscape_requir(1, 7, "CBN"))
        rhes_f_2 <- suppressMessages(
            generate_n_f_landscape_requir(1, 7, "HESBCN",
                                          hesbcn_relations = c("AND", "OR", "XOR")))

        ## Intervene on the fitness landscape
        fl_i_rcbn_f_2 <- intervene_fitness_landscape_every_gene(rcbn_f_2[[1]])
        fl_i_rhes_f_2 <- intervene_fitness_landscape_every_gene(rhes_f_2[[1]])

        ## Intervene via the DAG
        i_rcbn_f_2 <- intervene_cpm_every_gene(rcbn_f_2[[1]]$other, "CBN",
                                               verbose = TRUE)
        i_rhes_f_2 <- intervene_cpm_every_gene(rhes_f_2[[1]]$other, "HESBCN",
                                               verbose = TRUE)

        ## Check identical
        stop_unless_intervention_identical(fl_i_rcbn_f_2, i_rcbn_f_2)
        stop_unless_intervention_identical(fl_i_rhes_f_2, i_rhes_f_2)
    }
})


### CBN, H-ESBCN, MHN: remove, from transition rate matrix, genotypes with intervened gene
## Remember: DO NOT USE THIS in general as limited to a few methods

local({
    set.seed(NULL)
    ## Run a few times
    total_iters <- 10
    cat("CBN, H-ESBCN, MHN: remove, from transition rate matrix, genotypes with intervened gene")
    for (i in 1:total_iters) {
        cat("\n #################### Doing iteration ", i, "\n\n")
        ## If you set.seed(1) you get CBN to have all depend on A

        ## Generate models
        rcbn <- evamtools::random_evam(7, model = "CBN")
        rhes <- evamtools::random_evam(7, model = "HESBCN")
        rmhn <- evamtools::random_evam(7, model = "MHN")

        ## Intervention with standard procedure
        i_rcbn <- intervene_cpm_every_gene(rcbn, "CBN", verbose = TRUE)
        i_rhes <- intervene_cpm_every_gene(rhes, "HESBCN", verbose = TRUE)
        i_rmhn <- intervene_cpm_every_gene(rmhn, "MHN", verbose = TRUE)

        ## Intervention rm genotypes from trm
        i_rm_trm_rcbn <- intervene_cpm_trm_rm_every_gene(rcbn, "CBN")
        i_rm_trm_rhes <- intervene_cpm_trm_rm_every_gene(rhes, "HESBCN")
        i_rm_trm_rmhn <- intervene_cpm_trm_rm_every_gene(rmhn, "MHN")

        stop_unless_intervention_identical(i_rcbn, i_rm_trm_rcbn)
        stop_unless_intervention_identical(i_rhes, i_rm_trm_rhes)
        stop_unless_intervention_identical(i_rmhn, i_rm_trm_rmhn)
    }
})

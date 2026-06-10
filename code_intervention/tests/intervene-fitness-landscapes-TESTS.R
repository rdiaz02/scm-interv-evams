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


### What this is is
## Some tests that the code for intervention on fitness landscapes
## behaves sensibly.
## Correspondence with the other killing methods
## is also explored exhaustively in the kill-CBN-TESTS.R
## and kill-HESBCN-TESTS.R.



## options(intervention_every_gene_cores = parallel::detectCores())
library(testthat)
pwd <- getwd()
setwd("../")
source("intervention.R")
source("generate_all_fitness_landscape.R")
setwd(pwd)

set.seed(NULL)

### Intervene on fitness landscape: example with only two viable genots


test_that("Intervene on fitness landscape: example with only two viable genots", {
  #### Example 1
  ## Create a fitness landscape with only A and AB viables
  ## This should happen
  ##   Intervene on A: only WT remains
  ##   Intervene on B: only WT and A
  ##   Intervene on C to G: no effect.

  r1 <- rfitness(7)
  r1[, 8] <- 1e-9
  r1[c(1, 2, 9), 8] <- c(1, 1.1, 1.2)

  ## Prepare the object for input
  fit_land_letters <- no_evam_rfitness_to_letter(r1)
  tmp <- no_evam_genots_2_fgraph_and_trans_mat_rf(fit_land_letters)
  rel_diff <- tmp$relative_fitness_differences
  the_c <- 3
  r1o <- list(fitness_landscape = r1,
              c = the_c,
              trm_scaled = get_square_matrix(the_c * rel_diff, 0))

  ## Intervene
  r1i <- intervene_fitness_landscape_every_gene(r1o)

  ## Checks
  just_wt <- c(WT = 1)
  just_wt_hp <- c(WT = 1.0)

  expect_true(isTRUE(all.equal(r1i[["I:A"]]$genot_freqs, just_wt)))
  expect_true(isTRUE(all.equal(r1i[["I:A"]]$hitting_probs_from_WT, just_wt_hp)))
  for (g in c("C", "D", "E", "F", "G")) {
      expect_true(isTRUE(all.equal(r1i[["no_intervention"]]$genot_freqs,
                                   r1i[[paste0("I:", g)]]$genot_freqs)))
      expect_true(isTRUE(all.equal(r1i[["no_intervention"]]$hitting_probs_from_WT,
                                   r1i[[paste0("I:", g)]]$hitting_probs_from_WT)))
  }

  ## Frequency of A should be larger when B is killed
  expect_true(r1i[["I:B"]]$genot_freqs[2] > r1i[["no_intervention"]]$genot_freqs[2])
  ## Frequency of WT should be the same when B is killed and no interv.
  expect_equal(r1i[["I:B"]]$genot_freqs[1], r1i[["no_intervention"]]$genot_freqs[1])
})


test_that("Intervene on fitness landscape: example with only two viable genots, example 2", {
  #### Example 2
  ## Create a fitness landscape with only B, BC, BD viables
  ## This should happen
  ##   Intervene on B: only WT remains
  ##   Intervene on C: B and BD
  ##   Intervene on D: B and BC
  ##   Intervene on others: no effect

  r2 <- rfitness(7)
  r2[, 8] <- 1e-9
  r2[c(1, 3, 15, 16), 8] <- c(1, 1.02, 1.15, 1.1)

  ## Prepare the object for input
  fit_land_letters2 <- no_evam_rfitness_to_letter(r2)
  tmp2 <- no_evam_genots_2_fgraph_and_trans_mat_rf(fit_land_letters2)
  rel_diff2 <- tmp2$relative_fitness_differences
  the_c2 <- 2
  r2o <- list(fitness_landscape = r2,
              c = the_c2,
              trm_scaled = get_square_matrix(the_c2 * rel_diff2, 0))

  ## Intervene
  r2i <- intervene_fitness_landscape_every_gene(r2o)

  ## Checks
  just_wt <- c(WT = 1)
  just_wt_hp <- c(WT = 1.0)

  expect_true(isTRUE(all.equal(r2i[["I:B"]]$genot_freqs, just_wt)))
  expect_true(isTRUE(all.equal(r2i[["I:B"]]$hitting_probs_from_WT, just_wt_hp)))
  for (g in c("A", "E", "F", "G")) {
      expect_true(isTRUE(all.equal(r2i[["no_intervention"]]$genot_freqs,
                                   r2i[[paste0("I:", g)]]$genot_freqs)))
      expect_true(isTRUE(all.equal(r2i[["no_intervention"]]$hitting_probs_from_WT,
                                   r2i[[paste0("I:", g)]]$hitting_probs_from_WT)))
  }

  expect_true(isTRUE(all.equal(names(r2i[["I:C"]]$genot_freqs), c("WT", "B", "B, D"))))
  expect_true(isTRUE(all.equal(names(r2i[["I:D"]]$genot_freqs), c("WT", "B", "B, C"))))

  ## Frequency of B should be larger when C or D are killed
  expect_true(r2i[["I:C"]]$genot_freqs[2] > r2i[["no_intervention"]]$genot_freqs[2])
  expect_true(r2i[["I:D"]]$genot_freqs[2] > r2i[["no_intervention"]]$genot_freqs[2])

  ## Frequency of WT should be the same when C is killed and no interv.
  expect_equal(r2i[["I:C"]]$genot_freqs[1], r2i[["no_intervention"]]$genot_freqs[1])
  ## Ditto for D
  expect_equal(r2i[["I:D"]]$genot_freqs[1], r2i[["no_intervention"]]$genot_freqs[1])
})

test_that("Intervene on fitness landscape: example with only two viable genots, example 3", {
  #### Example 3
  ## Create a fitness landscape with only E, AE, ACE viables
  ## This should happen
  ##   Intervene on E: only WT remains
  ##   Intervene on A: E
  ##   Intervene on C: E and AE
  ##   Intervene on others: no effect

  r3 <- rfitness(7)
  r3[, 8] <- 1e-9
  r3[c(1, 6, 12, 36), 8] <- c(1, 1.02, 1.5, 2.1)

  ## Prepare the object for input
  fit_land_letters2 <- no_evam_rfitness_to_letter(r3)
  tmp2 <- no_evam_genots_2_fgraph_and_trans_mat_rf(fit_land_letters2)
  rel_diff3 <- tmp2$relative_fitness_differences
  the_c3 <- 7.5
  r3o <- list(fitness_landscape = r3,
              c = the_c3,
              trm_scaled = get_square_matrix(the_c3 * rel_diff3, 0))

  ## Intervene
  r3i <- intervene_fitness_landscape_every_gene(r3o)

  ## Checks
  just_wt <- c(WT = 1)
  just_wt_hp <- c(WT = 1.0)

  expect_true(isTRUE(all.equal(r3i[["I:E"]]$genot_freqs, just_wt)))
  expect_true(isTRUE(all.equal(r3i[["I:E"]]$hitting_probs_from_WT, just_wt_hp)))
  for (g in c("B", "D", "F", "G")) {
      expect_true(isTRUE(all.equal(r3i[["no_intervention"]]$genot_freqs,
                                   r3i[[paste0("I:", g)]]$genot_freqs)))
      expect_true(isTRUE(all.equal(r3i[["no_intervention"]]$hitting_probs_from_WT,
                                   r3i[[paste0("I:", g)]]$hitting_probs_from_WT)))
  }

  expect_true(isTRUE(all.equal(names(r3i[["I:A"]]$genot_freqs), c("WT", "E"))))
  expect_true(isTRUE(all.equal(names(r3i[["I:C"]]$genot_freqs), c("WT", "E", "A, E"))))

  ## Frequency of E should be larger when A is killed
  expect_true(r3i[["I:A"]]$genot_freqs[2] > r3i[["no_intervention"]]$genot_freqs[2])
  ## And that of A, E larger when C is killed
  expect_true(r3i[["I:C"]]$genot_freqs[3] > r3i[["no_intervention"]]$genot_freqs[3])

  ## Frequency of WT should be the same when A is killed and no interv.
  expect_equal(r3i[["I:A"]]$genot_freqs[1], r3i[["no_intervention"]]$genot_freqs[1])
  ## Ditto for C
  expect_equal(r3i[["I:C"]]$genot_freqs[1], r3i[["no_intervention"]]$genot_freqs[1])
})


### Intervene on fitness landscape: no intervention return same genots as original
##   with CPM-based landscape, also identical to DAG-based genotype frequencies
##   Note that this second is also tested in kill-gene-equivalences-TESTS.R
##     CBN and H-ESBCN: intervene by modifying the fitness landscape identical
##     to DAG intervention

test_that(paste0("F landscape: no intervention same genots as original;",
                 " with CPM-based landscape, also identical to DAG-based ",
                 " genotype frequencies"),
{

  local_edition(3)

  ## Helper: compute hitting_probs_from_WT from trm_scaled independently
  ## (same logic as in intervene_fitness_landscape_every_gene no-intervention case,
  ## including the WT-preserving filter via filter_hp_keep_wt).
  fl_hp_from_trm <- function(fl_obj) {
      trm <- fl_obj$trm_scaled
      rs <- rowSums(trm)
      embedded <- trm
      embedded[rs > 0, ] <- trm[rs > 0, ] / rs[rs > 0]
      hp <- hitting_probs_from_WT(embedded)
      return(filter_hp_keep_wt(hp))
  }

  ## NK 1
  fl_nk_1 <- generate_n_f_landscape_requir(1, 7, "NK", K = 1)
  fl_nk_1_pred <- genots_from_trm(fl_nk_1[[1]]$trm_scaled)
  fl_nk_1_pred <- fl_nk_1_pred[fl_nk_1_pred > 0]

  i_fl_nk_1 <- intervene_fitness_landscape_every_gene(fl_nk_1[[1]])
  expect_true(all.equal(fl_nk_1_pred, i_fl_nk_1[["no_intervention"]]$genot_freqs))
  expect_true(all.equal(fl_hp_from_trm(fl_nk_1[[1]]), i_fl_nk_1[["no_intervention"]]$hitting_probs_from_WT))

  ## NK 2
  fl_nk_2 <- generate_n_f_landscape_requir(1, 7, "NK", K = 2)
  fl_nk_2_pred <- genots_from_trm(fl_nk_2[[1]]$trm_scaled)
  fl_nk_2_pred <- fl_nk_2_pred[fl_nk_2_pred > 0]

  i_fl_nk_2 <- intervene_fitness_landscape_every_gene(fl_nk_2[[1]])
  expect_true(all.equal(fl_nk_2_pred, i_fl_nk_2[["no_intervention"]]$genot_freqs))
  expect_true(all.equal(fl_hp_from_trm(fl_nk_2[[1]]), i_fl_nk_2[["no_intervention"]]$hitting_probs_from_WT))

  ## NK 3
  fl_nk_3 <- generate_n_f_landscape_requir(1, 7, "NK", K = 3)
  fl_nk_3_pred <- genots_from_trm(fl_nk_3[[1]]$trm_scaled)
  fl_nk_3_pred <- fl_nk_3_pred[fl_nk_3_pred > 0]

  i_fl_nk_3 <- intervene_fitness_landscape_every_gene(fl_nk_3[[1]])
  expect_true(all.equal(fl_nk_3_pred, i_fl_nk_3[["no_intervention"]]$genot_freqs))
  expect_true(all.equal(fl_hp_from_trm(fl_nk_3[[1]]), i_fl_nk_3[["no_intervention"]]$hitting_probs_from_WT))

  ## RMF
  fl_rmf <- generate_n_f_landscape_requir(1, 7, "RMF")
  fl_rmf_pred <- genots_from_trm(fl_rmf[[1]]$trm_scaled)
  fl_rmf_pred <- fl_rmf_pred[fl_rmf_pred > 0]

  i_fl_rmf <- intervene_fitness_landscape_every_gene(fl_rmf[[1]])
  expect_true(all.equal(fl_rmf_pred, i_fl_rmf[["no_intervention"]]$genot_freqs))
  expect_true(all.equal(fl_hp_from_trm(fl_rmf[[1]]), i_fl_rmf[["no_intervention"]]$hitting_probs_from_WT))

  ## CBN
  fl_cbn <- generate_n_f_landscape_requir(1, 7, "CBN")
  fl_cbn_pred <- genots_from_trm(fl_cbn[[1]]$trm_scaled)
  fl_cbn_pred <- fl_cbn_pred[fl_cbn_pred > 0]

  i_fl_cbn <- intervene_fitness_landscape_every_gene(fl_cbn[[1]])
  expect_true(all.equal(fl_cbn_pred, i_fl_cbn[["no_intervention"]]$genot_freqs))
  expect_true(all.equal(fl_hp_from_trm(fl_cbn[[1]]), i_fl_cbn[["no_intervention"]]$hitting_probs_from_WT))

  cpm_pred_cbn <- fl_cbn[[1]]$other[["CBN_predicted_genotype_freqs"]]
  cpm_pred_cbn <- cpm_pred_cbn[cpm_pred_cbn > 0]
  expect_true(all.equal(fl_cbn_pred, cpm_pred_cbn))


  ## HESBCN, AND and XOR
  fl_hesbcn <- generate_n_f_landscape_requir(1, 7, "HESBCN",
                                             hesbcn_relations = c("AND", "XOR"))
  fl_hesbcn_pred <- genots_from_trm(fl_hesbcn[[1]]$trm_scaled)
  fl_hesbcn_pred <- fl_hesbcn_pred[fl_hesbcn_pred > 0]

  i_fl_hesbcn <- intervene_fitness_landscape_every_gene(fl_hesbcn[[1]])
  expect_true(all.equal(fl_hesbcn_pred, i_fl_hesbcn[["no_intervention"]]$genot_freqs))
  expect_true(all.equal(fl_hp_from_trm(fl_hesbcn[[1]]), i_fl_hesbcn[["no_intervention"]]$hitting_probs_from_WT))

  cpm_pred_hesbcn <- fl_hesbcn[[1]]$other[["HESBCN_predicted_genotype_freqs"]]
  cpm_pred_hesbcn <- cpm_pred_hesbcn[cpm_pred_hesbcn > 0]
  expect_true(all.equal(fl_hesbcn_pred, cpm_pred_hesbcn))


  ## HESBCN, OR and XOR
  fl_hesbcn2 <- generate_n_f_landscape_requir(1, 7, "HESBCN",
                                              hesbcn_relations = c("OR", "XOR"))
  fl_hesbcn2_pred <- genots_from_trm(fl_hesbcn2[[1]]$trm_scaled)
  fl_hesbcn2_pred <- fl_hesbcn2_pred[fl_hesbcn2_pred > 0]

  i_fl_hesbcn2 <- intervene_fitness_landscape_every_gene(fl_hesbcn2[[1]])
  expect_true(all.equal(fl_hesbcn2_pred, i_fl_hesbcn2[["no_intervention"]]$genot_freqs))
  expect_true(all.equal(fl_hp_from_trm(fl_hesbcn2[[1]]), i_fl_hesbcn2[["no_intervention"]]$hitting_probs_from_WT))

  cpm_pred_hesbcn2 <- fl_hesbcn2[[1]]$other[["HESBCN_predicted_genotype_freqs"]]
  cpm_pred_hesbcn2 <- cpm_pred_hesbcn2[cpm_pred_hesbcn2 > 0]
  expect_true(all.equal(fl_hesbcn2_pred, cpm_pred_hesbcn2))


  ## HESBCN, AND, OR,  and XOR
  fl_hesbcn3 <- generate_n_f_landscape_requir(1, 7, "HESBCN",
                                              hesbcn_relations = c("AND", "OR", "XOR"))
  fl_hesbcn3_pred <- genots_from_trm(fl_hesbcn3[[1]]$trm_scaled)
  fl_hesbcn3_pred <- fl_hesbcn3_pred[fl_hesbcn3_pred > 0]

  i_fl_hesbcn3 <- intervene_fitness_landscape_every_gene(fl_hesbcn3[[1]])
  expect_true(all.equal(fl_hesbcn3_pred, i_fl_hesbcn3[["no_intervention"]]$genot_freqs))
  expect_true(all.equal(fl_hp_from_trm(fl_hesbcn3[[1]]), i_fl_hesbcn3[["no_intervention"]]$hitting_probs_from_WT))

  cpm_pred_hesbcn3 <- fl_hesbcn3[[1]]$other[["HESBCN_predicted_genotype_freqs"]]
  cpm_pred_hesbcn3 <- cpm_pred_hesbcn3[cpm_pred_hesbcn3 > 0]
  expect_true(all.equal(fl_hesbcn3_pred, cpm_pred_hesbcn3))

  ## BEWARE: we are NOT using OT or OncoBN fitness landscapes
  ## in this paper (unclear interpretation). But the code does what
  ## it should.

  ## OT
  fl_ot <- suppressWarnings(generate_n_f_landscape_requir(1, 7, "OT"))
  fl_ot_pred <- genots_from_trm(fl_ot[[1]]$trm_scaled)
  fl_ot_pred <- fl_ot_pred[fl_ot_pred > 0]

  i_fl_ot <- intervene_fitness_landscape_every_gene(fl_ot[[1]])
  expect_true(all.equal(fl_ot_pred, i_fl_ot[["no_intervention"]]$genot_freqs))
  expect_true(all.equal(fl_hp_from_trm(fl_ot[[1]]), i_fl_ot[["no_intervention"]]$hitting_probs_from_WT))

  cpm_pred_ot <- fl_ot[[1]]$other[["OT_predicted_genotype_freqs"]]
  cpm_pred_ot <- cpm_pred_ot[cpm_pred_ot > 0]
  ## You should NOT expect these to be equal, as
  ## scaling and fitness landscape relationships need
  ## not hold with untimed models.
  ## expect_true(all.equal(fl_ot_pred, cpm_pred_ot))

  ## OncoBN
  fl_oncobn <- suppressWarnings(generate_n_f_landscape_requir(1, 7, "OncoBN"))
  fl_oncobn_pred <- genots_from_trm(fl_oncobn[[1]]$trm_scaled)
  fl_oncobn_pred <- fl_oncobn_pred[fl_oncobn_pred > 0]

  i_fl_oncobn <- intervene_fitness_landscape_every_gene(fl_oncobn[[1]])
  expect_true(all.equal(fl_oncobn_pred, i_fl_oncobn[["no_intervention"]]$genot_freqs))
  expect_true(all.equal(fl_hp_from_trm(fl_oncobn[[1]]), i_fl_oncobn[["no_intervention"]]$hitting_probs_from_WT))

  cpm_pred_oncobn <- fl_oncobn[[1]]$other[["OncoBN_predicted_genotype_freqs"]]
  cpm_pred_oncobn <- cpm_pred_oncobn[cpm_pred_oncobn > 0]
  ## You should NOT expect these to be equal, as
  ## scaling and fitness landscape relationships need
  ## not hold with untimed models.
  ## expect_true(all.equal(fl_oncobn_pred, cpm_pred_oncobn))
})


#### Killing in the fitness landscape, leaving only one genotype as viable
##  This does not affect OT or OncoBN or MHN directly, since
##  fitness landscapes are not generated from them.
##  And there are tests
##  for the "only two viable genots" that, after intervention,
##  leave only two, one, or no genotypes.
##  There is no real need to test this, but here we go anyway
##  Getting an H-ESBCN that has this pattern is very unlikely
##  and since this test is not really necessary, anyway, I stop trying.

test_that("Killing in the fitness landscape, leaving only one genotype as viable", {

  stop_unless_intervention_identical <- function(x, y) {
    expect_true(isTRUE(all(unlist(lapply(1:length(x),
                                         function(i)
                                           all.equal(x[[i]], y[[i]]))))))
  }

  stop_unless_intervention_identical_rm_0 <- function(x, y) {
    expect_true(isTRUE(all(unlist(lapply(1:length(x),
                                         function(i) {
                                           x1 <- x[[i]]
                                           y1 <- y[[i]]
                                           x1 <- x1[x1 > 0]
                                           y1 <- y1[y1 > 0]
                                           all.equal(x1, y1)}
                                         )))))
  }

  ## If we kill A, only G is viable
  set.seed(1)
  rcbn_f <- suppressMessages(generate_n_f_landscape_requir(1, 7, "CBN"))
  ## See how we kill A, only G is viable
  (rcbn_f[[1]][["other"]][["CBN_model"]])
  ##
  rcbn_i_l <- intervene_fitness_landscape_every_gene(rcbn_f[[1]])
  rcbn_i_s <- intervene_cpm_every_gene(rcbn_f[[1]]$other, "CBN", verbose = TRUE)
  rcbn_i_p <- intervene_cpm_every_gene(rcbn_f[[1]]$other, "CBN",
                                       kill_gene_funct = kill_gene_by_params_to_0,
                                       verbose = TRUE)
  ##
  gf_l <- lapply(rcbn_i_l, `[[`, "genot_freqs")
  gf_s <- lapply(rcbn_i_s, `[[`, "genot_freqs")
  gf_p <- lapply(rcbn_i_p, `[[`, "genot_freqs")
  hp_l <- lapply(rcbn_i_l, `[[`, "hitting_probs_from_WT")
  hp_s <- lapply(rcbn_i_s, `[[`, "hitting_probs_from_WT")
  hp_p <- lapply(rcbn_i_p, `[[`, "hitting_probs_from_WT")

  stop_unless_intervention_identical(gf_l, gf_s)
  stop_unless_intervention_identical(gf_l, gf_p)

  stop_unless_intervention_identical(hp_l, hp_s)

  stop_unless_intervention_identical_rm_0(hp_l, hp_p)


  set.seed(NULL)
})




set.seed(NULL)

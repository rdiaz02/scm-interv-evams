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



## TL;DR

## HyperHMM (because of how Baum-Welch works) can return matrices where
## more than a 1/3 of the entries are below 1e-12. For a few 20 fits I
## looked at, probabilities below 1e-12 were 22%–50%, median 38%, and the
## smallest positive entry seen: 4.8e-320. When this happens,
## markovchain::hittingProbabilities can run into trouble. The cause seems
## to be the huge dynamic range (not the multiple absorbing states) that
## can appear with HyperHMM: thresholding at 1e-12 makes the trouble go
## away while leaving all of those absorbing states in place (it even
## creates a few more). I reported this
## https://github.com/spedygiorgio/markovchain/issues/233
## on 2026-07-04. I reported a problem with the fix on 2026-09-17.
## More fixes on 2026-09-18.
## And then, code was written to avoid this problem altogether using
## the textbook approach to find hitting probs (and just from WT).
## But the most sensible solution is probably to threshold at
## HyperHMM's output at, say, 1e-12, and set those values to 0,
## and use markovchain::hittingProbabilities.

## Why thresholding, and why it is kept as the default, is explained in
## detail in the (long) comments before function
## threshold_transition_matrix, in file trm.R. Please refer to those.


## Just in case, we compute now hitting probs from WT using markovchain::
## and our custom code, and report the largest difference. But using
## markovchain:: has the virtue of using a well tested package by people
## who really know what they are doing.

## As a comparison, none of 2200 random fitness landscapes from
## a variety of models have shown any problems with
## markovchain::hittingProbabilities (where the transition probs.)
## are computed from fitness differences (without using the thresholding).

## The rest of this file and tests provides further details. Most
## were drafted pre-bug report of markovchain::hittingProbabilities,
## most of them.

## Tests for the fix to the slow / "system is singular"-flooding HyperHMM
## intervention hitting-probability computation. Two pieces:
##
##  A) kill_gene_HyperHMM_drop_unreachable: drop the killed genotypes
##     (inert, unreachable states) from the transition matrix. It is
##     lossless: it must reproduce the downstream predicted genotype
##     frequencies and WT hitting probabilities (on the surviving
##     genotypes; the dropped ones are 0).
##
##  B) hitting_probs_from_WT_direct: compute the WT hitting probabilities
##     directly (fundamental matrix), instead of via the full all-pairs
##     markovchain::hittingProbabilities.
##

## The ground truth here is a Monte-Carlo simulation of the chain, not
## markovchain::hittingProbabilities. markovchain builds the whole
## all-pairs matrix, solving one system per target; with entries as tiny as
## those above, those systems are singular and its approximate solves
## return wrong answers (off by up to ~0.1, and even negative
## "probabilities" ---what was reported on 2026-07-04).
## The direct method is exact and matches the simulation.
## So we validate B against Monte Carlo, and additionally against
## markovchain only in the single- absorbing-state regime (no
## intervention), where markovchain is reliable.

## Even without the numerical problem above, you might wonder:
## (and I started here, before finding the "negative probabilities")
## does this really matter? Yes, for simulations all of
## this can make a huge difference in speed. In one example with 9 genes,
## killing one of the genes (F in this case) lead to 45 absorbing states:
## these 45 genotypes had their only forward move, in the original fit, go
## to a genotype that gains F; they spaned 3-8 mutations (6 at 3 muts, 14
## at 4, 19 at 5, 5 at 6, 1 at 8), e.g. A, B, I → A, B, F, I and
## B, D, E → B, D, E, F. If we remove the unreachable genotypes,
## markovchain::hittingProbabilities still would still be looking
## at the hitting probabilities starting from 45 states that are
## absorbing, but those are not the problem we are trying to solve.
## If we don't remove the unreachable genotypes, and just use the
## hitting_probs_from_WT_direct we skip the problem completely
## (and the "negative probabilities" too :-) ).
## That specific 9-gene model is included as fixture_hyperhmm_9genes.rds
## and used in the regression test at the bottom.

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

set_and_print_seed()

## Compare two named numeric vectors on their common names, with a numeric
## tolerance, ignoring attributes/order.
same_on_common <- function(a, b, tol = 1e-6) {
    common <- intersect(names(a), names(b))
    isTRUE(all.equal(a[common], b[common],
                     tolerance = tol, check.attributes = FALSE))
}

## Independent *GROUND TRUTH*: Monte-Carlo estimate of P(ever reach each state
## from WT). We simulate the JUMP chain (skip self-loops: they never change
## which states are visited), so it is fast regardless of how much mass sits on
## the diagonal. NOTE: this is exact for every OFF-diagonal target, but returns
## 0 for WT itself (it skips WT's self-loop "return"). So all comparisons below
## use max_diff_off_wt(), which excludes the WT self-entry -- a self-loop
## return-probability convention, checked separately against markovchain in the
## single-absorbing no-intervention case where markovchain is reliable.
mc_hitting_from_WT <- function(m, nsim) {
  m <- as.matrix(m)
  states <- rownames(m)
  n <- length(states)
  wt <- which(states == "WT")
  off <- m
  diag(off) <- 0
  rs <- rowSums(off)
  terminal <- rs <= 1e-12
  jump <- off
  jump[!terminal, ] <- off[!terminal, ] / rs[!terminal]
  hitc <- integer(n)
  for (s in seq_len(nsim)) {
    cur <- wt
    vis <- logical(n)
    while (!terminal[cur]) {
      cur <- sample.int(n, size = 1, prob = jump[cur, ])
      vis[cur] <- TRUE
    }
    hitc <- hitc + vis
  }
  return(setNames(hitc / nsim, states))
}

## max |a - b| over common names, EXCLUDING the WT self-entry.
max_diff_off_wt <- function(a, b) {
    common <- setdiff(intersect(names(a), names(b)), "WT")
    return(max(abs(a[common] - b[common])))
}


test_that("A) drop_unreachable is lossless (genot freqs and hitting probs)", {
  for (ng in c(3, 4, 5)) {
    for (rep in 1:4) {
      rmhn <- random_evam(model = "MHN", ngenes = ng)
      sample_mhn <- sample_evam(rmhn, N = round(runif(1, 500, 3000)),
                                obs_noise = 0.03)
      dd <- sample_mhn$MHN_sampled_genotype_counts_as_data
      h1 <- evam_like_HyperHMM(dd)
      tm <- h1$HyperHMM_trans_mat

      for (gene in LETTERS[1:ng]) {
        full <- kill_gene_HyperHMM(tm, gene)
        drop <- kill_gene_HyperHMM_drop_unreachable(tm, gene)

        surv <- rownames(drop)
        expect_true(all(surv %in% rownames(full)))
        expect_true(nrow(drop) <= nrow(full))
        expect_true("WT" %in% surv)

        of <- get_full_output(full)
        od <- get_full_output(drop)

        ## (i) predicted genotype frequencies match on surviving
        ## genotypes; the dropped genotypes carried ~0 frequency.
        gf_full <- of$HyperHMM_predicted_genotype_freqs
        gf_drop <- od$HyperHMM_predicted_genotype_freqs
        expect_true(same_on_common(gf_full, gf_drop))
        dropped <- setdiff(names(gf_full), surv)
        if (length(dropped)) {
          expect_true(max(abs(gf_full[dropped])) < 1e-9)
        }

        ## (ii) WT hitting probabilities match on surviving genotypes
        ## (get_full_output uses the direct method for both); the
        ## dropped genotypes had ~0 hitting probability.
        hp_full <- of$HyperHMM_hitting_probs_from_WT
        hp_drop <- od$HyperHMM_hitting_probs_from_WT
        expect_true(same_on_common(hp_full, hp_drop))
        dropped_hp <- setdiff(names(hp_full), surv)
        if (length(dropped_hp)) {
          expect_true(max(abs(hp_full[dropped_hp])) < 1e-8)
        }
      }
    }
  }
})


test_that("B) direct hitting probs are valid and match Monte Carlo", {
  for (ng in c(3, 4, 5)) {
    for (rep in 1:3) {
      rmhn <- random_evam(model = "MHN", ngenes = ng)
      sample_mhn <- sample_evam(rmhn, N = round(runif(1, 500, 3000)),
                                obs_noise = 0.03)
      dd <- sample_mhn$MHN_sampled_genotype_counts_as_data
      h1 <- evam_like_HyperHMM(dd)
      tm <- h1$HyperHMM_trans_mat

      ## No intervention: single absorbing state, so markovchain is
      ## reliable here -- use it as a TIGHT exactness anchor for direct.
      expect_true(same_on_common(hitting_probs_from_WT_direct(tm),
                                 hitting_probs_from_WT(tm), tol = 1e-6))

      ## Every kill: direct must be a valid probability vector and match
      ## the Monte-Carlo ground truth.
      for (gene in LETTERS[1:ng]) {
        full <- kill_gene_HyperHMM(tm, gene)
        hp <- hitting_probs_from_WT_direct(full)
        expect_true(all(hp >= -1e-9 & hp <= 1 + 1e-9))
        set.seed(1)
        mc <- mc_hitting_from_WT(full, nsim = 40000)
        expect_true(max_diff_off_wt(hp, mc) < 0.03)
      }
    }
  }
})


test_that("dropping unreachable does not change direct hitting probs", {
  for (ng in c(3, 5)) {
    for (rep in 1:3) {
      rmhn <- random_evam(model = "MHN", ngenes = ng)
      sample_mhn <- sample_evam(rmhn, N = round(runif(1, 500, 3000)),
                                obs_noise = 0.03)
      dd <- sample_mhn$MHN_sampled_genotype_counts_as_data
      h1 <- evam_like_HyperHMM(dd)
      tm <- h1$HyperHMM_trans_mat
      for (gene in LETTERS[1:ng]) {
        full <- kill_gene_HyperHMM(tm, gene)
        drop <- kill_gene_HyperHMM_drop_unreachable(tm, gene)
        hp_full <- hitting_probs_from_WT_direct(full)
        hp_drop <- hitting_probs_from_WT_direct(drop)
        ## identical on the surviving (reachable) genotypes ...
        expect_true(same_on_common(hp_full, hp_drop))
        ## ... and 0 on the dropped (unreachable) ones.
        dropped <- setdiff(names(hp_full), rownames(drop))
        if (length(dropped)) {
          expect_true(max(abs(hp_full[dropped])) < 1e-10)
        }
      }
    }
  }
})


## Regression test on the included 9-gene example whose gene F kill produces 45
## absorbing states. This is the case that showed the markovchain bug: direct
## must match the Monte-Carlo ground truth and be a valid probability vector,
## while markovchain was demonstrably wrong here: far from the simulation.
test_that("45-absorbing fixture: direct correct (formerly also legacy markovchain wrong)", {
  tm <- readRDS("fixture_hyperhmm_9genes.rds")
  expect_equal(attributes(tm)$num_features, 9)

  killF <- kill_gene_HyperHMM(tm, "F")
  m <- as.matrix(killF)
  eps <- sqrt(.Machine$double.eps)
  n_absorbing <- sum((rowSums(m) <= eps) | (diag(m) >= 1 - eps))
  ## 256 inert killed-gene genotypes + the reachable dead-ends
  expect_true(n_absorbing >= 256 + 40)

  hp_direct <- hitting_probs_from_WT_direct(killF)
  ## legacy markovchain is unusable on this ill-conditioned chain. Its exact
  ## failure mode is NOT stable across machines even with the same markovchain
  ## version: the per-target systems are solved by Armadillo (arma::solve in
  ## markovchain's compiled Rcpp code), so the compiled RcppArmadillo/BLAS/
  ## LAPACK stack decides whether solve() returns an approximate (negative,
  ## mathematically impossible) solution or throws "solve(): solution not
  ## found". We accept EITHER as proof that markovchain is broken here.
  ## Well, no; as said below:
  ## this is all very fragile; markovchaing could fix this completely
  ## or introduce other changes (this is in a state of flux) and this
  ## would fail for reasons I have no control over and that would
  ## not affect our code (as we have the divergence check)
  ## hp_legacy <- tryCatch(suppressWarnings(hitting_probs_from_WT(killF)),
  ##                       error = function(e) e)

  set.seed(7)
  mc <- mc_hitting_from_WT(killF, nsim = 100000)

  ## direct: valid probabilities, and matches the simulation
  expect_true(all(hp_direct >= -1e-9 & hp_direct <= 1 + 1e-9))
  expect_true(max_diff_off_wt(hp_direct, mc) < 0.03)

  ## markovchain: was wrong here, up until 2026-09-18 -- either it errors
  ## on the ill-conditioned solve, or it returns a solution clearly off the
  ## simulation. We do not assert HOW it is wrong, as that depends on the
  ## version: before the fix of 2026-08-12 it returned negative
  ## "probabilities"; after the fix it returns values inside [0, 1] that
  ## are still far from the simulation.

  ## As of 2026-09-18 the fix/hittingprobs branch commit db73e689
  ## works fine. It differs from the Monte Carlo because we use
  ## 1-e12 thresholding in the Monte Carlo, but not in the matrix
  ## we give to markovchain. Actually, no:
  ## this is all very fragile; markovchaing could fix this completely
  ## or introduce other changes (this is in a state of flux) and this
  ## would fail for reasons I have no control over and that would
  ## not affect our code (as we have the divergence check)
  ## Commented for now.
  ## if (inherits(hp_legacy, "error")) {
  ##   succeed("markovchain errors on the ill-conditioned chain")
  ## } else {
  ##   expect_true(max_diff_off_wt(hp_legacy, mc) > 0.05)
  ## }
})

set.seed(NULL)

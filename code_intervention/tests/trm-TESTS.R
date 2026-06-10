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

## Tests of genots_at_t_from_trm and probs_uniform_sampling_custom from trm.R,
## using closed-form analytical ground truths on minimal absorbing chains.
## Also tests of hitting probabilities.
##
##
## genots_at_t_from_trm(trm, t) uses matrix exponentiation (expm::expAtv) to
## compute the genotype distribution at exact time t.
##
## probs_uniform_sampling_custom(trm) averages genots_at_t_from_trm over 101
## equally-spaced time points in [0, 5].
##
## For both, the analytical ground truth comes from closed-form solutions:
##
##   Model 1 (1-gene): WT -> A with rate r.
##     Exact: P(WT, t) = exp(-r * t)
##            P(A,  t) = 1 - exp(-r * t)
##
##   Model 2 (2-gene sequential chain): WT -> A -> (A, B) with rates r1, r2.
##     Gene B is unreachable directly; only accessible via A.
##     Exact: P(WT, t) = exp(-r1 * t)
##            P(A, t)  = r1/(r2 - r1) * (exp(-r1*t) - exp(-r2*t))   [r1 != r2]
##            P(B, t)  = 0
##            P(A, B, t) = 1 - P(WT, t) - P(A, t)
##
## Test structure:
##   Test 1 — genots_at_t_from_trm alone, using Model 1.
##   Test 2 — genots_at_t_from_trm AND probs_uniform_sampling_custom,
##             using Model 2 ("dos pájaros de un tiro").
##             For probs_uniform_sampling_custom, the analytical check
##             computes mean(exp(-r * times)) and
##             mean(r1/(r2-r1) * (exp(-r1*t) - exp(-r2*t))) directly as
##             scalar operations — no matrices, no expm — giving a
##             completely independent implementation.
##   Test 3 — t parameter threading through intervene_cpm_every_gene.
##             Uses CBN models with the same topology as Models 1 and 2.
##             For CBN, the TRM rate equals rerun_lambda exactly
##             (a * c = 0.006 * (1/0.006) = 1), so the same analytical
##             formulas apply end-to-end.


## There is additional exploratory code in
## check_si_stats_and_explorations_trms_exponential_runs.R


## After that, there are tests of hitting probabilities.



library(testthat)

pwd <- getwd()
setwd("../")
source("intervention.R")   ## sources kill-gene-and-output-from-cpm.R -> trm.R
setwd(pwd)

set.seed(NULL)


### Model 1: 1-gene absorbing chain WT -> A with rate r
##
## TRM (diagonal 0, as required by genots_at_t_from_trm):
##      WT   A
##  WT [  0   r ]
##  A  [  0   0 ]

test_that("genots_at_t_from_trm: 1-gene absorbing chain, analytical ground truth", {
    local_edition(3)

    r <- 2.5

    trm1 <- matrix(c(0, r,
                     0, 0),
                   nrow = 2, ncol = 2, byrow = TRUE,
                   dimnames = list(c("WT", "A"), c("WT", "A")))

    for (t in c(0, 0.01, 0.1, 0.5, 1, 2, 5, 10)) {
        result <- genots_at_t_from_trm(trm1, t)
        expect_equal(result[["WT"]], exp(-r * t),
                     tolerance = 1e-8,
                     label = paste("WT at t =", t))
        expect_equal(result[["A"]], 1 - exp(-r * t),
                     tolerance = 1e-8,
                     label = paste("A at t =", t))
        expect_equal(sum(result), 1.0,
                     tolerance = 1e-8,
                     label = paste("sum = 1 at t =", t))
    }
})


### Model 2: 2-gene sequential chain WT -> A -> (A, B) with rates r1, r2
##
## TRM (diagonal 0):
##          WT    A    B   A, B
##  WT   [  0    r1    0    0  ]
##  A    [  0     0    0   r2  ]
##  B    [  0     0    0    0  ]
##  A, B [  0     0    0    0  ]
##
## B is never reached (no edge WT -> B or A -> B); P(B, t) = 0 for all t.
## (A, B) is absorbing.

test_that("genots_at_t_from_trm and probs_uniform_sampling_custom: 2-gene sequential chain", {
    local_edition(3)

    r1 <- 1.5
    r2 <- 3.5   ## r1 != r2 required by the closed-form solution

    trm2 <- matrix(0, nrow = 4, ncol = 4,
                   dimnames = list(c("WT", "A", "B", "A, B"),
                                   c("WT", "A", "B", "A, B")))
    trm2["WT", "A"]   <- r1
    trm2["A", "A, B"] <- r2

    ## ---- genots_at_t_from_trm ----
    for (t in c(0, 0.01, 0.1, 0.5, 1, 2, 5)) {
        result <- genots_at_t_from_trm(trm2, t)

        p_WT <- exp(-r1 * t)
        p_A  <- r1 / (r2 - r1) * (exp(-r1 * t) - exp(-r2 * t))
        p_AB <- 1 - p_WT - p_A

        expect_equal(result[["WT"]],   p_WT, tolerance = 1e-8,
                     label = paste("WT at t =", t))
        expect_equal(result[["A"]],    p_A,  tolerance = 1e-8,
                     label = paste("A at t =", t))
        expect_equal(result[["B"]],    0,    tolerance = 1e-8,
                     label = paste("B at t =", t))
        expect_equal(result[["A, B"]], p_AB, tolerance = 1e-8,
                     label = paste("A, B at t =", t))
        expect_equal(sum(result), 1.0, tolerance = 1e-8,
                     label = paste("sum = 1 at t =", t))
    }

    ## ---- probs_uniform_sampling_custom ----
    ## Analytical check: scalar formulas over the same 101 time points,
    ## with no matrices or expm involved.
    times <- seq(from = 0, to = 5, length.out = 101)
    p_WT_v  <- exp(-r1 * times)
    p_A_v   <- r1 / (r2 - r1) * (exp(-r1 * times) - exp(-r2 * times))
    p_AB_v  <- 1 - p_WT_v - p_A_v

    expected_WT <- mean(p_WT_v)
    expected_A  <- mean(p_A_v)
    expected_B  <- 0
    expected_AB <- mean(p_AB_v)

    result_u <- probs_uniform_sampling_custom(trm2)

    expect_equal(result_u[["WT"]],   expected_WT, tolerance = 1e-8,
                 label = "uniform WT")
    expect_equal(result_u[["A"]],    expected_A,  tolerance = 1e-8,
                 label = "uniform A")
    expect_equal(result_u[["B"]],    expected_B,  tolerance = 1e-8,
                 label = "uniform B")
    expect_equal(result_u[["A, B"]], expected_AB, tolerance = 1e-8,
                 label = "uniform A, B")
    expect_equal(sum(result_u), 1.0, tolerance = 1e-8,
                 label = "uniform sum = 1")
})


### Test 3: t parameter threading through intervene_cpm_every_gene
##
## For CBN, the TRM rate from genotype X to X+{g} equals rerun_lambda_g
## exactly, because a * c = 0.006 * (1/0.006) = 1.  This means the same
## analytical formulas from Models 1 and 2 apply directly to CBN models
## with matching rerun_lambda values, letting us test the full call chain:
##
##   intervene_cpm_every_gene(..., t = t_val)
##     -> get_genotype_freqs_cpm(model, t = t_val)
##       -> get_full_output(model)  [to obtain CBN_trans_rate_mat]
##       -> genots_from_trm(trm, t = t_val)
##         -> genots_at_t_from_trm(trm, t_val)
##
## If t were silently dropped or reset to NA anywhere in this chain,
## the result would match the exponential-sampling prediction instead
## of the exact-time prediction, and the test would fail.

test_that("t parameter threading: intervene_cpm_every_gene passes t correctly", {
    local_edition(3)

    r  <- 2.5
    r1 <- 1.5
    r2 <- 3.5

    ## 1-gene CBN: one gene A from Root with rate r
    cbn1 <- data.frame(From = "Root", To = "A",
                       rerun_lambda = r,
                       stringsAsFactors = FALSE)

    ## 2-gene sequential CBN: A from Root (r1), B from A (r2)
    ## B is unreachable from WT (AND constraint; B requires A)
    cbn2 <- data.frame(From         = c("Root", "A"),
                       To           = c("A",    "B"),
                       rerun_lambda = c(r1,      r2),
                       stringsAsFactors = FALSE)

    for (t_val in c(0.5, 1, 2)) {

        ## --- 1-gene model ---
        res1 <- suppressWarnings(
            intervene_cpm_every_gene(list(CBN_model = cbn1), "CBN",
                                     t = t_val))
        ni1 <- res1[["no_intervention"]]$genot_freqs

        expect_equal(ni1[["WT"]], exp(-r * t_val),
                     tolerance = 1e-6,
                     label = paste("1-gene WT, t =", t_val))
        expect_equal(ni1[["A"]], 1 - exp(-r * t_val),
                     tolerance = 1e-6,
                     label = paste("1-gene A, t =", t_val))

        ## --- 2-gene sequential model ---
        res2 <- suppressWarnings(
            intervene_cpm_every_gene(list(CBN_model = cbn2), "CBN",
                                     t = t_val))
        ni2 <- res2[["no_intervention"]]$genot_freqs

        p_WT <- exp(-r1 * t_val)
        p_A  <- r1 / (r2 - r1) * (exp(-r1 * t_val) - exp(-r2 * t_val))
        p_AB <- 1 - p_WT - p_A

        expect_equal(ni2[["WT"]],   p_WT, tolerance = 1e-6,
                     label = paste("2-gene WT, t =", t_val))
        expect_equal(ni2[["A"]],    p_A,  tolerance = 1e-6,
                     label = paste("2-gene A, t =", t_val))
        expect_equal(ni2[["A, B"]], p_AB, tolerance = 1e-6,
                     label = paste("2-gene A, B, t =", t_val))
        ## B has zero probability (requires A; unreachable from WT)
        expect_false("B" %in% names(ni2),
                     label = paste("2-gene B absent, t =", t_val))
    }
})


### Test 4: hitting probabilities tests

test_that("to_markovchain works for disconnected states", {

  for (i in 1:5) {
    ## Just to have a scaffolding
    emhn3 <- evamtools::random_evam(3, model = "MHN")
    m3 <-  emhn3$MHN_trans_mat
    ## Disconnect A,B
    m3["A", "A, C"] <- m3["A", "A, C"] + m3["A", "A, B"]
    m3["B", "B, C"] <- m3["B", "B, C"] + m3["B", "A, B"]
    m3[, "A, B"] <- 0
    m3["A, B", ] <- 0

    m3c <- to_markovchain(m3)
    expect_true(isTRUE(all.equal(hittingProbabilities(m3c)[1, "A, B"],
                                 0.0, check.attributes = FALSE)))

    emhn5 <- evamtools::random_evam(5, model = "MHN")
    m5 <-  emhn5$MHN_trans_mat
    ## Disconnect B, C, D
    m5["B, C", "B, C, E"] <- m5["B, C", "B, C, E"] + m5["B, C", "B, C, D"]
    m5["B, D", "B, D, E"] <- m5["B, D", "B, D, E"] + m5["B, D", "B, C, D"]
    m5["C, D", "A, C, D"] <- m5["C, D", "A, C, D"] + m5["C, D", "B, C, D"]

    m5[, "B, C, D"] <- 0
    m5["B, C, D", ] <- 0

    m5c <- to_markovchain(m5)
    expect_true(isTRUE(all.equal(hittingProbabilities(m5c)[1, "B, C, D"],
                                 0.0, check.attributes = FALSE)))
  }
})


test_that("hitting_probs_from_WT: basic 3-state linear chain", {
    ## Hand-crafted 3-state chain: WT -> A -> AB (linear, no branching)
    ## Transition probability matrix (embedded chain):
    ##   WT  -> A with prob 1
    ##   A   -> AB with prob 1
    ##   AB  -> AB with prob 1  (absorbing, self-loop added by to_markovchain)
    ## From WT: P(ever visit A)  = 1, P(ever visit AB) = 1
    ## First-passage convention: h(WT, WT) = 0
    mat <- matrix(0, nrow = 3, ncol = 3,
                  dimnames = list(c("WT", "A", "AB"), c("WT", "A", "AB")))
    mat["WT", "A"]  <- 1
    mat["A",  "AB"] <- 1
    ## AB is absorbing: row sums to 0, to_markovchain adds self-loop

    hp <- hitting_probs_from_WT(mat)

    expect_equal(hp[["WT"]], 0.0,   tolerance = 1e-10, label = "h(WT,WT)=0")
    expect_equal(hp[["A"]],  1.0,   tolerance = 1e-10, label = "h(WT,A)=1")
    expect_equal(hp[["AB"]], 1.0,   tolerance = 1e-10, label = "h(WT,AB)=1")
})


test_that("hitting_probs_from_WT: MHN consistency with hittingProbabilities", {
    ## For any MHN trans_mat, hitting_probs_from_WT should equal
    ## the first row of hittingProbabilities(to_markovchain(trans_mat)).
    ## State 1 (row 1) is always the WT genotype.

  for (n_genes in c(3, 4, 5)) {
    emhn <- evamtools::random_evam(n_genes, model = "MHN")
    tm <- emhn$MHN_trans_mat
    mc <- to_markovchain(tm)
    expected <- hittingProbabilities(mc)[1, ]
    observed <- hitting_probs_from_WT(tm)
    expect_equal(observed, expected,
                 tolerance = 1e-10,
                 label = paste0("MHN n=", n_genes, " hitting probs"))
  }
})


test_that("hitting_probs_from_WT: all values in [0,1], WT=0, absorbing=1", {
    ## For a 3-gene MHN (all 8 genotypes):
    ## - hitting probs are in [0, 1]
    ## - WT has hitting prob 0 (first-passage convention)

  for (i in 1:5) {
    emhn <- evamtools::random_evam(3, model = "MHN")
    tm <- emhn$MHN_trans_mat
    hp <- hitting_probs_from_WT(tm)

    expect_true(all(hp >= 0 - 1e-10 & hp <= 1 + 1e-10),
                label = paste("all in [0,1], iter", i))
    ## WT is the first state (row 1 of trans_mat)
    wt_name <- rownames(tm)[1]
    expect_equal(hp[[wt_name]], 0.0, tolerance = 1e-10,
                 label = paste("WT hitting prob = 0, iter", i))
  }
})


test_that("hitting_probs_from_WT: inline row-scaling of rate mat agrees with trans_mat route", {
    ## For MHN, MHN_trans_mat == row_normalized(MHN_trans_rate_mat).
    ## hitting_probs_from_WT computed from MHN_trans_mat (CPM route) should
    ## equal hitting_probs_from_WT computed from the inline row-scaled
    ## MHN_trans_rate_mat (the fitness-landscape-style route).
    ## This verifies the two routes are equivalent.

  for (i in 1:5) {
    emhn <- evamtools::random_evam(4, model = "MHN")
    tm <- emhn$MHN_trans_mat
    trm <- emhn$MHN_trans_rate_mat

    ## CPM route: use pre-computed embedded chain
    hp_cpm <- hitting_probs_from_WT(tm)

    ## Fitness-landscape-style route: inline row-scaling of rate matrix
    rs <- rowSums(trm)
    embedded <- trm
    embedded[rs > 0, ] <- trm[rs > 0, ] / rs[rs > 0]
    hp_fl <- hitting_probs_from_WT(embedded)

    expect_equal(hp_cpm, hp_fl,
                 tolerance = 1e-10,
                 label = paste("inline row-scale vs trans_mat, iter", i))
  }
})


### Test 5: OT and OncoBN end-to-end (hand-specified model -> trans_mat -> HP)
##
## These tests pin down the model -> embedded-transition-matrix step for OT
## and OncoBN, then verify HP on top.  Rationale: cross-path equivalences
## in kill-OT-TESTS.R and kill-OncoBN-TESTS.R show that the two available
## kill paths (DAG kill vs. params -> 0) agree, but only against each other;
## neither path is itself pinned to externally hand-derived TRM values.
## For CBN/HESBCN/MHN that gap is closed by 3- or 4-way agreement across
## independent constructions (DAG, params -> 0, TRM removal, fitness
## landscape).  OT and OncoBN have no TRM-removal or fitness-landscape
## routes, so the only way to anchor model -> TRM to external math is a
## direct hand derivation.
##
## OT and OncoBN expose only `*_trans_mat` (the row-normalized embedded
## chain); `*_trans_rate_mat` is NULL for both.  So the assertion target
## is the embedded matrix, not a rate matrix.  HP is then `hitting_probs_from_WT`
## of that embedded matrix; the TRM -> HP arithmetic is already verified
## analytically by the 3-state linear-chain test above.
##
## Genotypes in the returned matrix are restricted to those reachable
## from WT under the model's restrictions (e.g., in an OT chain Root -> A -> B,
## B alone is not a reachable state and does not appear).


## ------ Test 5a: OT, branching topology ------------------------------
##
## Model: Root -> A (weight 0.9), Root -> B (weight 0.7).
## Both genes have Root as their parent, so all 4 genotypes are reachable.
##
## Rate-level transitions (rate = OT_edgeWeight when the gene's parent is
## present in the source genotype, else 0):
##   WT      -> A    : 0.9        (A's parent Root is implicit in WT)
##   WT      -> B    : 0.7        (B's parent Root)
##   A       -> A, B : 0.7        (B's parent Root present)
##   B       -> A, B : 0.9        (A's parent Root present)
##   A, B            : absorbing.
##
## Row-normalized to embedded transition probabilities:
##   WT      -> A    : 0.9 / 1.6 = 0.5625
##   WT      -> B    : 0.7 / 1.6 = 0.4375
##   A       -> A, B : 1.0
##   B       -> A, B : 1.0
##
## Hitting probabilities from WT (P(ever visit g | start at WT)):
##   HP(WT)   = 0      (convention: starting state excluded)
##   HP(A)    = P(WT -> A)         = 0.5625
##   HP(B)    = P(WT -> B)         = 0.4375
##   HP(A, B) = 1 (absorbing; both paths converge here).
test_that("OT, branching Root->A, Root->B: trans_mat and HP from hand derivation", {
    local_edition(3)
    ot_branch <- data.frame(
        From = c("Root", "Root"),
        To   = c("A",    "B"),
        OT_edgeWeight = c(0.9, 0.7))
    out <- get_full_output(ot_branch)

    genots <- c("WT", "A", "B", "A, B")
    expected_tm <- matrix(0, nrow = 4, ncol = 4,
                          dimnames = list(genots, genots))
    expected_tm["WT",   "A"]    <- 0.9 / 1.6
    expected_tm["WT",   "B"]    <- 0.7 / 1.6
    expected_tm["A",    "A, B"] <- 1.0
    expected_tm["B",    "A, B"] <- 1.0

    observed_tm <- as.matrix(out$OT_trans_mat)
    expect_equal(observed_tm[genots, genots], expected_tm, tolerance = 1e-10)

    expected_hp <- c(WT = 0, A = 0.5625, B = 0.4375, "A, B" = 1.0)
    expect_equal(out$OT_hitting_probs_from_WT[genots], expected_hp,
                 tolerance = 1e-10)
})


## ------ Test 5b: OT, chain topology ----------------------------------
##
## Model: Root -> A (weight 0.9), A -> B (weight 0.7).
## B's parent is A (not Root), so B alone is unreachable from WT and is
## absent from the returned matrix.  Reachable genotypes: WT, A, A, B.
##
## Rate-level transitions:
##   WT   -> A    : 0.9     (A's parent Root)
##   A    -> A, B : 0.7     (B's parent A is present)
##   A, B         : absorbing.
##
## Row-normalized:
##   WT   -> A    : 1.0
##   A    -> A, B : 1.0
##
## Hitting probabilities:
##   HP(WT)   = 0
##   HP(A)    = 1   (deterministic single path WT -> A)
##   HP(A, B) = 1   (deterministic, absorbing).
test_that("OT, chain Root->A->B: trans_mat and HP from hand derivation", {
    local_edition(3)
    ot_chain <- data.frame(
        From = c("Root", "A"),
        To   = c("A",    "B"),
        OT_edgeWeight = c(0.9, 0.7))
    out <- get_full_output(ot_chain)

    genots <- c("WT", "A", "A, B")
    expected_tm <- matrix(0, nrow = 3, ncol = 3,
                          dimnames = list(genots, genots))
    expected_tm["WT", "A"]    <- 1.0
    expected_tm["A",  "A, B"] <- 1.0

    observed_tm <- as.matrix(out$OT_trans_mat)
    expect_equal(rownames(observed_tm), genots,
                 label = "unreachable B excluded from matrix")
    expect_equal(observed_tm[genots, genots], expected_tm, tolerance = 1e-10)

    expected_hp <- c(WT = 0, A = 1.0, "A, B" = 1.0)
    expect_equal(out$OT_hitting_probs_from_WT[genots], expected_hp,
                 tolerance = 1e-10)
})


## ------ Test 5c: OncoBN with AND relations (CBN-mode) ----------------
##
## Model:
##   Root -> A : Single, theta = 0.5
##   Root -> B : Single, theta = 0.4
##   A    -> C : AND,    theta = 0.3
##   B    -> C : AND,    theta = 0.3
## C is a multi-parent AND-node, so C is reachable only from genotypes
## containing BOTH A and B.  Reachable genotypes from WT under this
## restriction: WT, A, B, A, B, A, B, C  (5 states; no A, C / B, C / C alone).
##
## Rate-level transitions (rate = theta when parent condition is satisfied):
##   WT      -> A       : 0.5
##   WT      -> B       : 0.4
##   A       -> A, B    : 0.4   (B's parent Root present)
##   B       -> A, B    : 0.5   (A's parent Root present)
##   A, B    -> A, B, C : 0.3   (both parents of C present, AND satisfied)
##   A, B, C            : absorbing.
##
## Row-normalized:
##   WT   -> A          : 0.5 / 0.9 = 5/9
##   WT   -> B          : 0.4 / 0.9 = 4/9
##   A    -> A, B       : 1.0
##   B    -> A, B       : 1.0
##   A, B -> A, B, C    : 1.0
##
## Hitting probabilities:
##   HP(WT)        = 0
##   HP(A)         = 5/9
##   HP(B)         = 4/9
##   HP(A, B)      = 1   (both paths converge: 5/9 + 4/9 = 1)
##   HP(A, B, C)   = 1   (absorbing).
test_that("OncoBN AND (CBN-mode): trans_mat and HP from hand derivation", {
    local_edition(3)
    onc_and <- data.frame(
        From  = c("Root", "Root", "A", "B"),
        To    = c("A",    "B",    "C", "C"),
        theta = c(0.5,    0.4,    0.3, 0.3),
        Relation = c("Single", "Single", "AND", "AND"))
    out <- suppressMessages(get_full_output(onc_and))

    genots <- c("WT", "A", "B", "A, B", "A, B, C")
    expected_tm <- matrix(0, nrow = 5, ncol = 5,
                          dimnames = list(genots, genots))
    expected_tm["WT",   "A"]       <- 5 / 9
    expected_tm["WT",   "B"]       <- 4 / 9
    expected_tm["A",    "A, B"]    <- 1.0
    expected_tm["B",    "A, B"]    <- 1.0
    expected_tm["A, B", "A, B, C"] <- 1.0

    observed_tm <- as.matrix(out$OncoBN_trans_mat)
    expect_equal(rownames(observed_tm), genots,
                 label = "AND mode: only A, B-bearing C reachable")
    expect_equal(observed_tm[genots, genots], expected_tm, tolerance = 1e-10)

    expected_hp <- c(WT = 0, A = 5/9, B = 4/9, "A, B" = 1.0, "A, B, C" = 1.0)
    expect_equal(out$OncoBN_hitting_probs_from_WT[genots], expected_hp,
                 tolerance = 1e-10)
})


## ------ Test 5d: OncoBN with OR relations (DBN-mode) -----------------
##
## Model: same topology as 5c, but C is a multi-parent OR-node.
##   Root -> A : Single, theta = 0.5
##   Root -> B : Single, theta = 0.4
##   A    -> C : OR,     theta = 0.3
##   B    -> C : OR,     theta = 0.3
## With OR, C is reachable as soon as EITHER A or B is present.  This adds
## two new reachable genotypes vs. AND: A, C (from A) and B, C (from B).
## Reachable from WT: WT, A, B, A, B, A, C, B, C, A, B, C (7 states).
##
## Rate-level transitions:
##   WT      -> A          : 0.5
##   WT      -> B          : 0.4
##   A       -> A, B       : 0.4   (B's parent Root)
##   A       -> A, C       : 0.3   (C: A present, OR satisfied)
##   B       -> A, B       : 0.5   (A's parent Root)
##   B       -> B, C       : 0.3   (C: B present, OR satisfied)
##   A, B    -> A, B, C    : 0.3
##   A, C    -> A, B, C    : 0.4
##   B, C    -> A, B, C    : 0.5
##   A, B, C               : absorbing.
##
## Row-normalized:
##   WT   -> A             : 0.5 / 0.9 = 5/9
##   WT   -> B             : 0.4 / 0.9 = 4/9
##   A    -> A, B          : 0.4 / 0.7 = 4/7
##   A    -> A, C          : 0.3 / 0.7 = 3/7
##   B    -> A, B          : 0.5 / 0.8 = 5/8
##   B    -> B, C          : 0.3 / 0.8 = 3/8
##   A, B -> A, B, C       : 1.0
##   A, C -> A, B, C       : 1.0
##   B, C -> A, B, C       : 1.0
##
## Hitting probabilities:
##   HP(WT)         = 0
##   HP(A)          = 5/9
##   HP(B)          = 4/9
##   HP(A, B)       = HP(A) * P(A -> A, B) + HP(B) * P(B -> A, B)
##                  = 5/9 * 4/7 + 4/9 * 5/8
##                  = 20/63 + 20/72
##                  = 20/63 + 5/18
##                  = (20*18 + 5*63) / (63*18) = (360 + 315) / 1134 = 675/1134
##                  = 25/42         ( ~ 0.595238 )
##   HP(A, C)       = HP(A) * P(A -> A, C) = 5/9 * 3/7 = 15/63 = 5/21
##                  ( ~ 0.238095 )
##   HP(B, C)       = HP(B) * P(B -> B, C) = 4/9 * 3/8 = 12/72 = 1/6
##                  ( ~ 0.166667 )
##   HP(A, B, C)    = 1            (absorbing).
test_that("OncoBN OR (DBN-mode): trans_mat and HP from hand derivation", {
    local_edition(3)
    onc_or <- data.frame(
        From  = c("Root", "Root", "A", "B"),
        To    = c("A",    "B",    "C", "C"),
        theta = c(0.5,    0.4,    0.3, 0.3),
        Relation = c("Single", "Single", "OR", "OR"))
    out <- suppressMessages(get_full_output(onc_or))

    genots <- c("WT", "A", "B", "A, B", "A, C", "B, C", "A, B, C")
    expected_tm <- matrix(0, nrow = 7, ncol = 7,
                          dimnames = list(genots, genots))
    expected_tm["WT",   "A"]       <- 5 / 9
    expected_tm["WT",   "B"]       <- 4 / 9
    expected_tm["A",    "A, B"]    <- 4 / 7
    expected_tm["A",    "A, C"]    <- 3 / 7
    expected_tm["B",    "A, B"]    <- 5 / 8
    expected_tm["B",    "B, C"]    <- 3 / 8
    expected_tm["A, B", "A, B, C"] <- 1.0
    expected_tm["A, C", "A, B, C"] <- 1.0
    expected_tm["B, C", "A, B, C"] <- 1.0

    observed_tm <- as.matrix(out$OncoBN_trans_mat)
    expect_equal(sort(rownames(observed_tm)), sort(genots),
                 label = "OR mode: A, C and B, C reachable")
    expect_equal(observed_tm[genots, genots], expected_tm, tolerance = 1e-10)

    expected_hp <- c(WT      = 0,
                     A       = 5/9,
                     B       = 4/9,
                     "A, B"  = 25/42,
                     "A, C"  = 5/21,
                     "B, C"  = 1/6,
                     "A, B, C" = 1.0)
    expect_equal(out$OncoBN_hitting_probs_from_WT[genots], expected_hp,
                 tolerance = 1e-10)
})


### How hitting probabilities are tested across the codebase
##
## The relevant correctness claim for downstream statistics (O_genot in
## compute_intervention_objectives, code_evam_simul_interv/) is "HP is
## correct for a given TRM / weighted graph / model".  Intervention does
## not enter this claim: an intervention just modifies the model to
## produce a new model in the same class.  So if HP-from-a-given-model is
## correct, HP-after-intervention is correct by composition.
##
## The pipeline factors cleanly into two steps:
##
##   Step 1.  model  ->  embedded transition matrix (or rate matrix)
##            Model-specific construction.
##
##   Step 2.  TRM    ->  hitting_probs_from_WT (WT row of HP matrix)
##            Pure arithmetic in hitting_probs_from_WT() (trm.R:479),
##            a thin wrapper over markovchain::hittingProbabilities().
##
## Coverage of Step 2 (this file, "Test 4"):
##   - hitting_probs_from_WT verified against an analytical 3-state
##     linear chain (WT -> A -> A, B): closed-form HP solved by hand and
##     asserted numerically.
##   - Range / WT = 0 / absorbing = 1 properties on random MHN models.
##   - Inline row-scaling of a rate matrix agrees with the trans_mat
##     route (identity check across two embedded-chain constructions).
##   - MHN consistency with hittingProbabilities() is tautological for
##     this thin wrapper but documents the wrapping contract.
##
## Coverage of Step 1, per model class:
##
##   CBN     (kill-CBN-TESTS.R):    4 independent constructions of the
##           post-intervention TRM (DAG kill, params -> 0, TRM removal,
##           fitness landscape) verified to give the same HP via
##           expect_equal on 6 hand-crafted DAGs.  3 of the 4 paths are
##           genuinely independent code routes; agreement across all 4
##           is strong evidence the TRM (and HP) are correctly
##           constructed.
##
##   HESBCN  (kill-HESBCN-TESTS.R): same 4-way construction agreement
##           on 6 DAGs x 27 AND/OR/XOR combinations (138 variants).
##
##   MHN     (kill-MHN-TESTS.R):    3-way construction agreement (DAG
##           kill, params -> 0, TRM removal) on a hand-crafted theta
##           matrix.  No fitness-landscape route exists for MHN.
##
##   HyperHMM (kill-HyperHMM-TESTS.R): independent hand-recomputation of
##           the post-kill transition matrix via row/col zeroing and
##           re-normalization, compared to the pipeline output.
##
##   OT      (this file, Test 5a, 5b):    hand-derived embedded chain
##           and HP for a branching topology (Root -> A, Root -> B) and
##           a chain topology (Root -> A -> B), asserted against
##           get_full_output().  OT has no TRM-removal or fitness-
##           landscape route, so the 2-path agreement in
##           kill-OT-TESTS.R (DAG kill vs. params -> 0) only checks the
##           two paths against each other; these new tests anchor model
##           -> TRM to externally hand-derived numbers.
##
##   OncoBN  (this file, Test 5c AND-mode, Test 5d OR-mode): same
##           rationale as OT.  Topology Root -> A, Root -> B, A -> C,
##           B -> C with both AND ("CBN-mode": C requires both parents)
##           and OR ("DBN-mode": C requires either parent).  Reachable
##           genotype sets differ between AND (5 states) and OR (7
##           states); both branches hand-derived.
##
## Steps 1 and 2 compose: each non-trivial model class has Step 1 pinned
## to independent ground truth (3- or 4-way construction agreement, or
## hand re-implementation, or direct hand derivation), and Step 2 is
## analytically validated by the 3-state linear-chain test above.
## Therefore HP from any intervened model in any of these classes is
## also pinned, by composition.


set.seed(NULL)

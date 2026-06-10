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

## Ground-truth tests of MHN interventions using a hand-crafted 4-gene model.
##
## MHN has no cascade: killing gene X means removing X's row and column from
## the theta matrix. Every other gene's behaviour is unaffected except through
## the direct interactions with X. There is no AND/OR/XOR logic and no
## propagation to downstream nodes.
##
## For each gene kill, three procedures are compared against each other and
## against an explicit ground truth:
##
##   Procedure 1: standard kill_gene (removes row/col from theta matrix)
##   Procedure 2: kill_gene_by_params_to_0 (zeroes parameters)
##   Procedure 3: intervene_cpm_trm_rm_every_gene (removes genotypes from TRM)
##
## Ground truth for killing gene g: run get_full_output on the sub-matrix
## theta[-g, -g], constructed by direct row/column indexing — entirely
## independent of kill_gene.
##
## Procedure 4 (fitness landscape) is not applicable to MHN.

## options(intervention_every_gene_cores = parallel::detectCores())
library(testthat)
pwd <- getwd()
setwd("../")
source("intervention.R")
setwd(pwd)

set.seed(NULL)


### Hand-crafted 4-gene MHN theta matrix
##
## Genes: A, B, C, D.
## Diagonal entries: baseline log-hazard rates (negative = slow accumulation).
## Off-diagonal theta[i,j]: effect of gene j on gene i's rate
##   (positive = promoting, negative = inhibiting).
##
## The values are chosen to give a non-degenerate model with a mix of
## promoting and inhibiting interactions, so all gene kills produce
## meaningfully different outputs.

theta <- matrix(
    c(-2.0,  0.5, -0.3,  0.0,
       0.3, -1.5,  0.4, -0.2,
      -0.1,  0.2, -1.0,  0.3,
       0.0, -0.3,  0.1, -3.0),
    nrow = 4, ncol = 4, byrow = TRUE,
    dimnames = list(c("A", "B", "C", "D"), c("A", "B", "C", "D"))
)

## Wrap in the list structure expected by intervene_cpm_every_gene
mhn_cpm <- list(MHN_theta = theta)



### Utility

## Predicted genotype frequencies and hitting probs from a theta (sub-)matrix.
preds_from_mhn_theta <- function(th) {
    tmp <- suppressWarnings(get_full_output(th))
    f <- tmp$MHN_predicted_genotype_freqs
    hp <- tmp$MHN_hitting_probs_from_WT
    list(genot_freqs = f[f > 0], hitting_probs_from_WT = filter_hp_keep_wt(hp))
}

## Extract a single named intervention from the full intervention output as a
## list(genot_freqs, hitting_probs_from_WT), for direct comparison with
## preds_from_mhn_theta.
get_interv <- function(res, name) {
    return(res[[name]])
}


### Tests

test_that("MHN 4-gene: procedures 1, 2, 3 agree on all gene kills", {
    local_edition(3)
    cat("\n MHN 4-gene: procedures 1, 2, 3 vs each other\n")

    ## Procedure 1
    i1 <- intervene_cpm_every_gene(mhn_cpm, "MHN")

    ## Procedure 2
    i2 <- suppressWarnings(
        intervene_cpm_every_gene(mhn_cpm, "MHN",
                                 kill_gene_funct = kill_gene_by_params_to_0))

    ## Procedure 3: TRM removal
    full_out <- suppressWarnings(get_full_output(theta))
    full_out$MHN_theta <- theta
    i3 <- intervene_cpm_trm_rm_every_gene(full_out, "MHN")

    expect_equal(i1, i2, label = "MHN i1 vs i2")
    expect_equal(i1, i3, label = "MHN i1 vs i3")
})


test_that("MHN 4-gene: ground-truth comparison for all gene kills", {
    local_edition(3)
    cat("\n MHN 4-gene: ground-truth comparison\n")

    ## Ground truth: for each gene g, run the 3x3 sub-matrix theta[-g, -g].
    ## Constructed by direct name-based indexing, independent of kill_gene.
    genes <- colnames(theta)

    i1 <- intervene_cpm_every_gene(mhn_cpm, "MHN")

    for (g in genes) {
        others <- setdiff(genes, g)
        expected <- preds_from_mhn_theta(theta[others, others])
        expect_equal(get_interv(i1, paste0("I:", g)), expected,
                     label = paste0("MHN ground truth kill ", g))
    }
})

### Another example of a 4x4
##  The above is better, though, as how we remove the genes
##  is completely different from what we do in the code

test_that("MHN, another 4x4", {
    t1 <- matrix(runif(36, -4, 4), ncol = 6)
    colnames(t1) <- rownames(t1) <- LETTERS[1:6]
    expect_equal(kill_gene(t1, "A"), t1[-1, -1])
    expect_equal(kill_gene(t1, "D"), t1[-4, -4])
})

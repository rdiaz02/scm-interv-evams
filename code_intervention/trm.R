## Copyright 2022 Ramon Diaz-Uriarte, Íñigo Ríos Arroyo

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

## Functions related to transition rate matrices:
##   - scaling
##   - verifying the matrices lead to samples that fulfill requirements
##   - obtaining samples
##   - obtaining a scaled trm from a fitness landscape

#### Caveat: gene names
# This code only works with genes named from A-Z
# If we wanted to use custom named genes (from a dataset for example)
# we need to rename the genes to letters and then execute the code



### Loading/sourcing dependencies

library(evamtools)
library(Matrix)
suppressMessages(library(OncoSimulR, quietly = TRUE, verbose = FALSE))
library(stringi)
library(expm)
library(markovchain) ## hitting probabilities

source("rfitness_to_trm.R")
source("utils.R")



### Code


## fitness landscape, as a matrix (as returned from rfitness or
##       generate_f_landscape) -> scaled transition rate matrix
## For scaling, r ~ Unif(1, 6) is now set in stone.
## If you pass a c, that is the one used.
## Why the "adaptive?" Because if you do not pass a "c" one
## is generated "adaptively". This naming is a relict
## from the initial explorations.
get_scaled_trm_adaptive <- function(fit_land, c = NA) {
  if (all(is.na(fit_land))) {
    message("Fitness landscape is NA.")
    return(list(trm_scaled = NA,
                c = NA,
                rand_target_mean_rate = NA,
                actual_mean_rate = NA,
                rel_diff = NA,
                si_stats = c(si_mean = NA,
                             si_median = NA,
                             si_max  = NA,
                             si_min = NA)
                ))
  }

  fit_land_letters <- no_evam_rfitness_to_letter(fit_land)
  tmp <- no_evam_genots_2_fgraph_and_trans_mat_rf(fit_land_letters)
  accessible_genotypes <- tmp$accessible_genotypes
  rel_diff <- tmp$relative_fitness_differences

  if (length(accessible_genotypes) == 0) {
    ## Paranoid checks
    if (length(rel_diff) > 1) stop("This should not happen")
    if (!all(is.na(rel_diff))) stop("This should not happen")

    message("No accesible genotypes: returning NA from get_scaled_trm")
    return(list(trm_scaled = NA,
                c = NA,
                rand_target_mean_rate = NA,
                actual_mean_rate = NA,
                rel_diff = NA,
                si_stats = c(si_mean = NA,
                             si_median = NA,
                             si_max  = NA,
                             si_min = NA)
                ))
  }

  if ((is.na(as.matrix(rel_diff)[1, 1])) ||
        any(is.na(rel_diff))) stop("This should not happen")

  rdtmp <- rel_diff[which(rel_diff > 0)]
  mean_si <- mean(rdtmp)

  ## rel_diff, from no_evam_genots_2_fgraph_and_trans_mat_rf,
  ## removes rows that are only destinations.
  ## Must be squared for genots_at_t_from_trm and probs_from_trm to work
  rel_diff_sq <- get_square_matrix(rel_diff, 0)

  ## This uses the relationship in Gillespie, 1984
  ## "Molecular evolution over the mutational landscape", Evolution.
  ## transition rate = s_i * scaling term

  ## rand_target_mean_rate <- rgamma(1, shape = 0.8, rate = 0.1) ## runif(1, .1, 100)
  ## log-uniform in 0.1, 100
  ## This interval is way too large. But we want to know how things behave
  ## We will decrease it later.
  ## rand_target_mean_rate <- 10^(runif(1, -1, 2))
  ## Yes, 1 to 6 are sensible limits when sampling time exponential rate 1
  rand_target_mean_rate <- runif(1, 1, 6)
  if (is.na(c))
    c <- rand_target_mean_rate/mean_si
  else
    message("c not computed, as given as ", c)

  message("    actual mean rate = ", c * mean_si,
          "  target = ", rand_target_mean_rate)

  trm_scaled <- c * rel_diff_sq

  ## Some of this output is silly for CBN/H-ESBCN
  ##      c is fixed (1/0.006 as of 2025)
  ##      rand_target_mean_rate is never used
  ##      trm_scaled is just c * original trm from the evam model

  ## And for NK and RMF, rand_target_mean_rate = actual_mean_rate
  ## by construction
  return(list(trm_scaled = trm_scaled,
              c = c,
              rand_target_mean_rate = rand_target_mean_rate,
              actual_mean_rate = c * mean_si,
              rel_diff = rel_diff,
              si_stats = c(si_mean = mean_si,
                           si_median = median(rdtmp),
                           si_max  = max(rdtmp),
                           si_min = min(rdtmp))
              ))
}

## Returns TRUE if trm meets the requirements of gene and genotype frequencies
##   This function will use BLAS and multiple threads.
check_frequencies_of_genotypes_in_trm <- function(trm,
                                                  genes,
                                                  min_gene_freq = 0.1,
                                                  max_gene_freq = 0.9,
                                                  min_WT_freq = 0.02,
                                                  max_WT_freq = 0.2,
                                                  min_genots_at_least_01 = 5,
                                                  min_average_num_muts = 2,
                                                  custom_sampling = FALSE,
                                                  custom_sampling_function = probs_uniform_sampling_custom) {
    if (all(is.na(trm))) {
        return(list(test_OK = NA,
                    min_gene_freq = NA,
                    max_gene_freq = NA,
                    freq_WT = NA,
                    genots_at_least_01 = NA,
                    average_num_muts = NA,
                    average_prop_diff_trm_1_2 = NA
                    ))
    }

    if (!custom_sampling) {
        genots_freq <- evamtools:::probs_from_trm(trm)
    } else {
        genots_freq <- custom_sampling_function(trm)
    }

    genots_freq_no_WT <- genots_freq[-which(names(genots_freq) == "WT")]
    genes_freq <- vapply(genes, function(gene){
        sum(genots_freq_no_WT[grep(gene, names(genots_freq_no_WT), fixed = TRUE)])
    }, FUN.VALUE = numeric(1))
    less_than_min_freq <- any(genes_freq < min_gene_freq)
    more_than_max_freq <- any(genes_freq > max_gene_freq)

    freq_WT <- genots_freq["WT"]
    WT_less_than_min <- (freq_WT < min_WT_freq)
    WT_more_than_max <- (freq_WT > max_WT_freq)

    no_num_genots_at_least_01 <-
        (sum(genots_freq > 0.01) < min_genots_at_least_01)

    num_muts <- stringi::stri_count_fixed(names(genots_freq_no_WT), ",") + 1
    average_num_muts <- as.vector(num_muts %*% genots_freq_no_WT)
    no_min_average_num_muts <- (average_num_muts < min_average_num_muts)


    ## This compares the change in frequency of genotypes if all the sample
    ## taken   at exactly time 1 vs. all at exactly time 2. If we are
    ## already in the   limiting distribution, there will be little
    ## change. This is just   descriptive for now, we  aren't testing
    ## on it.  abs on the freqs. because for tiny values can return negative.

    p1 <- try(genots_at_t_from_trm(trm, 1))
    p2 <- try(genots_at_t_from_trm(trm, 2))
    if (inherits(p1, "try-error") || inherits(p2, "try-error"))
    average_prop_diff_trm_1_2 <- NA
  else
    average_prop_diff_trm_1_2 <- sum(abs(abs(p1) - abs(p2)))/2


  test_OK <- (!less_than_min_freq) &&
    (!more_than_max_freq) &&
    (!WT_less_than_min) &&
    (!WT_more_than_max) &&
    (!no_num_genots_at_least_01) &&
    (!no_min_average_num_muts)


  if (less_than_min_freq)
    message("CFGT_Check failed: Some genes less than min freq")
  if (more_than_max_freq)
    message("CFGT_Check failed: Some genes more than max freq")
  if (WT_less_than_min) message("CFGT_Check failed: WT less than min")
  if (WT_more_than_max) message("CFGT_Check failed: WT more than max")
  if (no_num_genots_at_least_01)
    message("CFGT_Check failed: Num. genots at least 0.01 not reached")
  if (no_min_average_num_muts)
    message("CFGT_Check failed: Average number of mutations too small")


  return(list(test_OK = test_OK,
              min_gene_freq = min(genes_freq),
              max_gene_freq = max(genes_freq),
              freq_WT = unname(freq_WT),
              genots_at_least_01 = sum(genots_freq > 0.01),
              average_num_muts = average_num_muts,
              average_prop_diff_trm_1_2 = average_prop_diff_trm_1_2,
              ## This was not returned originally,
              ## before the uniform runs
              ## but we could
              ## have, since we've computed the probs of genots.
              ## In fact, silly that we did not.
              genots_freq_from_check = genots_freq
              ))
}


## transition rate matrix, time -> genotype composition at time t
## Use the exact solution
genots_at_t_from_trm <- function(trm, t) {
  ## Recall our trms have zero in diagonal
  stopifnot(isTRUE(all(diag(trm) == 0)))
  dim_trm <- dim(trm)
  if (dim_trm[1] != dim_trm[2]) stop("transition rate matrix not square")
  diag(trm) <- -1 * rowSums(trm)
  out <- expm::expAtv(A = t(trm), v = c(1, rep(0, ncol(trm) - 1)), t = t)$eAtv
  names(out) <- colnames(trm)
  ## Could this fail sometimes when it should not?
  ## stopifnot(isTRUE(all.equal(sum(out), 1.0)))
  if (!isTRUE(all.equal(sum(out), 1.0))) {
      ## Both warning and message, to simplify detecting it
      ## when using mclapply
      message("WARNING_all_equal_from_genots_at_t_from_trm.",
              " t = ", t,
              " Difference = ", all.equal(sum(out), 1.0))
      warning("WARNING_all_equal_from_genots_at_t_from_trm.",
              " t = ", t,
              " Difference = ", all.equal(sum(out), 1.0))
      ## This attributes hack allows me to carry through
      ## the message. Warnings not very useful when doing thousands
      ## and using mclapply
      attr(out, "all_equal_difference") <- list(t = t,
                                                all_equal_message = all.equal(sum(out), 1.0))
  }
  out
}

## transition rate matrix, optionally time -> predicted genotype composition
##  if time is NA, assume sampling with exponentially distributed
##  time of rate 1. Otherwise, sample at exactly t.
##  Recall MHN, CBN, H-ESBCN assume sampling with exponentially distributed
##  time of rate 1. OT and OncoBN have no such thing; neither does HyperHMM.
##
##  If custom_sampling = TRUE, use custom_sampling_function(trm) instead.
##  The default (probs_uniform_sampling_custom) averages the genotype
##  distribution over 101 equally-spaced time points in [0, 5]; this is
##  the "uniform sampling" regime. When custom_sampling = TRUE, t is
##  ignored (and we warn if it was set, to surface unintended mixing).
genots_from_trm <- function(trm, t = NA,
                            custom_sampling = FALSE,
                            custom_sampling_function = probs_uniform_sampling_custom) {
  if (custom_sampling) {
    if (!is.na(t))
      warning("custom_sampling = TRUE: argument t is ignored")
    return(custom_sampling_function(trm))
  }
  if (is.na(t)) {
    return(evamtools:::probs_from_trm(trm))
  } else {
    warning("Predicting genotype frequencies at specified t = ", t)
    return(genots_at_t_from_trm(trm = trm, t = t))
  }
}




## trans. rate matrix, number of genes in the trans. rate matrix,
## number of samples -> matrix of binary genotypes

## Outputs a matrix of binary genotypes. Each row corresponds to one
## sample obtained from the trm, when sampling with time
## distributed as exponential

sample_trm <- function(trm, n_genes, n_samples,
                       custom_sampling = FALSE,
                       custom_sampling_function = probs_uniform_sampling_custom) {
    if (n_samples < 1 || n_samples %% 1 != 0)
    stop("n_samples must be an integer greater than or equal to 1")
  if (n_genes < 1 || n_genes %% 1 != 0)
    stop("n_genes must be an integer greater than or equal to 1")

  ## Gets list of classes that the object inherits from
  ## check if it inherits from the matrix of basic R (in lower case)
  ## or if it is a Matrix from the Matrix package (upper case)
  if (!("matrix" %in% stringi::stri_trans_tolower(is(trm))))
    stop("trm must be a matrix")
  if (!custom_sampling) {
    genots_freq <- evamtools:::probs_from_trm(trm)
  } else {
    genots_freq <- custom_sampling_function(trm)
  }
  samp <- evamtools:::genot_probs_2_pD_ordered_sample(genots_freq,
                                                      n_genes,
                                                      LETTERS[1:n_genes],
                                                      n_samples,
                                                      out = "vector")
  ## repeat each genotype (names(samp)) the number of times it appears (samp)
  genots_to_bin(rep(names(samp), samp), n_genes)
}


## Fitness landscape, scaling -> scaled transition rate matrix, indicator
##   accessible genotypes
##  This function recomputes the trm after we intervene by manipulating
##  the fitness lanscape directly.
fitness_landscape_2_scaled_trm <- function(x, c) {
  fit_land_letters <- no_evam_rfitness_to_letter(x)
  rel_diff <-
    no_evam_genots_2_fgraph_and_trans_mat_rf(fit_land_letters)$relative_fitness_differences

  if (is.na(as.matrix(rel_diff)[1, 1])) {
    warning("No accesible genotypes: warning from get_scaled_trm")
    return(list(trm_scaled = NA, no_accessible_genotypes = TRUE))
  }

  ## rel_diff, from no_evam_genots_2_fgraph_and_trans_mat_rf,
  ## removes rows that are only destinations.
  ## Must be squared for genots_at_t_from_trm and probs_from_trm to work
  rel_diff_sq <- get_square_matrix(rel_diff, 0)
  trm_scaled <- c * rel_diff_sq
  return(list(trm_scaled = trm_scaled, no_accessible_genotypes = FALSE))
}


## ## For playing around: consequences of changing the value of "c"
## what_if_c <- function(x, c) {
##   p1 <- evamtools:::probs_from_trm(x[[1]]$trm)
##   fit_land_letters <- no_evam_rfitness_to_letter(x[[1]]$fitness_landscape)
##   rel_diff <-
##     no_evam_genots_2_fgraph_and_trans_mat_rf(fit_land_letters)$relative_fitness_differences
##   rel_diff_sq <- get_square_matrix(rel_diff, 0)
##   trm <- c * rel_diff_sq
##   return(list(p1 = p1, p2 = evamtools:::probs_from_trm(trm),
##               t1 = x[[1]]$trm, t2 = trm))
## }





## trm, vector of times, which are, now 100 points between 0 and 3 ->
##  probabilities of genotypes
probs_uniform_sampling_custom <- function(trm,
                                          times = seq(from = 0, to = 5,
                                                      length.out = 101)) {
    ## Does not preserve the error attributes
    ## pm <- vapply(times,
    ##              function(x) genots_at_t_from_trm(trm, x),
    ##              rep(0.0, nrow(trm))
    ##              )
    pm2 <- lapply(times, function(x) genots_at_t_from_trm(trm, x))
    pm <- do.call(cbind, pm2)
    p <- rowMeans(pm)

    ## Check if at any time we had a discrepancy with 1.0
    which_err_attr <- which(unlist(lapply(pm2, function(x) !(is.null(attr(x, "all_equal_difference"))))))
    if (length(which_err_attr)) {
        all_err_attr <- lapply(pm2[which_err_attr], function(x) attr(x, "all_equal_difference"))

    }

    ## Could this fail sometimes when it should not?
    ## stopifnot(isTRUE(all.equal(sum(p), 1)))

    if (!isTRUE(all.equal(sum(p), 1.0))) {
      ## Both warning and message, to simplify detecting it
      ## when using mclapply
      message("WARNING_all_equal_from_probs_uniform_sampling_custom.",
              " Difference = ", all.equal(sum(p), 1.0))
      warning("WARNING_all_equal_from_probs_uniform_sampling_custom.",
              " Difference = ", all.equal(sum(p), 1.0))
      sum_p_not_1 <- TRUE
    } else {
      sum_p_not_1 <- FALSE
    }

    ## Always return all genotypes, for consistency with probs_from_trm
    ## From probs_from_trm
    gene_names <- evamtools:::evam_string_sort(
      setdiff(unique(unlist(strsplit(colnames(trm),
                                     split = ", "))),
              "WT"))
    number_genes <- length(gene_names)
    num_genots <- 2^number_genes

    if (length(p) == num_genots) {
      if (length(which_err_attr))  attr(p, "all_err_attr") <- all_err_attr
      if (sum_p_not_1) attr(p, "sum_p_not_1") <- TRUE
      return(p)
    }

    ## These are two procedures of doing the same.
    ## I was once bitten by a similar issue when dealing with HyperTraPS
    ## So this is a paranoid procedure. rm one of the procedures eventually
    ## Procedure 1
    allGts <- evamtools:::genes_2_genotypes_standard_order(gene_names)
    p_all <- rep(0.0, length = length(allGts))
    names(p_all) <- allGts
    ## Next line should preclude any possible errors
    if (!(all(names(p) %in% allGts))) stop("mismatch in gene names")
    p_all[names(p)] <- p
    ## return(p_all)

    ## Procedure 2. Probably slightly slower?
    p <- evamtools:::reorder_to_standard_order(p)
    p[is.na(p)] <- 0
    stopifnot(identical(p, p_all))
    if (length(which_err_attr))  attr(p, "all_err_attr") <- all_err_attr
    if (sum_p_not_1) attr(p, "sum_p_not_1") <- TRUE
    return(p)
}


## Example usage custom sampling
## t1 <- random_evam(7, model = "CBN")
## trm1 <- t1$CBN_trans_mat
## probs_uniform_sampling_custom(trm1)





## When calling hittingProbabilities from markovchain
## all rows must sum to 1.
## So: anything without a 1 in the row is an absorbing state
## even if it were not, and it were disconnected, no
## problem
to_markovchain <- function(x) {
  no_sum_1 <- which(abs(rowSums(x) - 1) > sqrt(.Machine$double.eps))
  if (length(no_sum_1)) {
    x[cbind(no_sum_1, no_sum_1)] <- 1.0
  }
  return(new("markovchain", transitionMatrix = x))
}


## Transition probability matrix (embedded chain: rows sum to 1 for transient
## states, 0 for absorbing states — to_markovchain adds the self-loop) ->
## named vector of hitting probabilities from WT, using the
## first-passage convention: h(WT, WT) = 0.
## WT is located by name, not by position, since CPM trans_mats may not have
## WT as the first row.
hitting_probs_from_WT <- function(trans_mat) {
  mc <- to_markovchain(trans_mat)
  hp <- hittingProbabilities(mc)
  wt_row <- which(rownames(trans_mat) == "WT")
  if (length(wt_row) != 1) stop("WT not found exactly once in rownames(trans_mat)")
  result <- hp[wt_row, ]
  ## hittingProbabilities may drop names for small matrices; always restore them.
  names(result) <- colnames(trans_mat)
  return(result)
}


## Same output as hitting_probs_from_WT (named vector of hitting
## probabilities WT), but computed directly and just for WT as the origin.
## Used to check the output of markovchain:hittingProbabilities. Why? See
## comments before function threshold_transition_matrix.

hitting_probs_from_WT_direct <- function(trans_mat, absorb_tol = 1e-12) {

  ## What this code does it avoid the per-target (destination) linear
  ## solves that become singular whenever the chain has absorbing states
  ## other than the target; this can happen for kill-gene HyperHMM chains,
  ## which can have many absorbing genotypes. Our code avoids this, so it
  ## doesn't produce "system is singular" messages and is a lot cheaper
  ## (one solve on the transient block instead of one per target).

  ## Method (fundamental matrix). Divide the states into absorbing (no
  ## transition to a different state: off-diagonal row mass ~ 0) and
  ## transient (the rest). With Q = P[transient, transient] and R =
  ## P[transient, absorbing], the fundamental matrix N = (I - Q)^{-1} gives
  ## the expected number of visits.

  ## Then, starting from WT (which must be transient):

  ##   - for a transient target j != WT:
  ##                   P(ever reach j) = N[WT, j] / N[j, j];

  ##   - for an absorbing target j:
  ##                   P(reach j)  = (N %*% R)[WT, j];

  ##   - for WT itself: the RETURN probability 1 - 1 / N[WT, WT]. This is
  ##     the same as markovchain::hittingProbabilities, whose diagonal is
  ##     the return probability, not a forced 0: it is 0 only when WT has
  ##     no self-loop (e.g. no intervention), but killing puts probability
  ##     mass on WT's diagonal, so WT can re-hit itself and the return
  ##     probability equals that self-loop mass.
  ##

  ## States unreachable from WT get 0 automatically (N[WT, .] = 0 for
  ## them). (I - Q) is nonsingular because every transient state reaches an
  ## absorbing state with probability 1, so Q has spectral radius < 1.
  ## absorb_tol: a state is treated as absorbing when its total probability
  ## of moving to a different state is <= absorb_tol. This matters because
  ## HyperHMM's matrices after killing can have transition probabilities
  ## spanning a huge range (we have seen off-diagonal entries from ~1e-176
  ## to ~1). A state with a tiny-but-positive probability of moving out is,
  ## under the run-until-absorption (infinite time) semantics of
  ## markovchain::hittingProbabilities, a transient state (it escapes with
  ## probability 1). Thus, the WT hitting probabilities can change a lot
  ## (e.g., ~0.15 in one 9-gene example) depending on where we put the
  ## cutoff.

  ## What we do here is use a default of 1e-12, that considers any
  ## probability of moving out (of escaping) with (total) probability less
  ## than 1e-12 as zero (i.e., it does not move out). This also keeps (I -
  ## Q) well enough conditioned that this can be solved. See
  ## tests/hitting-probs-direct-and-drop-unreachable-TESTS.R (the MC ground
  ## truth uses the same cutoff).

  ## Note: this absorb_tol is not the same as the thresholding in function
  ## threshold_transition_matrix: in threshold_transition_matrix we set
  ## each single entry below 1e-12 to 0; here, we treat a genotype as
  ## absorbing if its whole row adds up to ≤ 1e-12 (i.e., we consider can't
  ## move away from this genotype unless the total movement, over all
  ## possible destination genotypes, is more than 1e-12 of probability
  ## mass).

  m <- as.matrix(trans_mat)
  n <- nrow(m)
  states <- rownames(m)
  wt <- which(states == "WT")
  if (length(wt) != 1) stop("WT not found exactly once in rownames(trans_mat)")

  result <- setNames(rep(0.0, n), states)

  ## Absorbing = no (non-negligible) probability of moving to a different
  ## state.
  off_diag_out <- rowSums(m) - diag(m)
  absorbing <- which(off_diag_out <= absorb_tol)
  transient <- setdiff(seq_len(n), absorbing)

  ## If WT itself is absorbing, nothing else is reachable: from WT you
  ## never move, so every other state has hitting probability 0, and WT
  ## itself is hit with probability 1 (recall you start there). WT = 1
  ## matches markovchain::hittingProbabilities convention for an absorbing
  ## state and the "everything killed -> only WT" hardcoded paths in
  ## get_genotype_freqs_cpm / the intervention functions.

  if (wt %in% absorbing) {
    result[wt] <- 1.0
    return(result)
  }

  Q <- m[transient, transient, drop = FALSE]
  N <- solve(diag(nrow(Q)) - Q)
  wt_t <- match(wt, transient) ## index of WT within the transient block

  ## Transient targets (j != WT): N[WT, j] / N[j, j].
  diagN <- diag(N)
  hp_trans <- N[wt_t, ] / diagN
  result[transient] <- hp_trans
  ## WT itself: return probability (0 if WT has no self-loop, else the
  ## self-loop mass); matches markovchain::hittingProbabilities' diagonal.
  result[wt] <- 1.0 - 1.0 / N[wt_t, wt_t]

  ## Absorbing targets: absorption probabilities (N %*% R)[WT, ].
  if (length(absorbing) > 0) {
    R <- m[transient, absorbing, drop = FALSE]
    absorb_probs <- N[wt_t, ] %*% R
    result[absorbing] <- as.vector(absorb_probs)
  }

  return(result)
}


## Set to zero the off-diagonal transition probabilities below `tol` and
## renormalize each row to sum to 1.

## Why? Fitted models --- HyperHMM especially --- can produce a large
## fraction of transition probabilities that are just numerical noise (in
## examples we checked, ~1/3 of a HyperHMM fit's entries are < 1e-12, and
## can be as small as ~1e-320). The model cannot represent an exact zero,
## so "never happens" can become ~1e-130. Those tiny entries can break
## markovchain::hittingProbabilities; see
## https://github.com/spedygiorgio/markovchain/issues/233.
##
##   - Before the fix of 2026-08-12, they made the per-target linear
##     solves ill-conditioned. markovchain returned negative
##     "probabilities", accompanied by a "system is singular" message
##     from Armadillo for every target.
##
##   - After the fix of 2026-08-12, things can still be very wrong, and
##     quietly (reported 2026-09-17)

##    - A further fix (branch fix/hittingprobs, commit db73e68, not on CRAN
##     as of 2026-09-18) gives, on our 9-gene example, the exact answer to
##     ~1e-14 with and without thresholding.

##  For now, we keep thresholding. Suppose a genotype (X) that could
##  transition with prob. = 1 - (1e-13) to a genotype with the killed gene,
##  and with prob. 1e-13 to another genotype (Y). If we do not threshold,
##  after killing, the hitting probability of X is affected by what other
##  tiny entries might have been left (e.g., genotypes that transition to
##  via tiny probs. to something else, instead of to X). Moreover, suppose
##  the only way to get to Y were via X; if we threshold, Y's hitting
##  probability will be 0, but if we don't, it won't (and it might actually
##  be relatively large if X's own hitting prob. is large); this behavior,
##  where X and Y can end up with the same hitting prob. is not sensible
##  when the transition prob. from X to Y is 1e-12 or less, given our
##  sample sizes.

##  And thresholding has another virtue: all three versions of markovchain
##  since we reported the bug (right before the report, the 2026-08-12 fix,
##  and the current branch fix/hittingprobs as of 2026-09-18) behave the
##  same and also the same as our own code for hitting probs and with
##  Monte-Carlo results. Given the state of flux of markovchain, and given
##  that the fix as of 2026-09-18 is not yet in CRAN, this seems the
##  sensible choice. (Other options, such as thresholding for our code but
##  not for markovchain, and using the latest markovchain version, put too
##  large a burden on users and developers, and might give alarming
##  messages of a large `divergence` from `hitting_probs_from_WT_both` that
##  have no practical relevance)

## VERY IMPORTANT: this *ONLY THRESHOLDS* tiny transition probabilities to
## zero. It does *NOT* drop or otherwise remove genotypes (e.g. killed
## genotypes with all-zero rows are left in place, just as they are) ---
## dropping is a separate issue (see function
## kill_gene_HyperHMM_drop_unreachable).

threshold_transition_matrix <- function(m, tol = 1e-12) {
  tm <- as.matrix(m)
  offdiag <- row(tm) != col(tm)
  tm[offdiag & (tm < tol)] <- 0
  rs <- rowSums(tm)
  pos <- rs > 0
  tm[pos, ] <- tm[pos, ] / rs[pos]
  return(tm)
}


## Compute the WT hitting probabilities both ways: with markovchain package
## (hitting_probs_from_WT) and our code above for the direct method
## (hitting_probs_from_WT_direct). Return both, with markovchain's result
## as the canonical, and record their disagreement (and give warning if
## large). See details is hitting-probs-direct-and-drop-unreachable-TESTS.R
##
## 'context' is a short string (e.g. the method, or method and gene) shown
## in the warning to help locate the computation with the problem.
##
## However, after using the thresholding on HyperHMM results (see function
## threshold_transition_matrix) we should not see any serious diveregences.

hitting_probs_from_WT_both <- function(trans_mat,
                                       context = "", tol = 1e-6) {

  ## The markovchain result is kept as the canonical value (hp) for now;
  ## the direct result (hp_direct) is also returned. Their disagreement is
  ## always recorded as divergence = max(abs(hp - hp_direct)) over the
  ## shared genotypes --- even when it is below tol. A warning (naming the
  ## worst genotypes) is emitted only when divergence > tol. We compare on
  ## max(abs(.)), not all.equal: all.equal averages the (relative)
  ## differences, so a localized spike, which is how the markovchain bug
  ## was found, could be diluted. max(abs(.)) catches the worst genotype
  ## directly and is the same quantity we return in the results.
  ##
  ## markovchain's own "system is computationally singular" warnings are no
  ## longer suppressed here: our explicit divergence is the difference
  ## between the two methods; this can catch wrong markovchain answers even
  ## when they are returned silently, but we want to see the warnings, if
  ## they appear. And it seems that the actual message could differ between
  ## R/Armadillo/BLAS versions.

  ##  hp <- suppressWarnings(hitting_probs_from_WT(trans_mat))
  hp <- hitting_probs_from_WT(trans_mat)
  hp_direct <- hitting_probs_from_WT_direct(trans_mat)
  nm <- names(hp)
  absdiff <- abs(hp[nm] - hp_direct[nm])
  divergence <- max(absdiff)
  if (is.finite(divergence) && (divergence > tol)) {
    worst <- nm[order(absdiff, decreasing = TRUE)]
    worst <- worst[seq_len(min(5, length(worst)))]
    warning("hitting_probs divergence (markovchain vs direct)",
            if (nzchar(context)) paste0(" [", context, "]") else "",
            ": max abs diff = ", format(divergence),
            "; worst genotypes: ", paste(worst, collapse = ", "))
  }
  return(list(hp = hp, hp_direct = hp_direct, divergence = divergence))
}


library(codetools)
checkUsageEnv(env = .GlobalEnv)

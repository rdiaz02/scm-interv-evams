### "Standard" HyperHMM code: before the HyperHMM package existed

## This is no longer necessary, because a HyperHMM package now
## exists: https://github.com/StochasticBiology/hyperhmm

if (FALSE) {
  ## Based on
  ## https://github.com/rdiaz02/HT_test_run/blob/igj-hyperhmm/K_1_example/example_run.R

  library(Rcpp)
  library(RcppArmadillo)
  library(igraph)
  library(ggraph)
  library(stringr)
  ## The next file must either exist or by a symlink to the
  ## existing one
  ## It comes from https://github.com/StochasticBiology/hypercube-hmm
  sourceCpp("./HyperHMM_cpp_and_other_code-obsolete/hyperhmm-r.cpp")
}


### Load the HyperHMM package. If needed, install

library(hyperhmm)
## If not installed, install as
## remotes::install_github("StochasticBiology/hyperhmm")

### Our code for running HyperHMM

library(evamtools)
library(Matrix)


## transition matrix, vector of weights of each "time", scalar (number features) ->
##        predicted genotype frequencies, and predicted genotype frequencies at
##        each time period

## Much of this is already available from the output of HyperHMM
## but this function computes everything from transition matrix
## and can be reused for interventions.
## There is testing of this function in intervention-TESTS.R
probs_from_HyperHMM <- function(trans_mat,
                                num_prob.set,
                                num_features) {

    v <- rep(0.0, length = nrow(trans_mat))
    names(v) <- rownames(trans_mat)
    v["WT"] <- 1.0

    genot_pred_t <- Matrix(0,
                           nrow = length(v),
                           ncol = num_features + 1,
                           sparse = TRUE)
    rownames(genot_pred_t) <- names(v)
    colnames(genot_pred_t) <- c(0, 1:num_features)
    stopifnot(isTRUE(identical(rownames(trans_mat), names(v))))
    stopifnot(isTRUE(identical(colnames(trans_mat), names(v))))
    genot_pred_t[, 1] <- v

    for (i in 1:num_features) {
        genot_pred_t[, (i + 1)] <- t(genot_pred_t[, i]) %*% trans_mat
    }

    stopifnot(isTRUE(all.equal(colSums(genot_pred_t),
                               rep(1.0, num_features + 1),
                               check.attributes = FALSE)))

    predicted_genotype_freqs <- genot_pred_t %*% matrix(num_prob.set, ncol = 1)
    pgf <- as.vector(predicted_genotype_freqs)
    names(pgf) <- rownames(predicted_genotype_freqs)
    stopifnot(isTRUE(all.equal(sum(pgf), 1.0)))
    return(list(predicted_genotype_freqs = pgf,
                predicted_genotype_freq_at_t = genot_pred_t,
                num_prob.set = num_prob.set))
}

## Denoising of the HyperHMM transition matrix.
##
## This is done once, and only once, in run_HyperHMM, right after the
## HyperHMM package returns its fit, and before anything else is
## computed from the transition matrix (predicted genotype frequencies,
## and, later, interventions and hitting probabilities). The raw,
## untouched, matrix is kept in the output (HyperHMM_trans_mat_raw, in
## evam_like_HyperHMM). HyperHMM ought to give this already, but it
## does not, so we sanitize it first. The threshold and the rule used
## are recorded in the output (HyperHMM_threshold,
## HyperHMM_threshold_rule).

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

##  An alternative (suggested by Iain Johnston), not implemented.

##  The EM algorithm of HyperHMM does not perfectly converge. Setting these
##  values to 0, as we do here, rules those transitions out, which the data
##  cannot support either. His alternative: set every entry below the
##  threshold to the threshold itself (e.g., 1e-12), and then renormalize
##  each row. Then no tiny exit is favoured over another because of noise,
##  and a genotype whose large exit is killed can still move on, with the
##  same probability to each of the remaining exits. Neither rule is
##  supported by the data: one says "never", the other says "all equally
##  likely".
##
##  Both rules need a threshold. Iain sometimes uses 0.01/N (N, the
##  number of samples), since a data set of that size cannot support
##  smaller probabilities; this is much larger than our 1e-12.
##


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


run_HyperHMM <- function(x, opts = list(seed = -1, prob.set = "observed",
                                        threshold = 1e-12,
                                        threshold_rule = "set_to_zero")) {
    ## The seed only matters for bootstrap and for the
    ## $stats component, that we do not use (for now)
    if (is.null(opts$seed) || is.na(opts$seed) || (opts$seed == -1))
        opts$seed <- round(runif(1, 1, 1e9))

    if (is.null(opts$prob.set)) opts$prob.set <- "observed"
    if (is.null(opts$threshold)) opts$threshold <- 1e-12
    if (is.null(opts$threshold_rule)) opts$threshold_rule <- "set_to_zero"

    time_out <- system.time({
        out <- HyperHMM(x, nboot = 0, seed = opts$seed)
    })["elapsed"]

    ## Next copies heavily from run_HyperTraPS
    ## file evam_main_utils_run_methods.R
    num_features <- ncol(x)
    feature_labels <- colnames(x)
    states <- unique(c(out$transitions$From, out$transitions$To))

    decoded_states <- vapply(states, evamtools:::decode_state,
                             character(1),
                             num_features = num_features,
                             feature_labels = feature_labels)

    trans_mat <- Matrix(0,
                        nrow = length(states),
                        ncol = length(states),
                        sparse = TRUE)

    rownames(trans_mat) <- colnames(trans_mat) <-
        evamtools:::reorder_genotypes_2_standard_order(decoded_states)

    flux_mat <- trans_mat

    for (i in 1:nrow(out$transitions)) {
        from_state <- out$transitions$From[i]
        to_state <- out$transitions$To[i]
        probability <- out$transitions$Probability[i]
        flux <- out$transitions$Flux[i]
        from_decoded <- evamtools:::decode_state(from_state,
                                                 num_features,
                                                 feature_labels)
        to_decoded <- evamtools:::decode_state(to_state,
                                               num_features,
                                               feature_labels)
        trans_mat[from_decoded, to_decoded] <- probability
        flux_mat[from_decoded, to_decoded] <- flux
    }

    ## Denoise the transition matrix: this is the only place where this
    ## is done. See the comments before threshold_transition_matrix.
    ## We keep the raw, untouched, matrix, and everything else below
    ## (predicted genotype frequencies) and later (interventions, hitting
    ## probabilities) uses the denoised matrix.
    ## Rule "make_all_tiny_equally_tiny" is the alternative suggested by
    ## Iain Johnston (see comments before threshold_transition_matrix);
    ## not implemented yet.
    trans_mat_raw <- trans_mat
    if (opts$threshold_rule == "set_to_zero") {
        trans_mat <- Matrix(threshold_transition_matrix(trans_mat_raw,
                                                        tol = opts$threshold),
                            sparse = TRUE)
    } else if (opts$threshold_rule == "make_all_tiny_equally_tiny") {
        stop("threshold_rule make_all_tiny_equally_tiny not implemented yet")
    } else {
        stop("Unrecognized threshold_rule option")
    }

    if ((length(opts$prob.set) == 1) &&
        (is.character(opts$prob.set))) {
        if (opts$prob.set == "uniform") {
            num_prob.set <- rep(1/(ncol(x) + 1), ncol(x) + 1)
            names(num_prob.set) <- 0:(ncol(x))
        } else if (opts$prob.set == "observed") {
            num_prob.set <- evamtools:::props_num_muts(x)
        } else {
            stop("Unrecognized prob.set option")
        }
    } else {
        num_prob.set <- opts$prob.set
    }

    message("HyperHMM prob.set option = ",
            paste(opts$prob.set, collapse = " "),
            ". Value passed as num_prob.set = ",
            paste(num_prob.set, collapse = " "))


    time2 <- system.time({
        probs_hyper_hmm <- probs_from_HyperHMM(trans_mat, num_prob.set, num_features)
    })["elapsed"]

    ## For interventions
    attr(trans_mat, "method_output") <- "HyperHMM_trans_mat"
    attr(trans_mat, "num_prob.set") <- num_prob.set
    attr(trans_mat, "num_features") <- num_features
    ## The same attributes for the raw matrix, so it can be used directly
    ## if ever needed
    attr(trans_mat_raw, "method_output") <- "HyperHMM_trans_mat"
    attr(trans_mat_raw, "num_prob.set") <- num_prob.set
    attr(trans_mat_raw, "num_features") <- num_features
    return(list(time_out = time_out + time2,
                out = c(primary_output = list(out),
                        trans_mat = list(trans_mat)
                      , trans_mat_raw = list(trans_mat_raw)
                      , threshold = list(opts$threshold)
                      , threshold_rule = list(opts$threshold_rule)
                      , predicted_genotype_freqs = list(probs_hyper_hmm$predicted_genotype_freqs)
                      , conditional_genotype_freqs = list(probs_hyper_hmm$predicted_genotype_freq_at_t)
                      , used_prob.set = list(num_prob.set)
                        ## , Prob_Cond_Prob_df = list(tmp$HyperTraps_Prob_Cond_Prob)
                        )))
}

## Copies heavily from evam_like_MHN_python
## in file MHN_python.R
evam_like_HyperHMM <- function(x,
                               opts = list(seed = -1, prob.set = "observed",
                                           threshold = 1e-12,
                                           threshold_rule = "set_to_zero"),
                               max_cols = 15) {
    cat("\n Starting a HyperHMM run\n")
    ## ########      Preprocessing: common to all methods
    x <- evamtools:::df_2_mat_integer(x)
    ## xoriginal <- x

    x <- evamtools:::add_pseudosamples(x)
    ## remove.constant makes no difference IFF we add pseudosamples, as
    ## there can be no constant column when we add pseudosamples
    x <- evamtools:::pre_process(x, remove.constant = FALSE,
                                 min.freq = 0, max.cols = max_cols)

    theout <- run_HyperHMM(x, opts)

    outlist <- list(
        ## The denoised transition matrix: the one used for everything
        HyperHMM_trans_mat                = theout$out$trans_mat,
        ## The raw, untouched, transition matrix from HyperHMM
        HyperHMM_trans_mat_raw            = theout$out$trans_mat_raw,
        HyperHMM_threshold                = theout$out$threshold,
        HyperHMM_threshold_rule           = theout$out$threshold_rule,
        HyperHMM_predicted_genotype_freqs = theout$out$predicted_genotype_freqs,
        HyperHMM_elapsed_time             = theout$time_out,
        HyperHMM_used_prob.set            = theout$out$used_prob.set,
        HyperHMM_conditional_genotype_freqs = theout$out$conditional_genotype_freqs,
        ## Oooops, I could have removed parts that are already above.
        ## Oh well. Beware component "viz" can be huge in size.
        ## FIXME: for evamtools, give option not to save, and make default
        HyperHMM_rest_stuff               = theout$out$primary_output,
        ## Makes intervention a lot simpler
        HyperHMM_gene_names               = colnames(x)
    )
    cat("\n      Finished a HyperHMM run. Elapsed = ",
        theout$time_out, "\n")
    return(outlist)
}




#### Example
if (FALSE) {
    local({
        rmhn <- random_evam(model = "MHN", ngenes = 5)
        sample_mhn <- sample_evam(rmhn, N = 1000, obs_noise = 0.05)
        dd <- sample_mhn$MHN_sampled_genotype_counts_as_data
        ## The next are silly reruns that show the seed does not change
        o1 <- evam_like_HyperHMM(dd)
        o1b <- evam_like_HyperHMM(dd)
        os2 <- evam_like_HyperHMM(dd, opts= list(seed = 2))
        os3 <- evam_like_HyperHMM(dd, opts= list(seed = 3))
    })
}


library(codetools)
checkUsageEnv(env = .GlobalEnv)

### "Standard" HyperHMM code
## From
## https://github.com/rdiaz02/HT_test_run/blob/igj-hyperhmm/K_1_example/example_run.R

library(Rcpp)
library(RcppArmadillo)
library(igraph)
library(ggraph)
library(stringr)
## The next file must either exist or by a symlink to the
## existing one
## It comes from https://github.com/StochasticBiology/hypercube-hmm
sourceCpp("./HyperHMM_cpp_and_other_code/hyperhmm-r.cpp")

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

run_HyperHMM <- function(x, opts = list(seed = -1, prob.set = "observed")) {
    ## The seed only matters for bootstrap and for the
    ## $stats component, that we do not use (for now)
    if (is.null(opts$seed) || is.na(opts$seed) || (opts$seed == -1))
        opts$seed <- round(runif(1, 1, 1e9))

    if (is.null(opts$prob.set)) opts$prob.set <- "observed"

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
    return(list(time_out = time_out + time2,
                out = c(primary_output = list(out),
                        trans_mat = list(trans_mat)
                      , predicted_genotype_freqs = list(probs_hyper_hmm$predicted_genotype_freqs)
                      , conditional_genotype_freqs = list(probs_hyper_hmm$predicted_genotype_freq_at_t)
                      , used_prob.set = list(num_prob.set)
                        ## , Prob_Cond_Prob_df = list(tmp$HyperTraps_Prob_Cond_Prob)
                        )))
}

## Copies heavily from evam_like_MHN_python
## in file MHN_python.R
evam_like_HyperHMM <- function(x,
                               opts = list(seed = -1, prob.set = "observed"),
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
        HyperHMM_trans_mat                = theout$out$trans_mat,
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

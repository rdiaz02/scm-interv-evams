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

##  The main intervention code: these are the upper-level functions,
##  - intervene on every gene in a CPM,
##  - intervene on every gene in a fitness landscape





### Loading/sourcing dependencies

source("kill-gene-and-output-from-cpm.R")
## source("trm.R") ## pulled from the above
library(parallel)

### Code


## Drop entries with value <= 0 from a hitting-probs vector, but always
## keep WT. WT's hitting prob from itself is 0 in any non-absorbing chain
## (by the convention used in hitting_probs_from_WT: starting state excluded
## from "hitting"), so the bare > 0 filter would drop WT. We keep WT
## explicitly so the HP vector parallels the genot_freqs vector in WT-
## presence, and downstream consumers (e.g., interv_genotype_list_2_matrix
## inside compute_intervention_objectives) can rely on WT being there.
filter_hp_keep_wt <- function(hp) {
    return(hp[(hp > 0) | (names(hp) == "WT")])
}


## a list with one of more CPM models, method, optionally time, filename,
##   procedure for killing, verbose -> for each intervention (no_intervention
##       and one per gene), a list(genot_freqs = ..., hitting_probs_from_WT = ...)
##       with the post-intervention predictions. Names are "no_intervention"
##       and "I:<gene>".
## time: if NA, assume sampling with exponentially distributed rate 1
##       for CBN, HESBCN, MHN.
## cpm_output is a list, with output from possibly several methods
## procedure is either the "standard" one, intervening in each DAG
## (kill_gene) or the "set parameters to 0" (kill_gene_by_params_to_0)
## verbose = TRUE shows the killing function used
intervene_cpm_every_gene <- function(cpm_output,
                                     method = c("OT", "OncoBN",
                                                "CBN", "HESBCN",
                                                "MHN", "HyperHMM"),
                                     t = NA,
                                     filename = NA,
                                     kill_gene_funct = kill_gene,
                                     verbose = FALSE,
                                     mc.cores = getOption("intervention_every_gene_cores",
                                                          detectCores())) {
    if (length(method) != 1) stop("length(method) != 1")
    method <- match.arg(method)
    if ((method == "HyperHMM")  &&
        (deparse(substitute(kill_gene_funct)) != "kill_gene"))
      stop("HyperHMM can only use 'kill_gene' for kill_gene_funct")
  if (method == "MHN") {
      model <- cpm_output$MHN_theta
  } else if (method == "HyperHMM") {
      model <- cpm_output$HyperHMM_trans_mat
  } else {
      model <- cpm_output[[paste0(method, "_", "model")]]
  }
  ## item <- ifelse(method == "MHN", "theta", "model")
  ## model <- cpm_output[[paste0(method, "_", item)]]
  if (length(model) == 0)
      stop("Input contains no CPM with requested method")

  ## Get all the gene names
  if (method == "MHN") {
    genes <- colnames(model)
  } else if (method == "HyperHMM") {
      genes <- sort(cpm_output$HyperHMM_gene_names)
  } else {
      genes <- unique(c(model[, "From"], model[, "To"]))
      genes <- genes[-which(genes == "Root")]
      genes <- sort(genes)
  }

  NO_INTERV_STR <- "no_intervention"
  interventions <- c(NO_INTERV_STR, genes)
  intervene_all_genes <- mclapply(interventions, function(gene) {
      if (gene == NO_INTERV_STR) {
          ## If we had passed a fitted model with predicted
          ## genot freqs. already there, this would recompute them
          ## (as get_genotype_freqs_cpm calls get_full_output).
          ## The virtue of doing it this way is that we can
          ## just pass any model, no need for all the rest of the stuff
          ## and this will work. Use more complex logic to
          ## avoid recomputing if existing?
          tmp <- get_genotype_freqs_cpm(model, t = t)
      } else {
      model_after_intervention <- kill_gene_funct(model, gene,
                                                  verbose = verbose)
          tmp <- get_genotype_freqs_cpm(model_after_intervention, t = t)
    }
      return(list(genot_freqs = tmp$genot_freqs[tmp$genot_freqs > 0],
                  hitting_probs_from_WT = filter_hp_keep_wt(tmp$hitting_probs_from_WT)))
  }, mc.cores = mc.cores)
  intervention_names <- c(NO_INTERV_STR, paste0("I:", interventions[-1]))
  out_list <- setNames(intervene_all_genes, intervention_names)
  if (!is.na(filename)) saveRDS(out_list, filename)
  return(out_list)
}






## object with fitness landscape (and c scaling and trans rate mat ) ->
##    genotype frequencies
##    after making each gene lethal.
##    If t = NA, sample with exponentially distributed time of 1,
##     o.w., sample at exactly t
##    If custom_sampling = TRUE, genotype frequencies are obtained via
##     custom_sampling_function(trm) instead, and t is ignored. The
##     default custom_sampling_function is probs_uniform_sampling_custom,
##     which averages the distribution over 101 equally-spaced time
##     points in [0, 5] (the "uniform sampling" regime).
##    Hitting probabilities are sampling-time independent and are not
##     affected by custom_sampling.
##    Optionally, write to rds if filename is not NA.
## The object is like the ones generated from generate_n_f_landscape_requir
## Used to be called intervene_fitness_every_gene
intervene_fitness_landscape_every_gene <- function(x,
                                                   t = NA,
                                                   filename = NA,
                                                   custom_sampling = FALSE,
                                                   custom_sampling_function = probs_uniform_sampling_custom,
                                                   mc.cores = getOption("intervention_every_gene_cores",
                                                                        detectCores())) {
    NO_INTERV_STR <- "no_intervention"
    genes <- colnames(x$fitness_landscape)[-ncol(x$fitness_landscape)]

    out <- mclapply(genes, function(gene) {
        intervened_fitness <-
            kill_gene_fitness_landscape(x$fitness_landscape, gene)
        intervened_trm_scaled <- fitness_landscape_2_scaled_trm(intervened_fitness, c = x$c)
        ## If not valid TRM because there are no accesible genotypes
        if (intervened_trm_scaled$no_accessible_genotypes) {
            return(list(genot_freqs = c(WT = 1),
                        hitting_probs_from_WT = c(WT = 1.0)))
        }
        trm <- intervened_trm_scaled$trm_scaled
        tmp_genot_freqs <- genots_from_trm(trm, t = t,
                                           custom_sampling = custom_sampling,
                                           custom_sampling_function = custom_sampling_function)
        ## Embedded chain via competing exponentials (row-scale rows with
        ## rowSums > 0; absorbing states keep all-zero rows for to_markovchain)
        rs <- rowSums(trm)
        embedded <- trm
        embedded[rs > 0, ] <- trm[rs > 0, ] / rs[rs > 0]
        tmp_hp <- hitting_probs_from_WT(embedded)
        return(list(genot_freqs = tmp_genot_freqs[tmp_genot_freqs > 0],
                    hitting_probs_from_WT = filter_hp_keep_wt(tmp_hp)))
    }, mc.cores = mc.cores)

    ## Prepend no intervention
    tmp_ni <- genots_from_trm(x$trm_scaled, t = t,
                              custom_sampling = custom_sampling,
                              custom_sampling_function = custom_sampling_function)
    trm_ni <- x$trm_scaled
    rs_ni <- rowSums(trm_ni)
    embedded_ni <- trm_ni
    embedded_ni[rs_ni > 0, ] <- trm_ni[rs_ni > 0, ] / rs_ni[rs_ni > 0]
    hp_ni <- hitting_probs_from_WT(embedded_ni)

    out <- c(list(list(genot_freqs = tmp_ni[tmp_ni > 0],
                       hitting_probs_from_WT = filter_hp_keep_wt(hp_ni))),
             out)

    names(out)[1] <- NO_INTERV_STR
    names(out)[-1] <- paste0("I:", genes)

    if (!is.na(filename)) saveRDS(out, filename)
    return(out)
}



## CPM model (CBN, HESBCN, MHN only), gene, method -> transition rate
##   matrix after removing all genotypes with gene "gene"
## DO NOT USE THIS in general as limited to a few methods.
rm_genots_trm <- function(x, gene, method) {
  trm <- x[[paste0(method, "_trans_rate_mat")]]
  rown <- rownames(trm)
  coln <- colnames(trm)
  rm_rown <- grep(gene, rown, fixed = TRUE)
  rm_coln <- grep(gene, coln, fixed = TRUE)

  if (length(rm_rown) != length(rm_coln))
      stop("length(rm_rown) != length(rm_coln)")
  if(length(rm_rown) == 0) {
      warning("No genotypes to remove?")
      trm2 <- trm
  } else {
      ## Should always be square
      ## If it ends up with 1 row and 1 column, and just a 0, it means only WT
      trm2 <- trm[-rm_rown, -rm_coln, drop = FALSE]
  }
  return(trm2)
}

## a list with one of more CPM models (only CBN, HESBCN, MHN), method,
##      optionally time -> genotype frequencies after intervening in each
##       gene
## time: if NA, assume sampling with exponentially distributed rate 1
## DO NOT USE THIS in general as limited to a few methods.
## This is used for testing the equivalence of the standard intervention
## (removal from the TRM)
intervene_cpm_trm_rm_every_gene <- function(cpm_output,
                                            method = c("CBN", "HESBCN",
                                                       "MHN"),
                                            t = NA,
                                            filename = NA,
                                            mc.cores = getOption("intervention_every_gene_cores",
                                                                  detectCores())) {
    if (length(method) != 1) stop("length(method) != 1")
    method <- match.arg(method)
    item <- ifelse(method == "MHN", "theta", "model")
    model <- cpm_output[[paste0(method, "_", item)]]
    if (length(model) == 0)
    stop("Input contains no CPM with requested method")

  ## Get all the gene names
  if (method == "MHN") {
    genes <- colnames(model)
  } else {
    genes <- unique(c(model[, "From"], model[, "To"]))
    genes <- genes[-which(genes == "Root")]
    genes <- sort(genes)
  }

  NO_INTERV_STR <- "no_intervention"
  interventions <- c(NO_INTERV_STR, genes)

  intervene_all_genes <- mclapply(interventions, function(gene) {
      if (gene == NO_INTERV_STR) {
          tmp <- get_genotype_freqs_cpm(model, t = t)
          tmp_genot_freqs <- tmp$genot_freqs
          tmp_hp <- tmp$hitting_probs_from_WT
      } else {
          trm_after_intervention <- rm_genots_trm(cpm_output, gene, method)
          if (nrow(trm_after_intervention) > 1) {
              tmp_genot_freqs <-
                  genots_from_trm(trm_after_intervention, t = t)
          } else { ## All genots killed
              ## Must be a 0.
              stopifnot(as.vector(trm_after_intervention) == 0)
              tmp_genot_freqs <- c(WT = 1)
          }
          ## Compute embedded chain from Q matrix (negative diagonal, rows sum to 0):
          ## zero out diagonal to get off-diagonal rates, then row-normalize.
          off_q <- trm_after_intervention
          diag(off_q) <- 0
          rs <- rowSums(off_q)
          embedded <- off_q
          embedded[rs > 0, ] <- off_q[rs > 0, ] / rs[rs > 0]
          tmp_hp <- hitting_probs_from_WT(embedded)
      }
      return(list(genot_freqs = tmp_genot_freqs[tmp_genot_freqs > 0],
                  hitting_probs_from_WT = filter_hp_keep_wt(tmp_hp)))
  }, mc.cores = mc.cores)
  names(intervene_all_genes) <- interventions
  names(intervene_all_genes)[-1] <- paste0("I:", interventions[-1])

  if (!is.na(filename)) saveRDS(intervene_all_genes, filename)
  return(intervene_all_genes)
}



library(codetools)
checkUsageEnv(env = .GlobalEnv)

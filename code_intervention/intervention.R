## Copyright 2022 Ramon Diaz-Uriarte

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


### What this is

##  The main intervention code: these are the upper-level functions,
##  - intervene on every gene in a CPM,
##  - intervene on every gene in a fitness landscape


### Notation, equivalences of killing, mapping to the notation in the manuscript


## The same intervention on gene g can be implemented in more than one way,
## and all of them give the same predictions. This is what the paper ("A
## structural causal framework for interventions on evolutionary
## accumulation models" https://arxiv.org/pdf/2606.12597) calls each of
## them, and where each one is in the code. The equivalence between them is
## checked in tests/kill-gene-equivalences-TESTS.R and in the per-model
## files tests/kill-CBN-TESTS.R, kill-HESBCN-TESTS.R, kill-MHN-TESTS.R,
## kill-OT-TESTS.R, kill-OncoBN-TESTS.R.

## The main entry point is intervene_cpm_every_gene (this file). Its
## argument kill_gene_funct defaults to kill_gene, the default
## procedures below. Passing kill_gene_by_params_to_0 gives the
## parameters-to-0 route, etc.

#### The canonical way: set parameter to 0: DAG_{lambda_g = 0}, Theta_{g, g} = 0

## The "canonical definition" in the paper uses the idea of "set the
## parameter to 0" which is a very direct mapping from the do operator.

## The common entry point is : kill_gene_by_params_to_0, with specific
## versions according to model.

## DAG_{lambda_g = 0}, for CBN and H-ESBCN, and the same idea for OT
##   (pi_{parent(g), g} = 0) and OncoBN (theta_g = 0): set to 0 the
##   parameter of the edges that lead into g. Function
##   kill_gene_DAG_param_0, called from kill_gene_by_params_to_0, in
##   file kill-gene-and-output-from-cpm.R. We use this to check that
##   the default gives the same predictions.

## Theta_{g, g} = 0, for MHN: since we store log-thetas, this means
##   theta_{g, g} = -Inf. Function kill_gene_MHN_theta_minus_Inf,
##   called from kill_gene_by_params_to_0, in file
##   kill-gene-and-output-from-cpm.R. Again, used for checking.

#### Removing entries "in the model": DAG_{-g}, Theta_{-g}

## For computational reasons, this is what the code uses by default.

## DAG_{-g}, for CBN, H-ESBCN, OT and OncoBN: remove, from the DAG of
##   restrictions, the node of gene g and any other node that only
##   appears, in accessible genotypes, together with g. This is the
##   default. Function kill_gene_DAG, called from kill_gene, in file
##   kill-gene-and-output-from-cpm.R.


## Theta_{-g}, for MHN: remove the row and the column of gene g from
##   the matrix of log-thetas. This is the default. Function
##   kill_gene_MHN, called from kill_gene, in file
##   kill-gene-and-output-from-cpm.R.

#### Removing entries from the transition rate matrix. Q_{-g}

## Limited to models with transition rate matrix, so CBN, H-ESBCN, MHN

## Q_{-g}, for CBN, H-ESBCN and MHN: remove, from the transition rate
##   matrix, the rows and columns of all genotypes that contain gene
##   g. Function rm_genots_trm, called from
##   intervene_cpm_trm_rm_every_gene, both in this file. This one is
##   only used in the tests. See file Q_g_paranoid_checks.org.


#### HyperHMM, R_{-g}

## R_{-g}, for HyperHMM: the conditional transition matrix without the
##   genotypes with gene g. Function
##   kill_gene_HyperHMM_drop_unreachable, called from kill_gene, in
##   file kill-gene-and-output-from-cpm.R. HyperHMM has no equivalent
##   parameter to set to 0.


#### Modifying the fitness landscape, for CBN and H-ESBCN under SSWM:

##   As said: this requires a fitness landscape, etc. We make lethal all
##   genotypes with gene g.

##   Function  kill_gene_fitness_landscape, called from
##   intervene_fitness_landscape_every_gene (this file).







### Loading/sourcing dependencies

source("kill-gene-and-output-from-cpm.R")
## source("trm.R") ## pulled from the above
library(parallel)

### Code

## Set a new random seed and print it. Used in the tests instead of
## set.seed(NULL), so that a failing run can be replayed: look for the
## last "Seed used was" printed before the failure, and replace the
## call to set_and_print_seed() there by set.seed(<that number>).
set_and_print_seed <- function() {
  ## First, set.seed(NULL): it re-seeds from the clock and the process
  ## ID. Without it, if an earlier line did, say, set.seed(1), the
  ## "random" seed drawn below would be the same number in every run.
  set.seed(NULL)
  ## sample.int gives a whole number between 1 and 1e9, a valid seed
  ## (the largest integer R allows is about 2.1e9).
  seed <- sample.int(1e9, 1)
  set.seed(seed)
  cat("\n Seed used was ", seed, "\n")
  return(invisible(seed))
}



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
                  hitting_probs_from_WT = filter_hp_keep_wt(tmp$hitting_probs_from_WT),
                  hitting_probs_from_WT_direct = filter_hp_keep_wt(tmp$hitting_probs_from_WT_direct),
                  divergence_hitting_prob_calculation = tmp$divergence_hitting_prob_calculation))
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
                  hitting_probs_from_WT = c(WT = 1.0),
                  hitting_probs_from_WT_direct = c(WT = 1.0),
                  divergence_hitting_prob_calculation = 0.0))
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
    both <- hitting_probs_from_WT_both(embedded, context = paste0("fitness_landscape I:", gene))
    return(list(genot_freqs = tmp_genot_freqs[tmp_genot_freqs > 0],
                hitting_probs_from_WT = filter_hp_keep_wt(both$hp),
                hitting_probs_from_WT_direct = filter_hp_keep_wt(both$hp_direct),
                divergence_hitting_prob_calculation = both$divergence))
  }, mc.cores = mc.cores)

  ## Prepend no intervention
  tmp_ni <- genots_from_trm(x$trm_scaled, t = t,
                            custom_sampling = custom_sampling,
                            custom_sampling_function = custom_sampling_function)
  trm_ni <- x$trm_scaled
  rs_ni <- rowSums(trm_ni)
  embedded_ni <- trm_ni
  embedded_ni[rs_ni > 0, ] <- trm_ni[rs_ni > 0, ] / rs_ni[rs_ni > 0]
  both_ni <- hitting_probs_from_WT_both(embedded_ni, context = "fitness_landscape no_intervention")

  out <- c(list(list(genot_freqs = tmp_ni[tmp_ni > 0],
                     hitting_probs_from_WT = filter_hp_keep_wt(both_ni$hp),
                     hitting_probs_from_WT_direct = filter_hp_keep_wt(both_ni$hp_direct),
                     divergence_hitting_prob_calculation = both_ni$divergence)),
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
  ## Row names and column names are the genotypes. They must be the
  ## same and in the same order, so that one index works for both.
  if (!identical(rownames(trm), colnames(trm)))
    stop("rownames and colnames of the transition rate matrix differ")
  ## Our transition rate matrices have a zero diagonal: they store only
  ## the rates between different genotypes, and the diagonal, when
  ## needed, is computed as minus the sum of the other entries of the
  ## row. This is crucial here.
  ## Suppose the new matrix kept the diagonal of the full matrix. That
  ## diagonal is minus the sum of all the rates in the row of the full
  ## matrix, including the rates to the genotypes we remove. After the
  ## removal those rates are gone from the row, but the old diagonal
  ## still counts them, and that is why probability is lost: we would
  ## have a chain that loses probability, not the intervention we want.
  ## The convention that is followed in the code is always the same:
  ## diagonals are 0 because they are -sum of the rest of the row,
  ## and thus they contain implicit info from the rest of the row; if
  ## the diagonal is not 0, it is easy to forget this fact; if it is 0 and
  ## we need the diagonal, it is easy to compute the sum on demand.
  ## So we test things are sane here (again ---this is tested elsewhere
  ## in the code too):
  stopifnot(isTRUE(all(diag(trm) == 0)))

  ## Genotype names are genes separated by ", " (e.g., "A, B, D"),
  ## and "WT" for the wildtype. We split each genotype name into its
  ## genes and look for an exact match of the gene. We do not search
  ## for the gene name inside the genotype name (e.g., with grep):
  ## gene "W" would then match "WT", and gene "G1" would match "G10".
  genes_in_genots <- stringi::stri_split_fixed(rownames(trm), ", ")
  rm_genots <- which(vapply(genes_in_genots,
                            function(z) gene %in% z,
                            logical(1)))
  if (length(rm_genots) == 0) {
    warning("No genotypes to remove?")
    trm2 <- trm
  } else {
    ## Should always be square
    ## If it ends up with 1 row and 1 column, and just a 0, it means only WT
    trm2 <- trm[-rm_genots, -rm_genots, drop = FALSE]
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
          tmp_hpd <- tmp$hitting_probs_from_WT_direct
          tmp_div <- tmp$divergence_hitting_prob_calculation
      } else {
          trm_after_intervention <- rm_genots_trm(cpm_output, gene, method)
          if (nrow(trm_after_intervention) > 1) {
              tmp_genot_freqs <-
                  genots_from_trm(trm_after_intervention, t = t)
              ## Compute the embedded chain from the transition rate
              ## matrix: divide each row by its sum. The diagonal
              ## should be zero, so each row sum is the sum of the
              ## rates to the other genotypes.
              off_q <- trm_after_intervention
              ## Check the diagonal is zero. This is dead code, as this
              ## would be caught earlier, in rm_genots_trm; left here
              ## to make it clear what we expect.
              stopifnot(isTRUE(all(diag(off_q) == 0)))
              rs <- rowSums(off_q)
              embedded <- off_q
              embedded[rs > 0, ] <- off_q[rs > 0, ] / rs[rs > 0]
              both <- hitting_probs_from_WT_both(embedded,
                                                 context = paste0(method,
                                                                  " trm_rm I:",
                                                                  gene))
              tmp_hp <- both$hp
              tmp_hpd <- both$hp_direct
              tmp_div <- both$divergence
          } else { ## All genots killed: only WT survives.
              ## Must be a 0.
              stopifnot(as.vector(trm_after_intervention) == 0)
              tmp_genot_freqs <- c(WT = 1)

              ## Hardcode WT = 1, to match the kill_gene path (the nrow ==
              ## 0 branch of get_genotype_freqs_cpm). On the trivial
              ## WT-only absorbing chain hitting_probs_from_WT_direct
              ## returns 0 while markovchain returns 1, so computing them
              ## here would both disagree with the kill_gene path AND
              ## record a spurious divergence of 1.
              tmp_hp <- c(WT = 1.0)
              tmp_hpd <- c(WT = 1.0)
              tmp_div <- 0.0
          }
      }
      return(list(genot_freqs = tmp_genot_freqs[tmp_genot_freqs > 0],
                  hitting_probs_from_WT = filter_hp_keep_wt(tmp_hp),
                  hitting_probs_from_WT_direct = filter_hp_keep_wt(tmp_hpd),
                  divergence_hitting_prob_calculation = tmp_div))
  }, mc.cores = mc.cores)
  names(intervene_all_genes) <- interventions
  names(intervene_all_genes)[-1] <- paste0("I:", interventions[-1])

  if (!is.na(filename)) saveRDS(intervene_all_genes, filename)
  return(intervene_all_genes)
}



library(codetools)
checkUsageEnv(env = .GlobalEnv)

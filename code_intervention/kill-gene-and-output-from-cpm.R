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


### What this code does:



## A. Given a CPM model, return that same model after making a gene lethal
##    Works for all of the models we use (OT, OncoBN, MHN, CBN, HESBCN)
##    This is done in several different ways (modifying DAG, setting
##    parameters to 0)
##
##    We also modify directly a fitness landscape by making
##    all genotypes with a given gene lethal.

## B. For a CPM model (original or modified after killing) of types MHN, CBN,
##   HESBCN, obtain the transition rates and, then, the predicted frequencies
##   of genotypes  at any arbitrary time.
##   For OT and OncoBN, there is no way to sample at
##    arbitrary times, so we return the predicted genot. freq. in
##    a "standard sample".

## The model is specified as:
##  - MHN: matrix of log-thetas
##  - remaining models: data frame as is standard in evamtools with
##     From, To, etc. columns

## See examples under "Examples of use" in file kill-gene-examples.R.

## The code has been tested (see section "Tests") but there might still
## be strange corner cases.


######################################################################
######################################################################

### Loading/sourcing dependencies


library(gtools)
library(igraph)
library(expm)
source("evam_v2.R")
source("trm.R")
source("HyperHMM-wrapper.R")

### Modifying evamtools code

## We need to modify a few evamtools functions for corner cases.
## This should probably eventually be included in evamtools.
## How we do it is from https://stackoverflow.com/a/58238931

## Give a named vector for the predicted freqs of genotypes
## changed here: add ", drop = FALSE". Oh man, I still trip on this ;-)
ev2_dist_oncotree_output_2_named_genotypes <- function(odt) {
  ## In the output of distribution.oncotree there
  ## is a column called Root, unlike OncoBN
  gpnroot <- which(colnames(odt) == "Root")
  gpnfr <- which(colnames(odt) == "Prob")
  gpn_names <- genot_matrix_2_vector(odt[, -c(gpnroot, gpnfr), drop = FALSE])
  odt <- as.vector(odt[, "Prob"])
  names(odt) <- gpn_names

  ## If with.errors = FALSE, there can be missing genotypes
  odt <- reorder_to_standard_order(odt)
  if (length(is.na(odt)))
    odt[is.na(odt)] <- 0
  return(odt)
}

environment(ev2_dist_oncotree_output_2_named_genotypes) <- asNamespace("evamtools")

assignInNamespace("dist_oncotree_output_2_named_genotypes",
                  ev2_dist_oncotree_output_2_named_genotypes,
                  ns = "evamtools")



## Give a named vector for the predicted freqs of genotypes
## changed here: add ", drop = FALSE".
ev2_DBN_est_genots_2_named_genotypes <- function(odt) {
  ## There is no column called Root, unlike OT
  gpnfr <- which(colnames(odt) == "Prob")
  gpn_names <- genot_matrix_2_vector(odt[, -gpnfr, drop = FALSE])
  odt <- as.vector(odt[, "Prob"])
  names(odt) <- gpn_names
  return(reorder_to_standard_order(odt))
}

environment(ev2_DBN_est_genots_2_named_genotypes) <- asNamespace("evamtools")

assignInNamespace("DBN_est_genots_2_named_genotypes",
                  ev2_DBN_est_genots_2_named_genotypes,
                  ns = "evamtools")



### Code



## fitness landscape, gene to make lethal -> intervened fitness landscape
##  I.e.: return new fitness landscape where we set to
##  "intervened_fitness" the fitness of all genotypes
##  that have the intervened gene mutated
##  Used to be called intervene_fitness_landscape_gene
kill_gene_fitness_landscape <- function(fitness_landscape, gene,
                                        intervened_fitness = 1e-9){
  if (!(gene %in% colnames(fitness_landscape)))
    stop("Gene not found in fitness landscape")
  genots_with_gene <- which(fitness_landscape[,gene] == 1)
  fitness_birth_column <- ifelse("Fitness" %in% colnames(fitness_landscape),
                                 "Fitness", "Birth")
  fitness_landscape[genots_with_gene, fitness_birth_column] <-
    intervened_fitness
  return(fitness_landscape)
}




## Given a model and a gene, return the model after
## making that gene lethal
## This is the modification of the DAG of restrictions
## for models with DAGs
## Model (in standard form: log-Theta matrix for MHN, data frame for the rest)
## and gene to be killed -> model after making the gene lethal.
## log-Theta matrix: in the output of evamtools this is MHN_theta (it is NOT
## MHN_exp_theta: MHN_exp_theta is the matrix that results from taking the exp of
## each element of MHN_theta). The values of MHN_theta can go from -infinity to
## +infinity. Those of exp_theta are (0, infinity) and are the actual hazards.

kill_gene <- function(x, gene, verbose = FALSE) {
  if (verbose) message("Standard killing: DAG/theta matrix removal")
  if (length(gene) != 1) stop("length(gene) != 1")
  gene <- gsub(" ", "", gene)

  ## For reasons I do not remember, I wrote all initial code
  ## without passing the method, so I inferred it from the type of output
  ## Becomes more complicated if we have both MHN and HyperHMM
  ## The first condition is an addition, a lot after the rest of the code
  ## was in place
  if (isTRUE(attributes(x)$method_output == "HyperHMM_trans_mat")) {
      kill_gene_HyperHMM(x, gene)
  } else if (is.matrix(x) &&
             (all(colnames(x) == rownames(x))) &&
      is.numeric(x)) {
      kill_gene_MHN(x, gene)
  } else if (is.data.frame(x) &&
             ("From" %in% colnames(x)) &&
             ("To" %in% colnames(x))) {
      kill_gene_DAG(x, gene)
  } else {
      stop("kill_gene called with unrecognized structure")
  }
}


## kill_gene for all models except MHN and HyperHMM.
## Actually, kill gene for OT, OncoBN, CBN, H-ESBCN
## Modify the DAG of restrictions removing the appropriate rows
## from the model
kill_gene_DAG <- function(x, gene) {
  if (length(gene) != 1) stop("length(gene) != 1")
  tmp <- suppressMessages(evamtools:::cpm2tm(list(edges = x))$weighted_fgraph)
  ## tmp <- suppressMessages(cpm2tm(list(edges = x))$weighted_fgraph)

  if (sum(colSums(tmp) == 0) > 1) { ## WT always has colSum = 0
    warning("weighted_fgraph contains unreachable destinations. ",
            "Removing those by decree")
    access_genots_before <- c("WT", colnames(tmp)[colSums(tmp > 0)])
  } else {
    access_genots_before <- colnames(tmp)
  }

  ## Genes that should disappear from models are those that only
  ## appear with the killed gene. Steps:

  ## 1. Find genotypes that disappear
  ## 2. Find genes that are always and only with those in 1.: other_genes
  ## 3. Remove edges that have any of c(gene, other_genes) as From or To

  ## 1. Find genotypes that disappear
  gk_st <- grep(paste0("^", gene, ", "), access_genots_before)
  gk_en <- grep(paste0(", ", gene, "$"), access_genots_before)
  gk_in <- grep(paste0(", ", gene, ","), access_genots_before, fixed = TRUE)
  gk_sg <- grep(paste0("^", gene, "$"), access_genots_before)
  gk <- unique(c(gk_st, gk_en, gk_in, gk_sg))

  if (length(gk) == 0) {
    warning("Removal of gene ", gene,
            "has no effect.")
    return(x)
  }

  genots_to_kill <- access_genots_before[gk]

  ## 2. Genes that are always and only with genotypes in 1: Genes that are
  ##    present in only removed genotypes.
  ##    These are genes that are not present in genotypes that remain after
  ##    killing, but were present originally.
  ##    Thus:
  ##      2.1 Obtain genes in genotypes that remain after killing
  ##      2.2 Obtain all genes in original set of genotypes
  ##      2.3 Diff of those genes


  ## 2.1 Genes that are present in at least one accessible genotype after
  access_genots_after_kill <- setdiff(access_genots_before, genots_to_kill)
  all_genes_after <- evamtools:::genes_in_genotypes(access_genots_after_kill)
  ## all_genes_after <- genes_in_genotypes(access_genots_after_kill)

  ## 2.2
  all_genes_original <- setdiff(unique(c(x$From, x$To)), "Root")

  ## 2.3
  genes_to_rm <- unique(setdiff(all_genes_original, all_genes_after))


  ## 3. rm genes that disappear
  if (length(genes_to_rm) == 0) stop("eh??!! No genes to rm??!!")
  rma <- which(x$From %in% genes_to_rm)
  rmb <- which(x$To %in% genes_to_rm)
  rm <- unique(c(rma, rmb))
  if (length(rm) == 0) stop("eh??!! No rows to rm??!!")

  xret <- x[-rm, , drop = FALSE]
  if (nrow(xret) == 0) warning("Model has 0 rows")
  return(xret)
}

## kill_gene for MHN
kill_gene_MHN <- function(x, gene) {
  if (length(gene) != 1) stop("length(gene) != 1")
  ## For MHN if you make a gene lethal the net result
  ## is making that gene disappear from the matrix of thetas
  ## as no genotype with that gene is viable.
  rc <- which(colnames(x) == gene)
  rr <- which(rownames(x) == gene)
  if (length(rc) == 0)
        stop("Gene ", gene, " is not part of the model.")
    if (length(rc) > 1)
        stop("Gene ", gene, " is repeated.")
    stopifnot(rc == rr)
  return(x[-rr, -rc, drop = FALSE])
}


## Model (in standard form: log-Theta matrix for MHN, data frame for the rest)
##     and gene to be killed -> model after making the gene lethal by setting
##     the parameter to 0
##     Generally only used for testing the equivalence of the
##     standard procedure above.
kill_gene_by_params_to_0 <- function(x, gene, verbose = TRUE) {
  if (verbose) message("Killing by setting parameters to 0")
  if (length(gene) != 1) stop("length(gene) != 1")
  gene <- gsub(" ", "", gene)
  if (is.matrix(x) &&
        (all(colnames(x) == rownames(x))) &&
        is.numeric(x)) {
    return(kill_gene_MHN_theta_minus_Inf(x, gene))
  } else if (is.data.frame(x) &&
               ("From" %in% colnames(x)) &&
               ("To" %in% colnames(x))) {
    return(kill_gene_DAG_param_0(x, gene))
  } else {
    stop("kill_gene called with unrecognized structure")
  }
}

## kill gene by setting the Theta of that gene = 0
##   Setting Theta directly has no effect, so we set theta to -Inf
kill_gene_MHN_theta_minus_Inf <- function(x, gene) {
  if (length(gene) != 1) stop("length(gene) != 1")
  x[gene, gene] <- -Inf
  x[, gene] <- -Inf
  ## This next one is not really necessary if x[gene, gene] is -Inf
  ## x[gene, ] <- -Inf
  return(x)
}


## kill gene by setting the parameter (lambda or cond. prob) of
##   that gene = 0
kill_gene_DAG_param_0 <- function(x, gene) {
  if (length(gene) != 1) stop("length(gene) != 1")

  ## Relation column: HESBCN or OncoBN; O.w. CBN or OT
  if ("Relation" %in% colnames(x)) {
    if ("theta" %in% colnames(x)) {
        ## method <- "OncoBN"
        set_to_zero <- which(x[, "To"] == gene)
        x[set_to_zero, "theta"] <- 0
    } else if ("Lambdas" %in% colnames(x))  {
        ## method <- "HESBCN"
        set_to_zero <- which(x[, "To"] == gene)
        x[set_to_zero, "Lambdas"] <- 0
    } else {
      stop("Model structure not recognized")
    }
  } else {
    if ("OT_edgeWeight" %in% colnames(x)) {
        ## method <- "OT"
        set_to_zero <- which(x[, "To"] == gene)
        x[set_to_zero, "OT_edgeWeight"] <- 0
    } else if ("rerun_lambda" %in% colnames(x)) {
        ## method <- "CBN"
        set_to_zero <- which(x[, "To"] == gene)
        x[set_to_zero, "rerun_lambda"] <- 0
    } else {
      stop("Model structure not recognized")
    }
  }
  return(x)
}


## Given a modified model provide the usual
## standard full output (transition matrices,
## predicted genotype frequencies, etc)
get_full_output <- function(x, epos = 0) {
    ## For OT and OncoBN, with epos = 0,
    ## we assume model is completely faithful.
    ## This is coherent with CBN and H-ESBCN.
    ## But we could use, for OncoBN, its estimated error
    ## which is deviations from the model, and for
    ## OT, epos (though this combines observation and model error).
    ## To remove the observation error,
    ## for OT we could set this simulation epos as
    ## epos <- max(0, epos - eneg)
    ## See suppl mat for EvAM-Tools,
    ## sections 5.2
    ## "Error models and obtaining finite samples (or
    ## sampled genotype  counts)" and section 2.3.2
    ## "OncoBN" and 2.3.6 "CPMs: Error models"
    if (isTRUE(attributes(x)$method_output == "HyperHMM_trans_mat")) {
        method <- "HyperHMM"
        out <- list()
        tmph <- probs_from_HyperHMM(x,
                                    attributes(x)$num_prob.set,
                                    attributes(x)$num_features)
        out$HyperHMM_trans_mat <- x
        out$HyperHMM_predicted_genotype_freqs <- tmph$predicted_genotype_freqs
        out$HyperHMM_conditional_genotype_freqs <- tmph$predicted_genotype_freq_at_t
    }  else if (is.matrix(x) &&
                (all(colnames(x) == rownames(x))) &&
                is.numeric(x)
                ) {
        method <- "MHN"
        out <- ev2_MHN_from_thetas_allow_neg_Inf(x)
    } else if (is.data.frame(x) &&
               ("From" %in% colnames(x)) &&
               ("To" %in% colnames(x))
               ) {
        ## Relation column: HESBCN or OncoBN; O.w. CBN or OT
        if ("Relation" %in% colnames(x)) {
            if ("theta" %in% colnames(x)) {
                method <- "OncoBN"
                out <- evamtools:::OncoBN_model_2_output(x, epos)
            } else {
                method <- "HESBCN"
                pset <- evamtools:::parent_set_from_edges(x)
                out <- evamtools:::HESBCN_model_2_output(x, pset)
            }
        } else {
            if ("OT_edgeWeight" %in% colnames(x)) {
                method <- "OT"
                out <- evamtools:::OT_model_2_output(x, epos)
            } else {
                method <- "CBN"
                out <- evamtools:::CBN_model_2_output(x)
            }
        }
    } else {
        stop("get_full_output called with unrecognized structure")
    }

    if (method %in% c("CBN", "MHN", "HESBCN")) {
        outname <- paste0(method, "_predicted_genotype_freqs")
        inname  <- paste0(method, "_trans_rate_mat")
        out[[outname]] <- evamtools:::probs_from_trm(out[[inname]])
    }
    trans_mat_key <- paste0(method, "_trans_mat")
    out[[paste0(method, "_hitting_probs_from_WT")]] <-
        hitting_probs_from_WT(out[[trans_mat_key]])
    return(out)
}


## model (possibly modified), optionally time ->
##   list(genot_freqs = named vector, hitting_probs_from_WT = named vector)
## For OT and OncoBN recall untimed.
## We deal with special cases (no model, single-row model)
## and for more general cases, call get_full_output once.
## We ASSUME epos = 0. This is done in the single-row and no-model cases.

## Except if you pass a t != NA, and use CBN/HESBCN/MHN, all the work
## is really done by get_full_output.
get_genotype_freqs_cpm <- function(model, t = NA) {
    if (nrow(model) == 0) {
        return(list(genot_freqs = c(WT = 1),
                    hitting_probs_from_WT = c(WT = 1.0)))
    }

    ## Find out the method.
    ## For reasons I do not remember, I wrote all initial code
    ## without passing the method, so I inferred it from the type of output.
    ## Becomes more complicated if we have both MHN and HyperHMM.
    ## The first condition is an addition, a lot after the rest of the code
    ## was in place.
    if (isTRUE(attributes(model)$method_output == "HyperHMM_trans_mat")) {
        method <- "HyperHMM"
    }  else if (is.matrix(model) &&
                (all(colnames(model) == rownames(model))) &&
               is.numeric(model)) {
        method <- "MHN"
    } else if (is.data.frame(model) &&
               ("From" %in% colnames(model)) &&
               ("To" %in% colnames(model))) {
        if ("Relation" %in% colnames(model)) {
            if ("theta" %in% colnames(model)) {
                method <- "OncoBN"
            } else if ("Lambdas" %in% colnames(model))  {
                method <- "HESBCN"
            } else {
                stop("Model structure not recognized")
            }
        } else {
            if ("OT_edgeWeight" %in% colnames(model)) {
                method <- "OT"
            } else if ("rerun_lambda" %in% colnames(model)) {
                method <- "CBN"
            } else {
                stop("Model structure not recognized")
            }
        }
    } else {
        stop("unrecognized structure")
    }

    if (!is.na(t) && (method %in% c("OT", "OncoBN", "HyperHMM"))) {
        warning("With methods ",
                "OT, OncoBN, and HyperHMM ",
                "get_genotype_freqs_cpm ignores the value of t.")
    }

    if (nrow(model) == 1 &&
        method %in% c("OT", "OncoBN")) {
        mut_gene <- model[1, 2]
        freq_mut <- model[1, ifelse(method == "OT", "OT_edgeWeight", "theta")]
        freqs <- c(1 - freq_mut, freq_mut)
        names(freqs) <- c("WT", mut_gene)
        ## From WT the only accessible state is mut_gene, so hitting prob = 1
        hp <- c(0.0, 1.0)
        names(hp) <- c("WT", mut_gene)
        return(list(genot_freqs = freqs, hitting_probs_from_WT = hp))
    }

    output <- get_full_output(model, epos = 0)
    hp <- output[[paste0(method, "_hitting_probs_from_WT")]]

    if (method %in% c("OT", "OncoBN")) {
        return(list(genot_freqs = output[[paste0(method, "_predicted_genotype_freqs")]],
                    hitting_probs_from_WT = hp))
    }

    if (method %in% c("HyperHMM")) {
        return(list(genot_freqs = output$HyperHMM_predicted_genotype_freqs,
                    hitting_probs_from_WT = hp))
    }

    ## CBN, HESBCN, MHN
    if (is.na(t)) {
        genot_freqs <- output[[paste0(method, "_predicted_genotype_freqs")]]
    } else {
        trans_name <- paste0(method, "_trans_rate_mat")
        trans_rate_mat <- output[[trans_name]]
        ## The predicted genotype frequencies are already part
        ## of the CBN, HESBCN, and MHN output object
        ## but, here, we allow ourselves to pass a different t.
        genot_freqs <- genots_from_trm(trans_rate_mat, t = t)
    }
    return(list(genot_freqs = genot_freqs, hitting_probs_from_WT = hp))
}


kill_gene_HyperHMM <- function(x, gene) {
    if (length(gene) != 1) stop("length(gene) != 1")
    ## 1. Find genotypes that are killed
    ## 2. Assign that probability to the diagonal
    ## 3. Zero the entries killed

    ## 1. Find genotypes
    genots_before <- colnames(x)
    gk_st <- grep(paste0("^", gene, ", "), genots_before)
    gk_en <- grep(paste0(", ", gene, "$"), genots_before)
    gk_in <- grep(paste0(", ", gene, ","), genots_before, fixed = TRUE)
    gk_sg <- grep(paste0("^", gene, "$"), genots_before)
    gk <- unique(c(gk_st, gk_en, gk_in, gk_sg))

    if (length(gk) == 0) {
        warning("Removal of gene ", gene,
                "has no effect.")
        return(x)
    }

    x1 <- x
    ## (part of 3. zeroing, and making life simpler for sanity check)
    x1[gk, ] <- 0
    ## x2 is for a different procedure, below, to check
    x2 <- x1

    ## Sanity check. Should only kill one per row,
    ## except when killed is the row (that was zeroed)
    to_zero <- apply(x1[, gk], 1, function(z) sum(z > 0))
    if (any(to_zero > 1)) {
        message("Killing more than one!")
        browser()
        stop("Killing more than one!")
    }

    ## Procedure A
    ## 2. Assign to diagonal
    prob_to_killed <- apply(x1[, gk], 1, function(z) sum(z))
    diag(x1) <- prob_to_killed
    ## 3. Zero the destination
    x1[, gk] <- 0

    paranoid_check <- TRUE
    if (paranoid_check) {
        ## Procedure B
        ## zero destination, get diagonal by difference
        x2[, gk] <- 0
        diag(x2) <- 1.0 - rowSums(x2)
        x2[gk, ] <- 0 ## they were assigned a 1 in diag
        stopifnot(isTRUE(all.equal(x1, x2)))
    }

    ## Set required attributes
    attr(x1, "method_output") <- "HyperHMM_trans_mat"
    attr(x1, "num_prob.set") <- attributes(x)$num_prob.set
    attr(x1, "num_features") <- attributes(x)$num_features

    return(x1)
}


library(codetools)
checkUsageEnv(env = .GlobalEnv)

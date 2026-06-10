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
## Functions to generate fitness landscapes
## (either specifying the fitness landscape model, like RMF or NK,
##  or using a CPM model ---CBN or H-ESBCN---) that fulfill
## requirements of gene and genotype frequencies.

#### Caveat: gene names
# This code only works with genes named from A-Z
# For example, if you call a gene "AB", this won't work.
# If we wanted to use custom named genes (from a dataset for example)
# we need to rename the genes to letters and then execute the code.


### Loading/sourcing dependencies

library(uuid)
library(evamtools)
library(OncoSimulR)
library(parallel)
source("trm.R")
source("evam_v2.R")



### Threads and OpenMP

Sys.getenv(c("OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "OPENBLAS_NUM_THREADS"))
## Not needed to control the fitness landscape generation
## Sys.setenv(OMP_NUM_THREADS = 1)
## Sys.setenv(OMP_THREAD_LIMIT = 1)
## Sys.setenv(OPENBLAS_NUM_THREADS = 1)
library(RhpcBLASctl)
message("omp_get_max_threads = ",
        RhpcBLASctl::omp_get_max_threads())

## RhpcBLASctl::blas_set_num_threads(1) ## Not needed
## The next is needed to prevent
## check_frequencies_of_genotypes_in_trm
## from using many threads
RhpcBLASctl::omp_set_num_threads(1)
Sys.getenv(c("OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "OPENBLAS_NUM_THREADS"))
message("omp_get_max_threads = ",
        RhpcBLASctl::omp_get_max_threads())


### Code

## number of genes, model, NK K order, hesbcn_relations ->
##                                             fitness landscape as matrix
##                                             (and, if applicable, model)
## Return fitness landscape of genotypes under different models
## with our default parameters making sure desired HESBCN relations
## are respected.  Calls, as needed,  rfitness or random_evam;
## in the second case, convert CPM model to fitness.
## For hesbcn_relations, two ways of specifying:
##    - as a vector like ("AND", "XOR")
##        will ensure at least some are AND and some are XOR
##    - as a vector like (AND = 3, XOR = 5)
##        will ensure that there are at least 3 edges (edges, not nodes)
##        involved in AND relationships and 5 in XOR.
## Used to be called get_fitness
generate_f_landscape <- function(n_genes,
                                 model = c("RMF", "NK",
                                           "CBN", "HESBCN",
                                           "OT", "OncoBN",
                                           "CBN-tree"),
                                 K = 2,
                                 hesbcn_relations = NA,
                                 hesbcn_max_iter = 5000,
                                 cbn_tree = TRUE) {
  model <- match.arg(model)
  if (model %in% c("OT", "OncoBN")) warning("We do not generate fitness ",
                                            "landscapes from OT/OncoBN ",
                                            "in the paper. ",
                                            "Be VERY CAREFUL if you use ",
                                            "them, as the scaling of ",
                                            "trans rate mat in CBN/HESBCN ",
                                            "need not hold here.")
  if (n_genes <= 1)
    stop("The number of genes should be greater than 1.")

  if (model %in% c("RMF", "NK")) {
    max_f_rmf_nk <- (1 + 0.01)^(round(n_genes/2) + 1)
    if (model == "RMF") {
      ## Numbers not identical to paper with Vasallo
      ## We use c centered around 0.5, sd = 1, mu = 1
      ## and then scale.

      ## For RMF: Paper with C.Vasallo
      ## 1.The reference genotype (i.e., the genotype with maximum fitness)
      ## was randomly chosen (setting reference = ’random’ in the rfitness
      ## function in OncoSimulR).
      ## 2.Sd, of the random normal variate was set to 0.2
      ## 3.The decrease in fitness (strictly, birth rate) of a genotype
      ## per each unit increase in Hamming distance from the reference
      ## genotype, c, was chosen from a uniform U(0, 0.2) distribution.
      ## BEWARE: the "c" above is not the same c of the scaling of the
      ## transition rate matrix.

      random_decrease_fitness <- runif(1, 0.1, 0.8)
      rfo <- rfitness(n_genes, model = "RMF",
                      c = random_decrease_fitness,
                      sd = 1, reference = "random",
                      scale = c(max_f_rmf_nk, 1e-9, 1))
    } else if (model == "NK") {
      s <- format(round(runif(1, 1, 2^40)), scientific = FALSE)
      rfo <- rfitness(n_genes, model = "NK", K = K,
                      scale = c(max_f_rmf_nk, 1e-9, 1),
                      seed_magellan = s)
    }
    colnames(rfo)[ncol(rfo)] <- "Fitness"
    ## The next is mostly useless here, since we do not have 0s
    rfo <- rfo[rfo[, "Fitness"] > 0, ]
    ## If we remove genotypes, Magellan_stats breaks
    ## if class is genotype_fitness_matrix
    ## as we skip a transformation to fitness landscape
    ## with all genotypes
    class(rfo) <- c("matrix", "array")
    return(list(fitness_landscape = rfo,
                n_genes = n_genes,
                model = model,
                K = ifelse(model == "NK", K, NA),
                other = NA,
                graph_density = NA,
                a_cpm_2_si = NA,
                hesbcn_relations = NA))
  }
  ## Get random model, then transform it into a OncoSimulR format
  ## so we can get the fitness of each genotype

  if (model %in% c("HESBCN", "CBN", "CBN-tree")) {
    ## An alternative would have been to do
    ## a <- 0.01/mean(lamba) so made indep of the U(low, up) range
    a_cpm_2_si <- 0.006 ## 1.7 [U(1/3, 3)] * a : approx. 0.01
  } else if (model %in% c("OT", "OncoBN") ) {
    a_cpm_2_si <- 0.02  ## 0.5 * 0.02 = 0.01; 0.5 mean of U(0, 1).
  } else if (model == "MHN") {
    stop("We do not know what to do with MHN here")
  }

  if (model == "HESBCN") {
    if (any(is.na(hesbcn_relations)))
      stop("Enter a valid value for hesbcn_relations")

    if (!is.numeric(hesbcn_relations)) {
      ## If hesbcn_relations are just strings
      ## we just want some of each, period.
      min_edges <- rep(0, 4)
      names(min_edges) <- c("AND", "OR", "XOR", "Single")
      min_edges[hesbcn_relations] <- 1
    } else {
      min_edges <- rep(0, 4)
      names(min_edges) <- c("AND", "OR", "XOR", "Single")
      min_edges[names(hesbcn_relations)] <- hesbcn_relations
    }
    ## This will just make it easier to get the required numbers
    hesbcn_relations_probs <-
      min_edges[c("AND", "OR", "XOR")]/sum(min_edges[c("AND", "OR", "XOR")])
    hesbcn_iter <- 0
    hesbcn_failed <- FALSE
    while (TRUE) {
      hesbcn_iter <- hesbcn_iter + 1
      graph_density <- runif(1, 0.1, 0.5)
      message(" ... generating suitable HESBCN model.  graph_density = ",
              graph_density)
      random_model <- random_evam(n_genes, model = model,
                                  hesbcn_probs = hesbcn_relations_probs,
                                  graph_density = graph_density,
                                  ## cbn_tree here is irrelevant
                                  ## as model = HESBCN
                                  cbn_tree = FALSE)
      ## Make sure at least those many edges
      tt <- table(random_model[["HESBCN_model"]]$Relation)
      tt <- pmax(tt[c("AND", "OR", "XOR", "Single")], 0, na.rm = TRUE)
      names(tt) <- c("AND", "OR", "XOR", "Single")
      if (all(tt >= min_edges)) {
        break
      } else {
        message("         edges condition not fulfilled. ",
                "hesbcn_iter = ", hesbcn_iter)
      }
      if (hesbcn_iter >= hesbcn_max_iter) {
        hesbcn_failed <- TRUE
        warning("Could not generate HESBCN with requested ",
                "relationships. Expect failure if running ",
                "generate_n_f_landscape_requir")
        return(list(fitness_landscape = NA,
                    n_genes = n_genes,
                    model = model,
                    K = ifelse(model == "NK", K, NA),
                    other = random_model,
                    graph_density = graph_density,
                    a_cpm_2_si = a_cpm_2_si,
                    hesbcn_relations = ifelse(model == "HESBCN",
                                              hesbcn_relations, NA),
                    hesbcn_iter = ifelse(model == "HESBCN",
                                         hesbcn_iter, NA),
                    hesbcn_failed = ifelse(model == "HESBCN",
                                           hesbcn_failed, NA)
                    )
               )
      }
    }
  } else { ## We only end up here if model was not H-ESBCN
      ## nor RMF nor NK
      graph_density <- runif(1, 0.1, 0.5)
      random_model <- random_evam(n_genes, model = model,
                                  graph_density = graph_density,
                                  ot_oncobn_epos = 0,
                                  cbn_tree = cbn_tree)
  }

  genot_fitness <-
    ev2_cpm_to_fitness_genots(random_model[[paste0(model, "_model")]],
                              a = a_cpm_2_si)
  rfo <- cbind(genots_to_bin(genot_fitness$Genotype, n_genes),
               Fitness = genot_fitness$Birth)
  rfo <- rfo[rfo[, "Fitness"] > 0, ]
  ## This shouldn't have been here, if we remove genotypes
  ## as Magellan_stats breaks
  ##  class(rfo) <- c("matrix", "array", "genotype_fitness_matrix")
  ## The exp fitness were created in a way that they keep
  ## the genotype_fitness_matrix class. So Magellan would break
  ## but it does not, because of our wrapper complete_fitness_landscape
  class(rfo) <- c("matrix", "array")
  return(list(fitness_landscape = rfo,
              n_genes = n_genes,
              model = model,
              K = ifelse(model == "NK", K, NA),
              other = random_model,
              graph_density = graph_density,
              a_cpm_2_si = a_cpm_2_si,
              hesbcn_relations = ifelse(model == "HESBCN",
                                        hesbcn_relations, NA),
              hesbcn_iter = ifelse(model == "HESBCN",
                                   hesbcn_iter, NA),
              hesbcn_failed = ifelse(model == "HESBCN",
                                     hesbcn_failed, NA),
              cbn_tree = ifelse(model == "CBN", cbn_tree, NA)
              )
         )
}


## Number of fitness landscapes, number of genes, model, maximum attempts,
##     NK K order -> fitness landscape, scaling, transition rate matrix

## Get multiple fitness landscapes (and value of c scaling of transition
## rate matrix) that fulfill the requirements of frequencies of genots.
## and genes.
## Used to be called get_n_fitness
generate_n_f_landscape_requir <- function(n, n_genes,
                                          model = c("RMF", "NK",
                                                    "CBN", "HESBCN",
                                                    "OT",
                                                    "OncoBN"),
                                          K = 2,
                                          hesbcn_relations = NA,
                                          max_iter = 5000,
                                          hesbcn_max_iter = max_iter,
                                          cbn_tree = TRUE,
                                          cores = detectCores(),
                                          mc.preschedule = FALSE,
                                          custom_sampling = FALSE,
                                          custom_label = NULL) {
    model <- match.arg(model)
    if (model %in% c("OT", "OncoBN")) warning("We do not generate fitness ",
                                              "landscapes from OT/OncoBN ",
                                              "in the paper.",
                                              "Be VERY CAREFUL if you use ",
                                              "them, as the scaling of ",
                                              "trans rate mat in CBN/HESBCN ",
                                              "need not hold here."
                                              )

    out <- mclapply(1:n, function(x) {

        ## generate an id string of the iteration
        text_id <- function() {
            tmp <- paste(model,
                         ifelse(model == "NK", paste0("K=", K), ""),
                         ifelse(model == "HESBCN",
                                paste0("relations_",
                                       paste(names(hesbcn_relations),
                                             hesbcn_relations,
                                             sep = "=", collapse = ",")),
                                ""),
                         paste0("n_genes=", n_genes),
                   paste0(custom_label),
                         paste0("rep.", x),
                         sep = "_"
                         )
            ## Yes, sometimes up to three
            tmp <- gsub("____", "_", tmp, fixed = TRUE)
            tmp <- gsub("___", "_", tmp, fixed = TRUE)
            return(gsub("__", "_", tmp, fixed = TRUE))
        }

        iter <- 0
        while (TRUE) {
            iter <- iter + 1
            message("  ... iter = ", iter)
            ogf <- generate_f_landscape(n_genes = n_genes,
                                        model = model,
                                        K = K,
                                        hesbcn_relations = hesbcn_relations,
                                        hesbcn_max_iter = hesbcn_max_iter)
            fit_landscape <- ogf$fitness_landscape
            if (model %in% c("CBN", "HESBCN", "CBN-tree")) {
                trm_and_c <-
                    suppressMessages(get_scaled_trm_adaptive(fit_landscape,
                                                             c = 1/ogf$a_cpm_2_si))
            } else {
                trm_and_c <-
                    suppressMessages(get_scaled_trm_adaptive(fit_landscape))
            }
            if (!all(is.na(trm_and_c$trm_scaled))) {
                gene_genots_checks <-
                    check_frequencies_of_genotypes_in_trm(trm_and_c$trm_scaled,
                                                          LETTERS[1:n_genes],
                                                          custom_sampling = custom_sampling)
                if (gene_genots_checks[["test_OK"]]) break
            } else {
                gene_genots_checks <- NA
            }
            if (iter >= max_iter) {
                warning(paste("Not able to generate the fitness landscape",
                              " under specified model and parameters"))
                return(list(
                    fitness_landscape = NA,
                    c = NA, ## trm_and_c$c,
                    trm_scaled = NA, ## trm_and_c$trm_scaled,
                    actual_mean_rate = trm_and_c$actual_mean_rate,
                    rand_target_mean_rate = trm_and_c$rand_target_mean_rate,
                    rel_diff = NA, ## trm_and_c$rel_diff,
                    si_stats = trm_and_c$si_stats,
                    gene_genot_checks = gene_genots_checks[-1],
                    iter = iter,
                    n_genes = ogf$n_genes,
                    model = ogf$model,
                    other = NA,
                    graph_density = NA,
                    K = ogf$K,
                    hesbcn_relations = ogf$hesbcn_relations,
                    hesbcn_iter = ogf$hesbcn_iter,
                    hesbcn_failed = ogf$hesbcn_failed,
                    failed_status = TRUE,
                    text_id = text_id(),
                    uuid = NA,
                    custom_sampling = custom_sampling,
                    genots_freq_from_gene_genots_check = NA
                ))
            }
        }

        return(list(fitness_landscape = fit_landscape,
                    c = trm_and_c$c,
                    trm_scaled = trm_and_c$trm_scaled,
                    actual_mean_rate = trm_and_c$actual_mean_rate,
                    rand_target_mean_rate = trm_and_c$rand_target_mean_rate,
                    rel_diff = trm_and_c$rel_diff,
                    si_stats = trm_and_c$si_stats,
                    gene_genot_checks = gene_genots_checks[-1],
                    iter = iter,
                    n_genes = ogf$n_genes,
                    model = ogf$model,
                    other = ogf$other,
                    graph_density = ogf$graph_density,
                    K = ogf$K,
                    hesbcn_relations = ogf$hesbcn_relations,
                    hesbcn_iter = ogf$hesbcn_iter,
                    hesbcn_failed = ogf$hesbcn_failed,
                    failed_status = FALSE,
                    text_id = text_id(),
                    uuid = UUIDgenerate(),
                    custom_sampling = custom_sampling,
                    ## Yes, already in gene_genots_check
                    genots_freq_from_gene_genots_check = gene_genots_checks$genots_freq_from_check))
    },
    mc.cores = cores,
    mc.preschedule = mc.preschedule)

    return(out)
}


## A simple way to see what is going on
## u <- get_n_fitness(1, 7, "RMF"); u[[1]]$c; u[[1]]$rand_target_mean_rate; plotFitnessLandscape(u[[1]]$fitness_landscape, only_accessible = TRUE)


## list of models and parameters, file name ->
##     list with all fitness landscapes, scaling, transition rate
##     written to file and returned
## Used to be called run_models_interv
write_f_landscapes <- function(model_list, file_name) {
  all_fitnesses <- lapply(model_list, function(model) {
    cat("\n Generating for model ",
        paste(unlist(model), collapse = ", "))
    do.call("generate_n_f_landscape_requir", model)
  })
  retobj <- list(all_fitnesses = all_fitnesses,
                 model_list = model_list,
                 file_name = file_name)
  saveRDS(retobj, file = file_name)
  return(retobj)
}


### example call:


## MODELS_INTERV_CPM <- list(
##     "NK_100-replic_7-genes_K1" = list(n = 100, n_genes = 7,
##                                       model = "NK", K = 1),
##     "NK_100-replic_7-genes_K2" = list(n = 100, n_genes = 7,
##                                       model = "NK", K = 2),
##     "NK_100-replic_7-genes_K3" = list(n = 100, n_genes = 7,
##                                       model = "NK", K = 3),
##     "RMF_200-replic_7-genes" = list(n = 200, n_genes = 7,
##                                     model = "RMF"),
##     "CBN_100-replic_7-genes" = list(n = 100, n_genes = 7,
##                                     model = "CBN"),
##     "HESBCN_50-replic_7-genes_AND-OR-XOR" = list(n = 50, n_genes = 7,
##                                                  model = "HESBCN",
##                                                  hesbcn_relations =
##                                                      c("AND", "OR", "XOR")),
##     "HESBCN_50-replic_7-genes_AND-OR" = list(n = 50, n_genes = 7,
##                                              model = "HESBCN",
##                                              hesbcn_relations =
##                                                  c("AND", "OR")),
##     "HESBCN_50-replic_7-genes_AND-XOR" = list(n = 50, n_genes = 7,
##                                               model = "HESBCN",
##                                               hesbcn_relations =
##                                                   c("AND", "XOR"))
## )

## write_f_landscapes(MODELS_INTERV_CPM,
##                    "fitness_landscapes.RDS")


library(codetools)
checkUsageEnv(env = .GlobalEnv)

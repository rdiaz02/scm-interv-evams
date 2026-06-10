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
## along with this program.  If not, see <http://www.gnu.org/licenses/>.


## These are functions we lift from evamtools but we modify the one
## we call directly, genots_2_fgraph_and_trans_mat_rf.
## To make sure we are calling what we want, we prepend "no_evam"
## to the function names.



library(OncoSimulR)
library(igraph)
library(Matrix)


######################################################################
######################################################################


## output of rfitness to data.frame or named vector with genotypes as
## letters
no_evam_rfitness_to_letter <- function(y, output = "vector") {
    col_fitness <- ncol(y)
    genotlabel <- function(u, nn = colnames(y)[-col_fitness]) {
        if (any(is.na(u)))
            return(NA)
        else {
            return(paste(sort(nn[as.logical(u)]), collapse = ", "))
        }
    }
    Genotype <- apply(y[, -col_fitness], 1, genotlabel)
    Genotype[Genotype == ""] <- "WT"
    if (output == "vector") {
        out <- y[, col_fitness]
        names(out) <- Genotype
    } else {
        out <- data.frame(Genotype = Genotype,
                     Fitness = y[, col_fitness])
    }
    return(out)
}


## From the function of the same name in evamtools (file evam-wrapper.R)

## list of accessible genotypes -> adjacency matrix of genotypes (fitness graph)
## return maximally connected fitness graph for a given set of accessible
## genotypes.   This list contains no WT (we add it)
## BEWARE! This is the maximally connected fitness graph. This works with
## CPMs. But this is wrong, if, say, fitnesses are: A = 2, B = 3, AB = 2.5.
## This will place an arrow between B and AB, but there should  be no such edge.
## function genots_2_fgraph_and_trans_mat
## only returns the truly accessible

## In gacc, each genotype is a vector of genes, and the vector is already
## sorted. This ensures that the naming of genotypes is consistent later.

no_evam_unrestricted_fitness_graph_sparseM <- function(gacc) {
  gs <- unlist(lapply(gacc, function(g) paste0(g, collapse = ", ")))
  gs <- c("WT", gs)
  nmut <- c(0, vapply(gacc, length, 1))
  gs <- gs[order(nmut, gs)]

  jj <- match(gs[which(nmut == 1)], gs)
  ii <- rep.int(match("WT", gs), length(jj))
  adjmat <- sparseMatrix(i = ii, j = jj, x = 1L,
                         dims = c(length(gs), length(gs)),
                         dimnames = list(gs, gs))

  for(m in 2:max(nmut)){
    g <- gs[which(nmut == m)]
    for (gn in g) {
      parents <- gacc[which(nmut == m-1)-1]
      gns <- unlist(strsplit(gn, ", "))
      parents <-
        parents[which(unlist(lapply(parents,
                                    function(p)
                                      length(setdiff(gns, p)))) == 1)]
      for (p in parents){
        ## Works but better via indices, I think
        ## adjmat[paste0(p, collapse = ", "), gn] <- 1L
        jjj <- match(gn, gs)
        iii <- rep.int(match(paste0(p, collapse = ", "), gs),
                       length(jjj))
        adjmat[iii, jjj] <- 1L
      }
    }
  }
  return(adjmat)
}


## from genots_2_fgraph_and_trans_mat in evamtools (file
## access_genots_from_oncosimul.R) but modified to give the relative fitness
## differences and, thus, the "s_i"

## named vector of genotype fitness -> fitness graph and transition matrix and
##                                     accessible genotypes
##                                   Only accessible genotypes shown in output.
##   We assume WT is 1 if not given explicitly
no_evam_genots_2_fgraph_and_trans_mat_rf <- function(x) {
  ## Logic:
  ##  - Construct a fully connected fitness graph
  ##    between the given genotypes. So any genotype connected to
  ##    genotypes with one extra mutation.
  ##    This matrix might contain genotypes that are not truly accessible.
  ##  - Construct matrix of fitness differences between ancestor and immediate
  ##    descendants. Likely slow if many genotypes.
  ##  - Set to non accessible if fitness difference <= 0

  ##  Could be done faster, by not creating the unrestricted fitness graph
  ##  and instead maybe using OncoSimulR's wrap_accessibleGenotypes
  ##  But would need to check that works with partial lists of genotypes
  ##  and use allGenotypes_to_matrix. And would need to change
  ##  OncoSimulR's wrap_accessibleGenotype and how it uses the th.

  ## We use this approach, to minimize the number of genotypes we call
  ## no_evam_unrestricted_fitness_graph_sparseM in cpm_to_trans_mat_oncosimul

  #Stop if any fitness is under or equal to 0
  if(any(x <= 0)) stop("All fitness values must be above 0.")

  which_wt <- which(names(x) == "WT")
  if (length(which_wt) == 1) {
    fit_wt <- x["WT"]
    if (fit_wt != 1.0) message("Your WT has fitness different from 1.",
                               " Using WT with the fitness you provided.")
    x <- x[-which_wt]
  } else {
    fit_wt <- c("WT" = 1.0)
  }


  ## Silly? Inside the next, we now put them together. FIXME?
  access_genots_as_list <- lapply(names(x),
                                  function(v) strsplit(v, ", ")[[1]])
  ## For ordered output
  no <-  order(unlist(lapply(access_genots_as_list, length)), names(x))
  access_genots_as_list <- access_genots_as_list[no]

  fgraph <- no_evam_unrestricted_fitness_graph_sparseM(access_genots_as_list)

  genots_fitness <- c(fit_wt, x)[colnames(fgraph)]
  mf <- matrix(rep(genots_fitness, nrow(fgraph)),
               nrow = nrow(fgraph), byrow = TRUE)
  stopifnot(identical(dim(mf), dim(fgraph)))

  fdiff <- mf - genots_fitness
  fdiff <- fgraph * fdiff

  fgraph[fdiff <= 0] <- 0

  tm <- fdiff
  tm[tm < 0] <- 0
  tm <- tm / ifelse(rowSums(tm) != 0, rowSums(tm), 1)
  ## This we know is always 0
  tm[nrow(tm), ] <- 0
  ## This we can set
  tm[rowSums(fgraph) == 0, ] <- 0

  ## First simple filtering for potentially expensive call to igraph
  accessible_genotypes_candidates <-
    genots_fitness[colnames(fgraph)[colSums(fgraph) >= 1]]

  ig_fgraph <- igraph::graph_from_adjacency_matrix(fgraph)

  num_paths_from_WT <-
    vapply(names(accessible_genotypes_candidates),
           function(x)
             length(igraph::all_simple_paths(ig_fgraph,
                                             from = "WT",
                                             to = x, mode = "out")),
           FUN.VALUE = 0)

  accessible_genotypes <-
    accessible_genotypes_candidates[num_paths_from_WT > 0]

  if (length(accessible_genotypes) == 0) {
    message("No accessible genotypes")
    return(
      list(fitness_graph = NA,
           transition_matrix = NA,
           fitness_differences = NA,
           relative_fitness_differences = NA,
           accessible_genotypes = accessible_genotypes))
  }

  ## Remove unaccessibe genotypes from matrices before returning
  col_ret <- c("WT", names(accessible_genotypes))
  ## And remove rows  that are only destinations
  row_ret <-
    setdiff(col_ret,
            names(which(rowSums(fgraph[col_ret, col_ret]) == 0)))
  fdiff[fdiff <= 0] <- 0

  return(
    list(fitness_graph = fgraph[row_ret, col_ret, drop = FALSE],
         transition_matrix = tm[row_ret, col_ret, drop = FALSE],
         fitness_differences = fdiff[row_ret, col_ret, drop = FALSE],
         relative_fitness_differences =
           fdiff[row_ret, col_ret, drop = FALSE] / genots_fitness[row_ret],
         ## sweep does not leave the 0 as sparse
         ## sweep(fdiff[row_ret, col_ret, drop = FALSE], 1, genots_fitness[row_ret], "/")
         accessible_genotypes = accessible_genotypes))
}

library(codetools)
checkUsageEnv(env = .GlobalEnv)

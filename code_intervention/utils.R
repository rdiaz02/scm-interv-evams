library(Matrix)
library(evamtools)

## vector of genotypes (strings), number of genes -> binary matrix of genotypes
## n must be at least as large as the index of the largest letter
## (i.e., if you pass "D, F, G", use at least n of 7)
genots_to_bin <- function(genots, n) {
    bin <- vapply(genots, evamtools:::str2binary, double(n), n = n)
    out <- if(is.vector(bin)) as.matrix(bin)
    else t(bin)
  out <- unname(out)
  colnames(out) <- LETTERS[1:n]
  out
}

## Converts matrix into a square.
## It does an outer join between the row-names and col-names.
## Fills missing data with fill_data
get_square_matrix <- function(mat, fill_data = 0){
    dim_names <- unique(c(colnames(mat), rownames(mat)))
    sq_matrix <- Matrix(fill_data, length(dim_names), length(dim_names),
                        dimnames = list(dim_names, dim_names), sparse = TRUE)
    sq_matrix[rownames(mat), colnames(mat)] <- mat
    sq_matrix
}


## All genotypes
## From code in theta_to_trans_rate_3 in evamtools
allGenotypesLetter <- function(g) {
    genots <- evamtools:::allGenotypes_3(g)
    geneNames <- LETTERS[1:g]
    genotNames <- unlist(
        lapply(genots$bin_genotype,
               function(x)
                   paste(geneNames[which(x == 1L)], sep = "", collapse = ", "))
    )
    genotNames[genotNames == ""] <- "WT"
    return(genotNames)
}


## Fitness-landscape statistics helpers (complete_fitness_landscape,
## magellan_stats_masked, fitness_landscape_stats, one_fl_stats) were
## moved to code_evam_simul_interv/utils-flandscape-stats.R during the
## 2026-05 split between code_intervention (pure intervention machinery)
## and code_evam_simul_interv (simulation, statistics, downstream analyses).


library(codetools)
checkUsageEnv(env = .GlobalEnv)

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


## This code will eventually end up in evamtools.
## It is already in branch future-3.3.5 and has been tested with
## the rest of evamtools.

## For now, this is self-contained. The only real change in is in
## cpm_out_to_oncosimul, where we have a scaling factor with default
## of 10⁻3 that multiplies lambdas and pi

suppressMessages(library(OncoSimulR, quietly = TRUE, verbose = FALSE))
library(evamtools)


## output from CPMs, sh if restrictions not satisfied ->
##                            input for OncoSimulR evalAllGenotypes
##  a: factor that multiplies lambda or p to give the s_i,
##     selection coefficient. Generally 0.006 for CBN and HESBCN,
##     0.02 for OT and OncoBN
ev2_cpm_out_to_oncosimul <- function(x, a = 6e-3,
                                     sh = -Inf) {

  if ("rerun_lambda" %in% names(x)) { ## CBN
    s <- x$rerun_lambda
    typeDep <- "AND"
  } else if ("lambda" %in% names(x)) { ## MCCBN-HCBN
    s <- x$lambda
    typeDep <- "AND"
  } else if ("Relation" %in% names(x)) { ## HESBCN (same thing as PMCE)
    ## Also using this for DBN, as it could return an AND
    if ("Lambdas" %in% names(x) ) ## HESBCN
      s <- x$Lambdas
    if ("theta" %in% names(x) ) ## DBN
      s <- x$theta
    typeDep <- x$Relation
    typeDep[typeDep == "Single"] <- "AND"
  } else if ("OT_edgeWeight" %in% names(x)) { ## OT
    s <- x$OT_edgeWeight
    typeDep <- "AND"
  } else {
    stop("Input not recognized")
  }
  ## To get the colors right
  ## Fix this in OncoSimulR
  typeDep[typeDep == "AND"] <- "MN"
  typeDep[typeDep == "OR"] <- "SM"
  typeDep[typeDep == "XOR"] <- "XMPN"
  x1 <- data.frame(parent = x$From,
                   child  = x$To,
                   s = a * s,
                   sh = sh,
                   typeDep = typeDep
                   )
  return(x1)
}


## fitness, target max fitness. WT fitness always 1.
ev2_scale_fitness_2 <- function(x, max_f) {
  max_x <- max(x)
  if (max_x > 1e10) {
      warning("Maximum fitness > 1e10. Expect numerical problems.")
  }
  return(1.0 +  (x - 1) * ((max_f - 1) / (max_x - 1)))
}


## output from CPMs, multiplying constant to obtain s_i,
##         max final fitness, sh when restrictions not stasified,
##         max num genots -> fitness of all genotypes
##    if max_f is NULL: no rescaling of fitness
##       max_f should have no effect in probs transition
##    max_genots: argument max of evalAllGenotypes
##    a: generally 0.006 for CBN and HESBCN, 0.02 for OT and OncoBN
##       (see )
ev2_cpm_to_fitness_genots <- function(x,
                                      a,
                                      max_f = NULL,
                                      sh = -Inf, max_genots = 2^15) {
  x1 <- ev2_cpm_out_to_oncosimul(x, a = a, sh = sh)
  x1 <- OncoSimulR::evalAllGenotypes(fitnessEffects = OncoSimulR::allFitnessEffects(rT = x1),
                                     addwt = TRUE, max = max_genots)

  ## In newer OncoSimulR, column names for Fitness can now be called Birth
  fitness_birth_column <- ifelse("Fitness" %in% colnames(x1),
                                 "Fitness", "Birth")
  if (!is.null(max_f)) {
    if (max_f < 1) stop("max_f must be larger than min_f")

    if (fitness_birth_column == "Fitness") {
      x1$Fitness[x1$Fitness > 0.0] <-
        ev2_scale_fitness_2(x1$Fitness[x1$Fitness > 0.0], max_f)
    } else if (fitness_birth_column == "Birth") {
      x1$Birth[x1$Birth > 0.0] <-
        ev2_scale_fitness_2(x1$Birth[x1$Birth > 0.0], max_f)
    } else {
      stop("The column should be called Birth or Fitness")
    }
  }
  return(x1)
}


cpm_to_trans_mat_oncosimul <- function(...) {
  message("This is a tripwire to ensure we are not calling a function ",
          "from evamtools that calls another cpm_to_fitness_genots")
  stop("DO NOT CALL THIS")
}

cpm2F2tm <- function(...) {
  message("This is a tripwire to ensure we are not calling a function ",
          "from evamtools that calls another cpm_to_fitness_genots")
  stop("DO NOT CALL THIS")
}


## Like evamtools::MHN_from_thetas, but skip the
## simple paranoid check of transition matrix if any theta is -Inf
## We could try a more sophisticated paranoid check
## but it would be circular (identify genotypes that never
## transition)?
ev2_MHN_from_thetas_allow_neg_Inf <- function(thetas) {
  ## A hack to make this function callable from others below
  inner_transitionRate_3_1 <- evamtools:::inner_transitionRate_3_1
  oindex <- evamtools:::evam_string_order(colnames(thetas))
  thetas <- thetas[oindex, oindex, drop = FALSE]
  output <- list()
  output[["MHN_theta"]] <- thetas
  paranoidCheck <- ifelse(any(thetas == -Inf), FALSE, TRUE)
  output[["MHN_trans_rate_mat"]] <-
    evamtools:::theta_to_trans_rate_3_SM(thetas,
                                         inner_transition = inner_transitionRate_3_1)
  output[["MHN_trans_mat"]] <-
    evamtools:::trans_rate_to_trans_mat(output[["MHN_trans_rate_mat"]],
                                        method = "competingExponentials",
                            paranoidCheck = paranoidCheck)
  output[["MHN_td_trans_mat"]] <-
    evamtools:::trans_rate_to_trans_mat(output[["MHN_trans_rate_mat"]],
                                        method = "uniformization",
                            paranoidCheck = paranoidCheck)
  output[["MHN_exp_theta"]] <- exp(thetas)
  return(output)
}




#### Modifying evamtools:::adjm_rm_no_access

##  adjm_rm_no_access is called from cpm2tm with OncoBN models
##  (cpm2tm is called from OncoBN_model_2_output, all in evamtools
##  and OncoBN_model_2_output is called from get_full_output).
##  adjm_rm_no_access does not expect to see a model where
##  no events can take place (because the first gene depends
##  on root, everything depends on the first gene, and
##  the prob. of that first gene is 0).
##  The original code would break in this case.
##  Note that this does NOT affect the standard killing procedure
##  where we manually deal with this case (models with 0 rows).
##  (see function get_genotype_freqs_cpm). And there are tests
##  for this case.

## From https://stackoverflow.com/a/58238931

ev2_adjm_rm_no_access <- function(x) {
    while (TRUE) {
        nacc <- which(colSums(x) == 0)
        wwt <- which(colnames(x) == "WT")
        nacc <- setdiff(nacc, wwt)
        if (length(nacc)) {
            x <- x[-nacc, -nacc, drop = FALSE]
        } else {
            break
        }
    }
    message("In ev2_adjm_rm_no_access")
    return(x)
}

environment(ev2_adjm_rm_no_access) <- asNamespace("evamtools")

assignInNamespace("adjm_rm_no_access",
                  ev2_adjm_rm_no_access,
                  ns = "evamtools")




## FIXME: more possible additions to evamtools in
## kill-gene-and-output-from-cpm.R
## the two ev2_ functions

## FIXME: possible additional changes in evamtools?
## Do, for all models, as we do for OncoBN in cpm2tm?
## So filter the weigthed_fgraph by adjm_rm_no_access?
## possibly using the new ev2_adjm_rm_no_access?
## Motivation: be consistent and do not have
## that annoying warning?
## See examples in
## "Explaining the warning about unreachable destinations"
## in kill-gene-equivalences-TESTS.R


## FIXME: TODO
## If we settle this stuff, then change the evamtools code.

library(codetools)
checkUsageEnv(env = .GlobalEnv)

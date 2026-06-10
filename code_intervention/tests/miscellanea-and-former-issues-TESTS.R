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

## options(intervention_every_gene_cores = parallel::detectCores())
library(testthat)
pwd <- getwd()
setwd("../")
source("intervention.R")
setwd(pwd)

set.seed(NULL)

### Examples of use
##  Though the tests provide many additional examples

local({
#### What are reasonable values of t?
    ##   If we assume sampling is exponential with rate 1,
    ##   t = 0.69 is the median
    ##   t = 1    is the mean
    ##   t = 1.39 is quantile 0.75
    ##   t = 1.9  is quantile 0.85


#### How the "standard, full output" looks like
    ## For all models, you can do this, to see examples
    ## This requires evamtools to work.

    if (require(evamtools)) {
        rcbn <- evamtools::random_evam(7, model = "CBN")
        rot <- evamtools::random_evam(7, model = "OT", ot_oncobn_epos = 0)
        rdbn_d <- evamtools::random_evam(7, model = "OncoBN", ot_oncobn_epos = 0,
                                         oncobn_model = "DBN")
        rdbn_c <- evamtools::random_evam(7, model = "OncoBN", ot_oncobn_epos = 0,
                                         oncobn_model = "CBN")
        rmhn <- evamtools::random_evam(7, model = "MHN")
        rhes <- evamtools::random_evam(7, model = "HESBCN")
    }
    ## But, below, I create the basic objects by hand.


#### MHN
    ## Define model: the log-Theta matrix.
    t1 <- matrix(runif(36, -4, 4), ncol = 6)
    colnames(t1) <- rownames(t1) <- LETTERS[1:6]

    ## Kill a gene
    t1b <- kill_gene(t1, "B")
    ## Obtain standard, full output
    t1bo <- get_full_output(t1b)

    ## Population at some fixed time, e.g., t = 2
    genots_at_t_from_trm(t1bo$MHN_trans_rate_mat, 2)



#### CBN
    ## Define model
    m1 <- data.frame(From = c("Root", "Root", "B", "B", "C", "C"),
                     To   = c("A",    "B",    "C", "D",  "E", "F"),
                     rerun_lambda = 1:6)
    ## Kill a gene
    m1c <- kill_gene(m1, "C")
    ## Obtain standard, full output
    m1co <- get_full_output(m1c)
    ## Sample at some other time, e.g., t = 2
    genots_at_t_from_trm(m1co$CBN_trans_rate_mat, 2)



#### HESBCN
    ## Define model
    m3x <- data.frame(
        From = c("Root", "Root", "Root", "A", "B", "D", "C", "E"),
        To =   c("A",    "B",    "D",    "C", "C", "E", "E", "F"),
        Lambdas = c(1, 2, 3, 4, 4, 5, 5, 6),
        Relation = c(rep("Single", 3), "AND", "AND", "XOR", "XOR", "Single"))
    ## Kill a gene
    m3xc <- kill_gene(m3x, "C")
    ## Obtain standard, full output
    m3xco <- get_full_output(m3xc)
    ## Sample at some other time, e.g., t = 2
    genots_at_t_from_trm(m3xco$HESBCN_trans_rate_mat, 2)


####  For OT and OncoBN we cannot get predicted at a given time.
    ##      And we need evamtools and other stuff
    if (require(evamtools)) {
###### OT
        require(Oncotree)
        motc <- kill_gene(rot$OT_model, "C")
        motco <- get_full_output(motc)
        motco$OT_predicted_genotype_freqs
        rot$OT_predicted_genotype_freqs

###### OncoBN

        mdbc <- kill_gene(rdbn_d$OncoBN_model, "C")
        mdbco <- get_full_output(mdbc)
        mdbco$OncoBN_predicted_genotype_freqs
        rdbn_d$OncoBN_predicted_genotype_freqs
    }

#### HyperHMM

    hm1 <- Matrix(0, nrow = 8, ncol = 8, sparse = TRUE)
    colnames(hm1) <- rownames(hm1) <- allGenotypesLetter(3)
    hm1[1, 2:4] <- c(0.1, 0.2, 0.7)
    hm1[2, 5:6] <- c(0.6, 0.4)
    hm1[3, c(5, 7)] <- c(0.2, 0.8)
    hm1[4, 6:7] <- c(0.1, 0.9)
    hm1[5:7, 8] <- 1
    attr(hm1, "method_output") <- "HyperHMM_trans_mat"
    attr(hm1, "num_prob.set") <- c(0.1, 0.2, 0.3, 0.4)
    attr(hm1, "num_features") <- 3

    ## Kill gene A
    hm1_A <- kill_gene(hm1, "A")
    ## Standard output
    hm1_bo <- get_full_output(hm1_A)

})



### Miscell tests, unreachable destinations, and fixing the first issue

test_that("We catch models that are just plainly wrong", {
    ## This is unnecessary. These models are just wrong.
    ## The sort of thing that someone who does not understand
    ## them could create, but that evamtools would never give.
    ## What is wrong? The same child node has different lambda

    m1 <- data.frame(From = c("Root", "Root", "A", "A", "B", "B"),
                     To   = c("A", "B", "C", "E", "E", "F"),
                     rerun_lambda = 1:6)

    m4 <- data.frame(From = c("Root", "Root", "E", "E", "B", "B", "C", "D"),
                     To = c("E", "B", "C", "D", "C", "D", "A", "A"),
                     Lambdas = 1:8,
                     Relation = c("Single", "Single", "And", "And", "And", "And", "OR", "OR"))

    for (gg in LETTERS[1:6]) {
        expect_error(kill_gene(m1, gg),
                     "Different lambda/weight for same destination gene")
    }
    ## Unless we use mc.cores= 1, the error is in mclapply and testthat
    ## complains
    expect_error(intervene_cpm_every_gene(list(CBN_model = m1), "CBN",
                                          mc.cores = 1),
                 "Different lambda/weight for same destination gene")

    for (gg in LETTERS[1:5]) {
        expect_error(kill_gene(m4, gg),
                     "Different lambda/weight for same destination gene")
    }

    expect_error(intervene_cpm_every_gene(list(HESBCN_model = m4), "HESBCN",
                                          mc.cores = 1),
                 "Different lambda/weight for same destination gene")

    expect_error(get_full_output(m1),
                 "Different lambda/weight for same destination gene.")

    expect_error(get_full_output(m4),
                 "Different lambda/weight for same destination gene.")
})

#### Explaining the warning about unreachable destinations

## What is the deal with the "weighted_fgraph contains unreachable destinations"?
## Function cpm2tm, called from function get_full_output,
## itself called from get_genotype_freqs_cpm, checks if there
## are, in the fitness graph, unreachable destinations,
## genotypes that have incoming connections, genotypes
## that are connected, as descendants, to other genotypes,
## but that are really not reachable.

## We get the warning from all models, except OncoBN,
## because, in cpm2tm, for OncoBN models function
## adjm_rm_no_access removes non-accessible genotypes.
## See the evamtools comments in code, but for OncoBN, in "regular use"
## sometimes some edges have weight 0, which leads to
## unreachable destinations.

## Here, I rewrite out get_genotype_freqs_cpm to produce
## the full output too, so we can see those transition rate matrices

## For CBN and H-ESBCN, note that we call, on the trans_rate_mat
## genots_from_trm; but the original freqs, stored
## in the output, are of course correct, as the trans. rate matrix
## computation properly accounts for the fact that there are
## rows that emit but nothing gets to them (so there is a 0)
## For OT and OncoBN we directly use predicted_genotype_freqs
## and this function does not use, for OT and OncoBN
## the transition/weighted_fraph matrices at all.

test_that("Explaining the warning about unreachable destinations", {
  ## Wrapped inside a "test_that" just to shut up testthat
  ## about code being run outside of test_that
  local({
    ## Modified for get_genotype_freqs_cpm, for demonstration
    gg2 <- function(model, t = NA) {
        if (nrow(model) == 0) return(c(WT = 1))

        ## Find out the method
        if (is.matrix(model) &&
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

        if (nrow(model) == 1 &&
            method %in% c("OT", "OncoBN")) {
            mut_gene <- model[1, 2]
            freq_mut <- model[1, ifelse(method == "OT", "OT_edgeWeight", "theta")]
            freqs <- c(1 - freq_mut, freq_mut)
            names(freqs) <- c("WT", mut_gene)
            return(freqs)
        }

        output <- get_full_output(model)
        if (method %in% c("OT", "OncoBN")) {
            return(list(
                all_out = output,
                pred_genots = output[[paste0(method, "_predicted_genotype_freqs")]]))
        }
        ## CBN, MHN, HESBCN
        trans_name <- paste0(method, "_trans_rate_mat")
        trans_rate_mat <- output[[trans_name]]
        return(list(all_out = output,
                    pred_genots = genots_from_trm(trans_rate_mat, t = t)))
    }


    ## Now, reuse those where everything depends on just one

    ## CBN
    m1 <- data.frame(From = c("Root", "A", "A", "A"),
                     To   = c("A",    "B", "C", "D"),
                     rerun_lambda = 1:4)

    ## OT
    m2 <- data.frame(From = c("Root", "B", "B", "B"),
                     To   = c("B",    "A", "C", "D"),
                     OT_edgeWeight = rep(0.3, 4))

    ## OncoBN
    m3 <- data.frame(From = c("Root", "Z", "Z", "Z"),
                     To   = c("Z",    "B", "C", "D"),
                     theta = rep(0.2, 4),
                     Relation = "Single")

    ## H-ESBCN
    m4 <- data.frame(From = c("Root", "E", "E", "E"),
                     To   = c("E",    "B", "C", "D"),
                     Lambdas = 1:4,
                     Relation = "Single")

    gg2m1 <- gg2(kill_gene_by_params_to_0(m1, "A", verbose = TRUE))
    gg2m2 <- gg2(kill_gene_by_params_to_0(m2, "B", verbose = TRUE))
    gg2m3 <- gg2(kill_gene_by_params_to_0(m3, "Z", verbose = TRUE))
    gg2m4 <- gg2(kill_gene_by_params_to_0(m4, "E", verbose = TRUE))


    #### CBN
    ## Notice there are non-zero entries in destinations like A, B  and A, C, etc
    gg2m1[["all_out"]][["CBN_trans_rate_mat"]]
    gg2m1[["all_out"]][["CBN_trans_mat"]]
    ## Notice the next one correctly says from WT we jumpt to WT with prob 1
    gg2m1[["all_out"]][["CBN_td_trans_mat"]]
    ## Only WT is non-zero
    gg2m1[["all_out"]][["CBN_predicted_genotype_freqs"]]
    gg2m1[["pred_genots"]]


    #### HESBCN
    ## Much like CBN
    gg2m4[["all_out"]][["HESBCN_trans_rate_mat"]]
    gg2m4[["all_out"]][["HESBCN_trans_mat"]]
    gg2m4[["all_out"]][["HESBCN_td_trans_mat"]]
    ## Only WT is non-zero
    gg2m4[["all_out"]][["HESBCN_predicted_genotype_freqs"]]
    gg2m4[["pred_genots"]]

    #### OT
    ## The first two show the inaccessible genotypes
    gg2m2[["all_out"]][["OT_f_graph"]]
    gg2m2[["all_out"]][["OT_trans_mat"]]
    ## Only WT
    gg2m2[["all_out"]][["OT_predicted_genotype_freqs"]]
    ## the next is identical to former, necessarily, from code
    gg2m2[["pred_genots"]]


    #### OncoBN
    ## Because we filter with adjm_rm_no_access
    ## here there are no inaccessible genotypes
    gg2m3[["all_out"]][["OncoBN_f_graph"]]
    gg2m3[["all_out"]][["OncoBN_trans_mat"]]
    ## Only WT
    gg2m3[["all_out"]][["OncoBN_predicted_genotype_freqs"]]
    ## the next is identical to former, necessarily, from code
    gg2m3[["pred_genots"]]



    #### More playing, if you want
    ## We can also play with this, but too cumbersome
    set.seed(1)
    rcbn <- evamtools::random_evam(7, model = "CBN")
    rhes <- evamtools::random_evam(7, model = "HESBCN")
    rot <- evamtools::random_evam(7, model = "OT", ot_oncobn_epos = 0)


    gg2cbn <- gg2(kill_gene_by_params_to_0(rcbn[[1]], "B", verbose = TRUE))
    gg2hes <- gg2(kill_gene_by_params_to_0(rhes[[2]], "B", verbose = TRUE))
    gg2ot <- gg2(kill_gene_by_params_to_0(rot[[1]], "B", verbose = TRUE))
    set.seed(NULL)
  })
})




#### Issue 1: https://github.com/rdiaz02/interv-CPM/issues/1

test_that("Issue 1 is solved", {
    ##  There are additional tests of issue in kill-gene-equivalences-TESTS.R

    ## Note that get_genotype_freqs_cpm already
    ## dealt with this case OK. It was get_full_output
    ## that didn't (get_genotype_freqs_cpm sidestepped
    ## calling get_full_output with single gene models and
    ## assumed epos = 0).
    ## And we also verify the "setting params to 0" works


  #### Generating models with a single gene fails. This is expected,
    ## and there is no reason to allow these things, at least not now.

    expect_error(random_evam(1, model = "CBN"))
    expect_error(random_evam(1, model = "HESBCN"))
    expect_error(random_evam(1, model = "MHN"))
    expect_error(random_evam(1, model =  "OncoBN"))
    expect_error(random_evam(1, model =  "OT"))

    set.seed(NULL)

  #### 2-gene models, where we kill one gene
    ## CBN and H-ESBCN work just fine
    (random_CBN <- random_evam(2, model = "CBN")$CBN_model)
  (one_gene_CBN <- suppressWarnings(kill_gene(random_CBN, "B")))
  o1 <- suppressWarnings(get_full_output(one_gene_CBN))
    o1$CBN_predicted_genotype_freqs

    (random_HESBCN <- random_evam(2, model = "HESBCN")$HESBCN_model)
  (one_gene_HESBCN <- suppressWarnings(kill_gene(random_HESBCN, "B")))
  o2 <- suppressWarnings(get_full_output(one_gene_HESBCN))
    o2$HESBCN_predicted_genotype_freqs

    ## Repeat a few times: this ain't a lucky thing
    ## Many generate a no gene model; fine.
    ##
    ## Note: structural kill reduces model dimensions (fewer genes), so
    ## get_genotype_freqs_cpm returns different-length vectors than
    ## params-to-0 kill (which keeps the full model structure).
    ## We compare after filtering > 0 to test equivalence of predictions.
    filter_preds <- function(x) list(
        genot_freqs = x$genot_freqs[x$genot_freqs > 0],
        hitting_probs_from_WT = x$hitting_probs_from_WT[x$hitting_probs_from_WT > 0])
    for (i in 1:10) {
        the_kg <- ifelse(i %% 2, "A", "B")

        (random_CBN <- random_evam(2, model = "CBN", graph_density = runif(1, 0.1, 0.5))$CBN_model)
    (one_gene_CBN <- suppressWarnings(kill_gene(random_CBN, the_kg)))
    o1 <- suppressWarnings(get_genotype_freqs_cpm(one_gene_CBN))
    o1_2 <- suppressWarnings(get_genotype_freqs_cpm(kill_gene_by_params_to_0(random_CBN, the_kg)))
        expect_equal(filter_preds(o1), filter_preds(o1_2))

        (random_HESBCN <- random_evam(2, model = "HESBCN", graph_density = runif(1, 0.1, 0.5))$HESBCN_model)
    (one_gene_HESBCN <- suppressWarnings(kill_gene(random_HESBCN, the_kg)))
    o2 <- suppressWarnings(get_genotype_freqs_cpm(one_gene_HESBCN))
    o2_2 <- suppressWarnings(get_genotype_freqs_cpm(kill_gene_by_params_to_0(random_HESBCN, the_kg)))
        expect_equal(filter_preds(o2), filter_preds(o2_2))
    }

    ## MHN used to fail , as it should as it leads to a single-gene model
    ## This is not relevant for us now, since we use many more than 2 genes
    ## But it is now solved too
    (random_MHN <- random_evam(2, model = "MHN")$MHN_theta)
    (one_gene_MHN <- kill_gene(random_MHN, "B"))
    mhn1_1 <- get_full_output(one_gene_MHN)
    ## This works, however, as expected
    mhn1_2 <- get_full_output(kill_gene_by_params_to_0(random_MHN, "B"))

    mhn1_1_pf <- mhn1_1$MHN_predicted_genotype_freqs
    mhn1_2_pf <- mhn1_2$MHN_predicted_genotype_freqs

    mhn1_1_pf <- mhn1_1_pf[mhn1_1_pf > 0]
    mhn1_2_pf <- mhn1_2_pf[mhn1_2_pf > 0]

    expect_equal(mhn1_1_pf, mhn1_2_pf)

    ## Nope, not a lucky thing
    set.seed(NULL)
    for (i in 1:10) {
        the_kg <- ifelse(i %% 2, "A", "B")
        (random_MHN <- random_evam(2, model = "MHN")$MHN_theta)
        (one_gene_MHN <- kill_gene(random_MHN, the_kg))
        mhn1_gfc <- get_genotype_freqs_cpm(one_gene_MHN)
        ## This works, however, as expected
        mhn1_p0 <- get_genotype_freqs_cpm(kill_gene_by_params_to_0(random_MHN, the_kg))
        ## Structural kill reduces model dimensions; params-to-0 keeps structure.
        ## Filter > 0 before comparing to test equivalence of predictions.
        expect_equal(filter_preds(mhn1_gfc), filter_preds(mhn1_p0))
        ## Refactor-resistance check (not an independent verification):
        ## for MHN, get_genotype_freqs_cpm internally calls
        ## get_full_output(model, epos = 0) and just renames the fields, so
        ## these two assertions are by construction true. They guard against
        ## accidental rename / removal of the field assignments in
        ## get_genotype_freqs_cpm (kill-gene-and-output-from-cpm.R). The
        ## meaningful equivalence check on the MHN predictions is the
        ## filter_preds(mhn1_gfc) == filter_preds(mhn1_p0) assertion above
        ## (structural kill vs params-to-0 kill).
        mhn1_full <- get_full_output(one_gene_MHN)
        expect_equal(mhn1_gfc$genot_freqs, mhn1_full$MHN_predicted_genotype_freqs)
        expect_equal(mhn1_gfc$hitting_probs_from_WT, mhn1_full$MHN_hitting_probs_from_WT)
    }


    (random_OT <- random_evam(2, model = "OT")$OT_model)
    (one_gene_OT <- kill_gene(random_OT, "B"))
    (o4 <- get_full_output(one_gene_OT))
    o4$OT_predicted_genotype_freqs

    (random_OncoBN <- random_evam(2, model = "OncoBN")$OncoBN_model)
    (one_gene_OncoBN <- kill_gene(random_OncoBN, "B"))
    (o5 <- get_full_output(one_gene_OncoBN))
    o5$OncoBN_predicted_genotype_freqs


  #### General models, where we end with a single gene model after killing

    ## The cases below all failed at some point. They were not failing before
    ## the fixes on 2024-08-01 when using get_genotype_freqs_cpm,
    ## as that function recognized models with only 1 (or 0)
    ## genes.

    set.seed(7)
    random_OT <- random_evam(7, model = "OT")$OT_model
    model_with_error <- kill_gene(random_OT, "A")
    get_full_output(model_with_error)

    set.seed(23)
    random_OT <- random_evam(7, model = "OT")$OT_model
    model_with_error <- kill_gene(random_OT, "A")
    get_full_output(model_with_error)

    set.seed(22)
    random_OT <- random_evam(7, model = "OT")$OT_model
    model_with_error <- kill_gene(random_OT, "B")
    get_full_output(model_with_error)

    set.seed(21)
    (random_7OncoBN <- random_evam(7, model = "OncoBN")$OncoBN_model)
    (seven_OncoBN <- kill_gene(random_7OncoBN, "B"))
    oo7 <- get_full_output(seven_OncoBN)

    set.seed(38)
    (random_7OncoBN <- random_evam(7, model = "OncoBN")$OncoBN_model)
    (seven_OncoBN <- kill_gene(random_7OncoBN, "A"))
    oo7 <- get_full_output(seven_OncoBN)

    set.seed(NULL)

})




### additional tests of issue 1 (https://github.com/rdiaz02/interv-CPM/issues/1)
#### 2-gene models, where we kill one gene: additional tests of issue 1

## For OT and OncoBN
## Note that for OncoBN when we have 0 genes, for kill_gene_by_params_to_0
## we must use the modified ev2_adjm_rm_no_access

## Repeat a few times: this ain't a lucky thing
## Recall that if we kill "A" we leave a 0 gene model
## and get_full_output will crash there. That is not a bug,
## as the no gene model is handled by get_genotype_freqs_cpm
## which is the function we use.
## In fact, it is questionable that we test, here,
## get_full_output, since that is never called directly.
test_that("Issue 1 is solved, additional", {
  local({
    set.seed(NULL)
    ## Structural kill reduces model dimensions; params-to-0 keeps structure.
    ## Filter > 0 before comparing to test equivalence of predictions.
    filter_preds <- function(x) list(
        genot_freqs = x$genot_freqs[x$genot_freqs > 0],
        hitting_probs_from_WT = x$hitting_probs_from_WT[x$hitting_probs_from_WT > 0])
    for (i in 1:10) {
        the_kg <- ifelse(i %% 2, "A", "B")
        (random_OT <- random_evam(2, model = "OT")$OT_model)
        (one_gene_OT <- kill_gene(random_OT, the_kg))
        o4 <- get_genotype_freqs_cpm(one_gene_OT)
        o42 <- get_genotype_freqs_cpm(kill_gene_by_params_to_0(random_OT, the_kg))
        expect_equal(filter_preds(o4), filter_preds(o42))

        if (the_kg == "B") {
            o43_full <- get_full_output(one_gene_OT)
            expect_equal(o4$genot_freqs, o43_full$OT_predicted_genotype_freqs)
            expect_equal(o4$hitting_probs_from_WT, o43_full$OT_hitting_probs_from_WT)
        }

        (random_OncoBN <- random_evam(2, model = "OncoBN")$OncoBN_model)
        (one_gene_OncoBN <- kill_gene(random_OncoBN, the_kg))
        o5 <- get_genotype_freqs_cpm(one_gene_OncoBN)
        o52 <- get_genotype_freqs_cpm(kill_gene_by_params_to_0(random_OncoBN, the_kg))
        expect_equal(filter_preds(o5), filter_preds(o52))
        if (the_kg == "B") {
            o53_full <- get_full_output(one_gene_OncoBN)
            expect_equal(o5$genot_freqs, o53_full$OncoBN_predicted_genotype_freqs)
            expect_equal(o5$hitting_probs_from_WT, o53_full$OncoBN_hitting_probs_from_WT)
        }
    }
  })
})

set.seed(NULL)

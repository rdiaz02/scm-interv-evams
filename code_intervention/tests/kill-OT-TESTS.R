## options(intervention_every_gene_cores = parallel::detectCores())
library(testthat)
pwd <- getwd()
setwd("../")
source("intervention.R")
setwd(pwd)

set.seed(NULL)

test_that("OT: compare predicted genotype freqs killing", {
    ## Kill with the two procedures for OT
    ## and compare with predicted genotype freqs of
    ## hand crafted structures
    ## and hand crafted structures with the 0 parameter


    full_output_just_preds <- function(x) {
        tmp <- get_full_output(x)
        f <- tmp$OT_predicted_genotype_freqs
        hp <- tmp$OT_hitting_probs_from_WT
        return(list(genot_freqs = f[f > 0], hitting_probs_from_WT = filter_hp_keep_wt(hp)))
    }

    get_interv <- function(res, name) {
        return(res[[name]])
    }

    set.seed(NULL)
    ot1 <- data.frame(
        From = c("Root", "Root", "C", "C", "C", "M", "M", "A", "E", "E", "X"),
        To   = c("C",  "M",  "A", "E", "B", "G", "D", "X", "H", "K", "I"),
        OT_edgeWeight = runif(11)
    )

    ## See this DAG of restrictions
    ## evamtools:::DAG_plot_graphAM(ot1, "ot1")

    ## This is slow. Oh well, 11 genes, which is much more than
    ## we use.
    i_ot1 <- intervene_cpm_every_gene(list(OT_model = ot1), "OT",
                                      verbose = TRUE)
    i_ot1_p0 <- intervene_cpm_every_gene(list(OT_model = ot1), "OT",
                                         verbose = TRUE,
                                         kill_gene_funct = kill_gene_by_params_to_0)


    ## The order in the model entry does not matter
    ot1r <- ot1[sample(1:nrow(ot1)), ]
    i_ot1r <- intervene_cpm_every_gene(list(OT_model = ot1r), "OT",
                                       verbose = TRUE)
    i_ot1r_p0 <- intervene_cpm_every_gene(list(OT_model = ot1r), "OT",
                                          verbose = TRUE,
                                          kill_gene_funct = kill_gene_by_params_to_0)

    ## The two killing procedures equal, and
    ## order in model makes no difference
    expect_equal(i_ot1, i_ot1_p0)
    expect_equal(i_ot1r, i_ot1r_p0)
    expect_equal(i_ot1, i_ot1r)

    ## DAG of restrictions removal of the entry
    ## Deliberately do not use code to remove rows, we give
    ## then by number, to minimize possibly repeating mistakes in code
    ot1_A <- ot1[-c(3, 8, 11), ]
    ot1_B <- ot1[-c(5), ]
    ot1_C <- ot1[-c(1, 3, 4, 5, 8, 9, 10, 11), ]
    ot1_D <- ot1[-c(7), ]
    ot1_E <- ot1[-c(4, 9, 10), ]
    ot1_G <- ot1[-c(6), ]
    ot1_H <- ot1[-c(9), ]
    ot1_I <- ot1[-c(11), ]
    ot1_K <- ot1[-c(10), ]
    ot1_M <- ot1[-c(2, 6, 7), ]
    ot1_X <- ot1[-c(8, 11), ]

    gene_names <- sort(setdiff(union(ot1$From, ot1$To), "Root"))

    ## No intervention. This is slow!
    expect_equal(get_interv(i_ot1, "no_intervention"), full_output_just_preds(ot1))

    ## Each intervention
    for (g in gene_names) {
        message("Doing gene ", g)
        expect_equal(get_interv(i_ot1, paste0("I:", g)),
                     full_output_just_preds(get(paste0("ot1_", g))))
    }

})


test_that("Simple OT testing", {
    local_edition(3)

    m1_ot <- data.frame(From = c("Root", "Root", "B", "B", "C", "C"),
                        To   = c("A",    "B",    "C", "D",  "E", "F"),
                        OT_edgeWeight = runif(6))
    kills_m1_ot <- lapply(LETTERS[1:6], function(v) kill_gene(m1_ot, v))

    names(kills_m1_ot) <- LETTERS[1:6]

    ## Compare structures
    expect_equal(kills_m1_ot[["A"]], m1_ot[-c(1), ])
    expect_equal(kills_m1_ot[["B"]], m1_ot[1, ])
    expect_equal(kills_m1_ot[["C"]], m1_ot[-c(3, 5, 6), ])
    expect_equal(kills_m1_ot[["D"]], m1_ot[-c(4), ])
    expect_equal(kills_m1_ot[["E"]], m1_ot[-c(5), ])
    expect_equal(kills_m1_ot[["F"]], m1_ot[-c(6), ])

    ## Compare predictions of both methods of killing
    m1_ot_i <- intervene_cpm_every_gene(list(OT_model = m1_ot), "OT",
                                        verbose = TRUE,
                                           kill_gene_funct = kill_gene)

    m1_ot_i_p0 <- intervene_cpm_every_gene(list(OT_model = m1_ot), "OT",
                                           verbose = TRUE,
                                        kill_gene_funct = kill_gene_by_params_to_0)

    expect_equal(m1_ot_i, m1_ot_i_p0)
})

set.seed(NULL)

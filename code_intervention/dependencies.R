## This file has two functions:
## a) tell users what they need
## b) lead to immediate test abortion if one is missing: this file is
##    sourced from ./run-all-tests.R

## Packages
library(parallel)
library(RhpcBLASctl)
library(OncoSimulR)
library(evamtools)
library(gtools)
library(igraph)
library(stringi)
library(Matrix)
library(expm)
library(hyperhmm)
library(uuid)
library(markovchain)
library(testthat)
## library(doesnotexist) ## tripwire.

## intervention  pulls:
##   kill-gene-and-output-from-cpm.R pulls:
##      HyperHMM-wrapper.R
##      evam_v2.R
##      trm.R
##        utils.R
##        rfitness_to_trm.R

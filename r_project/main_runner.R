##### MAIN RUNNER ######
if(sys.nframe() == 0L) rm(list = ls())
library(knitr)
library(here)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

# source(knitr::purl("empirical/prep_preqin.Rmd", output = tempfile(), quiet = TRUE))

data.out.folder <- "results/data_out_2026_02_18"

# OUTDATED: source("estim_model_empirical_runner.R")
source("run_empirical_study.R")

source("results/analyze_results_runner.R")

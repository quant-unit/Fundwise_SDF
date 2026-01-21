##### MAIN RUNNER ######
if(sys.nframe() == 0L) rm(list = ls())
library(knitr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

# source(knitr::purl("empirical/prep_preqin.Rmd", output = tempfile(), quiet = TRUE))

data.out.folder <- "results/data_out_2026-emp-max-vin-2019"

# source("estim_model_optimized.R")
source("estim_model_empirical_runner.R")

# source("results/analyze_result.R")
source("results/analyze_results_runner.R")

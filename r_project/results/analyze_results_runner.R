### analyze_results runner
# 0. Prologue ----
if(sys.nframe() == 0L) rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

source.internally <- FALSE

# set output folder
if(!exists("data.out.folder", envir = .GlobalEnv)) {
  data.out.folder <- "results/data_out_2026-XXX"
  data.out.folder <- "results/data_out_2026-emp-max-vin-2019"
  data.out.folder <- ""
}

# 1. analyze all results ----

prefix <- "q_factors_preqin_"
suffix <- "EW_VYP"

source(here("results/analyze_result.R"))

prefix <- "q_factors_preqin_"
suffix <- "FW_VYP"

source(here("results/analyze_result.R"))

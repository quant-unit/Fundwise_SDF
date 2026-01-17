### analyze_results runner
# 0. Prologue ----
if(sys.nframe() == 0L) rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

source.internally <- FALSE

data.out.folder <- "data_out_2026-emp"

# 1. analyze all results ----

prefix <- "q_factors_preqin_"
suffix <- "EW_VYP"

source("analyze_result.R")

prefix <- "q_factors_preqin_"
suffix <- "FW_VYP"

source("analyze_result.R")

if (FALSE) {
  prefix <- "q_factors_preqin_alpha_"
  suffix <- "EW_VYP"
  
  source("analyze_result.R")
  
  prefix <- "q_factors_preqin_alpha_"
  suffix <- "FW_VYP"
  
  source("analyze_result.R")
}

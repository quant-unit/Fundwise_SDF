### analyze_results runner
# 0. Prologue ----
if (sys.nframe() == 0L) rm(list = ls())

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

library(here)

source.internally <- FALSE

# set output folder
if (!exists("data.out.folder", envir = .GlobalEnv)) {
  data.out.folder <- "results/data_out_2026-XXX"
  data.out.folder <- "results/data_out_2026-emp-max-vin-2019"
  data.out.folder <- "results/data_out_2026_02_18"
  data.out.folder <- "results/data_out_2026_02_24"
  data.out.folder <- "results/data_out_2026_02_26"
}

# 1. analyze base results ----

prefix <- "q_factors_preqin_"
suffix <- "EW_VYP"

source(here("results/analyze_result.R"))

prefix <- "q_factors_preqin_"
suffix <- "FW_VYP"

source(here("results/analyze_result.R"))


# 2. Vintage Year Cutoffs ----

suffix <- paste0("EW_VYP_max_vin_", 2011)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_max_vin_", 2012)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_max_vin_", 2013)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_max_vin_", 2014)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_max_vin_", 2015)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_max_vin_", 2016)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_max_vin_", 2017)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_max_vin_", 2018)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_max_vin_", 2019)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_max_vin_", 2020)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_max_vin_", 2021)
source(here("results/analyze_result.R"))


suffix <- paste0("FW_VYP_max_vin_", 2011)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_max_vin_", 2012)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_max_vin_", 2013)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_max_vin_", 2014)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_max_vin_", 2015)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_max_vin_", 2016)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_max_vin_", 2017)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_max_vin_", 2018)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_max_vin_", 2019)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_max_vin_", 2020)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_max_vin_", 2021)
source(here("results/analyze_result.R"))

# 2b. Vintage Year Cutoffs (NC50) ----

prefix <- "q_factors_preqin_"

suffix <- paste0("EW_VYP_NC50_max_vin_", 2011)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_NC50_max_vin_", 2012)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_NC50_max_vin_", 2013)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_NC50_max_vin_", 2014)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_NC50_max_vin_", 2015)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_NC50_max_vin_", 2016)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_NC50_max_vin_", 2017)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_NC50_max_vin_", 2018)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_NC50_max_vin_", 2019)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_NC50_max_vin_", 2020)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_NC50_max_vin_", 2021)
source(here("results/analyze_result.R"))

suffix <- paste0("FW_VYP_NC50_max_vin_", 2011)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_NC50_max_vin_", 2012)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_NC50_max_vin_", 2013)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_NC50_max_vin_", 2014)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_NC50_max_vin_", 2015)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_NC50_max_vin_", 2016)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_NC50_max_vin_", 2017)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_NC50_max_vin_", 2018)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_NC50_max_vin_", 2019)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_NC50_max_vin_", 2020)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_NC50_max_vin_", 2021)
source(here("results/analyze_result.R"))

# 3. Fama French Factors ----

prefix <- "ff3_factors_preqin_alpha_ALL_"
suffix <- "EW_VYP"
source(here("results/analyze_result.R"))

suffix <- "FW_VYP"
source(here("results/analyze_result.R"))


# 4a. Pitchbook: Vintage Year Cutoffs ----

data.out.folder <- "results/data_out_2026_02_27_pitchbook"
prefix <- "q_factors_pitchbook_"

suffix <- paste0("EW_VYP_max_vin_", 2011)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_max_vin_", 2012)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_max_vin_", 2013)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_max_vin_", 2014)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_max_vin_", 2015)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_max_vin_", 2016)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_max_vin_", 2017)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_max_vin_", 2018)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_max_vin_", 2019)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_max_vin_", 2020)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_max_vin_", 2021)
source(here("results/analyze_result.R"))


suffix <- paste0("FW_VYP_max_vin_", 2011)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_max_vin_", 2012)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_max_vin_", 2013)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_max_vin_", 2014)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_max_vin_", 2015)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_max_vin_", 2016)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_max_vin_", 2017)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_max_vin_", 2018)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_max_vin_", 2019)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_max_vin_", 2020)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_max_vin_", 2021)
source(here("results/analyze_result.R"))

# 4b. Pitchbook: Vintage Year Cutoffs (NC50) ----


suffix <- paste0("EW_VYP_NC50_max_vin_", 2011)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_NC50_max_vin_", 2012)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_NC50_max_vin_", 2013)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_NC50_max_vin_", 2014)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_NC50_max_vin_", 2015)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_NC50_max_vin_", 2016)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_NC50_max_vin_", 2017)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_NC50_max_vin_", 2018)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_NC50_max_vin_", 2019)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_NC50_max_vin_", 2020)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_NC50_max_vin_", 2021)
source(here("results/analyze_result.R"))

suffix <- paste0("FW_VYP_NC50_max_vin_", 2011)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_NC50_max_vin_", 2012)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_NC50_max_vin_", 2013)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_NC50_max_vin_", 2014)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_NC50_max_vin_", 2015)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_NC50_max_vin_", 2016)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_NC50_max_vin_", 2017)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_NC50_max_vin_", 2018)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_NC50_max_vin_", 2019)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_NC50_max_vin_", 2020)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_NC50_max_vin_", 2021)
source(here("results/analyze_result.R"))

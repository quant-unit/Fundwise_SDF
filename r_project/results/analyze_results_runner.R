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
  # data.out.folder <- "results/data_out_2026_03_26"
  
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

# 5. Preqin North America USD results ----

data.out.folder <- "results/data_out_2026_03_26"

prefix <- "q_factors_preqin_"
suffix <- "EW_VYP_North America"

source(here("results/analyze_result.R"))

prefix <- "q_factors_preqin_"
suffix <- "FW_VYP_North America"

source(here("results/analyze_result.R"))

# 6. Preqin North America USD Vintage Year Cutoffs ----

data.out.folder <- "results/data_out_2026_03_26"

suffix <- paste0("EW_VYP_North America_max_vin_", 2011)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_North America_max_vin_", 2012)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_North America_max_vin_", 2013)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_North America_max_vin_", 2014)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_North America_max_vin_", 2015)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_North America_max_vin_", 2016)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_North America_max_vin_", 2017)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_North America_max_vin_", 2018)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_North America_max_vin_", 2019)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_North America_max_vin_", 2020)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_North America_max_vin_", 2021)
source(here("results/analyze_result.R"))


suffix <- paste0("FW_VYP_North America_max_vin_", 2011)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_North America_max_vin_", 2012)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_North America_max_vin_", 2013)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_North America_max_vin_", 2014)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_North America_max_vin_", 2015)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_North America_max_vin_", 2016)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_North America_max_vin_", 2017)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_North America_max_vin_", 2018)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_North America_max_vin_", 2019)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_North America_max_vin_", 2020)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_North America_max_vin_", 2021)
source(here("results/analyze_result.R"))

# 6b. Preqin North America USD Vintage Year Cutoffs (NC50) ----

prefix <- "q_factors_preqin_"

suffix <- paste0("EW_VYP_North America_NC50_max_vin_", 2011)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_North America_NC50_max_vin_", 2012)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_North America_NC50_max_vin_", 2013)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_North America_NC50_max_vin_", 2014)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_North America_NC50_max_vin_", 2015)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_North America_NC50_max_vin_", 2016)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_North America_NC50_max_vin_", 2017)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_North America_NC50_max_vin_", 2018)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_North America_NC50_max_vin_", 2019)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_North America_NC50_max_vin_", 2020)
source(here("results/analyze_result.R"))
suffix <- paste0("EW_VYP_North America_NC50_max_vin_", 2021)
source(here("results/analyze_result.R"))

suffix <- paste0("FW_VYP_North America_NC50_max_vin_", 2011)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_North America_NC50_max_vin_", 2012)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_North America_NC50_max_vin_", 2013)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_North America_NC50_max_vin_", 2014)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_North America_NC50_max_vin_", 2015)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_North America_NC50_max_vin_", 2016)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_North America_NC50_max_vin_", 2017)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_North America_NC50_max_vin_", 2018)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_North America_NC50_max_vin_", 2019)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_North America_NC50_max_vin_", 2020)
source(here("results/analyze_result.R"))
suffix <- paste0("FW_VYP_North America_NC50_max_vin_", 2021)
source(here("results/analyze_result.R"))

# 7. Preqin North America USD NAV Discounts ----

data.out.folder <- "results/data_out_2026_03_26"

prefix <- "q_factors_preqin_"
suffix <- "EW_VYP_North America_NC0"

source(here("results/analyze_result.R"))

prefix <- "q_factors_preqin_"
suffix <- "FW_VYP_North America_NC0"

source(here("results/analyze_result.R"))

prefix <- "q_factors_preqin_"
suffix <- "EW_VYP_North America_NC25"

source(here("results/analyze_result.R"))

prefix <- "q_factors_preqin_"
suffix <- "FW_VYP_North America_NC25"

source(here("results/analyze_result.R"))

prefix <- "q_factors_preqin_"
suffix <- "EW_VYP_North America_NC50"

source(here("results/analyze_result.R"))

prefix <- "q_factors_preqin_"
suffix <- "FW_VYP_North America_NC50"

source(here("results/analyze_result.R"))

prefix <- "q_factors_preqin_"
suffix <- "EW_VYP_North America_NC75"

source(here("results/analyze_result.R"))

prefix <- "q_factors_preqin_"
suffix <- "FW_VYP_North America_NC75"

source(here("results/analyze_result.R"))



# 8. Fama French Factors North America ----

data.out.folder <- "results/data_out_2026_03_26"


"cache_ff3_factors_preqin_alpha_ALL_EW_VYP_North America"
"cache_ff3_factors_preqin_alpha_Alpha_EW_VYP_North America"
"cache_ff3_factors_preqin_MKT_EW_VYP_North America"
"cache_ff3_factors_preqin_MKT_SMB_HML_EW_VYP_North America"

prefix <- "ff3_factors_preqin_alpha_ALL_"

prefix <- "ff3_factors_preqin_alpha_Alpha_"

suffix <- "EW_VYP_North America"
source(here("results/analyze_result.R"))
suffix <- "FW_VYP_North America"
source(here("results/analyze_result.R"))

prefix <- "ff3_factors_preqin_MKT_"

suffix <- "EW_VYP_North America"
source(here("results/analyze_result.R"))
suffix <- "FW_VYP_North America"
source(here("results/analyze_result.R"))

prefix <- "ff3_factors_preqin_MKT_SMB_HML_"

suffix <- "EW_VYP_North America"
source(here("results/analyze_result.R"))
suffix <- "FW_VYP_North America"
source(here("results/analyze_result.R"))
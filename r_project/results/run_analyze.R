library(here)
source.internally <- FALSE
prefix <- "ff3_factors_preqin_alpha_"
suffix <- "ALL_FW_VYP"
data.out.folder <- "results/data_out_2026_02_18"
max.months.to.keep <- c("0", "1", "30", "60", "120", "180", "240", "300", "360")

# Mock rstudioapi to prevent crash
rstudioapi <- new.env()
rstudioapi$getActiveDocumentContext <- function() {
  list(path = file.path(getwd(), "analyze_result.R"))
}
attach(rstudioapi, name = "rstudioapi")

source("analyze_result.R")

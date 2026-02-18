### estimate model runner: Emprical
# 0. Prologue ----
# if(sys.nframe() == 0L) rm(list = ls())
rm(list = setdiff(ls(), "data.out.folder"))

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

# 0. Set Constants for Empirical Runs (same for all) ----
source.internally <- FALSE
export.data <- FALSE
use.vintage.year.pfs <- TRUE
use.simulation <- FALSE
do.cache <- TRUE
do.parallel <- ifelse(.Platform$OS.type == "windows", FALSE, TRUE)
public.filename <- "q_factors"
max.months <- c(1, 30, 60, 120, 150, 180, 210, 240, 300, 360) # c(10, 20) * 12 # c(12.5, 15, 17.5) * 12
max.months <- c(1, 30, 60, 120, 150, 180)
lambdas <- 0
kernel.bandwidth <- 12
alpha.lower <- -0.01
alpha.upper <- 0.01
part.to.keep <- 1
no.partitions <- 1 # 10
max.vintage <- 2021

# set output folder
if (!exists("data.out.folder", envir = .GlobalEnv)) {
  data.out.folder <- "results/data_out_2029"
}
if (!dir.exists(data.out.folder)) dir.create(data.out.folder)

factors.to.use <- ""
cutoff <- ""

data.prepared.folder <- "data_prepared_2026"
private.source <- "preqin"
error.function <- "L2_Lasso"
sdf.model <- "linear"
include.alpha.term <- FALSE

# 1. FW & No CV ----

# CHOICES
do.cross.validation <- FALSE
weighting <- "FW"

# SET cache.folder.tag
if (use.vintage.year.pfs) weighting <- paste0(weighting, "_VYP")
cache.folder.tag <- paste0(private.source, ifelse(include.alpha.term, "_alpha_", "_"))
cache.folder.tag <- paste0(cache.folder.tag, weighting)
cache.folder.tag

# RUN
source("estim_model_optimized.R")

# 2. FW & CV ----

# CHOICES
do.cross.validation <- TRUE
weighting <- "FW"

# SET cache.folder.tag
if (use.vintage.year.pfs) weighting <- paste0(weighting, "_VYP")
cache.folder.tag <- paste0(private.source, ifelse(include.alpha.term, "_alpha_", "_"))
cache.folder.tag <- paste0(cache.folder.tag, weighting)
cache.folder.tag

# RUN
source("estim_model_optimized.R")

# 3. EW & No CV ----

# CHOICES
do.cross.validation <- FALSE
weighting <- "EW"

# SET cache.folder.tag
if (use.vintage.year.pfs) weighting <- paste0(weighting, "_VYP")
cache.folder.tag <- paste0(private.source, ifelse(include.alpha.term, "_alpha_", "_"))
cache.folder.tag <- paste0(cache.folder.tag, weighting)
cache.folder.tag

# RUN
source("estim_model_optimized.R")

# 4. EW & CV ----

# CHOICES
do.cross.validation <- TRUE
weighting <- "EW"

# SET cache.folder.tag
if (use.vintage.year.pfs) weighting <- paste0(weighting, "_VYP")
cache.folder.tag <- paste0(private.source, ifelse(include.alpha.term, "_alpha_", "_"))
cache.folder.tag <- paste0(cache.folder.tag, weighting)
cache.folder.tag

# RUN
source("estim_model_optimized.R")


# 5. FW + Max Vintage 2020 - 2010 Runs -----

# CHOICES
do.cross.validation <- FALSE
weighting <- "FW"
if (use.vintage.year.pfs) weighting <- paste0(weighting, "_VYP")
max.vintages <- 2011:2021

for (max.vintage in max.vintages) {
  # SET cache.folder.tag
  cache.folder.tag <- paste0(private.source, ifelse(include.alpha.term, "_alpha_", "_"))
  cache.folder.tag <- paste0(cache.folder.tag, weighting, "_max_vin_", max.vintage)
  cache.folder.tag

  # RUN
  source("estim_model_optimized.R")
}



# 6. EW + Max Vintage 2020 - 2010 Runs -----

# CHOICES
do.cross.validation <- FALSE
weighting <- "EW"
if (use.vintage.year.pfs) weighting <- paste0(weighting, "_VYP")
max.vintages <- 2011:2021

for (max.vintage in max.vintages) {
  # SET cache.folder.tag
  cache.folder.tag <- paste0(private.source, ifelse(include.alpha.term, "_alpha_", "_"))
  cache.folder.tag <- paste0(cache.folder.tag, weighting, "_max_vin_", max.vintage)
  cache.folder.tag

  # RUN
  source("estim_model_optimized.R")
}

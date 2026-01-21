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
part.to.keep <- 1
no.partitions <- 1 # 10

# set output folder
if(!exists("data.out.folder", envir = .GlobalEnv)) {
  data.out.folder <- "results/data_out_2026-emp-max-vin-2019"
}
if(!dir.exists(data.out.folder)) dir.create(data.out.folder)

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
if(use.vintage.year.pfs) weighting <- paste0(weighting, "_VYP")
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
if(use.vintage.year.pfs) weighting <- paste0(weighting, "_VYP")
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
if(use.vintage.year.pfs) weighting <- paste0(weighting, "_VYP")
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
if(use.vintage.year.pfs) weighting <- paste0(weighting, "_VYP")
cache.folder.tag <- paste0(private.source, ifelse(include.alpha.term, "_alpha_", "_"))
cache.folder.tag <- paste0(cache.folder.tag, weighting)
cache.folder.tag

# RUN
source("estim_model_optimized.R")


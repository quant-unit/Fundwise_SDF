### estimate model runner: Emprical
# 0. Prologue ----
if(sys.nframe() == 0L) rm(list = ls())

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
data.out.folder <- "results/data_out_2026-emp"
factors.to.use <- ""

# 1. Preqin, FW, No CV, L2_Lasso, linear, No Alpha ----

# CHOICES

do.cross.validation <- FALSE

private.source <- "pitchbook"
private.source <- "preqin"

cutoff <- "" # "_cutoff_2019" # only available for PitchBook (2024-09-28)
cutoff <- "_cutoff_2021"

weighting <- "EW"
weighting <- "FW"

#error.function <- "L1_Ridge"
error.function <- "L2_Lasso"

sdf.model <- "linear"
#sdf.model <- "exp.aff"

include.alpha.term <- FALSE

# SET cache.folder.tag
if(use.vintage.year.pfs) weighting <- paste0(weighting, "_VYP")
cache.folder.tag <- paste0(private.source, ifelse(include.alpha.term, "_alpha_", "_"))
cache.folder.tag <- paste0(cache.folder.tag, weighting)
cache.folder.tag

# RUN
source("estim_model_optimized.R")




# 2. Preqin, FW, CV, L2_Lasso, linear, No Alpha ----

# CHOICES

do.cross.validation <- TRUE

private.source <- "pitchbook"
private.source <- "preqin"

cutoff <- "" # "_cutoff_2019" # only available for PitchBook (2024-09-28)
cutoff <- "_cutoff_2021"

weighting <- "EW"
weighting <- "FW"

#error.function <- "L1_Ridge"
error.function <- "L2_Lasso"

sdf.model <- "linear"
#sdf.model <- "exp.aff"

include.alpha.term <- FALSE

# SET cache.folder.tag
if(use.vintage.year.pfs) weighting <- paste0(weighting, "_VYP")
cache.folder.tag <- paste0(private.source, ifelse(include.alpha.term, "_alpha_", "_"))
cache.folder.tag <- paste0(cache.folder.tag, weighting)
cache.folder.tag

# RUN
source("estim_model_optimized.R")


# 3. Preqin, EW, No CV, L2_Lasso, linear, No Alpha ----

# CHOICES

do.cross.validation <- FALSE

private.source <- "pitchbook"
private.source <- "preqin"

cutoff <- "" # "_cutoff_2019" # only available for PitchBook (2024-09-28)
cutoff <- "_cutoff_2021"

weighting <- "EW"
# weighting <- "FW"

#error.function <- "L1_Ridge"
error.function <- "L2_Lasso"

sdf.model <- "linear"
#sdf.model <- "exp.aff"

include.alpha.term <- FALSE

# SET cache.folder.tag
if(use.vintage.year.pfs) weighting <- paste0(weighting, "_VYP")
cache.folder.tag <- paste0(private.source, ifelse(include.alpha.term, "_alpha_", "_"))
cache.folder.tag <- paste0(cache.folder.tag, weighting)
cache.folder.tag

# RUN
source("estim_model_optimized.R")



# 4. Preqin, EW, CV, L2_Lasso, linear, No Alpha ----

# CHOICES

do.cross.validation <- TRUE

private.source <- "pitchbook"
private.source <- "preqin"

cutoff <- "" # "_cutoff_2019" # only available for PitchBook (2024-09-28)
cutoff <- "_cutoff_2021"

weighting <- "EW"
# weighting <- "FW"

#error.function <- "L1_Ridge"
error.function <- "L2_Lasso"

sdf.model <- "linear"
#sdf.model <- "exp.aff"

include.alpha.term <- FALSE

# SET cache.folder.tag
if(use.vintage.year.pfs) weighting <- paste0(weighting, "_VYP")
cache.folder.tag <- paste0(private.source, ifelse(include.alpha.term, "_alpha_", "_"))
cache.folder.tag <- paste0(cache.folder.tag, weighting)
cache.folder.tag

# RUN
source("estim_model_optimized.R")



# 5. Preqin, FW, No CV, L2_Lasso, linear, Alpha ----

# CHOICES

do.cross.validation <- FALSE

private.source <- "pitchbook"
private.source <- "preqin"

cutoff <- "" # "_cutoff_2019" # only available for PitchBook (2024-09-28)
cutoff <- "_cutoff_2021"

weighting <- "EW"
weighting <- "FW"

#error.function <- "L1_Ridge"
error.function <- "L2_Lasso"

sdf.model <- "linear"
#sdf.model <- "exp.aff"

include.alpha.term <- TRUE

# SET cache.folder.tag
if(use.vintage.year.pfs) weighting <- paste0(weighting, "_VYP")
cache.folder.tag <- paste0(private.source, ifelse(include.alpha.term, "_alpha_", "_"))
cache.folder.tag <- paste0(cache.folder.tag, weighting)
cache.folder.tag

# RUN
source("estim_model_optimized.R")




# 6. Preqin, FW, CV, L2_Lasso, linear, Alpha ----

# CHOICES

do.cross.validation <- TRUE

private.source <- "pitchbook"
private.source <- "preqin"

cutoff <- "" # "_cutoff_2019" # only available for PitchBook (2024-09-28)
cutoff <- "_cutoff_2021"

weighting <- "EW"
weighting <- "FW"

#error.function <- "L1_Ridge"
error.function <- "L2_Lasso"

sdf.model <- "linear"
#sdf.model <- "exp.aff"

include.alpha.term <- TRUE

# SET cache.folder.tag
if(use.vintage.year.pfs) weighting <- paste0(weighting, "_VYP")
cache.folder.tag <- paste0(private.source, ifelse(include.alpha.term, "_alpha_", "_"))
cache.folder.tag <- paste0(cache.folder.tag, weighting)
cache.folder.tag

# RUN
source("estim_model_optimized.R")


# 7. Preqin, EW, No CV, L2_Lasso, linear, Alpha ----

# CHOICES

do.cross.validation <- FALSE

private.source <- "pitchbook"
private.source <- "preqin"

cutoff <- "" # "_cutoff_2019" # only available for PitchBook (2024-09-28)
cutoff <- "_cutoff_2021"

weighting <- "EW"
# weighting <- "FW"

#error.function <- "L1_Ridge"
error.function <- "L2_Lasso"

sdf.model <- "linear"
#sdf.model <- "exp.aff"

include.alpha.term <- TRUE

# SET cache.folder.tag
if(use.vintage.year.pfs) weighting <- paste0(weighting, "_VYP")
cache.folder.tag <- paste0(private.source, ifelse(include.alpha.term, "_alpha_", "_"))
cache.folder.tag <- paste0(cache.folder.tag, weighting)
cache.folder.tag

# RUN
source("estim_model_optimized.R")



# 8. Preqin, EW, CV, L2_Lasso, linear, Alpha ----

# CHOICES

do.cross.validation <- TRUE

private.source <- "pitchbook"
private.source <- "preqin"

cutoff <- "" # "_cutoff_2019" # only available for PitchBook (2024-09-28)
cutoff <- "_cutoff_2021"

weighting <- "EW"
# weighting <- "FW"

#error.function <- "L1_Ridge"
error.function <- "L2_Lasso"

sdf.model <- "linear"
#sdf.model <- "exp.aff"

include.alpha.term <- TRUE

# SET cache.folder.tag
if(use.vintage.year.pfs) weighting <- paste0(weighting, "_VYP")
cache.folder.tag <- paste0(private.source, ifelse(include.alpha.term, "_alpha_", "_"))
cache.folder.tag <- paste0(cache.folder.tag, weighting)
cache.folder.tag

# RUN
source("estim_model_optimized.R")



# 9. Preqin, FW, No CV, L2_Lasso, linear, No Alpha, Single-FP ----

# CHOICES

use.vintage.year.pfs <- FALSE

do.cross.validation <- FALSE

private.source <- "pitchbook"
private.source <- "preqin"

cutoff <- "" # "_cutoff_2019" # only available for PitchBook (2024-09-28)
cutoff <- "_cutoff_2021"

weighting <- "EW"
weighting <- "FW"

#error.function <- "L1_Ridge"
error.function <- "L2_Lasso"

sdf.model <- "linear"
#sdf.model <- "exp.aff"

include.alpha.term <- FALSE

# SET cache.folder.tag
if(use.vintage.year.pfs) weighting <- paste0(weighting, "_VYP")
cache.folder.tag <- paste0(private.source, ifelse(include.alpha.term, "_alpha_", "_"))
cache.folder.tag <- paste0(cache.folder.tag, weighting)
cache.folder.tag

# RUN
source("estim_model_optimized.R")





# 10. Preqin, EW, No CV, L2_Lasso, linear, No Alpha, Single-FP ----

# CHOICES

use.vintage.year.pfs <- FALSE

do.cross.validation <- FALSE

private.source <- "pitchbook"
private.source <- "preqin"

cutoff <- "" # "_cutoff_2019" # only available for PitchBook (2024-09-28)
cutoff <- "_cutoff_2021"

weighting <- "EW"
# weighting <- "FW"

#error.function <- "L1_Ridge"
error.function <- "L2_Lasso"

sdf.model <- "linear"
#sdf.model <- "exp.aff"

include.alpha.term <- FALSE

# SET cache.folder.tag
if(use.vintage.year.pfs) weighting <- paste0(weighting, "_VYP")
cache.folder.tag <- paste0(private.source, ifelse(include.alpha.term, "_alpha_", "_"))
cache.folder.tag <- paste0(cache.folder.tag, weighting)
cache.folder.tag

# RUN
source("estim_model_optimized.R")





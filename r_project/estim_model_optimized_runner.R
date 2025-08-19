### estimate model runner
# 0. Prologue ----
if(sys.nframe() == 0L) rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

# 1. run: base case: VYP ----

## set parameters
export.data <- FALSE

use.vintage.year.pfs <- TRUE
use.simulation <- TRUE
do.cross.validation <- FALSE
do.cache <- TRUE
do.parallel <- ifelse(.Platform$OS.type == "windows", FALSE, TRUE)

private.source <- "pitchbook"
cutoff <- "_cutoff_2021"
public.filename <- "q_factors"

# CHOICES
weighting <- "EW"
error.function <- "L2_Lasso"

sdf.model <- "linear"

max.months <- c(1, 60, 120, 150, 180, 210, 240, 300, 360) # c(10, 20) * 12 # c(12.5, 15, 17.5) * 12

include.alpha.term <- FALSE
lambdas <- 0
kernel.bandwidth <- 12
if(use.vintage.year.pfs) weighting <- paste0(weighting, "_VYP")

cache.folder.tag <- "20250808_222540_simulated_cashflows_EW_VYP"
simulation.filename <- paste0("20250808_222540/", cache.folder.tag, ".csv")

part.to.keep <- 1
no.partitions <- 1

data.out.folder <- "data_out_2025-2"

factors.to.use <- ""

## trigger run

# source("estim_model_optimized.R")

# 1.a) r: base case: VYP - discount to duration date ----

## set parameters
cache.folder.tag <- "20250808_222540_simulated_cashflows_EW_VYP"
simulation.filename <- paste0("20250808_222540/", cache.folder.tag, ".csv")

max.months <- c(0, 10000)

sdf.model <- "linear.duration"

part.to.keep <- 1
no.partitions <- 1

## trigger run
# source("estim_model_optimized.R")

# 1.b) r: base case: VYP - discount to single date ----

## set parameters
cache.folder.tag <- "20250808_222540_simulated_cashflows_EW_VYP"
simulation.filename <- paste0("20250808_222540/", cache.folder.tag, ".csv")

max.months <- c(1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110,
                120, 150, 180, 210, 240, 300, 360)

sdf.model <- "linear.single.date"

part.to.keep <- 1
no.partitions <- 1

## trigger run
source("estim_model_optimized.R")

# 2. run: base case: cross-sectional unit ----

cache.folder.tag <- "20250808_222540_simulated_cashflows_EW"
simulation.filename <- paste0("20250808_222540/", cache.folder.tag, ".csv")

max.months <- c(1, 60, 120, 150, 180, 210, 240, 300, 360) # c(10, 20) * 12 # c(12.5, 15, 17.5) * 12

sdf.model <- "linear"

no.partitions <- 10

for (i in 1:10) {
  ## set parameters
  part.to.keep <- i
  ## trigger run
  # source("estim_model_optimized.R")
}

# 3. run: Shorter fund lifetime max(investing=5 + max.holding=5) ----

## set parameters
cache.folder.tag <- "20250808_223358_simulated_cashflows_EW_VYP"
simulation.filename <- paste0("20250808_223358/", cache.folder.tag, ".csv")

part.to.keep <- 1
no.partitions <- 1

## trigger run
# source("estim_model_optimized.R")

# 4. run: High Beta, Negative Alpha ----

## set parameters
cache.folder.tag <- "20250808_224228_simulated_cashflows_EW_VYP"
simulation.filename <- paste0("20250808_224228/", cache.folder.tag, ".csv")

part.to.keep <- 1
no.partitions <- 1

factors.to.use <- "Alpha"

## trigger run
# source("estim_model_optimized.R")

# 5. run: Exp aff models ----

## set parameters
cache.folder.tag <- "20250808_225102_simulated_cashflows_EW_VYP"
simulation.filename <- paste0("20250808_225102/", cache.folder.tag, ".csv")

sdf.model <- "exp.aff"

part.to.keep <- 1
no.partitions <- 1

factors.to.use <- "MKT"


## trigger run
# source("estim_model_optimized.R")

# 6. run: exp.aff, negative alpha, high beta ----

## set parameters
cache.folder.tag <- "20250808_225920_simulated_cashflows_EW_VYP"
simulation.filename <- paste0("20250808_225920/", cache.folder.tag, ".csv")

sdf.model <- "exp.aff"

part.to.keep <- 1
no.partitions <- 1

factors.to.use <- "Alpha"

## trigger run
# source("estim_model_optimized.R")

# 7. run: Big n/V: 40 Funds, vintages 1986-2005 ----

## set parameters
cache.folder.tag <- "20250808_231619_simulated_cashflows_EW_VYP"
simulation.filename <- paste0("20250808_231619/", cache.folder.tag, ".csv")

sdf.model <- "linear"

part.to.keep <- 1
no.partitions <- 1

factors.to.use <- "MKT"

## trigger run
# source("estim_model_optimized.R")

# 8. run: Big V: 10 Funds, vintages 1967-2005 ----

## set parameters
cache.folder.tag <- "20250808_232503_simulated_cashflows_EW_VYP"
simulation.filename <- paste0("20250808_232503/", cache.folder.tag, ".csv")

sdf.model <- "linear"
part.to.keep <- 1
no.partitions <- 1
factors.to.use <- ""

## trigger run
# source("estim_model_optimized.R")

# 9. run: Big V: 20 Funds, vintages 1967-2005 ----

## set parameters
cache.folder.tag <- "20250813_231430_simulated_cashflows_EW_VYP"
simulation.filename <- paste0("20250813_231430/", cache.folder.tag, ".csv")

sdf.model <- "linear"
part.to.keep <- 1
no.partitions <- 1
factors.to.use <- ""

## trigger run
# source("estim_model_optimized.R")

# 10. run: Small V: vintages 1986-1995 ----

## set parameters
cache.folder.tag <- "20250808_233803_simulated_cashflows_EW_VYP"
simulation.filename <- paste0("20250808_233803/", cache.folder.tag, ".csv")

sdf.model <- "linear"
part.to.keep <- 1
no.partitions <- 1
factors.to.use <- ""

## trigger run
# source("estim_model_optimized.R")

# 11. run: Small V: vintages 1996-2005 ----

## set parameters
cache.folder.tag <- "20250808_234201_simulated_cashflows_EW_VYP"
simulation.filename <- paste0("20250808_234201/", cache.folder.tag, ".csv")

sdf.model <- "linear"
part.to.keep <- 1
no.partitions <- 1
factors.to.use <- ""

## trigger run
# source("estim_model_optimized.R")

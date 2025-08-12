### estimate model runner
# 1. run ----

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

cache.folder.tag <- "20250808_222540_simulated_cashflows_EW"
simulation.filename <- paste0("20250808_222540/", cache.folder.tag, ".csv")

part.to.keep <- 1
no.partitions <- 10

data.out.folder <- "data_out_2025"

## trigger run

# source(estim_model_optimized.R)

# 2. run ----

for (i in 2:10) {
  ## set parameters
  part.to.keep <- i
  ## trigger run
  source(estim_model_optimized.R)
}




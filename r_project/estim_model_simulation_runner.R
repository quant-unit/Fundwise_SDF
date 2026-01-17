### estimate model runner: Simulation
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

data.out.folder <- "simulation/data_out_2026_new"

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
# source("estim_model_optimized.R")

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

max.months <- c(1, 30, 60, 90, 120, 150, 180, 210, 240, 300, 360) # c(10, 20) * 12 # c(12.5, 15, 17.5) * 12

part.to.keep <- 1
no.partitions <- 1

## trigger run
# source("estim_model_optimized.R")

# 3.a) r: Shorter fund lifetime max(investing=2 + max.holding=2) ----

## set parameters
cache.folder.tag <- "20250819_215442_simulated_cashflows_EW_VYP"
simulation.filename <- paste0("20250819_215442/", cache.folder.tag, ".csv")

part.to.keep <- 1
no.partitions <- 1

## trigger run
# source("estim_model_optimized.R")

# 3.b) r: Shorter fund lifetime max(investing=4 + max.holding=4) ----

## set parameters
cache.folder.tag <- "20250819_220526_simulated_cashflows_EW_VYP"
simulation.filename <- paste0("20250819_220526/", cache.folder.tag, ".csv")

part.to.keep <- 1
no.partitions <- 1

## trigger run
# source("estim_model_optimized.R")

# 3.c) r: Shorter fund lifetime max(investing=2 + max.holding=4) ----

## set parameters
cache.folder.tag <- "20250819_221615_simulated_cashflows_EW_VYP"
simulation.filename <- paste0("20250819_221615/", cache.folder.tag, ".csv")

part.to.keep <- 1
no.partitions <- 1

## trigger run
# source("estim_model_optimized.R")

# 3.d) r: Shorter fund lifetime max(investing=4 + max.holding=2) ----

## set parameters
cache.folder.tag <- "20250819_222708_simulated_cashflows_EW_VYP"
simulation.filename <- paste0("20250819_222708/", cache.folder.tag, ".csv")

part.to.keep <- 1
no.partitions <- 1

## trigger run
# source("estim_model_optimized.R")

# 3.e) r: Shorter fund lifetime max(investing=1 + max.holding=5) ----

## set parameters
cache.folder.tag <- "20250820_135255_simulated_cashflows_EW_VYP"
simulation.filename <- paste0("20250820_135255/", cache.folder.tag, ".csv")

part.to.keep <- 1
no.partitions <- 1

## trigger run
# source("estim_model_optimized.R")

# 3.f) r: Shorter fund lifetime max(investing=5 + max.holding=1) ----

## set parameters
cache.folder.tag <- "20250820_140352_simulated_cashflows_EW_VYP"
simulation.filename <- paste0("20250820_140352/", cache.folder.tag, ".csv")

part.to.keep <- 1
no.partitions <- 1

## trigger run
# source("estim_model_optimized.R")

# 3.g) r: Shorter fund lifetime max(investing=3 + max.holding=3) ----

## set parameters
cache.folder.tag <- "20250820_192412_simulated_cashflows_EW_VYP"
simulation.filename <- paste0("20250820_192412/", cache.folder.tag, ".csv")

part.to.keep <- 1
no.partitions <- 1

## trigger run
# source("estim_model_optimized.R")

# 3.h) r: Shorter fund lifetime max(investing=1 + max.holding=1) ----

## set parameters
cache.folder.tag <- "20250821_153327_simulated_cashflows_EW_VYP"
simulation.filename <- paste0("20250821_153327/", cache.folder.tag, ".csv")

part.to.keep <- 1
no.partitions <- 1

## trigger run
# source("estim_model_optimized.R")

# 3.i) r: Shorter fund lifetime max(investing=6 + max.holding=6) ----

## set parameters
cache.folder.tag <- "20250821_154353_simulated_cashflows_EW_VYP"
simulation.filename <- paste0("20250821_154353/", cache.folder.tag, ".csv")

part.to.keep <- 1
no.partitions <- 1

## trigger run
# source("estim_model_optimized.R")

# 3.j) r: Shorter fund lifetime max(investing=7 + max.holding=7) ----

## set parameters
cache.folder.tag <- "20250821_155431_simulated_cashflows_EW_VYP"
simulation.filename <- paste0("20250821_155431/", cache.folder.tag, ".csv")

part.to.keep <- 1
no.partitions <- 1

## trigger run
# source("estim_model_optimized.R")

# 3.k) r: Shorter fund lifetime max(investing=8 + max.holding=8) ----

## set parameters
cache.folder.tag <- "20250821_160511_simulated_cashflows_EW_VYP"
simulation.filename <- paste0("20250821_160511/", cache.folder.tag, ".csv")

part.to.keep <- 1
no.partitions <- 1

## trigger run
# source("estim_model_optimized.R")

# 3.l) r: Shorter fund lifetime max(investing=9 + max.holding=9) ----

## set parameters
cache.folder.tag <- "20250821_161553_simulated_cashflows_EW_VYP"
simulation.filename <- paste0("20250821_161553/", cache.folder.tag, ".csv")

part.to.keep <- 1
no.partitions <- 1

## trigger run
# source("estim_model_optimized.R")

# 12) No investment period (investing=0) -----

### 12.a) No investment period  (investing=0 + max.holding=1) 

## set parameters
cache.folder.tag <- "20250825_143009_simulated_cashflows_EW_VYP"
simulation.filename <- paste0("20250825_143009/", cache.folder.tag, ".csv")

part.to.keep <- 1
no.partitions <- 1

## trigger run
# source("estim_model_optimized.R")

### 12.b) No investment period  (investing=0 + max.holding=3) 

## set parameters
cache.folder.tag <- "20250825_143734_simulated_cashflows_EW_VYP"
simulation.filename <- paste0("20250825_143734/", cache.folder.tag, ".csv")

part.to.keep <- 1
no.partitions <- 1

## trigger run
# source("estim_model_optimized.R")

### 12.c) No investment period  (investing=0 + max.holding=5) 

## set parameters
cache.folder.tag <- "20250825_144508_simulated_cashflows_EW_VYP"
simulation.filename <- paste0("20250825_144508/", cache.folder.tag, ".csv")

part.to.keep <- 1
no.partitions <- 1

## trigger run
# source("estim_model_optimized.R")

### 12.d) No investment period  (investing=0 + max.holding=7) 

## set parameters
cache.folder.tag <- "20250825_145244_simulated_cashflows_EW_VYP"
simulation.filename <- paste0("20250825_145244/", cache.folder.tag, ".csv")

part.to.keep <- 1
no.partitions <- 1

## trigger run
# source("estim_model_optimized.R")

### 12.e) No investment period  (investing=0 + max.holding=10) 

## set parameters
cache.folder.tag <- "20250825_150023_simulated_cashflows_EW_VYP"
simulation.filename <- paste0("20250825_150023/", cache.folder.tag, ".csv")

part.to.keep <- 1
no.partitions <- 1

## trigger run
# source("estim_model_optimized.R")

# 13) No investment period (investing=0) - alpha & beta -----

### 13.a) No investment period  (investing=0 + max.holding=1) 

## set parameters
cache.folder.tag <- "20250825_143009_simulated_cashflows_EW_VYP"
simulation.filename <- paste0("20250825_143009/", cache.folder.tag, ".csv")

factors.to.use <- "Alpha"

part.to.keep <- 1
no.partitions <- 1

## trigger run
source("estim_model_optimized.R")

### 13.b) No investment period  (investing=0 + max.holding=3) 

## set parameters
cache.folder.tag <- "20250825_143734_simulated_cashflows_EW_VYP"
simulation.filename <- paste0("20250825_143734/", cache.folder.tag, ".csv")

part.to.keep <- 1
no.partitions <- 1

## trigger run
source("estim_model_optimized.R")

### 13.c) No investment period  (investing=0 + max.holding=5) 

## set parameters
cache.folder.tag <- "20250825_144508_simulated_cashflows_EW_VYP"
simulation.filename <- paste0("20250825_144508/", cache.folder.tag, ".csv")

part.to.keep <- 1
no.partitions <- 1

## trigger run
source("estim_model_optimized.R")

### 13.d) No investment period  (investing=0 + max.holding=7) 

## set parameters
cache.folder.tag <- "20250825_145244_simulated_cashflows_EW_VYP"
simulation.filename <- paste0("20250825_145244/", cache.folder.tag, ".csv")

part.to.keep <- 1
no.partitions <- 1

## trigger run
source("estim_model_optimized.R")

### 13.e) No investment period  (investing=0 + max.holding=10) 

## set parameters
cache.folder.tag <- "20250825_150023_simulated_cashflows_EW_VYP"
simulation.filename <- paste0("20250825_150023/", cache.folder.tag, ".csv")

part.to.keep <- 1
no.partitions <- 1

## trigger run
source("estim_model_optimized.R")

# 4. run: High Beta, Negative Alpha ----

## set parameters
cache.folder.tag <- "20250808_224228_simulated_cashflows_EW_VYP"
simulation.filename <- paste0("20250808_224228/", cache.folder.tag, ".csv")

max.months <- c(1, 60, 120, 150, 180, 210, 240, 300, 360) # c(10, 20) * 12 # c(12.5, 15, 17.5) * 12

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

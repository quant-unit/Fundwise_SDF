###################
### Compare to IRRs
###################
# 0) load data ------
if(sys.nframe() == 0L) rm(list = ls())

source("load_data.R")

# 1) calc IRRs ------

unique(df0$Fund.ID[df0$type == "BO"])

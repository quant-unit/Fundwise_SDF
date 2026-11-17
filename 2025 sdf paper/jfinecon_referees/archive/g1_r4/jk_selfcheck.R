# G1 rerun under the new inference (17 November 2026), see the execution log. Run from r_project/ with LANG=en_US.UTF-8:
#   Rscript "../2025 sdf paper/jfinecon_referees/archive/g1_r4/jk_selfcheck.R" 
suppressMessages(source("helper/result_cache.R")); invisible(capture.output(suppressMessages(source("helper/getNPVs.R"))))
root <- "results/data_out_2026-emp-c"
r <- Reduce(rbind.all.columns, lapply(list.dirs(root, recursive = FALSE), function(d) { x <- read.result.cache(d, verbose = FALSE); x$folder <- basename(d); x }))
j <- r[is.finite(r$jk.selfcheck), ]
j$npar <- ifelse(j$Factor %in% c("MKT"), 1, ifelse(j$Factor == "MKT+SMB+HML", 3, ifelse(j$Factor == "ALL", 5, 2)))
j$method <- sub("[0-9]+$", "", j$X)
cat("rows by method and self-check band:\n")
print(table(j$method, cut(j$jk.selfcheck, c(-1, 1e-6, 1e-4, 1e-3, 1e-2, Inf))))
cat("\nBO/VC only, self-check > 1e-3:\n")
b <- j[j$Type %in% c("BO", "VC") & j$jk.selfcheck > 1e-3, ]
print(table(b$folder, b$Factor))

# Calibration builds for the r5 pipeline (17 November 2026): conditional q = 1 Omega
# quantiles for the bandwidths of the S = ceil(1.3 sqrt V) rule, V = 2..80, and the
# joint-test quantile table with 1e6 draws (q = 1..4).
source("helper/parallel_workers.R"); source("helper/fixedb.R")
V <- 2:80; M <- sort(unique(vapply(V, function(v) fixedb.q1.M(hac.bandwidth.llsw(v) / v), 0L)))
cat(format(Sys.time()), "q1 precompute:", length(M), "bandwidths\n")
fixedb.q1.precompute(M, cores = n.parallel.workers(length(M), gb.per.worker = 1.5))
cat(format(Sys.time()), "q1 done; table rebuild 1e6\n")
tmp <- "helper/fixedb_quantiles.rds.tmp"
fixedb.build(reps = 1e6, file = tmp)
file.rename(tmp, FIXEDB.FILE)
cat(format(Sys.time()), "table done\n")

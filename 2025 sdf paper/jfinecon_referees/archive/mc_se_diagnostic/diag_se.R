# Step 4 SE diagnostic (plan log "Step 4, SE diagnostic", 19 November 2026; archived here,
# reproducible from the production draws in about 2 minutes): why does U_{v,H} (q_v = 1)
# over-reject in the equal-count designs but not in the Buyout core? Run from r_project/:
#   LANG=en_US.UTF-8 Rscript <this file> <out.rds>
# For every replication and H in {0, 15 years}, on the production draws
# (same seeds), at the true beta = 1:
#   s_v  = d l_v / d beta, the score of vintage v (l_v = w_v U_v, the
#          sandwich's unit loss), h_v = d2 l_v / d beta2;
# plus the production estimate and SE (rows). Nothing is written to the repo.
source("simulation/mc_driver.R")
library(parallel)
out.file <- commandArgs(TRUE)[1]
pub <- mc.public.factors(); D <- mc.load.designs()
est <- MC.ESTIMATORS$U.eq
designs <- list(core_BO = 1:1000, w_V25_n15 = 1:500, h_V25_n20 = 1:500)
H <- c(0, 180)

units.of <- function(fp, max.month) {
  w <- u.unit.weights(fp$group.n, est$vintage.weights)
  keep <- names(w)[!is.na(w)]
  units <- lapply(keep, function(k) fast.subset(fp, fp$key == k)); names(units) <- keep
  loss <- function(p, u) w[[u$key[1]]] * fast.u.criterion(fast.unit.errors(p, u, max.month), u, vintage.weights = est$vintage.weights)
  list(units = units, loss = loss, vint = vapply(units, function(u) u$vintage[1], 0))
}

one <- function(id, seed) {
  dr <- mc.draw(D[[id]]$design, seed, pub)
  fp <- mc.as.fast(dr, pub, "fund")
  do.call(rbind, lapply(H, function(m) {
    u <- units.of(fp, m); p0 <- c(MKT = 1)
    s <- vapply(u$units, function(x) num.grad(function(p) u$loss(p, x), p0, MC.SETTINGS$delta), 0)
    h <- vapply(u$units, function(x) num.hess(function(p) u$loss(p, x), p0, MC.SETTINGS$delta)[1, 1], 0)
    data.frame(design = id, rep = seed, max.month = m, vint = u$vint, s = s, h = h)
  }))
}

w <- n.parallel.workers(48)
res <- list()
for (id in names(designs)) {
  t0 <- Sys.time()
  res[[id]] <- do.call(rbind, mclapply(designs[[id]], function(k) one(id, k), mc.cores = w, mc.preschedule = FALSE))
  cat(format(Sys.time(), "%H:%M:%S"), id, nrow(res[[id]]), "rows,", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "s\n")
}
scores <- do.call(rbind, res); rownames(scores) <- NULL
rows <- do.call(rbind, lapply(names(designs), function(id) {
  r <- mc.read(id); r <- r[r$estimator == "U.eq" & r$model == "MKT" & r$max.month %in% H & r$rep %in% designs[[id]], ]
  data.frame(design = id, r[, c("rep", "max.month", "MKT", "SE.MKT", "hac.S", "hac.V", "hac.b", "at.bound")])
}))
# sanity: the production SE recomputed at the production estimate, first 3 replications per design and H
san <- do.call(rbind, lapply(names(designs), function(id) do.call(rbind, lapply(1:3, function(k) {
  dr <- mc.draw(D[[id]]$design, k, pub); fp <- mc.as.fast(dr, pub, "fund")
  do.call(rbind, lapply(H, function(m) { r <- rows[rows$design == id & rows$rep == k & rows$max.month == m, ]
    data.frame(design = id, rep = k, max.month = m, SE.row = r$SE.MKT,
               SE.again = unname(mc.sandwich(est, fp, c(MKT = r$MKT), m)$SE["SE.MKT"])) }))
}))))
saveRDS(list(scores = scores, rows = rows, sanity = san, commit = system("git rev-parse --short HEAD", intern = TRUE)), out.file)
print(san); cat("done\n")

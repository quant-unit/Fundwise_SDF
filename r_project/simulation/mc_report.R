# Step 4 Monte Carlo: summary of all production designs with the criterion
# targets (simulation/mc_analyze.R). Run from r_project/ after
# simulation/mc_production.R and simulation/mc_targets.R:
#   LANG=en_US.UTF-8 Rscript simulation/mc_report.R
# Writes simulation/data_out_mc/summary.rds and results/mc_step4/summary.csv
# (one row per design, estimator, model, parameter and horizon; simulated
# results only, no fund-level data).
source("simulation/mc_driver.R"); source("simulation/mc_analyze.R")
D <- mc.load.designs()
tg <- readRDS("simulation/data_out_mc/targets.rds")
targets <- tg$targets
s <- do.call(rbind, lapply(names(D), function(id) {
  r <- mc.read(id)
  x <- mc.summary(r, D[[id]]$design$theta0, targets[targets$design == id & targets$half == "all", ])
  cbind(design = id, x)
}))
# Monte Carlo precision of the core targets: half the difference between the
# targets on replications 1-100 and 101-200
h <- targets[targets$half != "all" & targets$method != "theta_0 (T7)", ]
if (nrow(h)) {
  w <- reshape(h[, c("design", "estimator", "model", "max.month", "half", "MKT.target", "Coef.target")],
               idvar = c("design", "estimator", "model", "max.month"), timevar = "half", direction = "wide")
  attr(s, "target.precision") <- data.frame(w[, c("design", "estimator", "model", "max.month")],
    MKT = abs(w$`MKT.target.1-100` - w$`MKT.target.101-200`) / 2, Coef = abs(w$`Coef.target.1-100` - w$`Coef.target.101-200`) / 2)
}
attr(s, "provenance") <- list(targets.commit = tg$commit, targets.dirty = tg$dirty,
  production = unique(unlist(lapply(names(D), function(id) vapply(attr(mc.read(id), "provenance"), `[[`, "", "git.commit")))))
saveRDS(s, "simulation/data_out_mc/summary.rds")
dir.create("results/mc_step4", showWarnings = FALSE, recursive = TRUE)
write.csv(s, "results/mc_step4/summary.csv", row.names = FALSE)
cat("summary:", nrow(s), "rows,", length(unique(s$design)), "designs; targets from commit", tg$commit, if (isTRUE(tg$dirty)) "(dirty)" else "(clean)", "\n")

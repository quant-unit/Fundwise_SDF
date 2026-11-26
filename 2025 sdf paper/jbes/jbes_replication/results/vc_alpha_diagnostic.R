# VC (alpha, beta) diagnostic: full sample vs without the 1993-1997 vintages.
# Regenerates results/figures_bovc_us_g1/vc_alpha_without_1993_1997.csv on the
# current data_prepared_2026_c pools (read by results/main_exhibits_u.R).
# Run from r_project/:  LANG=en_US.UTF-8 Rscript results/vc_alpha_diagnostic.R
source("simulation/mc_driver.R")
out.dir <- "results/figures_bovc_us_g1"
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
H.E2 <- c(0, 30, 60, 120, 180, 240, 360)

pub <- mc.public.factors()
fc <- unique(c("RF", setdiff(colnames(pub), "Date")))

load.vc <- function(w, vyp) {
  d <- read.csv(sprintf("empirical/data_prepared_2026_c/preqin_cashflows_%s%s_North America_NAV.csv", w, if (vyp) "_VYP" else ""), stringsAsFactors = FALSE)
  d$Date <- as.Date(d$Date)
  d$Fund.ID <- as.factor(paste(d$Fund.ID, d$type, sep = "_"))
  d <- pad.all.units(d, "first.cf", max(d$Date))
  d <- merge(d, pub, by = "Date", all.x = TRUE)
  d <- d[d$type == "VC" & d$Vintage <= 2010, ]
  d$Fund.ID <- droplevels(as.factor(as.character(d$Fund.ID)))
  prepare.fast(d, fc)
}

rows <- list()
for (est.name in c("L2", "U")) {
  for (w in c("EW", "FW")) {
    fp.full <- load.vc(w, vyp = (est.name == "L2"))
    fp.sub <- fast.subset(fp.full, !(fp.full$vintage %in% 1993:1997))
    est <- if (est.name == "L2") MC.ESTIMATORS$L2.vyp else if (w == "FW") MC.ESTIMATORS$U.eq.fw else MC.ESTIMATORS$U.eq
    # n.vintages header convention in main_exhibits_u.R uses L2 vintage counts (27 full, 22 without 1993-1997)
    n.vin.full <- length(unique(fp.full$vintage))
    n.vin.sub <- length(unique(fp.sub$vintage))
    for (smp in c("full", "without 1993-1997")) {
      fp <- if (smp == "full") fp.full else fp.sub
      nv <- if (smp == "full") n.vin.full else n.vin.sub
      for (m in H.E2) {
        r <- mc.fit(est, fp, "Alpha", m)
        rows[[length(rows) + 1L]] <- data.frame(
          type = "VC", weighting = w, estimator = est.name, sample = smp,
          n.vintages = nv, max.month = m, H = m / 12,
          beta = r$MKT, alpha = r$Coef, criterion = r$value,
          at.bound = r$at.bound, active = r$active,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}
out <- do.call(rbind, rows)
write.csv(out, file.path(out.dir, "vc_alpha_without_1993_1997.csv"), row.names = FALSE)
cat("wrote", file.path(out.dir, "vc_alpha_without_1993_1997.csv"), "with", nrow(out), "rows\n")

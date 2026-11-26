# Step 4 Monte Carlo: summary of all production designs with the criterion
# targets (simulation/mc_analyze.R). Run from r_project/ after
# simulation/mc_production.R and simulation/mc_targets.R:
#   LANG=en_US.UTF-8 Rscript simulation/mc_report.R [--partial]
# Without --partial, every design in config/mc_designs.yaml must be present
# with all requested replications and matching provenance (spec.md5,
# public.md5, MC.SETTINGS), and every cell must have a matching target. Target
# files are checked against the current specs only if they record provenance
# (mc_targets.R after 22 November 2026); older files raise a warning.
# Writes simulation/data_out_mc/summary.rds and results/mc_step4/summary.csv
# (one row per design, estimator, model, parameter and horizon; simulated
# results only, no fund-level data).
source("simulation/mc_driver.R"); source("simulation/mc_analyze.R")
args <- commandArgs(TRUE)
unknown <- setdiff(args, "--partial")
if (length(unknown)) stop("Unknown argument(s): ", paste(unknown, collapse = ", "))
partial.ok <- "--partial" %in% args

pub <- mc.public.factors()
pub.md5 <- value.md5(pub)
D <- mc.load.designs()

# targets of the first production (targets.rds), of the additions of
# 19 November (targets_additions.rds) and of later additions (other
# targets*.rds, simulation/mc_targets.R --out=); each design's targets must
# come from one file.
target.files <- sort(list.files("simulation/data_out_mc", pattern = "^targets.*\\.rds$", full.names = TRUE))
if (!length(target.files)) stop("No targets*.rds files found in simulation/data_out_mc")
tgs <- lapply(target.files, readRDS)
names(tgs) <- basename(target.files)
stopifnot(!anyDuplicated(unlist(lapply(tgs, function(x) unique(x$targets$design)))))
# Target files written before 22 November 2026 carry no provenance fields; they
# cannot be checked against the current specifications and are reported as such.
unvalidated <- character()
for (fn in names(tgs)) {
  x <- tgs[[fn]]
  if (is.null(x$public.md5) || is.null(x$settings) || is.null(x$spec.md5)) {
    unvalidated <- c(unvalidated, fn)
    warning("Target file ", fn, " has no provenance fields (spec.md5/public.md5/settings): only cell completeness is checked")
  }
  if (!is.null(x$public.md5) && !identical(x$public.md5, pub.md5)) {
    stop("Target file ", fn, " was computed with different public factors")
  }
  if (!is.null(x$settings) && !identical(x$settings, MC.SETTINGS)) {
    stop("Target file ", fn, " was computed with different MC.SETTINGS")
  }
  if (!is.null(x$spec.md5)) {
    for (id in intersect(names(x$spec.md5), names(D))) {
      if (!identical(unname(x$spec.md5[[id]]), mc.spec.md5(D[[id]]))) {
        stop("Target file ", fn, " spec.md5 mismatch for design ", id)
      }
    }
  }
}
tg <- list(commit = paste(unique(vapply(tgs, `[[`, "", "commit")), collapse = ","), dirty = any(vapply(tgs, function(x) isTRUE(x$dirty), TRUE)))
targets <- do.call(rbind, lapply(tgs, `[[`, "targets"))

prod.commits <- character()
s <- do.call(rbind, lapply(names(D), function(id) {
  dir.id <- file.path("simulation/data_out_mc", id)
  has.batches <- length(list.files(dir.id, "^batch_\\d+\\.rds$")) > 0L
  if (!has.batches) {
    if (partial.ok) return(NULL)
    stop("Design ", id, " has no finished batches (pass --partial to allow incomplete runs)")
  }
  exp.reps <- if (partial.ok) NULL else seq_len(D[[id]]$reps)
  r <- mc.read(id, expected.spec = D[[id]], expected.public = pub.md5,
               expected.settings = MC.SETTINGS, expected.reps = exp.reps)
  tg.id <- targets[targets$design == id & targets$half == "all", , drop = FALSE]
  exp.cells <- expand.grid(estimator = D[[id]]$estimators, model = D[[id]]$models,
                           max.month = D[[id]]$horizons, stringsAsFactors = FALSE)
  cell.key <- function(df) paste(df$estimator, df$model, df$max.month, sep = "|")
  missing.tg <- setdiff(cell.key(exp.cells), cell.key(tg.id))
  if (length(missing.tg)) {
    stop("Design ", id, ": missing target cell(s), e.g. ", missing.tg[1])
  }
  p <- attr(r, "provenance")
  if (!is.null(p)) prod.commits <<- unique(c(prod.commits, vapply(p, `[[`, "", "git.commit")))
  x <- mc.summary(r, D[[id]]$design$theta0, tg.id)
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
  production = prod.commits, partial = partial.ok, targets.unvalidated = unvalidated)
saveRDS(s, "simulation/data_out_mc/summary.rds")
dir.create("results/mc_step4", showWarnings = FALSE, recursive = TRUE)
write.csv(s, "results/mc_step4/summary.csv", row.names = FALSE)
cat("summary:", nrow(s), "rows,", length(unique(s$design)), "designs; targets from commit", tg$commit, if (isTRUE(tg$dirty)) "(dirty)" else "(clean)", "\n")

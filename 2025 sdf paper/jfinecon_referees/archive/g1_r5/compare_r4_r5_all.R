# G1 rerun under pipeline r5 (17 November 2026), see the execution log. Run from r_project/ with LANG=en_US.UTF-8:
#   Rscript "../2025 sdf paper/jfinecon_referees/archive/g1_r5/compare_r4_r5_all.R" <out.rds>
# r5 rows (admissible domain, constrained multi-start fit) vs r4 rows (pipeline
# label 2026-11-16-r3, kernel.bandwidth "llsw") in every G1 result folder.
# Expected: every r5 row matches an r4 row; a higher r5 objective only where
# the r4 estimate lies outside Theta (domain rebuilt from the row's own input
# and factor files); rows on a bound and jackknife flags listed.
suppressMessages(source("helper/result_cache.R")); invisible(capture.output(suppressMessages(source("helper/getNPVs.R"))))
source("helper/admissible.R")
root <- "results/data_out_2026-emp-c"
pc <- c("MKT", "ME", "IA", "ROE", "EG", "Alpha", "SMB", "HML")
k <- c("Type", "Factor", "max.month", "CV.key", "lambda")
months.cache <- list()
domain.of <- function(row, nm) {
  key <- paste(row$input.path, row$public.path, row$Type)
  if (is.null(months.cache[[key]])) {
    inp <- read.csv(row$input.path); inp$Date <- as.Date(inp$Date)
    pub <- read.csv(row$public.path); pub$Date <- as.Date(format(as.Date(pub$Date) + 1, "%Y-%m-01")) - 1; pub$Alpha <- 1
    d0 <- min(inp$Date[inp$type == row$Type]); d1 <- max(inp$Date)
    months.cache[[key]] <<- pub[pub$Date >= d0 & pub$Date <= d1, ]
  }
  stopifnot(length(row$gross.floor) == 1, is.finite(as.numeric(row$gross.floor)))
  admissible.constraints(months.cache[[key]], nm, as.numeric(row$gross.floor), c(-1, 4), c(-0.02, 0.02))
}
out <- list(); moved <- list()
for (d in list.dirs(root, recursive = FALSE)) {
  raw <- Reduce(rbind.all.columns, lapply(Sys.glob(file.path(d, "*_cached_res.csv")), read.csv))
  r4 <- dedupe.result.rows(raw[raw$pipeline.version %in% "2026-11-16-r3" & as.character(raw$kernel.bandwidth) == "llsw", ])
  r5 <- read.result.cache(d, verbose = FALSE)
  if (!nrow(r5)) { out[[basename(d)]] <- data.frame(folder = basename(d), r4 = nrow(r4), r5 = 0, matched = 0); next }
  m <- merge(r5, r4, by = k, suffixes = c(".n", ".o"))
  cols <- intersect(pc, sub("\\.n$", "", grep("\\.n$", names(m), value = TRUE)))
  m$dmax <- apply(sapply(cols, function(v) abs(m[[paste0(v, ".n")]] - m[[paste0(v, ".o")]])), 1, max, na.rm = TRUE)
  m$dval <- (m$value.n - m$value.o) / abs(m$value.o) # U_pair criteria can be negative
  m$r4.feasible <- NA
  for (i in which(m$dval > 1e-6)) {
    nm <- if (m$Factor[i] == "MKT") "MKT" else if (m$Factor[i] == "ALL") intersect(pc, cols)[is.finite(unlist(m[i, paste0(intersect(pc, cols), ".o")]))] else c("MKT", m$Factor[i])
    p <- setNames(vapply(nm, function(v) m[[paste0(v, ".o")]][i], 0), nm)
    row <- list(input.path = m$input.path.n[i], public.path = m$public.path.n[i], Type = m$Type[i], gross.floor = m$gross.floor.n[i])
    m$r4.feasible[i] <- in.domain(p, domain.of(row, nm))
  }
  m$folder <- basename(d)
  # r5-only fields carry the .n suffix when the r4 rows have the column (filled NA); absent in some folders
  for (v in c("active.constraints", "jk.n.bound", "jk.n.switch", "min.G")) m[[v]] <- if (!is.null(m[[paste0(v, ".n")]])) m[[paste0(v, ".n")]] else if (!is.null(m[[v]])) m[[v]] else NA
  moved[[basename(d)]] <- m[m$dmax > 1e-3 | m$dval > 1e-6 | m$at.bound.n %in% TRUE, c("folder", k, "MKT.o", "MKT.n", "value.o", "value.n", "dval", "dmax", "at.bound.n", "active.constraints", "min.G", "r4.feasible", "jk.n.bound", "jk.n.switch")]
  out[[basename(d)]] <- data.frame(folder = basename(d), r4 = nrow(r4), r5 = nrow(r5), matched = nrow(m),
    within.1e3 = sum(m$dmax <= 1e-3), lower = sum(m$dval < -1e-6), higher = sum(m$dval > 1e-6),
    higher.r4.feasible = sum(m$r4.feasible %in% TRUE), higher.unknown = sum(m$dval > 1e-6 & is.na(m$r4.feasible)),
    on.bound = sum(m$at.bound.n %in% TRUE), dirty = paste(unique(r5$git.dirty), collapse = ","), commit = paste(unique(r5$git.commit), collapse = ","))
}
o <- do.call(rbind, out); rownames(o) <- NULL
mv <- do.call(rbind, moved); rownames(mv) <- NULL
cat(sprintf("folders %d | rows r4 %d, r5 %d, matched %d | within 1e-3: %d | lower objective: %d | higher: %d (r4 inside Theta: %d, undetermined: %d) | on a bound: %d\n",
  nrow(o), sum(o$r4), sum(o$r5), sum(o$matched), sum(o$within.1e3), sum(o$lower), sum(o$higher), sum(o$higher.r4.feasible), sum(o$higher.unknown), sum(o$on.bound)))
cat("commits:", paste(unique(o$commit), collapse = " "), "| git.dirty:", paste(unique(o$dirty), collapse = " "), "\n")
print(o[o$matched != o$r5 | o$higher.r4.feasible > 0 | o$higher.unknown > 0, ], row.names = FALSE)
saveRDS(list(o = o, moved = mv), commandArgs(TRUE)[1])

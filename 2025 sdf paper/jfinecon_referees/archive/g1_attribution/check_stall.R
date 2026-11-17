# Is the paper's Buyout (alpha, beta) estimate a minimum? Objective at the old
# estimate (code d2bdf6e) vs. the current estimate, same pools (_b) and same
# window endpoint (own last cash flow). At H = 0 the criterion changed (B2), so
# only H >= 5 years compare like with like.
# Run from r_project/:  Rscript <this file> <OUT>
out <- commandArgs(TRUE)[1]
rd <- function(cell, w) do.call(rbind, lapply(Sys.glob(file.path(out, cell, sprintf("cache_q_factors_preqin_%s_VYP_North America", w), "*.csv")),
  function(f) read.csv(f)[, c("Type", "Factor", "max.month", "MKT", "Alpha", "value")]))
for (w in c("EW", "FW")) { o <- rd("old_b", w); n <- rd("ownend_b", w)
  for (h in c(60, 120, 180, 240)) { a <- o[o$Type == "BO" & o$Factor == "Alpha" & o$max.month == h, ]; b <- n[n$Type == "BO" & n$Factor == "Alpha" & n$max.month == h, ]
    cat(sprintf("BO %s H=%3d old (beta %.3f, alpha %.4f) obj %.5f | current (beta %.3f, alpha %.4f) obj %.5f | current lower by %.2f%%\n",
      w, h, a$MKT, a$Alpha, a$value, b$MKT, b$Alpha, b$value, 100 * (1 - b$value / a$value))) } }

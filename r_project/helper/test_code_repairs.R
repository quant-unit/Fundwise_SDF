######################
### Automated verification tests for the revision code repairs (revision plan, Section "Code and data repairs")
### Run from r_project/:  Rscript helper/test_code_repairs.R
######################
# Implemented so far:
#   B4  constant-composition VYP pooling  (helper/pool_vintage.R, empirical/prep_preqin.R)
#   B6  FW pools normalized by total vintage commitment
#   Sample rule: main sample = complete-history funds; all funds = robustness
#
# Each test prints PASS/FAIL; the script exits non-zero if any test fails.

suppressPackageStartupMessages(library(here))
source(here("helper/pool_vintage.R"))

n.fail <- 0
check <- function(name, ok, detail = "") {
  cat(sprintf("[%s] %s%s\n", if (isTRUE(ok)) "PASS" else "FAIL", name,
              if (nzchar(detail)) paste0("  (", detail, ")") else ""))
  if (!isTRUE(ok)) n.fail <<- n.fail + 1
}

# NPV of a dated cash-flow vector under a per-date discount function
npv <- function(dates, cf, disc) sum(cf * disc(dates))

# A family of discount curves: constant annual rates (incl. negative) and
# benchmark curves R^b(beta) = prod(1 + RF + beta * MKT) from the q5 factor file.
make.discount.curves <- function(origin) {
  curves <- lapply(c(-0.05, 0, 0.08, 0.25), function(r) {
    function(d) (1 + r)^(-as.numeric(as.Date(d) - origin) / 365.25)
  })
  q.file <- here("empirical/data_prepared_2026/q_factors.csv")
  if (file.exists(q.file)) {
    q <- read.csv(q.file)
    q$Date <- as.Date(q$Date)
    q <- q[q$Date >= origin, ]
    for (beta in c(0.5, 1, 1.5)) {
      cum <- cumprod(1 + q$RF + beta * q$MKT)
      ym <- format(q$Date, "%Y-%m")
      curves[[length(curves) + 1]] <- local({
        cum <- cum; ym <- ym
        function(d) {
          k <- match(format(as.Date(d), "%Y-%m"), ym)
          ifelse(is.na(k), 1, 1 / cum[k])
        }
      })
    }
  }
  curves
}

# ---------------------------------------------------------------------------
# B4/B6 unit tests on a synthetic vintage (no data needed)
# ---------------------------------------------------------------------------
cat("\n== B4/B6: synthetic vintage ==\n")
d <- as.Date(c("2000-03-31", "2000-06-30", "2000-09-30", "2000-12-31", "2001-03-31"))
toy <- rbind(
  # fund A reports on every date
  data.frame(Fund.ID = "A", Date = d, CF = c(-0.5, -0.3, 0.1, 0.4, 0.9), NAV = c(0.5, 0.8, 0.8, 0.5, 0.2), Size = 100),
  # fund B reports on three dates with a gap and liquidates early
  data.frame(Fund.ID = "B", Date = d[c(1, 3, 4)], CF = c(-0.8, 0.2, 1.1), NAV = c(0.8, 0.7, 0.3), Size = 300),
  # fund C starts late
  data.frame(Fund.ID = "C", Date = d[4:5], CF = c(-1.0, 1.6), NAV = c(1.0, 0.1), Size = 50)
)
toy$CON <- pmin(toy$CF, 0)
toy$DIS <- pmax(toy$CF, 0)
toy$Vintage <- 2000
toy$type <- "BO"
ids <- c("A", "B", "C")
sizes <- c(A = 100, B = 300, C = 50)

for (fsw in c(FALSE, TRUE)) {
  lab <- if (fsw) "FW" else "EW"
  f <- toy
  if (fsw) for (col in c("CF", "NAV", "CON", "DIS")) f[[col]] <- f[[col]] * f$Size
  p <- pool.vintage(f, fund.size.weighting = fsw)
  w <- if (fsw) sizes / sum(sizes) else rep(1 / 3, 3)
  names(w) <- ids

  # pooled CF on each date = sum over all funds / time-invariant denominator
  denom <- if (fsw) sum(sizes) else 3
  expect.cf <- sapply(d, function(t) sum(f$CF[f$Date == t])) / denom
  check(paste(lab, "pooled CF uses fixed denominator"), max(abs(p$CF - expect.cf)) < 1e-15)
  check(paste(lab, "N.funds and Commitment columns"),
        all(p$N.funds == 3) && all(p$Commitment == sum(sizes)))

  # exact NPV identity for every discount curve
  curves <- make.discount.curves(as.Date("1999-12-31"))
  err <- max(sapply(curves, function(disc) {
    pv.pool <- npv(p$Date, p$CF, disc)
    pv.funds <- sapply(ids, function(i) {
      fi <- f[f$Fund.ID == i, ]
      npv(fi$Date, fi$CF / if (fsw) sizes[[i]] else 1, disc)
    })
    abs(pv.pool - sum(w * pv.funds))
  }))
  check(paste(lab, "pooled NPV = weighted avg of fund NPVs (", length(curves), "curves)"),
        err < 1e-14, sprintf("max abs err %.1e", err))

  # NAV: carried forward within each fund's life, zero after liquidation
  nav.raw <- toy$NAV * if (fsw) toy$Size else 1
  nav.B.prev <- nav.raw[toy$Fund.ID == "B" & toy$Date == d[1]] # B has no 2000-06-30 report
  exp.nav.q2 <- (nav.raw[toy$Fund.ID == "A" & toy$Date == d[2]] + nav.B.prev) / denom
  exp.nav.q5 <- (nav.raw[toy$Fund.ID == "A" & toy$Date == d[5]] +
                 nav.raw[toy$Fund.ID == "C" & toy$Date == d[5]]) / denom # B liquidated
  check(paste(lab, "NAV carried forward / zero after liquidation"),
        abs(p$NAV[2] - exp.nav.q2) < 1e-15 && abs(p$NAV[5] - exp.nav.q5) < 1e-15)
}

# the old reporter-mean pooling violates the identity on the same data
old <- aggregate(CF ~ Date, data = toy, mean)
disc <- make.discount.curves(as.Date("1999-12-31"))[[3]]
err.old <- abs(npv(old$Date, old$CF, disc) -
               mean(sapply(ids, function(i) npv(toy$Date[toy$Fund.ID == i], toy$CF[toy$Fund.ID == i], disc))))
check("reporter-mean pooling fails the identity (sanity)", err.old > 1e-3, sprintf("err %.3f", err.old))

# ---------------------------------------------------------------------------
# B4/B6 tests on the rebuilt Preqin files (skipped if not yet generated)
# ---------------------------------------------------------------------------
main.dir <- here("empirical/data_prepared_2026_c")
regions <- c("North America", NA, "Europe")
tol.file <- 1e-12 # CSV round trip keeps 15 significant digits

for (data.dir in c(main.dir, paste0(main.dir, "_allfunds"))) {
is.main <- data.dir == main.dir
cat("\n== B4/B6: Preqin files in", data.dir, "==\n")
abs.sums <- list()
for (region in regions) {
  for (lab in c("EW", "FW")) {
    reg.tag <- if (is.na(region)) "" else paste0("_", region)
    f.fund <- file.path(data.dir, paste0("preqin_cashflows_", lab, reg.tag, "_NAV.csv"))
    f.vyp <- file.path(data.dir, paste0("preqin_cashflows_", lab, "_VYP", reg.tag, "_NAV.csv"))
    rlab <- paste(lab, if (is.na(region)) "ALL" else region)
    if (!file.exists(f.fund) || !file.exists(f.vyp)) {
      # a missing file means prep_preqin.R did not finish: fail, don't skip
      check(paste(rlab, "prepared files exist"), FALSE, basename(f.vyp))
      next
    }
    fund <- read.csv(f.fund)
    vyp <- read.csv(f.vyp)
    fund$Date <- as.Date(fund$Date)
    vyp$Date <- as.Date(vyp$Date)
    fsw <- lab == "FW"

    # sample rule: the main sample holds complete-history funds only, and its
    # pools match the funnel table; the all-funds sample adds truncated funds
    if (is.main) {
      check(paste(rlab, "main sample holds complete-history funds only"),
            all(fund$Complete) &&
              # live cutoff = data end (2022-03-31) minus live.window.months (12)
              all(fund$Age.Last > 5 | as.Date(fund$Last.Date) >= as.Date("2021-03-31")))
    }
    f.funnel <- file.path(main.dir, paste0("preqin_sample_funnel", reg.tag, ".csv"))
    if (file.exists(f.funnel)) {
      fn <- read.csv(f.funnel)
      nv <- vyp[!duplicated(paste(vyp$type, vyp$Vintage)), ]
      k <- match(paste(nv$type, nv$Vintage), paste(fn$type, fn$Vintage))
      expect.n <- if (is.main) fn$N.complete[k] else fn$N.complete[k] + fn$N.truncated[k]
      check(paste(rlab, "pool fund counts match the sample funnel"),
            !anyNA(k) && all(nv$N.funds == expect.n) &&
              all(fn$N.raw == fn$N.complete + fn$N.truncated + fn$N.late.start + fn$N.missing.nav))
    }

    max.err <- 0
    n.checked <- 0
    abs.sum <- c()
    for (key in unique(paste(vyp$type, vyp$Vintage))) {
      tp <- sub(" .*", "", key)
      vin <- as.integer(sub(".* ", "", key))
      pv <- vyp[vyp$type == tp & vyp$Vintage == vin, ]
      fv <- fund[fund$type == tp & fund$Vintage == vin, ]
      sz <- fv$Size[!duplicated(fv$Fund.ID)]
      names(sz) <- fv$Fund.ID[!duplicated(fv$Fund.ID)]
      w <- if (fsw) sz / sum(sz) else rep(1 / length(sz), length(sz))
      names(w) <- names(sz)
      abs.sum[key] <- sum(abs(pv$CF))
      curves <- make.discount.curves(as.Date(paste0(vin - 1, "-12-31")))
      for (disc in curves) {
        pv.pool <- npv(pv$Date, pv$CF, disc)
        pv.funds <- sapply(names(sz), function(i) {
          fi <- fv[fv$Fund.ID == i, ]
          npv(fi$Date, fi$CF / if (fsw) sz[[i]] else 1, disc)
        })
        scale <- max(1, abs(pv.pool))
        max.err <- max(max.err, abs(pv.pool - sum(w * pv.funds)) / scale)
        n.checked <- n.checked + 1
      }
      if (pv$N.funds[1] != length(sz)) max.err <- Inf
    }
    check(paste(rlab, "pooled NPV = weighted avg of fund NPVs"),
          max.err < tol.file,
          sprintf("%d vintage x curve checks, max rel err %.1e", n.checked, max.err))

    # B6: sum of |pooled CF| per vintage is O(1) for both EW and FW, on the
    # estimation sample (BO and VC, seasoned vintages <= 2010, n_v >= 2).
    # Unseasoned 2021 vintages and single-fund niche pools legitimately have
    # tiny totals and are excluded.
    nf <- vyp[!duplicated(paste(vyp$type, vyp$Vintage)), ]
    est <- paste(nf$type, nf$Vintage)[nf$type %in% c("BO", "VC") & nf$Vintage <= 2010 & nf$N.funds >= 2]
    est <- intersect(est, names(abs.sum))
    if (length(est) > 0 && is.main) {
      rng <- range(abs.sum[est])
      check(paste(rlab, "sum |pooled CF| per vintage is O(1) (BO/VC, V<=2010, n_v>=2)"),
            rng[1] > 0.3 && rng[2] < 20,
            sprintf("%d pools, range %.2f to %.2f", length(est), rng[1], rng[2]))
    }
    abs.sums[[rlab]] <- abs.sum[est]
  }
  # B6: FW and EW pools are on the same per-unit-of-commitment scale. Checked
  # on the main sample only: in the all-funds sample, near-empty truncated
  # funds make some small pools legitimately tiny, so the ratio is unstable.
  ew <- abs.sums[[paste("EW", if (is.na(region)) "ALL" else region)]]
  fw <- abs.sums[[paste("FW", if (is.na(region)) "ALL" else region)]]
  if (is.main && length(ew) > 0 && length(fw) > 0) {
    k <- intersect(names(ew), names(fw))
    ratio <- fw[k] / ew[k]
    check(paste(if (is.na(region)) "ALL" else region, "FW/EW scale ratio per vintage in [0.33, 3]"),
          all(ratio > 1 / 3 & ratio < 3),
          sprintf("range %.2f to %.2f", min(ratio), max(ratio)))
  }
}
}

cat(sprintf("\n%d test(s) failed\n", n.fail))
if (n.fail > 0) quit(status = 1)

#### estimate model
# 0) Prologue -----

if (!exists("source.internally", envir = .GlobalEnv)) {
  source.internally <- TRUE
}

if (source.internally) {
  if (FALSE) {
    files <- c("prep_preqin.R", "prep_public.R")
    for (file in files) {
      source(file)
      if (sys.nframe() == 0L) rm(list = ls())
    }
  }

  if (sys.nframe() == 0L) rm(list = ls())
  source.internally <- TRUE

  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  getwd()
}

# 1.1) PARAMETERS ----
if (source.internally) {
  export.data <- FALSE

  use.vintage.year.pfs <- TRUE
  use.simulation <- FALSE
  do.cross.validation <- FALSE
  do.cache <- TRUE
  do.parallel <- ifelse(.Platform$OS.type == "windows", FALSE, TRUE)

  private.source <- "pitchbook"
  # private.source <- "preqin"
  cutoff <- "" # "_cutoff_2019" # only available for PitchBook (2024-09-28)
  cutoff <- "_cutoff_2021"
  cutoff <- "_cutoff_2025"
  max.vintage <- 2011

  # public.filename <- "public_returns" # outdated
  # public.filename <- "msci_market_factors"
  public.filename <- "q_factors"
  # public.filename <- "DebtFactorsEURUSD" # outdated
  # public.filename <- "iBoxxFactorsMIX"
  # public.filename <- "iBoxxFactorsUSD"
  # public.filename <- "iBoxxFactorsEUR"


  # CHOICES
  weighting <- "EW"
  weighting <- "FW"
  # error.function <- "L1_Ridge"
  error.function <- "L2_Lasso"

  sdf.model <- "linear"
  # sdf.model <- "exp.aff"

  max.months <- c(1, 60, 120, 150, 180, 210, 240, 300, 360) # c(10, 20) * 12 # c(12.5, 15, 17.5) * 12
  max.months <- c(0, 180)

  include.alpha.term <- FALSE
  final.nav.discount <- 100
  lambdas <- 0
  kernel.bandwidth <- "llsw" # S = ceil(1.3 sqrt(V)), fixed-b tests (helper/fixedb.R)
  alpha.lower <- -0.01
  alpha.upper <- 0.01
  if (use.vintage.year.pfs) weighting <- paste0(weighting, "_VYP")
  cache.folder.tag <- paste0(private.source, "_")
  # cache.folder.tag <- paste0(private.source, ifelse(do.cross.validation, "_cv_", "_"))
  cache.folder.tag <- paste0(cache.folder.tag, ifelse(include.alpha.term, "alpha_", ""))
  cache.folder.tag <- paste0(cache.folder.tag, weighting)
  cache.folder.tag

  # cache.folder.tag <- "20250808_222540_simulated_cashflows_EW"
  simulation.filename <- paste0("20250808_222540/", cache.folder.tag, ".csv")

  part.to.keep <- 1
  no.partitions <- 1 # 10

  if (!exists("data.out.folder", envir = .GlobalEnv)) {
    data.out.folder <- "results/data_out_2028"
  }
  if (!dir.exists(data.out.folder)) dir.create(data.out.folder)

  factors.to.use <- ""
  # set by run_estimation(); cleared for direct sourcing
  input.md5 <- NULL
  result.identity.used <- NULL
  config.path <- NULL
  scenario.id <- NULL

  data.prepared.folder <- "data_prepared"

  # private.source <- "preqin"
  # data.prepared.folder <- "data_prepared_2026"
}

# Horizon origin and U-criterion settings (B2/B5, see helper/objective.R).
# run_estimation() always assigns them; these defaults apply when the script
# is sourced directly or from a legacy runner with source.internally = FALSE.
if (!exists("origin.convention")) origin.convention <- "first.cf" # or "predetermined"
if (!exists("vintage.weights")) vintage.weights <- "equal" # or "n.funds" (U_pair only)
if (!exists("config.path")) config.path <- NULL
if (!exists("scenario.id")) scenario.id <- NULL
source("helper/objective.R")
source("helper/parallel_workers.R")
# Admissible parameter domain (Assumption A7) and constrained multi-start fit
source("helper/admissible.R")
if (!exists("gross.floor")) gross.floor <- GROSS.FLOOR

# B9: Brent bracket for one-parameter (MKT) models, and provenance recorded in
# every result row (git commit, uncommitted changes, config file hash, scenario)
if (!exists("brent.interval")) brent.interval <- c(-1, 4)
provenance <- list(
  git.commit = tryCatch(system("git rev-parse --short HEAD", intern = TRUE, ignore.stderr = TRUE)[1], error = function(e) NA),
  git.dirty = tryCatch(length(system("git status --porcelain --untracked-files=no", intern = TRUE, ignore.stderr = TRUE)) > 0, error = function(e) NA),
  # hash of the config file actually used (run_estimation(config_path = ...));
  # the default file only when none was given
  config.md5 = unname(tools::md5sum(
    if (exists("config.path") && !is.null(config.path)) {
      config.path
    } else if (exists("use.simulation") && use.simulation) "config/simulation_scenarios.yaml" else "config/empirical_scenarios.yaml"
  )),
  config.path = if (exists("config.path") && !is.null(config.path)) config.path else NA,
  scenario.id = if (exists("scenario.id") && !is.null(scenario.id)) scenario.id else NA
)
source("helper/result_cache.R")

# 1.2) load data -----

# load public data
if (public.filename %in% c("q_factors", "ff3_factors")) {
  df.public <- read.csv(public.file.path(public.filename, data.prepared.folder))
} else {
  df.public <- read.csv2(paste0("empirical/", data.prepared.folder, "/", public.filename, ".csv"))
}
print(paste(public.filename, "max public date:", max(as.Date(df.public$Date))))

colnames(df.public) <- gsub("_World", "", colnames(df.public))
df.public$Date <- as.Date(df.public$Date)

# Normalize dates to last calendar day of month.
# FF3 data from Kenneth French uses last-business-day dates (e.g., Apr 29),
# while private cashflow data uses end-of-calendar-month (e.g., Apr 30).
# Without this, the merge by Date fails for ~30% of rows.
df.public$Date <- as.Date(format(df.public$Date, "%Y-%m-01")) # first of month
df.public$Date <- seq.Date(df.public$Date[1], by = "month", length.out = nrow(df.public))
# Convert to end-of-month: shift to first of next month, then subtract 1 day
df.public$Date <- as.Date(format(df.public$Date + 32, "%Y-%m-01")) - 1
bond.files <- c(
  "DebtFactorsUSD", "DebtFactorsEUR", "DebtFactorsEURUSD",
  "iBoxxFactorsEUR", "iBoxxFactorsUSD", "iBoxxFactorsMIX"
)
if (public.filename %in% bond.files) {
  df.public <- df.public[complete.cases(df.public[, 2:5]), ]
  df.public[is.na(df.public)] <- 0
  public.equity <- "msci_market_factors"
  df.pubeq <- read.csv2(paste0(data.prepared.folder, "/", public.equity, ".csv"))
  colnames(df.pubeq) <- gsub("_World", "", colnames(df.pubeq))
  df.pubeq$Date <- as.Date(df.pubeq$Date)

  df.public$RF <- NULL
  df.public <- merge(df.public, df.pubeq[, c("Date", "RF", "MKT")], by = "Date", all.x = TRUE)
  remove(df.pubeq)
  df.public <- df.public[df.public$Date > as.Date("2000-01-31"), ]
}

df.public$Alpha <- 1

if (export.data) {
  write.csv2(df.public, paste0("empirical/data_private_public/public_", public.filename, ".csv"))
}

# Load private data (path and hash: helper/result_cache.R, shared with the
# resume check in run_estimation())
private.file <- input.file.path(
  use.simulation, if (use.simulation) simulation.filename else NULL, private.source,
  data.prepared.folder, cutoff, weighting, if (exists("final.nav.discount")) final.nav.discount else 100
)
cat("Loading private data:", private.file, "\n")
df.private.cfs <- read.csv(private.file)
input.md5.used <- file.md5(private.file)
if (exists("input.md5") && !is.null(input.md5) && !identical(input.md5, input.md5.used)) {
  stop("Input file changed between the resume check and loading: ", private.file)
}

# split too large data.frames into partitions
if (no.partitions > 1) {
  data.table::setDT(df.private.cfs)
  lev <- unique(df.private.cfs$Fund.ID)
  part <- data.table::data.table(Fund.ID = lev, partition = as.integer(dplyr::ntile(seq_along(lev), no.partitions)))
  df.private.cfs <- as.data.frame(part[df.private.cfs, on = "Fund.ID"][partition == part.to.keep])
  # Only add partition suffix if not already present (prevents duplication when tag is pre-computed)
  partition_suffix <- paste0("_part", part.to.keep)
  if (!grepl(partition_suffix, cache.folder.tag, fixed = TRUE)) {
    cache.folder.tag <- paste0(cache.folder.tag, partition_suffix)
  }
  rm(lev)
  rm(part)
  df.private.cfs$partition <- NULL
}


df.private.cfs$Date <- as.Date(df.private.cfs$Date)
df.private.cfs$Fund.ID <- as.factor(paste(df.private.cfs$Fund.ID, df.private.cfs$type, sep = "_"))


# Monthly grid with zero cash flows (B5, assumption A5): every unit of a
# vintage starts at the common vintage origin t0 (origin.convention) and every
# unit ends at the common data end T.max, in the empirical and the simulation
# branch alike. Formerly each unit ran from its own first to its own last date
# (empirical) or to 2019-12-31 (simulation), so the horizon window depended on
# each vintage's realized liquidation date.
grid.t.max <- max(df.private.cfs$Date)
cat("Monthly grid: origin =", origin.convention, "| common end T.max =", format(grid.t.max), "\n")
system.time(df.private.cfs <- pad.all.units(df.private.cfs, origin.convention, grid.t.max))

stopifnot(sum(df.private.cfs$CF) > 0)


df.private.cfs$type <- as.factor(as.character(df.private.cfs$type))
df.private.cfs$Fund.ID <- as.factor(as.character(df.private.cfs$Fund.ID))
length(levels(df.private.cfs$type))



# merge private and public data
df0 <- base::merge(df.private.cfs, df.public, by = "Date", all.x = TRUE)
df0$Fund.ID <- as.factor(df0$Fund.ID)
rm(df.private.cfs)

# Override MKT (and RF) with simulated factors if available (simulate.mkt = TRUE in DGP)
if (use.simulation) {
  sim_factors_path <- sub(
    "_simulated_cashflows_.*\\.csv$", "_simulated_factors.csv",
    paste0("simulation/data_prepared_sim/", simulation.filename)
  )
  if (file.exists(sim_factors_path)) {
    cat("Loading simulated factors from:", basename(sim_factors_path), "\n")
    df.sim.factors <- read.csv(sim_factors_path)
    df.sim.factors$Date <- as.Date(df.sim.factors$Date)

    # Extract sample_id from type string (number after last ".")
    df0$sample_id <- as.integer(sub(".*\\.(\\d+)$", "\\1", as.character(df0$type)))

    # Determine which columns to override (MKT always; RF if available)
    sim_cols <- intersect(c("MKT", "RF"), colnames(df.sim.factors))
    merge_cols <- c("Date", "sample_id", sim_cols)

    # Merge simulated factors by (Date, sample_id)
    df0 <- merge(df0, df.sim.factors[, merge_cols],
      by = c("Date", "sample_id"), all.x = TRUE, suffixes = c(".hist", ".sim")
    )

    # Override MKT: use simulated where available, fall back to historical
    if ("MKT.sim" %in% colnames(df0)) {
      df0$MKT <- ifelse(is.na(df0$MKT.sim), df0$MKT.hist, df0$MKT.sim)
      df0$MKT.sim <- NULL
      df0$MKT.hist <- NULL
    }

    # Override RF: use simulated where available, fall back to historical
    if ("RF.sim" %in% colnames(df0)) {
      df0$RF <- ifelse(is.na(df0$RF.sim), df0$RF.hist, df0$RF.sim)
      df0$RF.sim <- NULL
      df0$RF.hist <- NULL
      cat("  -> RF column overridden with constant simulated value\n")
    }

    df0$sample_id <- NULL
    rm(df.sim.factors)
    cat("  -> MKT column overridden with simulated values per sample\n")
  }
}

# Recheck all loaded dependencies and resolved settings, not only cash flows.
loaded.identity <- make.estimation.identity(private.file, public.filename, data.prepared.folder,
  mget(RESULT.SETTINGS, envir = .GlobalEnv), config.path, use.simulation)
if (exists("result.identity.used") && !is.null(result.identity.used) &&
    !identical(loaded.identity$cache.identity, result.identity.used$cache.identity)) {
  stop("Inputs or estimator settings changed between the resume check and loading")
}
result.identity.used <- loaded.identity

table(df0$Date)

# check if we have enough public data
min(df0$Vintage)
min(df0$Date)

if (public.filename %in% c("DebtFactorsUSD", "DebtFactorsEUR", "DebtFactorsEURUSD")) {
  df0 <- df0[df0$Vintage > 2000, ]
}
df0 <- df0[df0$Vintage <= max.vintage, ]


# Export df0
if (export.data) {
  write.csv2(df0, paste0(
    "data_private_public/df0_",
    private.source, "_",
    public.filename, "_",
    weighting,
    cutoff, ".csv"
  ),
  row.names = FALSE
  )
}


# Summary statistics for df0
aggregate(Vintage ~ type, df0, min)
number.of <- function(x) length(table(x))
aggregate(Vintage ~ type, df0, number.of)
aggregate(as.character(Fund.ID) ~ type, df0, number.of)
table(df0$Vintage)

if (!use.simulation) {
  # print summary: funds per vintage
  df1 <- df0
  df1 <- df1[!duplicated(df1$Fund.ID), ]
  df1 <- as.data.frame.matrix(table(df1$Vintage, df1$type))
  # df1 <- df1[, colnames(df1) %in% c("BO", "DD", "INF", "MEZZ", "NATRES", "PD", "RE", "VC")]
  df1["Total", ] <- as.integer(colSums(df1))
  print(xtable::xtable(df1,
    caption = "Number of funds per vintage year.",
    label = paste0("tab:", private.source, "_data")
  ), include.rownames = TRUE)
  rm(df1)
}

# 2.1) getNPVs function ----
source("helper/getNPVs.R")

# 2.2) err.sqr.calc function ----

# B2: max.month = 0 is an alias of max.month = 1 (inception NPV at the origin
# t0, origin date not discounted); see helper/objective.R for the convention.
if (sdf.model == "linear") {
  f1 <- function(df.ss, max.month, par0) {
    return(horizon.error(
      df.ss$CF,
      exp(cumsum(log(1 + (as.matrix(df.ss[, names(par0)]) %*% par0)))),
      max.month
    ))
  }
}

if (sdf.model == "linear.duration") {
  f1 <- function(df.ss, max.month, par0) {
    return(getNPVsDuration(
      df.ss$CF,
      exp(cumsum(log(1 + (as.matrix(df.ss[, names(par0)]) %*% par0)))),
      max.month
    ))
  }
}

if (sdf.model == "linear.single.date") {
  f1 <- function(df.ss, max.month, par0) {
    return(getNPVsSingleDate(
      df.ss$CF,
      exp(cumsum(log(1 + (as.matrix(df.ss[, names(par0)]) %*% par0)))),
      max.month
    ))
  }
}

if (sdf.model == "exp.aff") {
  f1 <- function(df.ss, max.month, par0) {
    df.ss$Alpha <- exp(1) - 1

    return(horizon.error(
      df.ss$CF,
      exp(cumsum(as.matrix(log(1 + df.ss[, names(par0)])) %*% par0)),
      max.month
    ))
  }
}

# Speed: the objective works on prepared units (helper/objective.R): cash
# flows and factor columns of every unit are extracted once per data set
# instead of re-splitting the data frame at every evaluation. f1.prep() is the
# same computation as f1() above on the prepared matrices (identical results,
# tested in helper/test_code_repairs.R); f1() is kept as the reference.
factor.cols.all <- unique(c("RF", setdiff(colnames(df.public), "Date")))
as.prepared <- function(df) prepare.units(df, factor.cols.all)

# Admissible domain (Assumption A7, helper/admissible.R): the months of each
# type's full-sample window with their factor values. The gross-return floor
# applies in each of these months in the full-sample, CV and jackknife fits
# alike. Taken from df0, so simulated factors (simulation branch) count.
admissible.months <- lapply(split(df0[, c("Date", factor.cols.all)], as.character(df0$type)), function(d) {
  d <- d[!duplicated(d$Date), , drop = FALSE]
  d[order(d$Date), , drop = FALSE]
})
theta.domain <- function(type, par.names) {
  admissible.constraints(admissible.months[[type]], par.names, gross.floor, brent.interval, c(alpha.lower, alpha.upper))
}
f1.prep <- function(u, max.month, par0) unit.horizon.error(u, max.month, par0, sdf.model)
unit.errors <- function(par, df, max.month) {
  prep <- as.prepared(df)
  vapply(prep$units, f1.prep, 0, max.month = max.month, par0 = c("RF" = 1, par))
}

if (error.function == "L2_Lasso") {
  err.sqr.calc <- function(par, max.month, lambda, df) {
    npvs <- unit.errors(par, df, max.month)
    # Guard against NA/NaN/Inf in NPV computation
    # NA arises from missing factor data after merge; NaN from log(1+x) when x < -1
    if (any(!is.finite(npvs))) {
      return(1e10) # Large but finite penalty so optimizer can continue
    }
    p <- if (length(par) == 1) par else par[-1] # penalized coefficients (ifelse() kept only the first)
    return(
      sqrt(sum(npvs^2)) / length(npvs) + lambda / length(p) * sum(abs(p))
    )
  }
}

if (error.function == "L1_Ridge") {
  err.sqr.calc <- function(par, max.month, lambda, df) {
    npvs <- unit.errors(par, df, max.month)
    # Guard against NA/NaN/Inf in NPV computation
    if (any(!is.finite(npvs))) {
      return(1e10) # Large but finite penalty so optimizer can continue
    }
    p <- if (length(par) == 1) par else par[-1] # penalized coefficients (ifelse() kept only the first)
    return(
      sum(abs(npvs)) / length(npvs) + lambda / length(p) * sum(p^2)
    )
  }
}

# Cross-fund pair-product criterion U_{v,H} (B5, theorem T7): needs fund-level
# data. Funds of a vintage share the grid (pad.all.units), so their horizon
# errors share A_{v,H}(theta). Vintage weights q_v: vintage.weights.
if (error.function == "U_pair") {
  if (use.vintage.year.pfs) stop("error.function = 'U_pair' needs fund-level data (use_vintage_year_pfs: false)")
  use.size <- grepl("FW", weighting) && "Size" %in% colnames(df0)
  err.sqr.calc <- function(par, max.month, lambda, df) {
    prep <- as.prepared(df)
    npvs <- unit.errors(par, prep, max.month)
    if (any(!is.finite(npvs))) {
      return(1e10)
    }
    key <- vapply(prep$units, `[[`, "", "key")
    s <- if (use.size) vapply(prep$units, `[[`, 0, "size") else rep(1, length(npvs))
    u <- u.criterion(npvs, key, s, vintage.weights)
    if (!is.finite(u)) {
      return(1e10)
    }
    p <- if (length(par) == 1) par else par[-1] # penalized coefficients (ifelse() kept only the first)
    return(u + lambda / length(p) * sum(abs(p)))
  }
}

# Cross-sectional unit for standard errors (score / Hessian contributions):
# the fund or vintage pool for L2/L1 losses, the vintage for U_pair.
unit.split <- function(df) {
  if (error.function == "U_pair") split(df, paste(df$type, df$Vintage)) else split(df, df$Fund.ID)
}

max.month <- 180
lambda <- 0
par <- c(MKT = 1, Alpha = 0)
type <- levels(df0$type)[1]
df.in <- df0[df0$type == type, ]

# df.test <- df.in[df.in$Fund.ID == df.in$Fund.ID[1], ]
# plot(df.test$Date, cumsum(df.test$CF))
# plot(df.test$Date, exp(cumsum(log(1 + (as.matrix(df.test[, names(par)]) %*% par)))))
# f1(df.test, 180, par)


df.in$Fund.ID <- (as.character(df.in$Fund.ID))
system.time(
  print(err.sqr.calc(par, df = df.in, lambda = lambda, max.month = max.month))
)


# 2.3) Asymptotic gradient hessian -----
if (FALSE) {
  # Build bounds & choose method for Alpha at max.month == 1
  use.bounded <- ("Alpha" %in% names(par))
  if (use.bounded) {
    lower_bounds <- ifelse(names(par) == "Alpha", alpha.lower, -Inf)
    upper_bounds <- ifelse(names(par) == "Alpha", alpha.upper, Inf)
    optim_method <- "L-BFGS-B"
  } else {
    lower_bounds <- -Inf
    upper_bounds <- Inf
    optim_method <- "Nelder-Mead"
  }
  res <- tryCatch(
    {
      res_tmp <- optimx::optimx(par, err.sqr.calc,
        lambda = lambda,
        max.month = max.month,
        df = df.in,
        lower = lower_bounds,
        upper = upper_bounds,
        method = optim_method
      )
      if (nrow(res_tmp) == 0) stop("Optimx returned 0 rows in Block 1")
      res_tmp
    },
    error = function(e) {
      warning("optimx failed in Block 1 (", e$message, "). Falling back to native optimizers.")
      if (use.bounded) {
        ans <- stats::nlminb(par, err.sqr.calc,
          lower = lower_bounds, upper = upper_bounds,
          lambda = lambda, max.month = max.month, df = df.in
        )
        res_df <- as.data.frame(as.list(ans$par))
        res_df$value <- ans$objective
        res_df$fevals <- ans$evaluations[1]
        res_df$gevals <- ans$evaluations[2]
        res_df$niter <- ans$iterations
      } else {
        ans <- stats::optim(par, err.sqr.calc,
          method = "Nelder-Mead",
          lambda = lambda, max.month = max.month, df = df.in
        )
        res_df <- as.data.frame(as.list(ans$par))
        res_df$value <- ans$value
        res_df$fevals <- ans$counts[1]
        res_df$gevals <- ans$counts[2]
        res_df$niter <- NA
      }
      res_df$convcode <- ans$convergence
      res_df$kkt1 <- NA
      res_df$kkt2 <- NA
      res_df$xtime <- NA
      return(res_df)
    }
  )

  print("--- Block 1 Optimx Result ---")
  print(Class = class(res))
  print(paste("Dim:", dim(res)))
  print(res)

  names.par <- names(par)
  par <- as.numeric(res[1, names.par])
  names(par) <- names.par
}

# Test Data
dfx <- unit.split(df.in)
delta <- 0.0001

# Inference (B3, helper/inference.R). Unit loss l_i(theta) entering
# Q_n = mean_i l_i: squared horizon error (L2), absolute error (L1), or the
# vintage's U_{v,H} (U_pair). The horizon max.month is passed explicitly:
# the former SE functions read max.month as a free variable and therefore
# used the global test value 180 for every estimated horizon.
source("helper/inference.R")
# For U_pair the unit contribution is l_v = w_v U_v with the normalized
# vintage weights w_v of the estimation criterion (u.unit.weights, shared with
# u.criterion), and vintages with fewer than two funds are not units; before
# this fix the sandwich used unweighted U_v, wrong for vintage_weights = "n.funds".
# Units arrive prepared (as.prepared, speed): one prepared object per unit,
# holding one pool/fund (L2, L1) or the funds of one vintage (U_pair).
unit.loss.at <- function(max.month, unit.weights = NULL) {
  function(par, unit) {
    if (error.function == "U_pair") {
      w <- unit.weights[[unit$units[[1]]$key]]
      return(w * err.sqr.calc(par, max.month, 0, unit))
    }
    e <- f1.prep(unit$units[[1]], max.month, c("RF" = 1, par))
    if (error.function == "L1_Ridge") return(abs(e))
    e^2
  }
}

estim.CoVM <- function(dfx, par, indep = FALSE, max.month) {
  unit.weights <- NULL
  if (error.function == "U_pair") {
    n.funds <- vapply(dfx, function(d) length(unique(d$Fund.ID)), 0)
    names(n.funds) <- vapply(dfx, function(d) paste(d$type[1], d$Vintage[1]), "")
    unit.weights <- u.unit.weights(n.funds, vintage.weights)
    dfx <- dfx[!is.na(unit.weights)]
    unit.weights <- unit.weights[!is.na(unit.weights)]
  }
  sandwich.covariance(
    units = lapply(dfx, as.prepared), par = par, unit.loss = unit.loss.at(max.month, unit.weights),
    vintages = vapply(dfx, function(d) as.numeric(d$Vintage[1]), 0),
    bw = kernel.bandwidth, indep = indep, delta = delta
  )
}
# Leave-one-unit-out jackknife (VYP: one unit per vintage). Every refit uses
# the estimator's domain Theta (con) and fitting rule: Brent on the feasible
# interval for one parameter; otherwise starts at the full-sample estimate and
# at the distinct local optima of the full-sample search (endpoints), so that a
# refit moving to another basin is detected (refit.admissible). Returns the SEs,
# the number of refits on a bound of Theta and of refits that left the basin of
# the full-sample estimate, and selfcheck = largest change when the full sample
# is refitted this way.
jackknife.se <- function(par.hat, prep, max.month, lambda, con, endpoints, parscale, smooth = TRUE) {
  nm <- names(par.hat)
  na.out <- list(se = as.list(setNames(rep(NA_real_, length(par.hat)), nm)), selfcheck = NA_real_, n = length(prep$units), n.bound = NA, n.switch = NA)
  refit <- function(units) {
    d <- structure(list(units = units), class = "prepared.units")
    f <- function(p) err.sqr.calc(p, max.month, lambda, d)
    if (length(par.hat) == 1) {
      br <- feasible.interval(con)
      o <- stats::optim(par.hat, function(p) f(setNames(p, nm)), method = "Brent", lower = br[1], upper = br[2])
      p <- setNames(o$par, nm)
      return(list(par = p, at.bound = length(active.constraints(p, con)) > 0, switched = FALSE))
    }
    refit.admissible(f, par.hat, endpoints, con, parscale, smooth)
  }
  n <- length(prep$units)
  if (n < 3) return(na.out)
  full <- refit(prep$units)
  fits <- lapply(seq_len(n), function(i) refit(prep$units[-i]))
  if (is.null(full) || any(vapply(fits, is.null, TRUE))) return(na.out)
  est <- t(vapply(fits, function(r) as.numeric(r$par), numeric(length(par.hat))))
  if (length(par.hat) == 1) est <- matrix(est, ncol = 1)
  colnames(est) <- nm
  se <- sqrt((n - 1) / n * colSums(sweep(est, 2, colMeans(est))^2))
  list(se = as.list(se), selfcheck = max(abs(full$par - par.hat)), n = n,
       n.bound = sum(vapply(fits, `[[`, TRUE, "at.bound")), n.switch = sum(vapply(fits, `[[`, TRUE, "switched")))
}

# interactive check only (speed: skipped in runner calls, result unused)
if (source.internally) {
  system.time(
    out <- estim.CoVM(dfx, par, indep = FALSE, max.month = max.month)
  )
  system.time(
    out <- estim.CoVM(dfx, par, indep = TRUE, max.month = max.month)
  )
  out
}



# 3.1) cross validation -------

vintage.blocks <- list()

if (do.cross.validation) {
  # hv-block folds: 3-year validation windows buffered by 3-year blocks on each side.
  # For max_vintage <= 2017 (including the main max_vintage = 2010 sample), active
  # folds and training sets match the 1988..2015 / 1980:2020 schedule identically;
  # using df0's actual vintages avoids dropping post-2020 vintages if max_vintage > 2020.
  all.vintages <- sort(unique(df0$Vintage))
  cv.max.year <- max(2015L, max(all.vintages) - 2L)
  for (year in seq(1988L, cv.max.year, by = 3L)) {
    validate <- year:(year + 2L)
    block <- c((year - 3L):(year - 1L), (year + 3L):(year + 5L))

    estimate <- setdiff(all.vintages, c(validate, block))

    vintage.blocks[[paste(year)]] <- list(validate = validate, block = block, estimate = estimate)
  }

  for (year in names(vintage.blocks)) {
    for (key in c("estimate", "validate")) {
      vintage.blocks[[year]][[paste0("df.", key)]] <- df0[df0$Vintage %in% vintage.blocks[[year]][[key]], ]
      vintage.blocks[[year]][["cv.year"]] <- year
    }
  }

  viblo.list <- list()
  for (block in vintage.blocks) {
    viblo.list[[paste(block$validate[1])]] <- data.frame(
      training.before = paste0("start-", block$block[1] - 1),
      h.block.before = paste(block$block[1:3], collapse = ","),
      validation = paste(block$validate, collapse = ","),
      h.block.after = paste(block$block[4:6], collapse = ","),
      training.after = paste0(block$block[6] + 1, "-end")
    )
  }
  print(xtable::xtable(data.frame(do.call(rbind, viblo.list))), include.rownames = FALSE)
}

vintage.blocks[["ALL"]] <- list(df.estimate = df0, df.validate = df0, cv.year = "ALL")

names(vintage.blocks)

# 3.2) iterative run ----
library(optimx)
library(parallel)
library(doParallel)
library(foreach)

iter.run <- function(input.list) {
  # 1) INITIALIZE
  df.estimate <- input.list$df.estimate
  df.validate <- input.list$df.validate

  if (public.filename == "msci_market_factors") {
    factors <- c("HML", "SMB", "HDY", "QLT", "MOM", "ESG", "LOV")
    types <- levels(df0$type)
  }
  if (public.filename == "q_factors") {
    factors <- colnames(df.public)[2:6] # q_factors
    # factors <- "MKT"
    # factors <- "Alpha"
    factors <- c(factors, "Alpha")

    types <- c("PE", "VC", "PD", "RE", "NATRES", "INF", "BO", "GroBO") # asset classes
    if (use.simulation) types <- levels(df0$type)
  }
  if (public.filename == "ff3_factors") {
    factors <- c("SMB", "HML", "Alpha")
    types <- c("PE", "VC", "PD", "RE", "NATRES", "INF", "BO", "GroBO") # asset classes
    if (use.simulation) types <- levels(df0$type)
  }
  if (public.filename %in% bond.files) {
    factors <- c("TERM", "CORP", "HY", "LIQ")
    types <- levels(df0$type)
    # types <- c("PD", "DD", "MEZZ")
  }

  # Save all available factors before override (needed for "ALL" combined model)
  all_available_factors <- factors

  if (length(factors.to.use) == 1 && factors.to.use == "ALL") {
    # "ALL" = single combined model with all available factors
    factors <- "ALL"
  } else if (length(factors.to.use) > 1) {
    # Multiple explicit factors means we want one model combining them
    factors <- paste(factors.to.use, collapse = "+")
  } else if (!identical(factors.to.use, "")) {
    factors <- factors.to.use
  }

  # 2) RUNNER
  # Speed: types (asset classes; simulation samples) are estimated in
  # parallel on forked workers. The cluster used to be created but the
  # loop ran with %do%, i.e. serially. Results are identical (the
  # estimation is deterministic; .combine keeps the order of types).
  # One cluster per call (fold), not per horizon: stopping a fork cluster
  # and forking a new one at once fails intermittently with "error writing
  # to connection" (G1: 21 of 40 cycles in a standalone test). foreach
  # exports this function's small variables (max.month, lambda, ...) on
  # every call; nothing here assigns globals, so workers see the same state.
  # Memory: each task receives only its own type's rows (est.by.type,
  # val.by.type); the full data sets are not exported to every worker, and
  # the number of workers is capped by free memory (helper/parallel_workers.R).
  cv.key <- input.list$cv.year
  split.types <- function(df) split(df, factor(as.character(df$type), levels = types))[types]
  est.by.type <- split.types(df.estimate)
  # full-sample fit: validation data = estimation data, split only once
  val.by.type <- if (identical(df.validate, df.estimate)) est.by.type else split.types(df.validate)
  n.workers <- if (isTRUE(do.parallel)) n.parallel.workers(length(types)) else 1L
  use.parallel <- isTRUE(do.parallel) && n.workers > 1L
  if (use.parallel) {
    cl <- parallel::makeForkCluster(n.workers)
    on.exit(parallel::stopCluster(cl), add = TRUE) # also after an error
    doParallel::registerDoParallel(cl)
  }
  `%op%` <- if (use.parallel) foreach::`%dopar%` else foreach::`%do%`

  l <- list()
  for (lambda in lambdas) {
    for (max.month in max.months) {
      print(paste("Max.M", max.month, "Lambda", lambda, " - ", Sys.time()))

      output <- foreach(type = types, df.optim.in = est.by.type, df.val = val.by.type, .combine = "rbind",
                        .noexport = c("input.list", "df.estimate", "df.validate", "est.by.type", "val.by.type"),
                        .verbose = FALSE) %op% {
        df.optim.in$Fund.ID <- (as.character(df.optim.in$Fund.ID))
        df.val$Fund.ID <- (as.character(df.val$Fund.ID))
        # prepared once per type and fold, reused by every objective call
        if (nrow(df.optim.in) > 0) df.optim.prep <- as.prepared(df.optim.in)
        if (nrow(df.val) > 0) df.val.prep <- as.prepared(df.val)

        factor.list <- list()
        for (factor in factors) {
          if (nrow(df.optim.in) == 0 | nrow(df.val) == 0) next

          if (!include.alpha.term) {
            # two factor model
            par <- c("MKT" = 1)
          } else {
            # three factor model
            par <- c("MKT" = 1, "Alpha" = 0)
          }

          if (factor == "ALL") {
            for (fac in setdiff(all_available_factors, "ALL")) {
              par[fac] <- 0
            }
          } else if (length(factors.to.use) > 1) {
            # Add all requested multiple factors to the parameter vector
            for (fac in factors.to.use) {
              if (fac != "MKT" && fac != "Alpha" && !(fac %in% names(par))) {
                par[fac] <- 0
              }
            }
          } else if (factor == "Alpha") {
            par <- c("MKT" = 1, "Alpha" = 0)
          } else {
            if ("ALL" %in% factors) break
            if (!(factor %in% names(par))) par[factor] <- 0
          }
          # print(type)
          # print(par)

          if (length(par) == 1) {
            # B9: one-parameter models use Brent's method on the feasible
            # interval of Theta (brent.interval and the gross-return floor,
            # assumption A7, helper/admissible.R)
            con <- theta.domain(type, names(par))
            br <- feasible.interval(con)
            # optim(method = "Brent") passes an unnamed scalar: restore the name
            par.name <- names(par)
            ans <- stats::optim(par, function(p, ...) err.sqr.calc(setNames(p, par.name), ...),
              method = "Brent", lower = br[1], upper = br[2],
              lambda = lambda, max.month = max.month, df = df.optim.prep
            )
            res <- as.data.frame(as.list(setNames(ans$par, names(par))))
            res$value <- ans$value
            res$fevals <- ans$counts[1]
            res$gevals <- NA
            res$niter <- NA
            res$convcode <- ans$convergence
            res$kkt1 <- NA
            res$kkt2 <- NA
            res$xtime <- NA
            rownames(res) <- "Brent"
          } else {
            # Constrained multi-start fit on Theta (assumption A7; decision of
            # 17 November 2026, helper/admissible.R): fixed start set around the
            # pipeline start, Nelder-Mead, polished along the boundary when a
            # constraint is within reach. Replaces optimx (Nelder-Mead from one
            # start; L-BFGS-B when Alpha is estimated), which stopped in a local
            # minimum and had no domain beyond the alpha bounds. Parameter scales:
            # optim.parscale() (alpha ~ 0.001, loadings ~ 1).
            con <- theta.domain(type, names(par))
            par.scale <- optim.parscale(par)
            fit <- fit.admissible(function(p) err.sqr.calc(p, max.month, lambda, df.optim.prep),
              admissible.starts(par, con, par.scale), con, par.scale, smooth = error.function != "L1_Ridge")
            res <- as.data.frame(as.list(fit$par))
            res$value <- fit$value
            res$fevals <- NA
            res$gevals <- NA
            res$niter <- NA
            res$convcode <- fit$convergence
            res$kkt1 <- NA
            res$kkt2 <- NA
            res$xtime <- NA
            res$fit.method <- fit$method
            res$n.starts <- fit$n.starts
            rownames(res) <- fit$method
          }

          # Standard Error calc
          for (v in names(par)) par[v] <- res[1, v]
          # estimate on a bound of Theta (beta bracket, alpha bounds or the
          # gross-return floor): SEs and tests not valid; tables print a dagger
          act <- active.constraints(par, con)
          res$at.bound <- length(act) > 0
          res$active.constraints <- paste(act, collapse = ",")
          mg <- min.gross.return(par, con)
          res$min.G <- mg$value
          res$min.G.date <- format(mg$date)
          res$gross.floor <- gross.floor
          dfx <- unit.split(df.optim.in)

          # Inference (decisions of 17 November 2026). L1: none, because no
          # valid inference is established for the L1 criterion (its estimate
          # sits at a kink; State of evidence, sections E and I); the SE columns
          # are NA. Otherwise the Bartlett sandwich with S = ceil(1.3 sqrt V):
          # the rows store the Wald statistic of beta = 1 (others 0), q and
          # b = S / V; the fixed-b p-value is computed at analysis time
          # (results/analyze_result.R), so a recalibration of the fixed-b
          # distribution never invalidates estimation results. The chi-square
          # p-value does not depend on that calibration and stays in the row.
          # neff.<par>: effective number of units carrying each coefficient
          # (influence contributions, helper/inference.R).
          if (error.function == "L1_Ridge") {
            for (v in names(par)) res[, paste0("SE.", v)] <- res[, paste0("SE.", v, ".indep")] <- NA_real_
            res$inference <- "none (L1)"
          } else {
            cov <- estim.CoVM(dfx, par = par, indep = FALSE, max.month = max.month)
            for (v in names(cov$SE)) res[, v] <- cov$SE[v]
            res$Wald.stat.MKT_1 <- cov$Wald.stat
            res$Wald.q <- length(par)
            res$Wald.p.value.MKT_1.chisq <- cov$Wald.p.value.MKT_1.chisq
            res$hac.S <- cov$hac.S
            res$hac.V <- cov$hac.V
            res$hac.b <- cov$hac.b
            res$neff.units <- cov$neff.units # = neff.MKT
            for (v in names(par)) res[, paste0("neff.", v)] <- cov$neff[[v]]
            res$se.singular <- isTRUE(cov$se.singular) # Hessian not invertible: SEs NA

            # independent units (no kernel)
            cov.indep <- estim.CoVM(dfx, par = par, indep = TRUE, max.month = max.month)
            for (v in names(cov.indep$SE)) res[, paste0(v, ".indep")] <- cov.indep$SE[v]
            res$inference <- "sandwich"

            # Leave-one-vintage-out jackknife SE (robustness): full-sample fits
            # on vintage-year portfolios (one unit per vintage). Same domain
            # and refit rule as the estimator (jackknife.se). jk.n.bound: refits
            # on a bound of Theta; jk.n.switch: refits that left the basin of
            # the full-sample estimate; tables flag SE_JK when either is > 0.
            if (use.vintage.year.pfs && error.function == "L2_Lasso" && cv.key == "ALL") {
              ends <- if (length(par) > 1) fit$endpoints else list()
              jk <- jackknife.se(par, df.optim.prep, max.month, lambda, con, ends, optim.parscale(par))
              for (v in names(par)) res[, paste0("SE.", v, ".jk")] <- jk$se[[v]]
              res$jk.selfcheck <- jk$selfcheck
              res$jk.units <- jk$n
              res$jk.n.bound <- jk$n.bound
              res$jk.n.switch <- jk$n.switch
            }
          }

          # print(res)
          res$Factor <- factor
          res$Type <- type
          res$max.month <- max.month
          res$lambda <- lambda
          res$kernel.bandwidth <- kernel.bandwidth
          res$CV.key <- cv.key
          res$weighting <- weighting
          res$error.fun <- error.function
          res$origin.convention <- origin.convention
          res$vintage.weights <- if (error.function == "U_pair") vintage.weights else NA
          res$sdf.model <- sdf.model
          res$datetime <- Sys.time() + 2 * 60 * 60
          # B9 provenance
          res$git.commit <- provenance$git.commit
          res$git.dirty <- provenance$git.dirty
          res$scenario.id <- provenance$scenario.id
          for (field in names(result.identity.used)) res[[field]] <- result.identity.used[[field]]
          res$data.prepared.folder <- if (use.simulation) NA else data.prepared.folder

          # validation error on the held-out vintages: at the estimation
          # horizon (horizon-specific scale) and, for comparisons across
          # horizons, at the common H = 0 inception loss of the same
          # criterion (B8; for U_pair this is U_{v,0} of the held-out vintages)
          res$validation.error <- err.sqr.calc(par,
            max.month = max.month,
            lambda = 0, # set to zero
            df = df.val.prep
          )
          res$validation.error.H0 <- err.sqr.calc(par, max.month = 1, lambda = 0, df = df.val.prep)

          factor.list[[factor]] <- res
        }
        return(data.frame(Reduce(rbind.all.columns, factor.list)))
      }

      if (do.cache & nrow(output) > 0) {
        cache.path <- paste0(data.out.folder, "/cache_", paste0(public.filename, "_", cache.folder.tag), "/")
        if (!dir.exists(cache.path)) dir.create(cache.path)
        # file name carries fold, horizon, lambda and milliseconds: with a
        # seconds-only timestamp, folds/horizons finishing within the same
        # second overwrote each other's cache file (results silently lost)
        path.cache <- paste0(
          cache.path,
          format(Sys.time() + 2 * 60 * 60, "%Y-%m-%d_%H%M%OS3"),
          "_cv", input.list$cv.year, "_m", max.month, "_l", lambda,
          "_cached_res.csv"
        )
        write.csv(output, path.cache)
      }

      # print(output)
      l[[paste0(max.month, lambda)]] <- output
    }
  }
  df.res <- data.frame(Reduce(rbind.all.columns, l))
  return(df.res)
}

if (do.cross.validation) {
  for (v in names(vintage.blocks)) df.res <- iter.run(vintage.blocks[[v]])
  # All folds and the full-sample fit finished: mark these horizons complete
  # for this result identity (the resume check in run_estimation() requires
  # the marker; fold rows alone do not show that every fold ran)
  if (do.cache) {
    write.cv.complete.markers(
      paste0(data.out.folder, "/cache_", public.filename, "_", cache.folder.tag),
      result.identity.used$cache.identity, max.months, lambdas
    )
  }
} else {
  system.time(df.res <- iter.run(vintage.blocks$ALL))
}

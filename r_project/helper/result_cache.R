######################
### Result-cache identity (review fix 2026-11-16)
### Sourced by helper/estim_model_core.R, estim_model_optimized.R,
### results/analyze_result.R and simulation/analyze_simulation_bias.R
######################
#
# Resume and analysis use the same identity: version, resolved estimator
# settings, and paths/hashes of the data files (cash flows, factors).
# Analysis verifies those files; it never guesses freshness from timestamps.
# The configuration file is provenance only (config.path, config.md5;
# decision of 18 November 2026): what it resolves to is in settings.md5, so
# editing the YAML (comments, other scenarios) does not invalidate results.
#
# Data paths are stored relative to the project root (r_project/, the folder
# that holds estim_model_optimized.R), so result files carry no local
# absolute paths. cache.identity is recomputed from the stored columns
# whenever rows are read (current.result.rows), so rows written before
# 18 November 2026 (absolute paths, configuration in the hash) keep their
# identity.
PIPELINE.VERSION <- "2026-11-17-r5" # r5: admissible domain (A7, gross.floor) and constrained multi-start fit; r3: optimizer parameter scaling

RESULT.DEPENDENCIES <- c("input", "public", "simulation.factors", "equity")
RESULT.DEPENDENCY.COLUMNS <- as.vector(rbind(paste0(RESULT.DEPENDENCIES, ".path"), paste0(RESULT.DEPENDENCIES, ".md5")))
RESULT.IDENTITY.COLUMNS <- c("pipeline.version", "cache.identity", "settings.md5", RESULT.DEPENDENCY.COLUMNS)
RESULT.PROVENANCE.COLUMNS <- c("git.commit", "git.dirty", "scenario.id", "data.prepared.folder",
  "config.path", "config.md5", RESULT.IDENTITY.COLUMNS)
RESULT.SETTINGS <- c("sdf.model", "factors.to.use", "error.function", "include.alpha.term",
  "weighting", "max.vintage", "kernel.bandwidth", "alpha.lower", "alpha.upper",
  "origin.convention", "vintage.weights", "part.to.keep", "no.partitions", "brent.interval", "gross.floor")

# Input data file of an estimation run, relative to r_project/ (single
# definition, used by run_estimation() and estim_model_optimized.R)
input.file.path <- function(use.simulation, simulation.filename = NULL, private.source = "preqin",
                            data.prepared.folder = "", cutoff = "", weighting = "", final.nav.discount = 100) {
  if (use.simulation) return(file.path("simulation/data_prepared_sim", simulation.filename))
  nc.tag <- if (!is.null(final.nav.discount) && final.nav.discount != 100) paste0("_NC", final.nav.discount) else ""
  if (private.source == "preqin") {
    return(paste0("empirical/", data.prepared.folder, "/preqin_cashflows_", weighting, nc.tag, "_NAV.csv"))
  }
  if (private.source == "pitchbook") {
    return(paste0("empirical/", data.prepared.folder, cutoff, "/pitchbook_cashflows_", weighting, nc.tag, "_2026.csv"))
  }
  stop("unknown private.source: ", private.source)
}

file.md5 <- function(path) unname(tools::md5sum(path))

value.md5 <- function(value) {
  path <- tempfile("result-identity-")
  on.exit(unlink(path))
  writeBin(serialize(value, NULL, version = 2), path)
  file.md5(path)
}

# Project root: the nearest folder at or above the working directory that
# holds estim_model_optimized.R (r_project/; scripts also run from results/).
project.root <- function() {
  d <- normalizePath(getwd(), winslash = "/")
  repeat {
    if (file.exists(file.path(d, "estim_model_optimized.R"))) return(d)
    parent <- dirname(d)
    if (parent == d) return(normalizePath(getwd(), winslash = "/"))
    d <- parent
  }
}

is.absolute.path <- function(path) grepl("^(/|[A-Za-z]:)", path)

# Stored form of a data path: relative to the project root. Absolute paths of
# rows written in another checkout are cut after their r_project/ folder.
relative.data.path <- function(path, root = project.root()) {
  path <- gsub("\\\\", "/", as.character(path))
  abs <- !is.na(path) & nzchar(path) & is.absolute.path(path)
  inside <- abs & startsWith(path, paste0(root, "/"))
  path[inside] <- substring(path[inside], nchar(root) + 2L)
  other <- abs & !inside & grepl("/r_project/", path, fixed = TRUE)
  path[other] <- sub("^.*/r_project/", "", path[other])
  path
}

# Path to open: relative paths are taken from the project root.
resolve.data.path <- function(path, root = project.root()) {
  path <- as.character(path)
  rel <- !is.na(path) & nzchar(path) & !is.absolute.path(path)
  path[rel] <- file.path(root, path[rel])
  path
}

# cache.identity of each row: hash of the pipeline version, the settings
# hash and the data paths (stored form) and hashes; one hash per distinct row.
result.identity.hash <- function(df) {
  cols <- c("pipeline.version", "settings.md5", RESULT.DEPENDENCY.COLUMNS)
  vals <- lapply(df[cols], as.character)
  path.cols <- paste0(RESULT.DEPENDENCIES, ".path")
  vals[path.cols] <- lapply(vals[path.cols], relative.data.path)
  key <- do.call(paste, c(vals, sep = "\r"))
  first <- which(!duplicated(key))
  hashes <- vapply(first, function(i) value.md5(lapply(vals, `[`, i)), "")
  hashes[match(key, key[first])]
}

public.file.path <- function(public.filename, data.prepared.folder) {
  folder <- if (public.filename %in% c("q_factors", "ff3_factors")) "data_prepared_2026" else data.prepared.folder
  file.path("empirical", folder, paste0(public.filename, ".csv"))
}

# One namespace builder for estimation, analysis-only lookup and partitions.
result.cache.tag <- function(base.tag, error.function = "L2_Lasso",
                             origin.convention = "first.cf", vintage.weights = "equal", partition = NULL) {
  suffix <- if (error.function == "L2_Lasso") "" else paste0("_", error.function)
  if (error.function == "U_pair" && vintage.weights != "equal") suffix <- paste0(suffix, "_q", vintage.weights)
  if (origin.convention != "first.cf") suffix <- paste0(suffix, "_origin_", origin.convention)
  if (nzchar(suffix) && endsWith(base.tag, suffix)) base.tag <- substr(base.tag, 1L, nchar(base.tag) - nchar(suffix))
  if (!is.null(partition) && !endsWith(base.tag, paste0("_part", partition))) base.tag <- paste0(base.tag, "_part", partition)
  paste0(base.tag, suffix)
}

simulation.cache.tags <- function(params, no.partitions = 1L) {
  vapply(seq_len(no.partitions), function(i) result.cache.tag(params$cache_folder_tag,
    params$error_function, params$origin_convention, params$vintage_weights,
    partition = if (no.partitions > 1) i else NULL), "", USE.NAMES = FALSE)
}

# Optional simulated-factor paths are recorded even when absent, so adding
# a factor override later also invalidates an existing result.
make.result.identity <- function(input.path, public.path, settings, config.path = NULL,
                                 simulation.factors.path = NULL, equity.path = NULL) {
  paths <- list(input = input.path, public = public.path, simulation.factors = simulation.factors.path,
    equity = equity.path, config = config.path)
  fields <- list(pipeline.version = PIPELINE.VERSION)
  for (name in names(paths)) {
    path <- paths[[name]]
    present <- length(path) == 1L && !is.na(path) && nzchar(path)
    if (present) {
      path <- normalizePath(path.expand(path), winslash = "/", mustWork = FALSE)
      if (!is.absolute.path(path)) path <- file.path(getwd(), path)
      if (name != "simulation.factors" && !file.exists(path)) stop("Input file not found: ", path)
      hash <- if (file.exists(path)) file.md5(path) else NA_character_
      if (file.exists(path) && is.na(hash)) stop("Cannot hash input file: ", path)
      path <- relative.data.path(path)
    } else {
      path <- hash <- NA_character_
    }
    fields[[paste0(name, ".path")]] <- path
    fields[[paste0(name, ".md5")]] <- hash
  }
  # YAML integers and direct numeric arguments describe the same settings.
  settings <- settings[order(names(settings))]
  if ("factors.to.use" %in% names(settings)) settings$factors.to.use <- as.character(unlist(settings$factors.to.use))
  settings <- lapply(settings, function(x) if (is.numeric(x)) as.numeric(x) else x)
  fields$settings.md5 <- value.md5(settings)
  out <- as.data.frame(fields, stringsAsFactors = FALSE)
  out$cache.identity <- result.identity.hash(out)
  out
}

make.estimation.identity <- function(input.path, public.filename, data.prepared.folder,
                                     settings, config.path = NULL, use.simulation = FALSE) {
  if (is.null(config.path)) config.path <- file.path("config",
    if (use.simulation) "simulation_scenarios.yaml" else "empirical_scenarios.yaml")
  bond.files <- c("DebtFactorsUSD", "DebtFactorsEUR", "DebtFactorsEURUSD",
    "iBoxxFactorsEUR", "iBoxxFactorsUSD", "iBoxxFactorsMIX")
  make.result.identity(input.path, public.file.path(public.filename, data.prepared.folder),
    settings, config.path,
    simulation.factors.path = if (use.simulation) sub("_simulated_cashflows_.*\\.csv$", "_simulated_factors.csv", input.path) else NULL,
    equity.path = if (public.filename %in% bond.files) file.path(data.prepared.folder, "msci_market_factors.csv") else NULL)
}

result.identities.current <- function(df) {
  if (!all(RESULT.IDENTITY.COLUMNS %in% names(df))) return(rep(FALSE, nrow(df)))
  ok <- !is.na(df$pipeline.version) & df$pipeline.version == PIPELINE.VERSION &
    !is.na(df$cache.identity) & nzchar(df$cache.identity) & !is.na(df$settings.md5)
  # Hash each distinct dependency once, even for large simulation caches.
  hashes <- new.env(parent = emptyenv())
  current.hash <- function(path) {
    if (!exists(path, hashes, inherits = FALSE)) {
      assign(path, if (file.exists(path)) file.md5(path) else NA_character_, hashes)
    }
    get(path, hashes, inherits = FALSE)
  }
  for (name in RESULT.DEPENDENCIES) {
    paths <- as.character(df[[paste0(name, ".path")]])
    stored <- as.character(df[[paste0(name, ".md5")]])
    absent <- is.na(paths) | !nzchar(paths)
    if (name %in% c("input", "public")) ok[absent] <- FALSE
    ok[absent & !is.na(stored)] <- FALSE
    for (path in unique(paths[!absent])) {
      i <- which(!absent & paths == path)
      file <- resolve.data.path(path)
      actual <- current.hash(file)
      matches <- !is.na(stored[i]) & !is.na(actual) & stored[i] == actual
      if (name == "simulation.factors" && !file.exists(file)) matches <- is.na(stored[i])
      ok[i] <- ok[i] & matches
    }
  }
  ok
}

current.result.rows <- function(df, expected.identity = NULL, expected.input.md5 = NULL,
                                label = "results", verbose = TRUE) {
  if (nrow(df) == 0) return(df)
  if (!all(RESULT.IDENTITY.COLUMNS %in% names(df))) {
    if (verbose) message(label, ": ignored ", nrow(df), " result row(s) without a complete cache identity")
    return(df[FALSE, , drop = FALSE])
  }
  # stored form of the paths, and the identity under the current definition
  for (col in intersect(c(paste0(RESULT.DEPENDENCIES, ".path"), "config.path"), names(df))) df[[col]] <- relative.data.path(df[[col]])
  df$cache.identity <- result.identity.hash(df)
  ok <- result.identities.current(df)
  if (!is.null(expected.identity)) ok <- ok & !is.na(df$cache.identity) & df$cache.identity == expected.identity
  if (!is.null(expected.input.md5)) ok <- ok & !is.na(df$input.md5) & df$input.md5 == expected.input.md5
  if (verbose && any(!ok)) message(label, ": ignored ", sum(!ok), " stale or incompatible result row(s)")
  df <- df[ok, , drop = FALSE]
  if (length(unique(df$cache.identity)) > 1L) {
    stop(label, ": multiple current estimator settings; select an explicit cache identity or use separate result folders")
  }
  df
}

# Read cached result rows of one cache folder, keeping only current ones.
#   expected.identity selects the requested settings during resume. Analysis
#   without that argument still verifies every persisted input/config path.
read.result.csv <- function(path) {
  columns <- names(read.csv(path, nrows = 0, stringsAsFactors = FALSE))
  classes <- setNames(rep(NA_character_, length(columns)), columns)
  # Hashes and numeric-only git SHAs are identifiers, not numbers. In
  # particular, read.csv's inference must not discard their leading zeros.
  text <- intersect(setdiff(RESULT.PROVENANCE.COLUMNS, "git.dirty"), columns)
  classes[text] <- "character"
  read.csv(path, stringsAsFactors = FALSE, colClasses = classes)
}

read.result.cache <- function(cache.dir, expected.input.md5 = NULL, verbose = TRUE,
                              expected.identity = NULL) {
  files <- list.files(cache.dir, pattern = "_cached_res\\.csv$", full.names = TRUE)
  rows <- lapply(files, function(f) tryCatch(read.result.csv(f), error = function(e) NULL))
  rows <- rows[!vapply(rows, is.null, TRUE)]
  if (length(rows) == 0) return(data.frame())
  df <- as.data.frame(data.table::rbindlist(rows, use.names = TRUE, fill = TRUE))
  dedupe.result.rows(current.result.rows(df, expected.identity, expected.input.md5, basename(cache.dir), verbose))
}

# Summaries are results too: plotting an old CSV must not bypass validation.
read.result.summary <- function(path) {
  df <- read.result.csv(path)
  if (nrow(df) == 0) return(df)
  out <- current.result.rows(df, label = path)
  if (nrow(out) != nrow(df)) stop("Stale result summary; rerun estimation and analysis: ", path)
  out
}

# Cross-validation completion markers: one empty file per (identity,
# horizon, lambda), written after all folds and the full-sample fit finished.
cv.marker.path <- function(cache.dir, identity, max.month, lambda) {
  file.path(cache.dir, paste0("cv_complete_", identity, "_m", max.month, "_l", lambda, ".marker"))
}
write.cv.complete.markers <- function(cache.dir, identity, max.months, lambdas) {
  if (!dir.exists(cache.dir)) dir.create(cache.dir, recursive = TRUE)
  for (m in max.months) for (l in lambdas) file.create(cv.marker.path(cache.dir, identity, m, l))
  invisible(TRUE)
}

# Duplicate rows of the same estimate (same identity, fold, type, factor,
# horizon and lambda; e.g. a no-CV and a CV run sharing a folder, or a rerun
# after a failure) are reduced to the most recent one, so fold statistics are
# never double counted.
dedupe.result.rows <- function(df) {
  keys <- intersect(c("cache.identity", "CV.key", "Type", "Factor", "max.month", "lambda"), names(df))
  if (nrow(df) == 0 || length(keys) == 0) return(df)
  if ("datetime" %in% names(df)) df <- df[order(as.character(df$datetime)), , drop = FALSE]
  key <- do.call(paste, c(lapply(keys, function(k) as.character(df[[k]])), sep = "\r"))
  df[!duplicated(key, fromLast = TRUE), , drop = FALSE]
}

# Reporting policy shared by horizon, cutoff and direct analysis tables.
boundary.mask <- function(df) {
  if (!"at.bound" %in% names(df)) return(rep(FALSE, nrow(df)))
  flag <- as.logical(df$at.bound)
  !is.na(flag) & flag
}

boundary.marker <- function(df) ifelse(boundary.mask(df), "$^{\\dagger}$", "")

result.table.provenance <- function(dat) {
  commit <- tryCatch(system("git rev-parse --short HEAD", intern = TRUE, ignore.stderr = TRUE)[1], error = function(e) NA)
  vals <- function(col) {
    if (!col %in% names(dat)) return("unknown")
    v <- unique(as.character(dat[[col]]))
    v[is.na(v)] <- "unknown"
    paste(sort(v), collapse = ", ")
  }
  c(paste0("% Generated ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " at commit ", commit),
    paste0("% Estimation results: commit(s) ", vals("git.commit"), "; uncommitted changes: ", vals("git.dirty")),
    paste0("% Scenario(s): ", vals("scenario.id"), "; data folder(s): ", vals("data.prepared.folder")))
}

# Numeric CV means/SDs are computed separately. Identity fields must remain
# scalar, while commits and scenario IDs retain every contributing source.
aggregate.result.provenance <- function(df, key) {
  groups <- split(seq_len(nrow(df)), key)
  rows <- lapply(names(groups), function(id) {
    d <- df[groups[[id]], , drop = FALSE]
    out <- list(id = id)
    for (col in intersect(RESULT.PROVENANCE.COLUMNS, names(d))) {
      values <- unique(d[[col]])
      if (col %in% RESULT.IDENTITY.COLUMNS) {
        if (length(values) != 1L) stop("Mixed result identities in CV summary: ", id)
        out[[col]] <- values[1]
      } else if (col == "git.dirty") {
        out[[col]] <- if (any(values %in% TRUE)) TRUE else if (anyNA(values)) NA else FALSE
      } else {
        values <- as.character(values)
        values[is.na(values)] <- "unknown"
        out[[col]] <- paste(sort(values), collapse = ", ")
      }
    }
    out$at.bound <- any(boundary.mask(d)) # any fold hit a bound
    out$boundary.folds <- sum(boundary.mask(d))
    as.data.frame(out, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

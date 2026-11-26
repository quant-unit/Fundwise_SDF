######################
### Number of parallel workers (estimation and simulation)
### Sourced by estim_model_optimized.R and simulation/simulation_multifactor.R
######################
#
# Workers are capped by CPU cores and by the memory that is free right now,
# so a large run slows down instead of pushing the machine into swap (G1: a
# laptop froze with 31 R processes and 19 GB of swap). Overrides:
#   SDF_MAX_WORKERS=4          environment variable, hard upper limit
#   options(sdf.gb.per.worker) memory budget per worker in GB (default 1)
#   options(sdf.reserve.gb)    memory left for the system in GB (default 6)

# Memory available to new processes in GB (NA if unknown)
available.memory.gb <- function() {
  if (Sys.info()[["sysname"]] == "Darwin") {
    vm <- tryCatch(system("vm_stat", intern = TRUE), error = function(e) character())
    if (!length(vm)) return(NA_real_)
    page <- as.numeric(sub(".*page size of ([0-9]+) bytes.*", "\\1", vm[1]))
    pages <- function(key) {
      x <- grep(paste0("^", key, ":"), vm, value = TRUE)
      if (length(x)) as.numeric(gsub("[^0-9]", "", sub(".*:", "", x[1]))) else 0
    }
    # free + inactive + speculative pages: what macOS can hand out without swapping
    return((pages("Pages free") + pages("Pages inactive") + pages("Pages speculative")) * page / 1024^3)
  }
  if (file.exists("/proc/meminfo")) {
    x <- grep("^MemAvailable:", readLines("/proc/meminfo"), value = TRUE)
    if (length(x)) return(as.numeric(gsub("[^0-9]", "", x)) / 1024^2)
  }
  NA_real_
}

n.parallel.workers <- function(n.tasks, gb.per.worker = getOption("sdf.gb.per.worker", 1),
                               reserve.gb = getOption("sdf.reserve.gb", 6), verbose = TRUE,
                               cores = parallel::detectCores()) {
  by.cpu <- if (length(cores) == 1L && is.finite(cores)) {
    min(n.tasks, max(1L, as.integer(cores) - 2L), 10L)
  } else {
    1L
  }
  cap <- suppressWarnings(as.integer(Sys.getenv("SDF_MAX_WORKERS", "")))
  n <- if (!is.na(cap)) min(by.cpu, cap) else by.cpu
  avail <- available.memory.gb()
  if (is.finite(avail)) {
    by.memory <- floor((avail - reserve.gb) / gb.per.worker)
    if (by.memory < n && verbose) {
      message(sprintf("Parallel workers limited by memory: %d instead of %d (%.1f GB available, %.1f GB reserved, %.1f GB per worker)",
                      max(1, by.memory), n, avail, reserve.gb, gb.per.worker))
    }
    n <- min(n, by.memory)
  }
  max(1L, as.integer(n))
}

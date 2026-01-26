# estim_model_core.R
# ==================
# Clean wrapper around estim_model_optimized.R with explicit parameters.
# Provides a function-based interface instead of relying on global variables.
#
# Usage:
#   source("helper/estim_model_core.R")
#   result <- run_estimation(
#     simulation_file = "20250808_222540/20250808_222540_simulated_cashflows_EW_VYP.csv",
#     sdf_model = "linear",
#     factors_to_use = "ME",
#     weighting = "EW_VYP",
#     max_months = c(1, 60, 120, 180)
#   )

# ============================================================================
# run_estimation() - Main Entry Point
# ============================================================================

#' Run SDF estimation with explicit parameters
#'
#' This function provides a clean interface to the estimation engine by
#' setting up the required global variables and sourcing estim_model_optimized.R.
#' It wraps the legacy interface while providing a modern function-based API.
#'
#' @param simulation_file Path to simulated data CSV (relative to simulation/),
#'   or NULL for empirical data
#' @param private_source Data source for empirical runs: "preqin" or "pitchbook"
#' @param public_filename Public factors file name (default: "q_factors")
#' @param data_prepared_folder Folder containing prepared data
#' @param cutoff Date cutoff string (e.g., "_cutoff_2021")
#' @param sdf_model SDF specification: "linear", "exp.aff", "linear.duration",
#'   "linear.single.date"
#' @param factors_to_use Second factor: "" for MKT-only, or "Alpha", "ME",
#'   "IA", "ROE", "EG"
#' @param error_function Error function: "L2_Lasso", "L1_Ridge", etc.
#' @param include_alpha_term Include alpha term in model (legacy parameter)
#' @param weighting Weighting scheme: "EW", "FW", "EW_VYP", "FW_VYP"
#' @param max_months Vector of maximum compounding months
#' @param lambdas Regularization parameter(s)
#' @param kernel_bandwidth HAC bandwidth (default: 12)
#' @param do_cross_validation Enable cross-validation
#' @param part_to_keep Partition index to keep (for partitioned runs)
#' @param no_partitions Number of partitions
#' @param data_out_folder Output folder for results
#' @param cache_folder_tag Cache folder name (auto-generated if NULL)
#' @param do_cache Enable caching
#' @param do_parallel Enable parallel processing
#' @param verbose Print progress messages
#' @param scenario_id Optional scenario identifier for result tracking
#'
#' @return List with estimation results (or NULL on error)
run_estimation <- function(
    # Data paths
    simulation_file = NULL,
    private_source = "preqin",
    public_filename = "q_factors",
    data_prepared_folder = "data_prepared_2026",
    cutoff = "",
    # Model specification
    sdf_model = "linear",
    factors_to_use = "",
    error_function = "L2_Lasso",
    include_alpha_term = FALSE,
    # Estimation settings
    weighting = "EW_VYP",
    max_months = c(1, 60, 120, 150, 180, 210, 240, 300, 360),
    lambdas = 0,
    kernel_bandwidth = 12,
    # Cross-validation
    do_cross_validation = FALSE,
    part_to_keep = 1,
    no_partitions = 1,
    # Output
    data_out_folder = "simulation/data_out_2026_new",
    cache_folder_tag = NULL,
    # Control
    do_cache = TRUE,
    do_parallel = TRUE,
    verbose = TRUE,
    # Tracking
    scenario_id = NULL) {
    # -------------------------------------------------------------------------
    # Store current working directory to restore later
    # -------------------------------------------------------------------------
    orig_wd <- getwd()

    # Find the r_project root directory
    script_dir <- orig_wd
    if (!file.exists(file.path(script_dir, "estim_model_optimized.R"))) {
        # Try to find it
        if (file.exists(file.path(dirname(script_dir), "estim_model_optimized.R"))) {
            script_dir <- dirname(script_dir)
        } else {
            stop("Cannot find estim_model_optimized.R. Please run from r_project directory.")
        }
    }
    setwd(script_dir)

    on.exit(
        {
            setwd(orig_wd)
        },
        add = TRUE
    )

    # -------------------------------------------------------------------------
    # Print scenario info
    # -------------------------------------------------------------------------
    if (verbose) {
        cat("\n", paste(rep("=", 70), collapse = ""), "\n")
        if (!is.null(scenario_id)) {
            cat("Running scenario: ", scenario_id, "\n")
        }
        cat("SDF model: ", sdf_model, "\n")
        cat("Factors: MKT", if (factors_to_use != "") paste0(" + ", factors_to_use) else "", "\n")
        cat("Weighting: ", weighting, "\n")
        cat("Max months: ", paste(max_months, collapse = ", "), "\n")
        if (!is.null(simulation_file)) {
            cat("Simulation file: ", simulation_file, "\n")
        } else {
            cat("Data source: ", private_source, "\n")
        }
        cat(paste(rep("=", 70), collapse = ""), "\n\n")
    }

    # -------------------------------------------------------------------------
    # Determine run mode
    # -------------------------------------------------------------------------
    use_simulation <- !is.null(simulation_file)

    # -------------------------------------------------------------------------
    # Parse weighting to extract VYP flag
    # -------------------------------------------------------------------------
    use_vintage_year_pfs <- grepl("VYP", weighting)
    weighting_base <- gsub("_VYP", "", weighting)

    # -------------------------------------------------------------------------
    # Generate cache folder tag if not provided
    # -------------------------------------------------------------------------
    if (is.null(cache_folder_tag)) {
        if (use_simulation) {
            # Extract folder tag from simulation file path
            cache_folder_tag <- gsub("/.*", "", simulation_file)
            cache_folder_tag <- paste0(cache_folder_tag, "_simulated_cashflows_", weighting)
        } else {
            # Build tag from parameters
            alpha_str <- if (include_alpha_term) "_alpha_" else "_"
            cache_folder_tag <- paste0(private_source, alpha_str, weighting)
        }
    }

    # -------------------------------------------------------------------------
    # Set global variables for estim_model_optimized.R
    # These are the variables that the legacy script expects to find
    # -------------------------------------------------------------------------

    # Prevent estim_model_optimized.R from using its own defaults
    assign("source.internally", FALSE, envir = .GlobalEnv)

    # Core settings
    assign("export.data", FALSE, envir = .GlobalEnv)
    assign("use.vintage.year.pfs", use_vintage_year_pfs, envir = .GlobalEnv)
    assign("use.simulation", use_simulation, envir = .GlobalEnv)
    assign("do.cross.validation", do_cross_validation, envir = .GlobalEnv)
    assign("do.cache", do_cache, envir = .GlobalEnv)
    assign("do.parallel", do_parallel, envir = .GlobalEnv)

    # Data sources
    assign("private.source", private_source, envir = .GlobalEnv)
    assign("cutoff", cutoff, envir = .GlobalEnv)
    assign("public.filename", public_filename, envir = .GlobalEnv)
    assign("data.prepared.folder", data_prepared_folder, envir = .GlobalEnv)

    # Model settings
    assign("weighting", weighting, envir = .GlobalEnv)
    assign("error.function", error_function, envir = .GlobalEnv)
    assign("sdf.model", sdf_model, envir = .GlobalEnv)
    assign("include.alpha.term", include_alpha_term, envir = .GlobalEnv)
    assign("factors.to.use", factors_to_use, envir = .GlobalEnv)

    # Estimation parameters
    assign("max.months", max_months, envir = .GlobalEnv)
    assign("lambdas", lambdas, envir = .GlobalEnv)
    assign("kernel.bandwidth", kernel_bandwidth, envir = .GlobalEnv)

    # Partitioning
    assign("part.to.keep", part_to_keep, envir = .GlobalEnv)
    assign("no.partitions", no_partitions, envir = .GlobalEnv)

    # Output paths
    assign("data.out.folder", data_out_folder, envir = .GlobalEnv)
    assign("cache.folder.tag", cache_folder_tag, envir = .GlobalEnv)

    # Simulation-specific
    if (use_simulation) {
        assign("simulation.filename", simulation_file, envir = .GlobalEnv)
    }

    # Create output folder if needed
    if (!dir.exists(data_out_folder)) {
        dir.create(data_out_folder, recursive = TRUE)
    }

    # -------------------------------------------------------------------------
    # Run the estimation
    # -------------------------------------------------------------------------
    start_time <- Sys.time()

    tryCatch(
        {
            source("estim_model_optimized.R")

            end_time <- Sys.time()
            elapsed <- difftime(end_time, start_time, units = "secs")

            if (verbose) {
                cat("\nEstimation completed in", round(as.numeric(elapsed), 1), "seconds\n")
            }

            # -------------------------------------------------------------------------
            # Collect results
            # -------------------------------------------------------------------------
            result <- list(
                scenario_id = scenario_id,
                simulation_file = simulation_file,
                sdf_model = sdf_model,
                factors_to_use = factors_to_use,
                weighting = weighting,
                cache_folder_tag = cache_folder_tag,
                data_out_folder = data_out_folder,
                elapsed_time = elapsed,
                success = TRUE
            )

            # Try to extract estimation results from global environment
            if (exists("df.res", envir = .GlobalEnv)) {
                result$df_res <- get("df.res", envir = .GlobalEnv)
            }

            return(result)
        },
        error = function(e) {
            if (verbose) {
                cat("\nError during estimation:\n")
                cat(conditionMessage(e), "\n")
            }

            return(list(
                scenario_id = scenario_id,
                simulation_file = simulation_file,
                success = FALSE,
                error = conditionMessage(e)
            ))
        }
    )
}

# ============================================================================
# run_estimation_batch() - Run Multiple Scenarios
# ============================================================================

#' Run estimation for a batch of scenarios
#'
#' @param scenarios Named list of scenarios (from get_active_scenarios)
#' @param type "simulation" or "empirical"
#' @param verbose Print progress
#' @return Named list of results keyed by scenario_id
run_estimation_batch <- function(scenarios, type = "simulation", verbose = TRUE) {
    results <- list()
    n <- length(scenarios)

    for (i in seq_along(scenarios)) {
        scenario <- scenarios[[i]]
        id <- names(scenarios)[i]

        if (verbose) {
            cat("\n[", i, "/", n, "] Processing scenario: ", id, "\n")
        }

        if (type == "simulation") {
            source("helper/config_loader.R")
            params <- scenario_to_estimation_params(scenario)
        } else {
            source("helper/config_loader.R")
            params <- empirical_scenario_to_params(scenario)
        }

        result <- do.call(run_estimation, params)
        results[[id]] <- result
    }

    # Summary
    if (verbose) {
        cat("\n", paste(rep("=", 70), collapse = ""), "\n")
        cat(
            "Batch complete: ", sum(sapply(results, function(r) r$success)),
            "/", n, " succeeded\n"
        )
        cat(paste(rep("=", 70), collapse = ""), "\n")
    }

    return(results)
}

# ============================================================================
# Self-Test (if run directly)
# ============================================================================

if (sys.nframe() == 0L) {
    cat("estim_model_core.R loaded successfully.\n")
    cat("Use run_estimation() to run with explicit parameters.\n")
    cat("Use run_estimation_batch() to run multiple scenarios.\n")

    cat("\nExample usage:\n")
    cat("  result <- run_estimation(\n")
    cat('    simulation_file = "20250808_222540/20250808_222540_simulated_cashflows_EW_VYP.csv",\n')
    cat('    sdf_model = "linear",\n')
    cat('    factors_to_use = "ME",\n')
    cat('    weighting = "EW_VYP",\n')
    cat("    max_months = c(1, 60, 120, 180)\n")
    cat("  )\n")
}

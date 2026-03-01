# run_empirical_study.R
# =====================
# New empirical study runner with YAML configuration support.
# Provides clean interface for running empirical estimation scenarios.
#
# Usage:
#   source("run_empirical_study.R")
#
#   # Run all active scenarios
#   results <- run_empirical_study()
#
#   # Run specific scenarios
#   results <- run_empirical_study(scenario_ids = c("fw_no_cv", "ew_cv"))
#
#   # Run with different output folder
#   results <- run_empirical_study(
#     data_out_folder = "results/data_out_2026-my-experiment"
#   )

# ============================================================================
# Setup
# ============================================================================

# Ensure we're in the right directory
if (sys.nframe() == 0L) rm(list = ls())

# Try to set working directory (works in RStudio)
tryCatch(
    {
        setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
    },
    error = function(e) {
        # Not running in RStudio, assume already in r_project
    }
)

# Source dependencies
source("helper/config_loader.R")
source("helper/estim_model_core.R")

# ============================================================================
# Main Function: run_empirical_study()
# ============================================================================

#' Run empirical study with YAML configuration
#'
#' @param scenario_ids NULL for all active scenarios, or vector of scenario IDs
#' @param config_path Path to empirical_scenarios.yaml
#' @param data_out_folder Override output folder (default: from config)
#' @param verbose Print progress messages
#'
#' @return List of results keyed by scenario_id
run_empirical_study <- function(
    scenario_ids = NULL,
    config_path = "config/empirical_scenarios.yaml",
    data_out_folder = NULL,
    verbose = TRUE) {
    # -------------------------------------------------------------------------
    # Load configuration
    # -------------------------------------------------------------------------
    if (verbose) {
        cat("\n")
        cat(paste(rep("=", 70), collapse = ""), "\n")
        cat("EMPIRICAL STUDY RUNNER\n")
        cat(paste(rep("=", 70), collapse = ""), "\n")
        cat("Config: ", config_path, "\n")
    }

    config <- load_empirical_config(config_path)
    scenarios <- get_active_scenarios(config, scenario_ids)

    if (length(scenarios) == 0) {
        stop("No scenarios found. Check scenario_ids or 'active' flags in config.")
    }

    if (verbose) {
        cat("Scenarios to run: ", length(scenarios), "\n")
        cat("  ", paste(names(scenarios), collapse = ", "), "\n")
        cat("\n")
    }

    results <- list()

    # -------------------------------------------------------------------------
    # Run estimation for each scenario
    # -------------------------------------------------------------------------
    for (i in seq_along(scenarios)) {
        id <- names(scenarios)[i]
        scenario <- scenarios[[id]]

        if (verbose) {
            cat("[", i, "/", length(scenarios), "] ", id, "\n")
            cat("  Description: ", scenario$description, "\n")
        }

        # Validate scenario
        tryCatch(
            {
                validate_scenario(scenario, type = "empirical")
            },
            error = function(e) {
                cat("  -> INVALID:", conditionMessage(e), "\n\n")
                results[[id]] <- list(success = FALSE, error = conditionMessage(e))
                next
            }
        )

        # Convert scenario to estimation parameters
        params <- empirical_scenario_to_params(scenario)
        params$simulation_file <- NULL # Empirical mode

        # Override data_out_folder if provided
        if (!is.null(data_out_folder)) {
            params$data_out_folder <- data_out_folder
        }

        # ---------------------------------------------------------------------
        # Check for vintage sweep (max_vintages list)
        # ---------------------------------------------------------------------
        if (!is.null(params$max_vintages) && length(params$max_vintages) > 1) {
            max_vintages <- unlist(params$max_vintages)
            params$max_vintages <- NULL # Remove from params (not a run_estimation arg)

            if (verbose) {
                cat("  -> Vintage sweep:", paste(max_vintages, collapse = ", "), "\n")
            }

            for (mv in max_vintages) {
                if (verbose) {
                    cat("  [max_vintage = ", mv, "] ", sep = "")
                }

                # Set max_vintage for this iteration
                params$max_vintage <- mv

                # Build cache_folder_tag with _max_vin_{year} suffix
                # (matches legacy estim_model_empirical_runner.R convention)
                weighting <- params$weighting
                alpha_str <- if (isTRUE(params$include_alpha_term)) "_alpha_" else "_"
                params$cache_folder_tag <- paste0(
                    params$private_source, alpha_str, weighting,
                    if (!is.null(params$final_nav_discount) && params$final_nav_discount != 100) {
                        paste0("_NC", params$final_nav_discount)
                    } else {
                        ""
                    },
                    "_max_vin_", mv
                )

                sub_result <- do.call(run_estimation, params)

                if (verbose && isTRUE(sub_result$success)) {
                    cat("completed in", round(as.numeric(sub_result$elapsed_time), 1), "sec\n")
                } else if (verbose) {
                    cat("FAILED:", sub_result$error, "\n")
                }

                # Store result keyed by scenario_id + vintage
                results[[paste0(id, "_vin_", mv)]] <- sub_result
            }

            if (verbose) cat("\n")
        } else {
            # -----------------------------------------------------------------
            # Standard single-vintage run
            # -----------------------------------------------------------------
            params$max_vintages <- NULL # Remove from params (not a run_estimation arg)

            # Run estimation
            result <- do.call(run_estimation, params)
            results[[id]] <- result

            if (verbose && result$success) {
                cat("  -> Completed in", round(as.numeric(result$elapsed_time), 1), "sec\n\n")
            } else if (verbose) {
                cat("  -> FAILED:", result$error, "\n\n")
            }
        }
    }

    # -------------------------------------------------------------------------
    # Summary
    # -------------------------------------------------------------------------
    if (verbose) {
        cat("\n")
        cat(paste(rep("=", 70), collapse = ""), "\n")
        cat("SUMMARY\n")
        cat(paste(rep("-", 70), collapse = ""), "\n")

        success_count <- sum(sapply(results, function(r) isTRUE(r$success)))
        cat("Completed: ", success_count, "/", length(scenarios), " scenarios\n")

        if (success_count < length(scenarios)) {
            failed <- names(results)[sapply(results, function(r) !isTRUE(r$success))]
            cat("Failed: ", paste(failed, collapse = ", "), "\n")
        }

        cat(paste(rep("=", 70), collapse = ""), "\n")
    }

    return(results)
}

# ============================================================================
# Convenience Functions
# ============================================================================

#' List available empirical scenarios
#'
#' @param config_path Path to config file
#' @param active_only Show only active scenarios
list_empirical_scenarios <- function(
    config_path = "config/empirical_scenarios.yaml",
    active_only = FALSE) {
    config <- load_empirical_config(config_path)
    list_scenarios(config, show_inactive = !active_only)
}

#' Quick run: equivalent to the 4 runs in legacy estim_model_empirical_runner.R
#'
#' Runs: fw_no_cv, fw_cv, ew_no_cv, ew_cv (all one-factor models)
run_default_empirical_study <- function(data_out_folder = NULL, verbose = TRUE) {
    run_empirical_study(
        scenario_ids = c("fw_no_cv", "fw_cv", "ew_no_cv", "ew_cv"),
        data_out_folder = data_out_folder,
        verbose = verbose
    )
}

#' Run with a specific second factor for all weighting schemes
#'
#' @param factor Second factor: "Alpha", "ME", "IA", "ROE", or "EG"
run_two_factor_study <- function(factor = "ME", data_out_folder = NULL, verbose = TRUE) {
    valid_factors <- c("Alpha", "ME", "IA", "ROE", "EG")
    if (!factor %in% valid_factors) {
        stop("Invalid factor. Choose from: ", paste(valid_factors, collapse = ", "))
    }

    # Construct scenario IDs for FW and EW with this factor
    scenario_ids <- c(
        paste0("fw_no_cv_", factor),
        paste0("ew_no_cv_", factor)
    )

    run_empirical_study(
        scenario_ids = scenario_ids,
        data_out_folder = data_out_folder,
        verbose = verbose
    )
}

# ============================================================================
# Self-Test (if run directly)
# ============================================================================

if (sys.nframe() == 0L) {
    cat("\n=== run_empirical_study.R ===\n\n")
    cat("Available functions:\n")
    cat("  run_empirical_study()          - Main runner with YAML config\n")
    cat("  list_empirical_scenarios()     - List available scenarios\n")
    cat("  run_default_empirical_study()  - Quick run (FW/EW x CV/noCV)\n")
    cat("  run_two_factor_study()         - Run with a specific second factor\n")
    cat("\n")

    cat("Available scenarios:\n\n")
    list_empirical_scenarios(active_only = TRUE)

    data.out.folder <- "results/data_out_2026_02_26"
}

scenarios <- c(
    "fw_no_cv", "fw_cv", "ew_no_cv", "ew_cv",
    "fw_vintage_sweep", "ew_vintage_sweep",
    "fw_vintage_sweep_NC50", "ew_vintage_sweep_NC50",
    # Fama French
    "ff3_fw_no_cv_alpha", "ff3_ew_no_cv_alpha",
    "ff3_fw_cv_alpha", "ff3_ew_cv_alpha",
    "ff3_fw_no_cv_ALL", "ff3_ew_no_cv_ALL",
    "ff3_fw_cv_ALL", "ff3_ew_cv_ALL",
    "ff3_fw_cv_ALL_3f", "ff3_ew_cv_ALL_3f",
    "ff3_fw_no_cv_ALL_3f", "ff3_ew_no_cv_ALL_3f",
    "ff3_fw_no_cv", "ff3_ew_no_cv",
    "ff3_fw_cv", "ff3_ew_cv"
)
scenarios <- c("ff3_fw_no_cv", "ff3_ew_no_cv", "ff3_fw_cv", "ff3_ew_cv")

#data.out.folder <- "results/data_out_2026_02_27_pitchbook"
#scenarios <- c("pitchbook_fw_vintage_sweep_NC50", "pitchbook_ew_vintage_sweep_NC50", "pitchbook_fw_vintage_sweep", "pitchbook_ew_vintage_sweep")

results <- run_empirical_study(
    scenario_ids = scenarios,
    data_out_folder = data.out.folder
)

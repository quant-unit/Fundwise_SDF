# run_simulation_study.R
# ======================
# New simulation study runner with YAML configuration support.
# Provides clean interface for running simulation scenarios with
# optional data generation, estimation, and bias analysis.
#
# Usage:
#   source("run_simulation_study.R")
#
#   # Run all active scenarios
#   results <- run_simulation_study()
#
#   # Run specific scenarios
#   results <- run_simulation_study(scenario_ids = c("base_case_vyp", "high_beta_negative_alpha"))
#
#   # Generate data and run estimation
#   results <- run_simulation_study(scenario_ids = "base_case_ME", generate_data = TRUE)
#
#   # Run with bias analysis
#   bias <- run_simulation_study(scenario_ids = "base_case_vyp", analyze = TRUE)

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
# Main Function: run_simulation_study()
# ============================================================================

#' Run simulation study with YAML configuration
#'
#' @param scenario_ids NULL for all active scenarios, or vector of scenario IDs
#' @param config_path Path to simulation_scenarios.yaml
#' @param generate_data Whether to run data generation (default: FALSE)
#' @param estimate Whether to run estimation (default: TRUE)
#' @param analyze Whether to analyze estimation bias (default: FALSE)
#' @param data_out_folder Override output folder (default: from config)
#' @param verbose Print progress messages
#'
#' @return List of results keyed by scenario_id
run_simulation_study <- function(
    scenario_ids = NULL,
    config_path = "config/simulation_scenarios.yaml",
    generate_data = FALSE,
    estimate = TRUE,
    analyze = FALSE,
    data_out_folder = NULL,
    verbose = TRUE) {
    # -------------------------------------------------------------------------
    # Load configuration
    # -------------------------------------------------------------------------
    if (verbose) {
        cat("\n")
        cat(paste(rep("=", 70), collapse = ""), "\n")
        cat("SIMULATION STUDY RUNNER\n")
        cat(paste(rep("=", 70), collapse = ""), "\n")
        cat("Config: ", config_path, "\n")
    }

    config <- load_simulation_config(config_path)
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
    # Generate data (if requested)
    # -------------------------------------------------------------------------
    if (generate_data) {
        if (verbose) cat("--- DATA GENERATION PHASE ---\n\n")

        for (id in names(scenarios)) {
            scenario <- scenarios[[id]]

            if (verbose) cat("Generating data for:", id, "\n")

            # Source the simulation generator (multi-factor version)
            source("simulation/simulation_multifactor.R")

            # Extract DGP parameters
            dgp <- scenario$dgp

            # Call create.simulation.multifactor with all factor parameters
            tryCatch(
                {
                    create.simulation.multifactor(
                        no.deals = dgp$no_deals,
                        investment.period = dgp$investment_period,
                        max.holding.period = dgp$max_holding_period,
                        alpha = dgp$alpha,
                        beta_MKT = dgp$beta_MKT,
                        beta_ME = dgp$beta_ME,
                        beta_IA = dgp$beta_IA,
                        beta_ROE = dgp$beta_ROE,
                        beta_EG = dgp$beta_EG,
                        no.samples = dgp$no_samples,
                        no.funds = dgp$no_funds,
                        min.vin = dgp$min_vin,
                        max.vin = dgp$max_vin,
                        stdvs = dgp$stdvs,
                        exp.aff.sdf = dgp$exp_aff_sdf,
                        use.shifted.lognormal = if (is.null(dgp$use_shifted_lognormal)) FALSE else dgp$use_shifted_lognormal,
                        simulate.mkt = if (is.null(dgp$simulate_mkt)) FALSE else dgp$simulate_mkt,
                        deterministic.inv.timing = if (is.null(dgp$deterministic_inv_timing)) FALSE else dgp$deterministic_inv_timing,
                        scenario_id = id
                    )

                    if (verbose) cat("  -> Data generated successfully\n\n")
                },
                error = function(e) {
                    cat("  -> ERROR:", conditionMessage(e), "\n\n")
                }
            )
        }
    }

    # -------------------------------------------------------------------------
    # Run estimation
    # -------------------------------------------------------------------------
    if (estimate) {
        if (verbose) cat("\n--- ESTIMATION PHASE ---\n\n")

        for (i in seq_along(scenarios)) {
            id <- names(scenarios)[i]
            scenario <- scenarios[[id]]

            if (verbose) {
                cat("[", i, "/", length(scenarios), "] ", id, "\n")
            }

            # Note: data_folder will be auto-derived from scenario_id if null or "auto"

            # Validate scenario
            tryCatch(
                {
                    validate_scenario(scenario, type = "simulation")
                },
                error = function(e) {
                    cat("  -> INVALID:", conditionMessage(e), "\n\n")
                    results[[id]] <- list(success = FALSE, error = conditionMessage(e))
                    next
                }
            )

            # Convert scenario to estimation parameters
            params <- scenario_to_estimation_params(scenario)

            # Override data_out_folder if provided
            if (!is.null(data_out_folder)) {
                params$data_out_folder <- data_out_folder
            } else {
                params$data_out_folder <- "simulation/data_out_2026_new"
            }

            # Determine partition range
            n_partitions <- if (!is.null(scenario$estimation$no_partitions)) {
                scenario$estimation$no_partitions
            } else {
                1
            }

            # Store base cache tag (without _partN suffix)
            base_cache_tag <- params$cache_folder_tag

            # Run estimation (loop over partitions if needed)
            result <- NULL
            all_partitions_ok <- TRUE
            total_elapsed <- 0

            if (n_partitions > 1) {
                for (p in seq_len(n_partitions)) {
                    # Override partition-specific params
                    params$part_to_keep <- p
                    params$cache_folder_tag <- paste0(base_cache_tag, "_part", p)

                    if (verbose) cat("  Partition", p, "/", n_partitions, "\n")

                    result <- do.call(run_estimation, params)

                    if (!isTRUE(result$success)) {
                        all_partitions_ok <- FALSE
                        if (verbose) cat("  -> Partition", p, "FAILED:", result$error, "\n")
                        break
                    } else {
                        total_elapsed <- total_elapsed + as.numeric(result$elapsed_time)
                        if (verbose) {
                            cat(
                                "  -> Partition", p, "completed in",
                                round(as.numeric(result$elapsed_time), 1), "sec\n"
                            )
                        }
                    }
                }

                # Store combined result info for bias analysis
                result$cache_folder_tag <- base_cache_tag
                result$no_partitions <- n_partitions
                result$elapsed_time <- total_elapsed
                result$success <- all_partitions_ok
            } else {
                # Single partition (or no partitions) - run directly
                result <- do.call(run_estimation, params)
            }

            result$dgp <- scenario$dgp # Store DGP for bias analysis
            results[[id]] <- result

            if (verbose && isTRUE(result$success)) {
                cat("  -> Completed in", round(as.numeric(result$elapsed_time), 1), "sec\n\n")
            } else if (verbose && !isTRUE(result$success)) {
                cat("  -> FAILED:", result$error, "\n\n")
            }
        }
    }
    # -------------------------------------------------------------------------
    # Construct results for analysis if estimation was skipped
    # -------------------------------------------------------------------------
    if (analyze && !estimate) {
        if (verbose) cat("\n--- LOCATING EXISTING RESULTS ---\n\n")

        for (id in names(scenarios)) {
            scenario <- scenarios[[id]]

            # Convert to params to get cache tag and paths
            # (suppress warnings about missing files since we just want paths)
            params <- scenario_to_estimation_params(scenario)

            # Override data_out_folder if provided
            if (!is.null(data_out_folder)) {
                params$data_out_folder <- data_out_folder
            } else {
                params$data_out_folder <- "simulation/data_out_2026_new"
            }

            # Determine partition count
            n_partitions <- if (!is.null(scenario$estimation$no_partitions)) {
                scenario$estimation$no_partitions
            } else {
                1
            }

            # Construct a result object that mimics a successful estimation result
            # but only contains what analyze_simulation_bias() needs
            results[[id]] <- list(
                scenario_id = id,
                dgp = scenario$dgp,
                cache_folder_tag = params$cache_folder_tag,
                data_out_folder = params$data_out_folder,
                public_filename = "q_factors", # Assumed default since we don't have the original call
                no_partitions = n_partitions,
                success = TRUE # Assume success to allow analysis of whatever files exist
            )
        }

        if (verbose) cat("Ready for analysis.\n")
    }

    # -------------------------------------------------------------------------
    # Analyze bias
    # -------------------------------------------------------------------------
    if (analyze) {
        if (verbose) cat("\n--- BIAS ANALYSIS PHASE ---\n\n")

        # Source bias analysis script
        if (file.exists("simulation/analyze_simulation_bias.R")) {
            source("simulation/analyze_simulation_bias.R")
            bias_summary <- analyze_simulation_bias(results, scenarios)

            attr(results, "bias_summary") <- bias_summary

            # Save bias results to CSV for aggregation
            if (nrow(bias_summary) > 0) {
                # Get output folder from first successful result, or use default
                output_folder <- if (!is.null(data_out_folder)) {
                    data_out_folder
                } else {
                    # Try to get from first result
                    first_result <- results[[1]]
                    if (!is.null(first_result$data_out_folder)) {
                        first_result$data_out_folder
                    } else {
                        "simulation/data_out_2026_new" # Default
                    }
                }

                bias_output_dir <- file.path(output_folder, "bias_analysis")
                if (!dir.exists(bias_output_dir)) dir.create(bias_output_dir, recursive = TRUE)

                timestamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")

                # 1. Full detailed results (one row per scenario x horizon x factor)
                full_path <- file.path(bias_output_dir, paste0(timestamp, "_bias_full.csv"))
                write.csv(bias_summary, full_path, row.names = FALSE)

                # 2. Aggregated by scenario + horizon (for tables and charts)
                # Aggregate means for bias/estimate columns and sd for estimate columns
                # Use na.action = na.pass to handle columns with NA values
                # Wrap in tryCatch for robustness when certain columns have all NAs
                agg_by_scenario_horizon <- tryCatch(
                    {
                        # Aggregate mean columns
                        agg_mean <- aggregate(
                            cbind(
                                bias_MKT, bias_second,
                                est_beta_MKT, est_second,
                                true_beta_MKT, true_second
                            ) ~ scenario_id + max_month,
                            data = bias_summary,
                            FUN = function(x) mean(x, na.rm = TRUE),
                            na.action = na.pass
                        )

                        # Aggregate sd columns for MKT and second factor estimates
                        agg_sd <- aggregate(
                            cbind(est_beta_MKT, est_second) ~ scenario_id + max_month,
                            data = bias_summary,
                            FUN = function(x) sd(x, na.rm = TRUE),
                            na.action = na.pass
                        )
                        names(agg_sd)[names(agg_sd) == "est_beta_MKT"] <- "sd_MKT"
                        names(agg_sd)[names(agg_sd) == "est_second"] <- "sd_second"

                        # Merge mean and sd results
                        merge(agg_mean, agg_sd, by = c("scenario_id", "max_month"))
                    },
                    error = function(e) {
                        # Fallback: aggregate columns individually to handle all-NA cases
                        if (verbose) {
                            cat("Note: Using fallback aggregation due to NA handling\n")
                        }
                        # Get unique scenario/horizon combinations
                        groups <- unique(bias_summary[, c("scenario_id", "max_month")])

                        # Columns to aggregate with mean
                        cols_mean <- c(
                            "bias_MKT", "bias_second",
                            "est_beta_MKT", "est_second",
                            "true_beta_MKT", "true_second"
                        )

                        result <- groups
                        for (col in cols_mean) {
                            if (col %in% names(bias_summary)) {
                                agg_col <- aggregate(
                                    bias_summary[[col]],
                                    by = list(
                                        scenario_id = bias_summary$scenario_id,
                                        max_month = bias_summary$max_month
                                    ),
                                    FUN = function(x) mean(x, na.rm = TRUE)
                                )
                                names(agg_col)[3] <- col
                                result <- merge(result, agg_col[, c("scenario_id", "max_month", col)],
                                    by = c("scenario_id", "max_month"), all.x = TRUE
                                )
                            }
                        }

                        # Columns to aggregate with sd
                        cols_sd <- list(
                            list(src = "est_beta_MKT", dst = "sd_MKT"),
                            list(src = "est_second", dst = "sd_second")
                        )

                        for (item in cols_sd) {
                            if (item$src %in% names(bias_summary)) {
                                agg_col <- aggregate(
                                    bias_summary[[item$src]],
                                    by = list(
                                        scenario_id = bias_summary$scenario_id,
                                        max_month = bias_summary$max_month
                                    ),
                                    FUN = function(x) sd(x, na.rm = TRUE)
                                )
                                names(agg_col)[3] <- item$dst
                                result <- merge(result, agg_col[, c("scenario_id", "max_month", item$dst)],
                                    by = c("scenario_id", "max_month"), all.x = TRUE
                                )
                            }
                        }

                        result
                    }
                )

                # Sort by scenario then horizon
                agg_by_scenario_horizon <- agg_by_scenario_horizon[
                    order(agg_by_scenario_horizon$scenario_id, agg_by_scenario_horizon$max_month),
                ]

                # Add factor name columns for clarity
                # First factor is always MKT
                agg_by_scenario_horizon$first_factor_name <- "MKT"

                # Second factor name based on scenario configuration
                agg_by_scenario_horizon$second_factor_name <- sapply(
                    agg_by_scenario_horizon$scenario_id,
                    function(sid) {
                        if (sid %in% names(scenarios)) {
                            factors_cfg <- scenarios[[sid]]$estimation$factors_to_use
                            if (is.null(factors_cfg) || factors_cfg == "" || factors_cfg == "MKT") {
                                return(NA_character_) # One-factor model, no second factor
                            } else {
                                return(factors_cfg)
                            }
                        } else {
                            return(NA_character_)
                        }
                    }
                )

                agg_path <- file.path(bias_output_dir, paste0(timestamp, "_bias_by_scenario_horizon.csv"))
                write.csv(agg_by_scenario_horizon, agg_path, row.names = FALSE)

                if (verbose) {
                    cat("\nBias results saved to:\n")
                    cat("  ", full_path, "\n")
                    cat("  ", agg_path, "\n")
                }
            }

            if (verbose) {
                cat("\nBias analysis complete. Access via: attr(results, 'bias_summary')\n")
            }
        } else {
            warning("Bias analysis script not found: simulation/analyze_simulation_bias.R")
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

#' List available simulation scenarios
#'
#' @param config_path Path to config file
#' @param active_only Show only active scenarios
list_simulation_scenarios <- function(
    config_path = "config/simulation_scenarios.yaml",
    active_only = FALSE) {
    config <- load_simulation_config(config_path)
    list_scenarios(config, show_inactive = !active_only)
}

#' Get DGP info for a scenario
#'
#' @param scenario_id Scenario ID
#' @param config_path Path to config file
get_dgp_info <- function(
    scenario_id,
    config_path = "config/simulation_scenarios.yaml") {
    config <- load_simulation_config(config_path)

    if (!scenario_id %in% names(config$scenarios)) {
        stop("Unknown scenario: ", scenario_id)
    }

    return(config$scenarios[[scenario_id]]$dgp)
}

# ============================================================================
# Self-Test (if run directly)
# ============================================================================

if (sys.nframe() == 0L) {
    cat("\n=== run_simulation_study.R ===\n\n")
    cat("Available functions:\n")
    cat("  run_simulation_study()        - Main runner\n")
    cat("  list_simulation_scenarios()   - List available scenarios\n")
    cat("  get_dgp_info()                - Get DGP parameters for a scenario\n")
    cat("\n")

    cat("Available scenarios:\n\n")
    list_simulation_scenarios(active_only = FALSE)
    scenarios <- c(
        "base_case_vyp",
        "base_case_zero_alpha", "base_case_positive_alpha", "base_case_negative_alpha",
        "big_n_v_60funds", "small_n_10funds", "exp_aff_alpha",
        "big_n_v_40funds", "big_v_10funds_1967", "big_v_20funds_1967", "small_v_1986_1995", "small_v_1996_2005",
        "big_n_v_40funds_alpha", "big_v_10funds_1967_alpha", "big_v_20funds_1967_alpha", "small_v_1986_1995_alpha", "small_v_1996_2005_alpha",
        "exp_aff_base", "exp_aff_high_beta_alpha", "high_beta_alpha_two_factor",
        "base_case_ME", "base_case_IA", "base_case_ROE", "base_case_EG",
        "big_n_v_50funds_alpha_stdv30_shifted", "big_n_v_50funds_alpha_stdv30", "big_n_v_50funds_alpha", "big_n_v_50funds_alpha_stdv30_shifted_mkt", 
        "big_n_v_50funds_alpha_stdv30_shifted_mkt2", "big_n_v_50funds_alpha_stdv30_shifted_mkt3",
        "big_n_v_50funds_alpha_stdv30_shifted_mkt4", "big_n_v_50funds_alpha_stdv30_shifted_mkt6", 
        "big_n_v_50funds_alpha_stdv30_shifted_mkt7",
        "base_case_cross_sectional_zero_alpha", "base_case_cross_sectional"
    )

    # scenarios <- c("big_n_v_50funds_alpha_stdv30_shifted_mkt7")
    scenarios <- c("base_case_cross_sectional_zero_alpha", "base_case_cross_sectional")
    results <- run_simulation_study(
        scenario_ids = scenarios,
        generate_data = FALSE, estimate = TRUE, analyze = FALSE
    )
}

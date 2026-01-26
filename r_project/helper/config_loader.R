# config_loader.R
# ================
# Utility functions for loading and parsing YAML configuration files.
# Provides functions for simulation and empirical scenario configs.
#
# Usage:
#   source("helper/config_loader.R")
#   sim_config <- load_simulation_config("config/simulation_scenarios.yaml")
#   emp_config <- load_empirical_config("config/empirical_scenarios.yaml")

# ============================================================================
# Dependencies
# ============================================================================

if (!requireNamespace("yaml", quietly = TRUE)) {
    message("Installing 'yaml' package...")
    install.packages("yaml")
}
library(yaml)

# ============================================================================
# Load Simulation Configuration
# ============================================================================

#' Load and parse simulation scenarios from YAML config
#'
#' @param yaml_path Path to the simulation_scenarios.yaml file
#' @return List with 'defaults' and 'scenarios' elements
load_simulation_config <- function(yaml_path = "config/simulation_scenarios.yaml") {
    if (!file.exists(yaml_path)) {
        stop("Configuration file not found: ", yaml_path)
    }

    config <- yaml::read_yaml(yaml_path)

    # Merge defaults into each scenario
    defaults <- config$defaults
    scenarios <- config$scenarios

    for (id in names(scenarios)) {
        scenarios[[id]] <- merge_with_defaults(scenarios[[id]], defaults)
        scenarios[[id]]$id <- id # Store scenario ID within the scenario
    }

    return(list(
        defaults = defaults,
        scenarios = scenarios
    ))
}

# ============================================================================
# Load Empirical Configuration
# ============================================================================

#' Load and parse empirical scenarios from YAML config
#'
#' @param yaml_path Path to the empirical_scenarios.yaml file
#' @return List with 'defaults' and 'scenarios' elements
load_empirical_config <- function(yaml_path = "config/empirical_scenarios.yaml") {
    if (!file.exists(yaml_path)) {
        stop("Configuration file not found: ", yaml_path)
    }

    config <- yaml::read_yaml(yaml_path)

    # Merge defaults into each scenario
    defaults <- config$defaults
    scenarios <- config$scenarios

    for (id in names(scenarios)) {
        scenarios[[id]] <- merge_with_defaults_flat(scenarios[[id]], defaults)
        scenarios[[id]]$id <- id
    }

    return(list(
        defaults = defaults,
        scenarios = scenarios
    ))
}

# ============================================================================
# Get Active Scenarios
# ============================================================================

#' Get scenarios to run based on IDs or active flag
#'
#' @param config Config object from load_*_config()
#' @param scenario_ids NULL for all active, or vector of scenario IDs
#' @return Named list of scenarios to run
get_active_scenarios <- function(config, scenario_ids = NULL) {
    scenarios <- config$scenarios

    if (is.null(scenario_ids)) {
        # Return all active scenarios
        active <- sapply(scenarios, function(s) isTRUE(s$active))
        return(scenarios[active])
    } else {
        # Return specified scenarios (regardless of active flag)
        missing <- setdiff(scenario_ids, names(scenarios))
        if (length(missing) > 0) {
            warning("Unknown scenario IDs: ", paste(missing, collapse = ", "))
        }
        found <- intersect(scenario_ids, names(scenarios))
        return(scenarios[found])
    }
}

# ============================================================================
# List Available Scenarios
# ============================================================================

#' Print summary of available scenarios
#'
#' @param config Config object from load_*_config()
#' @param show_inactive Whether to include inactive scenarios
list_scenarios <- function(config, show_inactive = TRUE) {
    scenarios <- config$scenarios

    cat("Available scenarios:\n")
    cat(sprintf("%-30s %-8s %s\n", "ID", "Active", "Description"))
    cat(paste(rep("-", 80), collapse = ""), "\n")

    for (id in names(scenarios)) {
        s <- scenarios[[id]]
        active_str <- if (isTRUE(s$active)) "YES" else "no"
        desc <- if (!is.null(s$description)) s$description else ""

        if (show_inactive || isTRUE(s$active)) {
            cat(sprintf("%-30s %-8s %s\n", id, active_str, desc))
        }
    }
}

# ============================================================================
# Validate Scenario
# ============================================================================

#' Validate a scenario has required fields
#'
#' @param scenario Scenario object
#' @param type "simulation" or "empirical"
#' @return TRUE if valid, stops with error otherwise
validate_scenario <- function(scenario, type = "simulation") {
    if (type == "simulation") {
        required_dgp <- c("beta_MKT", "no_samples", "no_funds", "min_vin", "max_vin")
        required_est <- c("sdf_model", "weighting", "max_months")

        if (is.null(scenario$dgp)) {
            stop("Scenario missing 'dgp' section: ", scenario$id)
        }

        for (field in required_dgp) {
            if (is.null(scenario$dgp[[field]])) {
                stop("Scenario '", scenario$id, "' missing dgp.", field)
            }
        }

        if (is.null(scenario$estimation)) {
            stop("Scenario missing 'estimation' section: ", scenario$id)
        }

        for (field in required_est) {
            if (is.null(scenario$estimation[[field]])) {
                stop("Scenario '", scenario$id, "' missing estimation.", field)
            }
        }
    } else if (type == "empirical") {
        required <- c("weighting", "do_cross_validation")

        for (field in required) {
            if (is.null(scenario[[field]])) {
                stop("Scenario '", scenario$id, "' missing field: ", field)
            }
        }
    }

    return(TRUE)
}

# ============================================================================
# Helper: Merge with Defaults (Nested for Simulation)
# ============================================================================

#' Recursively merge scenario with defaults
#'
#' @param scenario Scenario from YAML
#' @param defaults Default values
#' @return Merged scenario with all defaults applied
merge_with_defaults <- function(scenario, defaults) {
    # Merge DGP defaults
    if (!is.null(defaults$dgp)) {
        if (is.null(scenario$dgp)) scenario$dgp <- list()
        for (key in names(defaults$dgp)) {
            if (is.null(scenario$dgp[[key]])) {
                scenario$dgp[[key]] <- defaults$dgp[[key]]
            }
        }
    }

    # Merge estimation defaults
    if (!is.null(defaults$estimation)) {
        if (is.null(scenario$estimation)) scenario$estimation <- list()
        for (key in names(defaults$estimation)) {
            if (is.null(scenario$estimation[[key]])) {
                scenario$estimation[[key]] <- defaults$estimation[[key]]
            }
        }
    }

    return(scenario)
}

# ============================================================================
# Helper: Merge with Defaults (Flat for Empirical)
# ============================================================================

#' Merge flat scenario with defaults (no nested structure)
#'
#' @param scenario Scenario from YAML
#' @param defaults Default values
#' @return Merged scenario with all defaults applied
merge_with_defaults_flat <- function(scenario, defaults) {
    for (key in names(defaults)) {
        if (is.null(scenario[[key]])) {
            scenario[[key]] <- defaults[[key]]
        }
    }
    return(scenario)
}

# ============================================================================
# Convert Scenario to Estimation Parameters
# ============================================================================

#' Convert simulation scenario to parameters for run_estimation()
#'
#' @param scenario Scenario object with merged defaults
#' @param base_path Base path for data files
#' @return Named list of parameters
scenario_to_estimation_params <- function(scenario, base_path = "simulation") {
    est <- scenario$estimation

    # Build simulation file path
    sim_file <- NULL
    if (!is.null(scenario$data_folder)) {
        folder_tag <- paste0(scenario$data_folder, "_simulated_cashflows_EW_VYP")
        sim_file <- paste0(scenario$data_folder, "/", folder_tag, ".csv")
    }

    # Build weighting string
    weighting <- est$weighting

    params <- list(
        simulation_file = sim_file,
        sdf_model = est$sdf_model,
        factors_to_use = if (!is.null(est$factors_to_use)) est$factors_to_use else "",
        error_function = est$error_function,
        weighting = weighting,
        max_months = est$max_months,
        lambdas = est$lambdas,
        kernel_bandwidth = est$kernel_bandwidth,
        do_cross_validation = isTRUE(est$do_cross_validation),
        part_to_keep = est$part_to_keep,
        no_partitions = est$no_partitions,
        scenario_id = scenario$id
    )

    return(params)
}

#' Convert empirical scenario to parameters for run_estimation()
#'
#' @param scenario Scenario object with merged defaults
#' @return Named list of parameters
empirical_scenario_to_params <- function(scenario) {
    # Build weighting string
    weighting <- scenario$weighting
    if (isTRUE(scenario$use_vintage_year_pfs)) {
        weighting <- paste0(weighting, "_VYP")
    }

    params <- list(
        private_source = scenario$private_source,
        public_filename = scenario$public_filename,
        data_prepared_folder = scenario$data_prepared_folder,
        cutoff = scenario$cutoff,
        sdf_model = scenario$sdf_model,
        error_function = scenario$error_function,
        include_alpha_term = isTRUE(scenario$include_alpha_term),
        weighting = weighting,
        max_months = scenario$max_months,
        lambdas = scenario$lambdas,
        kernel_bandwidth = scenario$kernel_bandwidth,
        do_cross_validation = isTRUE(scenario$do_cross_validation),
        part_to_keep = scenario$part_to_keep,
        no_partitions = scenario$no_partitions,
        data_out_folder = scenario$data_out_folder,
        factors_to_use = if (!is.null(scenario$factors_to_use)) scenario$factors_to_use else "",
        scenario_id = scenario$id
    )

    return(params)
}

# ============================================================================
# Self-Test (if run directly)
# ============================================================================

if (sys.nframe() == 0L) {
    cat("Testing config_loader.R...\n\n")

    # Test simulation config
    if (file.exists("config/simulation_scenarios.yaml")) {
        sim_cfg <- load_simulation_config()
        cat("Loaded", length(sim_cfg$scenarios), "simulation scenarios\n")
        list_scenarios(sim_cfg, show_inactive = FALSE)
        cat("\n")
    }

    # Test empirical config
    if (file.exists("config/empirical_scenarios.yaml")) {
        emp_cfg <- load_empirical_config()
        cat("Loaded", length(emp_cfg$scenarios), "empirical scenarios\n")
        list_scenarios(emp_cfg, show_inactive = FALSE)
    }
}

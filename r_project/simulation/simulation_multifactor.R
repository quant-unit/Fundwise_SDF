# simulation_multifactor.R
# ========================
# Extended simulation data generator with multi-factor DGP support.
# This file provides a new function that supports two-factor models
# while maintaining backward compatibility with the original simulation_improved.R.
#
# Supported factors: MKT + one of (Alpha, ME, IA, ROE, EG)
#
# Usage:
#   source("simulation/simulation_multifactor.R")
#   create.simulation.multifactor(
#     beta_MKT = 1,
#     beta_ME = 0.5,   # Second factor loading
#     ...
#   )

# ============================================================================
# Dependencies
# ============================================================================

if (!requireNamespace("data.table", quietly = TRUE)) {
    install.packages("data.table")
}
library(data.table)

# ============================================================================
# Load Q-Factors Data
# ============================================================================

# Handle both direct sourcing and sourcing from run_simulation_study.R
if (exists("df.q5", envir = .GlobalEnv) && is.data.frame(get("df.q5", envir = .GlobalEnv))) {
    # Already loaded, use it
} else {
    # Load from file
    current_dir <- getwd()

    # Try different paths
    possible_paths <- c(
        file.path(current_dir, "empirical", "data_prepared", "q_factors.csv"),
        file.path(current_dir, "empirical", "data_prepared_2026", "q_factors.csv"),
        file.path(dirname(current_dir), "empirical", "data_prepared", "q_factors.csv"),
        file.path(dirname(current_dir), "empirical", "data_prepared_2026", "q_factors.csv")
    )

    loaded <- FALSE
    for (path in possible_paths) {
        if (file.exists(path)) {
            df.q5 <- read.csv(path)
            df.q5$Date <- as.Date(df.q5$Date)
            loaded <- TRUE
            break
        }
    }

    if (!loaded) {
        stop("Cannot find q_factors.csv. Please run from r_project directory.")
    }
}

# Ensure Alpha column exists
if (!"Alpha" %in% colnames(df.q5)) {
    df.q5$Alpha <- 1
}

# ============================================================================
# create.simulation.multifactor() - Multi-Factor DGP Support
# ============================================================================

#' Create simulated cash flows with multi-factor DGP
#'
#' Extended version of create.simulation() that supports two-factor models.
#' The DGP is: R = alpha + RF + beta_MKT*MKT + beta_F*Factor + noise
#' where Factor is one of ME, IA, ROE, EG.
#'
#' @param no.deals Number of deals per fund
#' @param investment.period Investment period in years
#' @param max.holding.period Maximum holding period in years
#' @param alpha Monthly alpha (intercept)
#' @param beta_MKT Loading on market factor (MKT)
#' @param beta_ME Loading on ME factor (default: 0)
#' @param beta_IA Loading on IA factor (default: 0)
#' @param beta_ROE Loading on ROE factor (default: 0)
#' @param beta_EG Loading on EG factor (default: 0)
#' @param no.samples Number of Monte Carlo samples
#' @param no.funds Number of funds per vintage year
#' @param min.vin Earliest vintage year
#' @param max.vin Latest vintage year
#' @param stdvs Idiosyncratic volatility
#' @param exp.aff.sdf Use exponential-affine SDF (default: FALSE)
#' @param scenario_id Optional scenario identifier for metadata
#' @param output_folder Output folder for generated data
#'
#' @return Invisibly returns list with file paths
create.simulation.multifactor <- function(
    no.deals,
    investment.period,
    max.holding.period,
    alpha,
    beta_MKT,
    beta_ME = 0,
    beta_IA = 0,
    beta_ROE = 0,
    beta_EG = 0,
    no.samples,
    no.funds,
    min.vin,
    max.vin,
    stdvs,
    exp.aff.sdf = FALSE,
    use.shifted.lognormal = FALSE,
    scenario_id = NULL,
    output_folder = "simulation/data_prepared_sim") {
    set.seed(100) # For reproducibility
    months <- 12

    # -------------------------------------------------------------------------
    # Identify active second factor
    # -------------------------------------------------------------------------
    second_factors <- c(ME = beta_ME, IA = beta_IA, ROE = beta_ROE, EG = beta_EG)
    active_second_factor <- names(second_factors[second_factors != 0])

    if (length(active_second_factor) > 1) {
        warning(
            "Multiple non-zero second factors detected. Using only the first: ",
            active_second_factor[1]
        )
        active_second_factor <- active_second_factor[1]
    } else if (length(active_second_factor) == 0) {
        active_second_factor <- NULL
    }

    cat("DGP Configuration:\n")
    cat("  alpha =", alpha, "\n")
    cat("  beta_MKT =", beta_MKT, "\n")
    if (!is.null(active_second_factor)) {
        cat("  beta_", active_second_factor, " = ",
            second_factors[active_second_factor], "\n",
            sep = ""
        )
    }
    cat(
        "  Samples:", no.samples, "| Funds:", no.funds, "| Vintages:",
        min.vin, "-", max.vin, "\n\n"
    )

    # -------------------------------------------------------------------------
    # Inner function: Generate a single fund
    # -------------------------------------------------------------------------
    make.fund <- function(no = 0, vintage = 1990, stdv = 0.01, exp.aff = exp.aff.sdf,
                          use.shifted = use.shifted.lognormal,
                          df = df.q5) {
        # Filter data to start from vintage year
        df <- df[as.integer(format(df$Date, "%Y")) >= vintage, ]
        df$Vintage <- vintage
        df$CF <- 0

        # Build type string including factor info
        sdf_str <- ifelse(exp.aff, "ea", "sl")
        factor_str <- ifelse(!is.null(active_second_factor),
            paste0(" +", active_second_factor),
            ""
        )
        df$type <- paste0("sdf:", sdf_str, factor_str, " stdv:", stdv, " #d:", no.deals)
        df$Fund.ID <- paste(vintage, no, sep = "_")

        # Shifted lognormal parameters (for noise distribution)
        min.mkt <- 0.25
        foo <- function(sigma) {
            mu <- log((1 - min.mkt) / exp(0.5 * sigma^2))
            y <- exp(mu + 0.5 * sigma^2) * sqrt(exp(sigma^2) - 1)
            return((stdv - y)^2)
        }
        sigma <- optimize(foo, c(0, 10))$minimum
        mu <- log((1 - min.mkt) / exp(0.5 * sigma^2))

        # Simulate N deals
        for (i in 1:no.deals) {
            if (investment.period == 0) {
                start <- 1
            } else {
                start <- sample(1:(investment.period * months), 1)
            }
            end <- start + sample(1:(max.holding.period * months), 1)
            end <- min(end, nrow(df))

            if (exp.aff) {
                # Exponential-affine SDF
                df$deal <- alpha + log(1 + df$RF) +
                    beta_MKT * log(1 + df$MKT) +
                    beta_ME * log(1 + df$ME) +
                    beta_IA * log(1 + df$IA) +
                    beta_ROE * log(1 + df$ROE) +
                    beta_EG * log(1 + df$EG) +
                    rnorm(nrow(df), 0, stdv) - stdv^2 / 2
                m <- exp(sum(df$deal[start:end]))
            } else {
                # Simple linear SDF
                base_deal <- alpha + df$RF +
                    beta_MKT * df$MKT +
                    beta_ME * df$ME +
                    beta_IA * df$IA +
                    beta_ROE * df$ROE +
                    beta_EG * df$EG

                if (use.shifted) {
                    lower.bound <- -1
                    sigma <- sqrt(log(1 + stdv^2 / lower.bound^2))
                    mu <- log(-lower.bound) - sigma^2 / 2
                    df$deal <- base_deal + exp(rnorm(nrow(df), mu, sigma)) + lower.bound
                } else {
                    df$deal <- base_deal + rnorm(nrow(df), 0, stdv)
                }

                path <- cumprod(1 + df$deal[start:end])
                if (any(path < 0)) {
                    m <- 0 # Default
                } else {
                    m <- tail(path, 1)
                }
            }

            df$CF[start] <- df$CF[start] - 1
            df$CF[end] <- df$CF[end] + m
        }

        df <- df[df$Date <= max(df$Date[df$CF != 0]), ]
        df <- df[df$CF != 0, ]
        cols <- c("Date", "CF", "type", "Vintage", "Fund.ID")
        df <- df[, cols]
        return(df)
    }

    # -------------------------------------------------------------------------
    # Generate all samples
    # -------------------------------------------------------------------------
    cat("Generating", no.samples, "samples...\n")
    pb_interval <- max(1, no.samples %/% 10)

    output <- list()
    for (x in 1:no.samples) {
        if (x %% pb_interval == 0) {
            cat("  Progress:", x, "/", no.samples, "\n")
        }
        set.seed(x)
        l <- list()
        for (i in 1:no.funds) {
            for (v in min.vin:max.vin) {
                for (s in stdvs) {
                    l[[paste(i, v, s)]] <- make.fund(i, v, s)
                }
            }
        }
        df.out <- as.data.frame(data.table::rbindlist(l))
        df.out$type <- paste0(df.out$type, " #f:", no.funds, " #y:", max.vin - min.vin, " .", x)
        df.out$Fund.ID <- paste0(df.out$type, "+", df.out$Fund.ID)
        output[[paste0("i", x)]] <- df.out
    }

    df.out <- data.table::rbindlist(output)
    rm(output)

    # -------------------------------------------------------------------------
    # Aggregate to Vintage-Year Portfolios
    # -------------------------------------------------------------------------
    cat("Aggregating to vintage-year portfolios...\n")
    df.vyp <- df.out[, .(CF = sum(CF)), by = .(Date, type, Vintage)]
    df.vyp$Fund.ID <- paste0(df.vyp$type, " v:", df.vyp$Vintage)

    # -------------------------------------------------------------------------
    # Save output files
    # -------------------------------------------------------------------------
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

    # Use scenario_id for folder/file naming if provided, otherwise fall back to timestamp
    if (!is.null(scenario_id)) {
        folder_name <- scenario_id
    } else {
        folder_name <- timestamp
    }

    DATA_PREPARED_SIM_PATH <- file.path(output_folder, folder_name, "")
    dir.create(DATA_PREPARED_SIM_PATH, recursive = TRUE, showWarnings = FALSE)

    filename.vyp <- paste0(DATA_PREPARED_SIM_PATH, folder_name, "_simulated_cashflows_EW_VYP.csv")
    filename.fund <- paste0(DATA_PREPARED_SIM_PATH, folder_name, "_simulated_cashflows_EW.csv")
    filename.meta <- paste0(DATA_PREPARED_SIM_PATH, folder_name, "_simulated_cashflows_EW_meta.csv")

    # Build metadata with all factor loadings
    df.meta <- data.frame(
        scenario_id = ifelse(is.null(scenario_id), "", scenario_id),
        timestamp = timestamp,
        investment.period = investment.period,
        max.holding.period = max.holding.period,
        no.deals = no.deals,
        alpha = alpha,
        beta_MKT = beta_MKT,
        beta_ME = beta_ME,
        beta_IA = beta_IA,
        beta_ROE = beta_ROE,
        beta_EG = beta_EG,
        no.samples = no.samples,
        no.funds = no.funds,
        min.vin = min.vin,
        max.vin = max.vin,
        stdvs = stdvs,
        exp.aff.sdf = exp.aff.sdf,
        use.shifted.lognormal = use.shifted.lognormal
    )

    cat("Saving files to:", DATA_PREPARED_SIM_PATH, "\n")
    write.csv(df.vyp, filename.vyp, row.names = FALSE)
    write.csv(df.out, filename.fund, row.names = FALSE)
    write.csv(df.meta, filename.meta, row.names = FALSE)

    cat("Done! Generated files:\n")
    cat("  VYP data:", basename(filename.vyp), "\n")
    cat("  Fund data:", basename(filename.fund), "\n")
    cat("  Metadata:", basename(filename.meta), "\n\n")

    return(invisible(list(
        vyp_path = filename.vyp,
        fund_path = filename.fund,
        meta_path = filename.meta,
        folder = DATA_PREPARED_SIM_PATH,
        timestamp = timestamp,
        metadata = df.meta
    )))
}

# ============================================================================
# Backward-Compatible Wrapper
# ============================================================================

#' Backward-compatible wrapper for create.simulation.multifactor
#'
#' Maps old parameter names to new function. Can be used as drop-in replacement.
create.simulation <- function(
    no.deals,
    investment.period,
    max.holding.period,
    alpha,
    beta, # Old name for beta_MKT
    beta_ME = 0,
    beta_IA = 0,
    beta_ROE = 0,
    beta_EG = 0,
    no.samples,
    no.funds,
    min.vin,
    max.vin,
    stdvs,
    exp.aff.sdf,
    use.shifted.lognormal = FALSE,
    scenario_id = NULL) {
    create.simulation.multifactor(
        no.deals = no.deals,
        investment.period = investment.period,
        max.holding.period = max.holding.period,
        alpha = alpha,
        beta_MKT = beta,
        beta_ME = beta_ME,
        beta_IA = beta_IA,
        beta_ROE = beta_ROE,
        beta_EG = beta_EG,
        no.samples = no.samples,
        no.funds = no.funds,
        min.vin = min.vin,
        max.vin = max.vin,
        stdvs = stdvs,
        exp.aff.sdf = exp.aff.sdf,
        use.shifted.lognormal = use.shifted.lognormal,
        scenario_id = scenario_id
    )
}

# ============================================================================
# Self-Test (if run directly)
# ============================================================================

if (sys.nframe() == 0L) {
    cat("\n=== simulation_multifactor.R ===\n\n")
    cat("This script provides multi-factor DGP support for simulation studies.\n\n")
    cat("Available functions:\n")
    cat("  create.simulation.multifactor()  - Full multi-factor support\n")
    cat("  create.simulation()              - Backward-compatible wrapper\n")
    cat("\n")
    cat("Example:\n")
    cat("  create.simulation.multifactor(\n")
    cat("    no.deals = 15, investment.period = 5, max.holding.period = 10,\n")
    cat("    alpha = 0, beta_MKT = 1, beta_ME = 0.5,\n")
    cat("    no.samples = 100, no.funds = 5, min.vin = 1990, max.vin = 1995,\n")
    cat("    stdvs = 0.2, exp.aff.sdf = FALSE\n")
    cat("  )\n")
}

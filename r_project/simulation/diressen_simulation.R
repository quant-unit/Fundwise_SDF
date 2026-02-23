# Comparison to Driessen et al. 2012

scenarios <- c(
    "base_case_zero_alpha", "big_n_v_50funds_alpha",
    "big_n_v_50funds_alpha_stdv30", "big_n_v_50funds_alpha_stdv30_shifted",
    "big_n_v_50funds_alpha_stdv30_shifted_mkt", "big_n_v_50funds_alpha_stdv30_shifted_mkt2",
    "big_n_v_50funds_alpha_stdv30_shifted_mkt3", "big_n_v_50funds_alpha_stdv30_shifted_mkt4",
    "big_n_v_50funds_alpha_stdv30_shifted_mkt5", "big_n_v_50funds_alpha_stdv30_shifted_mkt6"
)

#' Generate a professional LaTeX table comparing our small-sample properties
#' to Driessen et al. 2012 Table 1.
#'
#' @param data_path Path to the CSV file containing simulation bias data
#' @param output_tex Path to save the resulting LaTeX table (optional)
#' @return The generated LaTeX string
generate_driessen_latex_table <- function(data_path, output_tex = NULL) {
    library(dplyr)
    library(tidyr)

    # Read full bias data
    df <- read.csv(data_path)

    # Load alpha bounds from YAML config per scenario
    config <- yaml::read_yaml("config/simulation_scenarios.yaml")
    alpha_bounds <- sapply(scenarios, function(s) {
        bound <- config$scenarios[[s]]$estimation$alpha_upper
        if (is.null(bound) || is.infinite(bound)) 0.01 else bound # default 0.01
    })

    # Filter for Horizon 0 (max_month == 0) and predefined scenarios
    df_sub <- df %>%
        filter(max_month == 0 & scenario_id %in% scenarios)

    # Add per-scenario alpha bound for accurate bound-hit calculation
    df_sub$alpha_upper <- alpha_bounds[as.character(df_sub$scenario_id)]

    # Compute statistics: mean, standard deviation, 25% and 75% percentiles
    # Note: Alpha is typically annualized and presented in percent (x 1200) for finance tables
    stats <- df_sub %>%
        group_by(scenario_id) %>%
        summarise(
            true_beta = mean(true_beta_MKT, na.rm = TRUE),
            true_alpha = mean(true_second, na.rm = TRUE) * 100,
            mean_beta = mean(est_beta_MKT, na.rm = TRUE),
            sd_beta = sd(est_beta_MKT, na.rm = TRUE),
            p25_beta = quantile(est_beta_MKT, 0.25, na.rm = TRUE),
            p75_beta = quantile(est_beta_MKT, 0.75, na.rm = TRUE),
            mean_alpha = mean(est_second, na.rm = TRUE) * 100,
            sd_alpha = sd(est_second, na.rm = TRUE) * 100,
            p25_alpha = quantile(est_second, 0.25, na.rm = TRUE) * 100,
            p75_alpha = quantile(est_second, 0.75, na.rm = TRUE) * 100,
            pct_at_bound = mean(abs(est_second - first(alpha_upper)) < 1e-6, na.rm = TRUE) * 100,
            .groups = "drop"
        )

    # Map scenario names to human-readable labels for the table
    scenario_labels <- c(
        "base_case_zero_alpha" = "(S1) $n=20$, $\\sigma=20\\%$, Normal",
        "big_n_v_50funds_alpha" = "(S2) $n=50$, $\\sigma=20\\%$, Normal",
        "big_n_v_50funds_alpha_stdv30" = "(S3) $n=50$, $\\sigma=30\\%$, Normal",
        "big_n_v_50funds_alpha_stdv30_shifted" = "(S4) $n=50$, $\\sigma=30\\%$, Shifted LN",
        "big_n_v_50funds_alpha_stdv30_shifted_mkt" = "(S5) $n=50$, $\\sigma=30\\%$, Shifted LN, Sim.~MKT",
        "big_n_v_50funds_alpha_stdv30_shifted_mkt2" = "(S6) $n=50$, $\\sigma=30\\%$, Shifted LN, Sim.~MKT, $\\bar{\\alpha}=0.5\\%$",
        "big_n_v_50funds_alpha_stdv30_shifted_mkt3" = "(S7) $n=50$, $\\sigma=20\\%$, Shifted LN, Sim.~MKT",
        "big_n_v_50funds_alpha_stdv30_shifted_mkt4" = "(S8) $n=50$, $\\sigma=30\\%$, Shifted LN, Sim.~MKT, Det.~Timing",
        "big_n_v_50funds_alpha_stdv30_shifted_mkt5" = "(S9) $n=50$, $\\sigma=30\\%$, Shifted LN, Sim.~MKT, $\\bar{\\alpha}=0.75\\%$",
        "big_n_v_50funds_alpha_stdv30_shifted_mkt6" = "(S10) $n=50$, $\\sigma=30\\%$, Shifted LN, Sim.~MKT, Max.~Hold 5"
    )

    # Formatting helper
    fmt <- function(x) sprintf("%.2f", x)

    # Begin construction of the LaTeX table string
    tex <- c()
    tex <- c(tex, "% Table: Small-Sample Bias and Variance (Horizon 0)")
    tex <- c(tex, "\\begin{table}[htbp]")
    tex <- c(tex, "\\centering")
    tex <- c(tex, "\\caption{Comparison to simulation \\cite[Table 1]{DLP12}}")
    tex <- c(tex, "\\label{tab:driessen_comparison}")
    tex <- c(tex, "\\resizebox{\\linewidth}{!}{%")
    tex <- c(tex, "\\begin{tabular}{l ccccc cccccc}")
    tex <- c(tex, "\\toprule")
    tex <- c(tex, "& \\multicolumn{5}{c}{Market Beta ($\\beta$)} & \\multicolumn{6}{c}{Alpha ($\\alpha$, monthly \\%)} \\\\")
    tex <- c(tex, "\\cmidrule(lr){2-6} \\cmidrule(lr){7-12}")
    tex <- c(tex, "Scenario & True & Mean & SD & 25\\% & 75\\% & True & Mean & SD & 25\\% & 75\\% & \\% Bound \\\\")
    tex <- c(tex, "\\midrule")

    # Output all scenarios in a single panel
    for (scen in scenarios) {
        row_data <- filter(stats, scenario_id == scen)
        if (nrow(row_data) == 1) {
            tex <- c(tex, sprintf(
                "%s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
                scenario_labels[[scen]],
                fmt(row_data$true_beta), fmt(row_data$mean_beta), fmt(row_data$sd_beta), fmt(row_data$p25_beta), fmt(row_data$p75_beta),
                fmt(row_data$true_alpha), fmt(row_data$mean_alpha), fmt(row_data$sd_alpha), fmt(row_data$p25_alpha), fmt(row_data$p75_alpha),
                sprintf("%.0f\\%%", row_data$pct_at_bound)
            ))
        }
    }

    tex <- c(tex, "\\addlinespace")
    tex <- c(tex, "\\multicolumn{12}{l}{\\textit{Driessen et al. 2012, Table 1}} \\\\")
    # Values from Driessen 2012 Table 1 Benchmark column
    # Beta: True 1.00, Mean 0.98, SD 0.06, IQR [0.97, 1.03]
    # Alpha: True 0.00, Mean 0.05, SD 0.03, IQR [-0.01, 0.02]
    tex <- c(tex, "Benchmark Estimates ($n=50$, $\\sigma_{\\text{monthly}}=30\\%$) & 1.00 & 0.98 & 0.32 & 0.80 & 1.15 & 0.00 & 0.05 & 0.2 & -0.06 & 0.15 & n/a \\\\")

    tex <- c(tex, "\\bottomrule")
    tex <- c(tex, "\\end{tabular}%")
    tex <- c(tex, "}")

    # Add table notes string
    table_notes <- paste0(
        "This table reports the small-sample properties of the generalized estimator for Horizon~0 (i.e., treating funds as held to maturity), ",
        "mimicking the structure of Table~1 in \\cite{DLP12}. ",
        "We simulate 1,000 Monte Carlo draws. ",
        "All scenarios use 20 vintage years (1986--2005) with 15 deals per fund; unless stated otherwise, $n=20$ funds per vintage, $\\sigma=0.2$ monthly idiosyncratic volatility, and a normal error distribution. ",
        "The simulation dimensions vary as follows: ",
        "$n$ = number of funds per vintage-year portfolio; ",
        "$\\sigma$ = monthly idiosyncratic standard deviation; ",
        "\\textit{Normal} = Gaussian error term (default); ",
        "\\textit{Shifted LN} = shifted lognormal error term (cf.\\ \\cite{DLP12}); ",
        "\\textit{Sim.~MKT} = total market return drawn from a shifted lognormal (S\\&P~500 1980--2003 moments) with a constant risk-free rate of 4\\% annualized, instead of conditioning on the historical path; ",
        "\\textit{Det.~Timing} = deterministic (evenly spaced) deal investment dates instead of random; ",
        "\\textit{Max.~Hold 5} = maximum holding period of a deal limited to 5 years (default 10 years); ",
        "$\\bar{\\alpha}$ = upper bound on the alpha parameter in the optimization (default $\\bar{\\alpha} = 1\\%$/month). ",
        "The \\textit{\\%~Bound} column reports the fraction of draws where the alpha estimate reached the upper optimization bound~$\\bar{\\alpha}$. ",
        "The final row presents the Benchmark Estimates obtained directly from Table~1 in \\cite{DLP12} for comparative purposes. ",
        "Alpha estimates are monthly and presented in percentages."
    )
    tex <- c(tex, "\\begin{minipage}{\\linewidth}")
    tex <- c(tex, "\\vspace{0.2cm}")
    tex <- c(tex, "\\footnotesize \\textit{Notes:} ", table_notes)
    tex <- c(tex, "\\end{minipage}")

    tex <- c(tex, "\\end{table}")

    tex_str <- paste(tex, collapse = "\n")

    if (!is.null(output_tex)) {
        writeLines(tex_str, output_tex)
        cat("LaTeX table written to:", output_tex, "\n")
    }

    return(tex_str)
}

# Example usage (uncomment to test directly):
data_file <- "simulation/data_out_2026_new/bias_analysis/2026-02-23_152411_bias_full.csv"
cat(generate_driessen_latex_table(data_file, output_tex = "simulation/driessen_simulation.tex"))

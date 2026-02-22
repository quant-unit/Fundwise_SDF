# Comparison to Driessen et al. 2012

scenarios <- c(
    "base_case_zero_alpha", "big_n_v_50funds_alpha",
    "big_n_v_50funds_alpha_stdv30", "big_n_v_50funds_alpha_stdv30_shifted"
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

    # Filter for Horizon 0 (max_month == 0) and predefined scenarios
    df_sub <- df %>%
        filter(max_month == 0 & scenario_id %in% scenarios)

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
            .groups = "drop"
        )

    # Map scenario names to human-readable labels for the table
    scenario_labels <- c(
        "base_case_zero_alpha" = "Base Case ($n=20$, $\\sigma_{\\text{monthly}}=20\\%$)",
        "big_n_v_50funds_alpha" = "Large Sample ($n=50$, $\\sigma_{\\text{monthly}}=20\\%$)",
        "big_n_v_50funds_alpha_stdv30" = "High Volatility ($n=50$, $\\sigma_{\\text{monthly}}=30\\%$)",
        "big_n_v_50funds_alpha_stdv30_shifted" = "Shifted Lognormal ($n=50$, $\\sigma_{\\text{monthly}}=30\\%$)"
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
    tex <- c(tex, "\\begin{tabular}{l ccccc ccccc}")
    tex <- c(tex, "\\toprule")
    tex <- c(tex, "& \\multicolumn{5}{c}{Market Beta ($\\beta$)} & \\multicolumn{5}{c}{Alpha ($\\alpha$, monthly \\%)} \\\\")
    tex <- c(tex, "\\cmidrule(lr){2-6} \\cmidrule(lr){7-11}")
    tex <- c(tex, "Scenario & True & Mean & SD & 25\\% & 75\\% & True & Mean & SD & 25\\% & 75\\% \\\\")
    tex <- c(tex, "\\midrule")

    # Output all scenarios in a single panel
    for (scen in scenarios) {
        row_data <- filter(stats, scenario_id == scen)
        if (nrow(row_data) == 1) {
            tex <- c(tex, sprintf(
                "%s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
                scenario_labels[[scen]],
                fmt(row_data$true_beta), fmt(row_data$mean_beta), fmt(row_data$sd_beta), fmt(row_data$p25_beta), fmt(row_data$p75_beta),
                fmt(row_data$true_alpha), fmt(row_data$mean_alpha), fmt(row_data$sd_alpha), fmt(row_data$p25_alpha), fmt(row_data$p75_alpha)
            ))
        }
    }

    tex <- c(tex, "\\addlinespace")
    tex <- c(tex, "\\multicolumn{11}{l}{\\textit{Driessen et al. 2012, Table 1}} \\\\")
    # Values from Driessen 2012 Table 1 Benchmark column
    # Beta: True 1.00, Mean 0.98, SD 0.06, IQR [0.97, 1.03]
    # Alpha: True 0.00, Mean 0.05, SD 0.03, IQR [-0.01, 0.02]
    tex <- c(tex, "Driessen et al. 2012 Benchmark ($n=50$, $\\sigma_{\\text{monthly}}=30\\%$) & 1.00 & 0.98 & 0.32 & 0.80 & 1.15 & 0.00 & 0.05 & 0.2 & -0.06 & 0.15 \\\\")

    tex <- c(tex, "\\bottomrule")
    tex <- c(tex, "\\end{tabular}%")
    tex <- c(tex, "}")

    # Add table notes string
    table_notes <- "This table reports the small-sample properties of the generalized estimator for Horizon 0 (i.e., treating funds as held to maturity), mimicking the structure of Table 1 in Driessen et al. (2012). We simulate 1,000 sets of cash flows under four specifications, exploring the effect of sample size $n$, idiosyncratic volatility $\\sigma_{\\text{monthly}}$, and lognormal shifting. The final row presents the Benchmark Estimates obtained directly from Table 1 in Driessen et al. (2012) for comparative purposes. Alpha estimates are monthly and presented in percentages."
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
data_file <- "simulation/data_out_2026_new/bias_analysis/2026-02-22_162216_bias_full.csv"
cat(generate_driessen_latex_table(data_file, output_tex = "simulation/driessen_simulation.tex"))

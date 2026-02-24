# results/driessen_comparison.R
# Generates a professional LaTeX table comparing 3-factor and 4-factor Fama-French models
# for horizon=0 and CV.key=ALL across EW and FW weightings.

library(dplyr)
library(tidyr)

# File paths
f_4f_ew <- "results/data_out_2026_02_24/cache_ff3_factors_preqin_alpha_ALL_EW_VYP/2026-02-24_165128_cached_res.csv"
f_4f_fw <- "results/data_out_2026_02_24/cache_ff3_factors_preqin_alpha_ALL_FW_VYP/2026-02-24_164657_cached_res.csv"
f_3f_ew <- "results/data_out_2026_02_24/cache_ff3_factors_preqin_MKT_SMB_HML_EW_VYP/2026-02-24_202316_cached_res.csv"
f_3f_fw <- "results/data_out_2026_02_24/cache_ff3_factors_preqin_MKT_SMB_HML_FW_VYP/2026-02-24_202131_cached_res.csv"

# Load and filter data function
load_and_prep <- function(path, model_name, w_val) {
    if (!file.exists(path)) {
        warning(paste("File not found:", path))
        return(NULL)
    }
    df <- read.csv(path, stringsAsFactors = FALSE)

    # Ensure necessary columns are present even if filled with NA
    cols_needed <- c("MKT", "SMB", "HML", "SE.MKT", "SE.SMB", "SE.HML", "Type", "max.month", "CV.key")
    for (col in cols_needed) {
        if (!col %in% colnames(df)) df[[col]] <- NA
    }
    if (!"Alpha" %in% colnames(df)) {
        df$Alpha <- NA
        df$SE.Alpha <- NA
    }

    df %>%
        filter(max.month == 0, CV.key == "ALL") %>%
        mutate(Model = model_name, Weighting = w_val) %>%
        select(Type, Model, Weighting, MKT, SMB, HML, Alpha, SE.MKT, SE.SMB, SE.HML, SE.Alpha)
}

# Load all 4 files
data_all <- bind_rows(
    load_and_prep(f_3f_ew, "3-Factor", "EW"),
    load_and_prep(f_4f_ew, "4-Factor", "EW"),
    load_and_prep(f_3f_fw, "3-Factor", "FW"),
    load_and_prep(f_4f_fw, "4-Factor", "FW")
)

if (is.null(data_all) || nrow(data_all) == 0) {
    stop("No data loaded. Check file paths.")
}

# Limit to primary asset classes for a clean fit
target_types <- c("BO", "VC", "PE", "RE")
data_filtered <- data_all %>% filter(Type %in% target_types)

# If any type is completely missing, adjust the target types
target_types <- intersect(target_types, unique(data_filtered$Type))

# Format a single value/SE pair
fmt_est <- function(val, se) {
    if (is.na(val)) {
        return(c("", ""))
    }
    c(sprintf("%.3f", val), sprintf("(%.3f)", se))
}

# Make LaTeX table
latex_out <- c()
latex_out <- c(latex_out, "\\begin{table}[ht]")
latex_out <- c(latex_out, "\\centering")
latex_out <- c(latex_out, "\\caption{Estimation Results: 3-Factor vs 4-Factor Models}")
latex_out <- c(latex_out, "\\label{tab:driessen_comparison_empirical}")
latex_out <- c(latex_out, "\\small")
latex_out <- c(latex_out, paste0("\\begin{tabular}{ll", paste(rep("c", length(target_types)), collapse = ""), "}"))
latex_out <- c(latex_out, "\\hline\\hline")

# Header
latex_out <- c(latex_out, paste0(" & & \\multicolumn{", length(target_types), "}{c}{Asset Class} \\\\ \\cline{3-", 2 + length(target_types), "}"))
latex_out <- c(latex_out, paste0("Model & Factor & ", paste(target_types, collapse = " & "), " \\\\ \\hline"))

# Helper for a Panel
build_panel <- function(weighting_label, w_filter) {
    lines <- c()
    lines <- c(lines, paste0("\\multicolumn{", 2 + length(target_types), "}{l}{\\textit{Panel ", ifelse(w_filter == "EW", "A", "B"), ": ", weighting_label, "}} \\\\"))

    for (mod in c("3-Factor", "4-Factor")) {
        lines <- c(lines, paste0("\\multirow{8}{*}{", mod, "}"))

        factors_to_show <- if (mod == "4-Factor") c("Alpha", "MKT", "SMB", "HML") else c("MKT", "SMB", "HML")
        for (fac in factors_to_show) {
            est_row <- paste0(" & ", fac)
            se_row <- " & "

            for (ty in target_types) {
                row_data <- data_filtered %>% filter(Weighting == w_filter, Model == mod, Type == ty)
                if (nrow(row_data) > 0) {
                    val <- row_data[[fac]][1]
                    se_col <- paste0("SE.", fac)
                    se_val <- row_data[[se_col]][1]

                    # Convert Alpha scale if needed. If it's pure decimal, leave as is or multiply by 100 for %.
                    # Driessen reports regular values. Leave as is, just 3 decimals.

                    fmt <- fmt_est(val, se_val)
                    est_row <- paste0(est_row, " & ", fmt[1])
                    se_row <- paste0(se_row, " & ", fmt[2])
                } else {
                    est_row <- paste0(est_row, " & ")
                    se_row <- paste0(se_row, " & ")
                }
            }
            lines <- c(lines, paste0(est_row, " \\\\"))
            lines <- c(lines, paste0(se_row, " \\\\"))
            if (fac != tail(factors_to_show, 1)) {
                lines <- c(lines, paste0(" & \\multicolumn{", length(target_types), "}{c}{} \\vspace{-0.2cm} \\\\")) # Small spacer
            }
        }
        lines <- c(lines, "\\hline")
    }
    return(lines)
}

latex_out <- c(latex_out, build_panel("Equal-Weighted (EW)", "EW"))
latex_out <- c(latex_out, build_panel("Fund-Weighted (FW)", "FW"))

latex_out <- c(latex_out, "\\end{tabular}")
latex_out <- c(latex_out, "\\par")
latex_out <- c(latex_out, "\\vspace{1ex}")
latex_out <- c(latex_out, "\\raggedright \\footnotesize \\textit{Note:} This table presents the estimation results of the Fama-French 3-factor (MKT, SMB, HML) and 4-factor (including $\\alpha$) models for different asset classes. Cash flows are discounted to horizon 0. Standard errors are reported in parentheses below the estimates. The models are estimated using the Least-Mean-Distance method. Panel A uses equal-weighted (EW) cash flows, while Panel B uses fund-size-weighted (FW) cash flows. The upper bound for $\\alpha$ was constrained to 0.02.")
latex_out <- c(latex_out, "\\end{table}")

# Write to file
out_path <- "results/driessen_comparison.tex"
writeLines(latex_out, out_path)
cat("Table successfully written to", out_path, "\n")

# results/driessen_comparison.R
# Generates a professional LaTeX table comparing 3-factor and 4-factor Fama-French models
# for horizon=0 and CV.key=ALL across EW and FW weightings.

library(dplyr)
library(tidyr)

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
}
getwd()

# File paths

# two-factor model: Alpha + MKT
f_2f_ew <- "results/data_out_2026_02_26/cache_ff3_factors_preqin_alpha_Alpha_EW_VYP/2026-02-26_142159_cached_res.csv"
f_2f_fw <- "results/data_out_2026_02_26/cache_ff3_factors_preqin_alpha_Alpha_FW_VYP/2026-02-26_141951_cached_res.csv"
# three-factor model: MKT + SMB + HML
f_3f_ew <- "results/data_out_2026_02_26/cache_ff3_factors_preqin_MKT_SMB_HML_EW_VYP/2026-02-26_053048_cached_res.csv"
f_3f_fw <- "results/data_out_2026_02_26/cache_ff3_factors_preqin_MKT_SMB_HML_FW_VYP/2026-02-26_052834_cached_res.csv"
# four-factor model: Alpha + MKT + SMB + HML
f_4f_ew <- "results/data_out_2026_02_26/cache_ff3_factors_preqin_alpha_ALL_EW_VYP/2026-02-26_041833_cached_res.csv"
f_4f_fw <- "results/data_out_2026_02_26/cache_ff3_factors_preqin_alpha_ALL_FW_VYP/2026-02-26_041328_cached_res.csv"

# folder for cross-validation
base_folder <- "results/data_out_2026_02_26/"
f_2f_ew_folder <- "cache_ff3_factors_preqin_alpha_Alpha_EW_VYP"
f_2f_fw_folder <- "cache_ff3_factors_preqin_alpha_Alpha_FW_VYP"
# three-factor model: MKT + SMB + HML
f_3f_ew_folder <- "cache_ff3_factors_preqin_MKT_SMB_HML_EW_VYP"
f_3f_fw_folder <- "cache_ff3_factors_preqin_MKT_SMB_HML_FW_VYP"
# four-factor model: Alpha + MKT + SMB + HML
f_4f_ew_folder <- "cache_ff3_factors_preqin_alpha_ALL_EW_VYP"
f_4f_fw_folder <- "cache_ff3_factors_preqin_alpha_ALL_FW_VYP"

# Goal: 1. read all files per folder, 2. filter appropriately, 3. determine CV mean and SD, 4. add CV results to table
# Filters for asymptotic results: max.month = 0, CV.key == "ALL" (then "use as is")
# Filters for CV results: max.month = 0, CV.key != "ALL" (then determine MEAN and SD)

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

# Load all files
data_all <- bind_rows(
    load_and_prep(f_2f_ew, "2-Factor", "EW"),
    load_and_prep(f_2f_fw, "2-Factor", "FW"),
    load_and_prep(f_3f_ew, "3-Factor", "EW"),
    load_and_prep(f_3f_fw, "3-Factor", "FW"),
    load_and_prep(f_4f_ew, "4-Factor", "EW"),
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
latex_out <- c(latex_out, "\\caption{Estimation Results: 2-, 3-, and 4-Factor Models}")
latex_out <- c(latex_out, "\\label{tab:driessen_comparison_empirical}")
latex_out <- c(latex_out, "\\resizebox{\\textwidth}{!}{")
latex_out <- c(latex_out, "\\begin{tabular}{l cccc c cccc}")
latex_out <- c(latex_out, "\\hline\\hline")

# Header
latex_out <- c(latex_out, "& \\multicolumn{4}{c}{Equal-Weighted (EW)} && \\multicolumn{4}{c}{Fund-Weighted (FW)} \\\\ \\cline{2-5} \\cline{7-10}")
latex_out <- c(latex_out, "Model & $\\alpha$ & MKT & SMB & HML && $\\alpha$ & MKT & SMB & HML \\\\ \\hline")

# Body
for (ty in target_types) {
    ty_label <- switch(ty,
        "BO" = "Buyout (BO)",
        "VC" = "Venture Capital (VC)",
        "PE" = "Private Equity (PE)",
        "RE" = "Real Estate (RE)",
        ty
    )
    latex_out <- c(latex_out, paste0("\\multicolumn{10}{l}{\\textbf{", ty_label, "}} \\\\"))

    for (mod in c("2-Factor", "3-Factor", "4-Factor")) {
        row_ew <- data_filtered %>% filter(Type == ty, Model == mod, Weighting == "EW")
        row_fw <- data_filtered %>% filter(Type == ty, Model == mod, Weighting == "FW")

        facs <- c("Alpha", "MKT", "SMB", "HML")
        est_cols <- c(mod)
        se_cols <- c("")

        for (w_filter in c("EW", "FW")) {
            if (w_filter == "FW") {
                est_cols <- c(est_cols, "")
                se_cols <- c(se_cols, "")
            }

            row_data <- if (w_filter == "EW") row_ew else row_fw
            for (fac in facs) {
                if (mod == "2-Factor" && fac %in% c("SMB", "HML")) {
                    est_cols <- c(est_cols, "")
                    se_cols <- c(se_cols, "")
                } else if (mod == "3-Factor" && fac == "Alpha") {
                    est_cols <- c(est_cols, "")
                    se_cols <- c(se_cols, "")
                } else if (nrow(row_data) > 0) {
                    val <- row_data[[fac]][1]
                    se_val <- row_data[[paste0("SE.", fac)]][1]
                    fmt <- fmt_est(val, se_val)
                    est_cols <- c(est_cols, fmt[1])
                    se_cols <- c(se_cols, fmt[2])
                } else {
                    est_cols <- c(est_cols, "")
                    se_cols <- c(se_cols, "")
                }
            }
        }

        latex_out <- c(latex_out, paste0(paste(est_cols, collapse = " & "), " \\\\"))
        if (mod == "4-Factor" && ty != tail(target_types, 1)) {
            latex_out <- c(latex_out, paste0(paste(se_cols, collapse = " & "), " \\\\[0.15cm]"))
        } else {
            latex_out <- c(latex_out, paste0(paste(se_cols, collapse = " & "), " \\\\"))
        }
    }
}

latex_out <- c(latex_out, "\\hline")
latex_out <- c(latex_out, "\\end{tabular}")
latex_out <- c(latex_out, "}")
latex_out <- c(latex_out, "\\par")
latex_out <- c(latex_out, "\\vspace{1ex}")
latex_out <- c(latex_out, "\\raggedright \\footnotesize \\textit{Note:} This table presents the estimation results of the 2-factor ($\\alpha$, MKT), Fama-French 3-factor (MKT, SMB, HML), and 4-factor ($\\alpha$, MKT, SMB, HML) models for different asset classes. Cash flows are discounted to horizon 0. Standard errors are reported in parentheses below the estimates. The estimation uses Equal-Weighted (EW) and Fund-Size-Weighted (FW) cash flows. The upper bound for $\\alpha$ was constrained to 0.02.")
latex_out <- c(latex_out, "\\end{table}")

# Write to file
out_path <- "results/driessen_comparison.tex"
writeLines(latex_out, out_path)
cat("Table successfully written to", out_path, "\n")

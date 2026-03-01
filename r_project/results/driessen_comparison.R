# results/driessen_comparison.R
# Generates a professional LaTeX table comparing 3-factor and 4-factor Fama-French models
# for horizon=0 and CV.key=ALL across EW and FW weightings.

library(dplyr)
library(tidyr)

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
}
getwd()

# folder for cross-validation
base_folder <- "results/data_out_2026_02_26/"
f_1f_ew_folder <- file.path(base_folder, "cache_ff3_factors_preqin_MKT_EW_VYP")
f_1f_fw_folder <- file.path(base_folder, "cache_ff3_factors_preqin_MKT_FW_VYP")
f_2f_ew_folder <- file.path(base_folder, "cache_ff3_factors_preqin_alpha_Alpha_EW_VYP")
f_2f_fw_folder <- file.path(base_folder, "cache_ff3_factors_preqin_alpha_Alpha_FW_VYP")
# three-factor model: MKT + SMB + HML
f_3f_ew_folder <- file.path(base_folder, "cache_ff3_factors_preqin_MKT_SMB_HML_EW_VYP")
f_3f_fw_folder <- file.path(base_folder, "cache_ff3_factors_preqin_MKT_SMB_HML_FW_VYP")
# four-factor model: Alpha + MKT + SMB + HML
f_4f_ew_folder <- file.path(base_folder, "cache_ff3_factors_preqin_alpha_ALL_EW_VYP")
f_4f_fw_folder <- file.path(base_folder, "cache_ff3_factors_preqin_alpha_ALL_FW_VYP")

# Goal: 1. read all files per folder (and rbind to one data.frame), 2. filter appropriately, 3. determine CV mean and SD, 4. add CV results to LaTeX table

# Load and filter data function
load_and_prep <- function(folder_path, model_name, w_val) {
    if (!dir.exists(folder_path)) {
        warning(paste("Directory not found:", folder_path))
        return(NULL)
    }

    files <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE)
    if (length(files) == 0) {
        warning(paste("No CSV files found in:", folder_path))
        return(NULL)
    }

    df_list <- lapply(files, function(pf) {
        d <- read.csv(pf, stringsAsFactors = FALSE)
        # Ensure necessary columns are present even if filled with NA
        cols_needed <- c("MKT", "SMB", "HML", "SE.MKT", "SE.SMB", "SE.HML", "Type", "max.month", "CV.key")
        for (col in cols_needed) {
            if (!col %in% colnames(d)) d[[col]] <- NA
        }
        if (!"Alpha" %in% colnames(d)) {
            d$Alpha <- NA
            d$SE.Alpha <- NA
        }
        d$CV.key <- as.character(d$CV.key)
        return(d)
    })

    df <- bind_rows(df_list)

    df %>%
        filter(max.month == 0) %>%
        mutate(Model = model_name, Weighting = w_val) %>%
        select(Type, Model, Weighting, CV.key, MKT, SMB, HML, Alpha, SE.MKT, SE.SMB, SE.HML, SE.Alpha)
}

# Load all files
data_all <- bind_rows(
    load_and_prep(f_1f_ew_folder, "1-Factor", "EW"),
    load_and_prep(f_1f_fw_folder, "1-Factor", "FW"),
    load_and_prep(f_2f_ew_folder, "2-Factor", "EW"),
    load_and_prep(f_2f_fw_folder, "2-Factor", "FW"),
    load_and_prep(f_3f_ew_folder, "3-Factor", "EW"),
    load_and_prep(f_3f_fw_folder, "3-Factor", "FW"),
    load_and_prep(f_4f_ew_folder, "4-Factor", "EW"),
    load_and_prep(f_4f_fw_folder, "4-Factor", "FW")
)

if (is.null(data_all) || nrow(data_all) == 0) {
    stop("No data loaded. Check file paths.")
}

# Limit to primary asset classes for a clean fit
target_types <- c("BO", "VC", "PE")
data_filtered <- data_all %>% filter(Type %in% target_types)

# If any type is completely missing, adjust the target types
target_types <- intersect(target_types, unique(data_filtered$Type))

# Separate into Asymptotic and CV results
data_asymp <- data_filtered %>%
    filter(CV.key == "ALL") %>%
    mutate(ResultType = "Asymp.") %>%
    select(Type, Model, Weighting, ResultType, MKT, SMB, HML, Alpha, SE.MKT, SE.SMB, SE.HML, SE.Alpha)

data_cv <- data_filtered %>%
    filter(CV.key != "ALL") %>%
    mutate(ResultType = "CV") %>%
    group_by(Type, Model, Weighting, ResultType) %>%
    summarize(
        MKT_mean = mean(MKT, na.rm = TRUE),
        SE.MKT = sd(MKT, na.rm = TRUE),
        SMB_mean = mean(SMB, na.rm = TRUE),
        SE.SMB = sd(SMB, na.rm = TRUE),
        HML_mean = mean(HML, na.rm = TRUE),
        SE.HML = sd(HML, na.rm = TRUE),
        Alpha_mean = mean(Alpha, na.rm = TRUE),
        SE.Alpha = sd(Alpha, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    rename(
        MKT = MKT_mean,
        SMB = SMB_mean,
        HML = HML_mean,
        Alpha = Alpha_mean
    ) %>%
    select(Type, Model, Weighting, ResultType, MKT, SMB, HML, Alpha, SE.MKT, SE.SMB, SE.HML, SE.Alpha)

data_final <- bind_rows(data_asymp, data_cv)

# Format a single value/SE pair
fmt_est <- function(val, se) {
    if (is.na(val)) {
        return(c("", ""))
    }
    c(sprintf("%.3f", val), sprintf("(%.3f)", se))
}

# Make LaTeX table generation a reusable function to support split tables
generate_table <- function(types_to_include, out_path, caption, label) {
    latex_out <- c()
    latex_out <- c(latex_out, "\\begin{table}[ht]")
    latex_out <- c(latex_out, "\\centering")
    latex_out <- c(latex_out, paste0("\\caption{", caption, "}"))
    latex_out <- c(latex_out, paste0("\\label{", label, "}"))
    latex_out <- c(latex_out, "\\resizebox{\\textwidth}{!}{")
    latex_out <- c(latex_out, "\\begin{tabular}{l cccc c cccc}")
    latex_out <- c(latex_out, "\\hline\\hline")

    # Header
    latex_out <- c(latex_out, "& \\multicolumn{4}{c}{Equal-Weighted (EW)} && \\multicolumn{4}{c}{Fund-Weighted (FW)} \\\\ \\cline{2-5} \\cline{7-10}")
    latex_out <- c(latex_out, "Model & $\\alpha$ & MKT & SMB & HML && $\\alpha$ & MKT & SMB & HML \\\\ \\hline")

    # Body
    for (ty in types_to_include) {
        ty_label <- switch(ty,
            "BO" = "Buyout (BO)",
            "VC" = "Venture Capital (VC)",
            "PE" = "Private Equity (PE)",
            "RE" = "Real Estate (RE)",
            ty
        )
        latex_out <- c(latex_out, paste0("\\multicolumn{10}{l}{\\textbf{", ty_label, "}} \\\\"))

        for (mod in c("1-Factor", "2-Factor", "3-Factor", "4-Factor")) {
            for (res_type in c("Asymp.", "CV")) {
                row_ew <- data_final %>% filter(Type == ty, Model == mod, ResultType == res_type, Weighting == "EW")
                row_fw <- data_final %>% filter(Type == ty, Model == mod, ResultType == res_type, Weighting == "FW")

                if (nrow(row_ew) == 0 && nrow(row_fw) == 0) {
                    next
                }

                facs <- c("Alpha", "MKT", "SMB", "HML")
                est_cols <- c(paste(mod, paste0("(", res_type, ")")))
                se_cols <- c("")

                for (w_filter in c("EW", "FW")) {
                    if (w_filter == "FW") {
                        est_cols <- c(est_cols, "")
                        se_cols <- c(se_cols, "")
                    }

                    row_data <- if (w_filter == "EW") row_ew else row_fw
                    for (fac in facs) {
                        if (mod == "1-Factor" && fac %in% c("Alpha", "SMB", "HML")) {
                            est_cols <- c(est_cols, "")
                            se_cols <- c(se_cols, "")
                        } else if (mod == "2-Factor" && fac %in% c("SMB", "HML")) {
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
                if (mod == "4-Factor" && res_type == "CV" && ty != tail(types_to_include, 1)) {
                    latex_out <- c(latex_out, paste0(paste(se_cols, collapse = " & "), " \\\\[0.15cm]"))
                } else {
                    latex_out <- c(latex_out, paste0(paste(se_cols, collapse = " & "), " \\\\"))
                }
            }
        }
    }

    latex_out <- c(latex_out, "\\hline")
    latex_out <- c(latex_out, "\\end{tabular}")
    latex_out <- c(latex_out, "}")
    latex_out <- c(latex_out, "\\par")
    latex_out <- c(latex_out, "\\vspace{1ex}")
    latex_out <- c(latex_out, "\\raggedright \\footnotesize \\textit{Note:} This table presents the estimation results of the single-factor (MKT), 2-factor ($\\alpha$, MKT), Fama-French 3-factor (MKT, SMB, HML), and 4-factor ($\\alpha$, MKT, SMB, HML) models. Cash flows are discounted to horizon 0. Both asymptotic (Asymp.) full-sample estimates and cross-validation (CV) means across hold-out sets are reported. Standard errors (or standard deviations of the estimates for CV) are reported in parentheses below the estimates. The estimation uses Equal-Weighted (EW) and Fund-Size-Weighted (FW) cash flows. The upper bound for $\\alpha$ was constrained to 0.02.")
    latex_out <- c(latex_out, "\\end{table}")

    writeLines(latex_out, out_path)
    cat("Table successfully written to", out_path, "\n")
}

# Generate separated tables
generate_table(c("BO", "VC"), "results/driessen_comparison_BO_VC.tex", "Estimation Results: 1-, 2-, 3-, and 4-Factor Models for Buyout and Venture Capital", "tab:driessen_comparison_bo_vc")
generate_table(c("PE"), "results/driessen_comparison_PE.tex", "Estimation Results: 1-, 2-, 3-, and 4-Factor Models for Private Equity", "tab:driessen_comparison_pe")

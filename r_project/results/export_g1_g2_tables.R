# Regenerates all G1 (L2, L1, NAV-discount, vintage-cutoff) and G2 (U_{v,H})
# CSVs, LaTeX tables, and PDFs in results/figures_bovc_us_g1/ from
# results/data_out_2026-emp-c/ using the exact functions defined in
# results/plot_empirical_estimates.R.
# Run from r_project/:  LANG=en_US.UTF-8 Rscript results/export_g1_g2_tables.R [--core-only]
suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(grid)
})
source(here("helper/result_cache.R"))
source(here("helper/fixedb.R"))

# Load only function definitions and helper constants from plot_empirical_estimates.R
# (avoid running its legacy bottom script block pointing at data_out_2026_03_27)
exprs <- parse(here("results/plot_empirical_estimates.R"))
for (ex in exprs) {
  if (is.call(ex) && identical(ex[[1]], as.name("<-"))) {
    lhs <- as.character(ex[[2]])[1]
    if (lhs %in% c("file.folder", "out.folder")) break
    eval(ex, .GlobalEnv)
  }
}

args <- commandArgs(TRUE)
core.only <- "--core-only" %in% args

file.folder <- "results/data_out_2026-emp-c"
out.folder <- "results/figures_bovc_us_g1"
dir.create(out.folder, showWarnings = FALSE, recursive = TRUE)

# Core G1 L2 tables/figures (MKT, Alpha, Multi) for BO and VC
for (ft in c("BO", "VC")) {
  plot_empirical_estimates(
    data_dir = file.folder, fund_type = ft, factors = c("MKT"),
    export_pdf = TRUE, export_csv = TRUE, export_latex = TRUE,
    height = 3.5, width = 6, y.max.mkt = 2.2, y.min.mkt = 0,
    main.linewidth = 1, cex = 0.8, weighting_filter = c("FW", "EW"),
    region = "North America",
    output_file = file.path(out.folder, sprintf("empirical_%s_US_MKT", ft))
  )
  plot_empirical_estimates(
    data_dir = file.folder, fund_type = ft, factors = c("Alpha"),
    export_pdf = TRUE, export_csv = TRUE, export_latex = TRUE,
    height = 5, width = 5, y.max.mkt = 2.5, y.min.mkt = -0.5,
    y.lim.second = list(Alpha = c(-0.02, 0.02)),
    main.linewidth = 1, cex = 0.8, exclude_cv = FALSE,
    weighting_filter = c("FW", "EW"), region = "North America",
    output_file = file.path(out.folder, sprintf("empirical_%s_US_Alpha", ft))
  )
  plot_empirical_estimates(
    data_dir = file.folder, fund_type = ft, factors = c("EG", "IA", "ME", "ROE"),
    export_pdf = TRUE, export_csv = TRUE, export_latex = TRUE,
    height = 5, width = 14, y.max.mkt = 2.5, y.min.mkt = -0.5,
    y.lim.second = list(Alpha = c(-0.02, 0.02)),
    main.linewidth = 1, cex = 0.8, exclude_cv = FALSE,
    weighting_filter = c("FW", "EW"), region = "North America",
    output_file = file.path(out.folder, sprintf("empirical_%s_US_Multi", ft))
  )
}

# Core G2 U_{v,H} tables
export_u_pair_tables(data_dir = file.folder, fund_type = "BO", region = "North America",
  output_file = file.path(out.folder, "upair_BO_US"))
export_u_pair_tables(data_dir = file.folder, fund_type = "VC", region = "North America",
  output_file = file.path(out.folder, "upair_VC_US"))

# NAV discount sensitivity plots (NC0, NC25, NC50, NC75, full NAV)
for (ft in c("BO", "VC")) {
  plot_nav_discounts_combined_mkt(
    data_dir = file.folder, fund_type = ft, vintages = 2010,
    nc_tags = c("", "_NC75", "_NC50", "_NC25", "_NC0"),
    nc_labels = c("NAV as Cashflow", "25% NAV Disc.", "50% NAV Disc.", "75% NAV Disc.", "Full Write-Off"),
    export_pdf = TRUE, export_svg = TRUE, export_csv = TRUE,
    width = 6, height = 3, y.max.mkt = 2.1, y.min.mkt = 0, main.linewidth = 0.8,
    weighting_filter = c("FW", "EW"), region = "North America",
    output_file = file.path(out.folder, sprintf("nav_disc_%s_US_MKT_combined_NC", ft))
  )
}

if (!core.only) {
  # Max-vintage sweeps (100% NAV and 50% NAV discount)
  for (ft in c("BO", "VC")) {
    plot_max_vintage_cutoff_combined_mkt(
      data_dir = file.folder, fund_type = ft, vintages = 2011:2021,
      nc_tags = c("", "_NC50"), nc_labels = c("100% NAV as Cashflow", "50% NAV-discount sample"),
      export_pdf = TRUE, export_svg = TRUE, export_csv = TRUE,
      width = 7, height = 3, y.max.mkt = 3, y.min.mkt = 0, main.linewidth = 0.8,
      weighting_filter = c("FW"), region = "North America",
      output_file = file.path(out.folder, sprintf("max_vintage_%s_US_MKT_combined_NC", ft))
    )
    plot_max_vintage_cutoff_combined_alpha_mkt(
      data_dir = file.folder, fund_type = ft, vintages = 2011:2021,
      nc_tags = c("", "_NC50"), nc_labels = c("100% NAV as Cashflow", "50% NAV-discount sample"),
      export_pdf = TRUE, export_svg = TRUE, export_csv = TRUE,
      width = 7, height = 6, y.max.mkt = 3, y.min.mkt = 0, y.max.alpha = 0.02, y.min.alpha = -0.02,
      main.linewidth = 0.8, weighting_filter = c("FW"), region = "North America",
      output_file = file.path(out.folder, sprintf("max_vintage_%s_US_Alpha_MKT_combined_NC", ft))
    )
    plot_empirical_estimates(
      data_dir = file.folder, fund_type = ft, factors = c("MKT"),
      export_pdf = TRUE, export_csv = TRUE, export_latex = TRUE,
      height = 3, width = 6, y.max.mkt = 2.5, y.min.mkt = -0.5,
      main.linewidth = 1, cex = 0.8, weighting_filter = c("FW", "EW"),
      region = "North America", error_function = "L1_Ridge",
      output_file = file.path(out.folder, sprintf("empirical_%s_US_MKT_L1", ft))
    )
    plot_empirical_estimates(
      data_dir = file.folder, fund_type = ft, factors = c("Alpha"),
      export_pdf = TRUE, export_csv = TRUE, export_latex = TRUE,
      height = 4, width = 5, y.max.mkt = 2.5, y.min.mkt = -0.5,
      y.lim.second = list(Alpha = c(-0.02, 0.02)),
      main.linewidth = 1, cex = 0.8, exclude_cv = FALSE,
      weighting_filter = c("FW", "EW"), region = "North America",
      error_function = "L1_Ridge",
      output_file = file.path(out.folder, sprintf("empirical_%s_US_Alpha_L1", ft))
    )
  }
}
cat("Done exporting G1/G2 tables to", out.folder, "\n")

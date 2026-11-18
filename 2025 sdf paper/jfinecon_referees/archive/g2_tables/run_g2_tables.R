# Gate G2 driver (18 November 2026). Run from r_project/ with LANG=en_US.UTF-8:
#   Rscript "../2025 sdf paper/jfinecon_referees/archive/g2_tables/run_g2_tables.R"
# Builds the U_{v,H} tables (export_u_pair_tables in results/plot_empirical_estimates.R)
# for Buyout and VC, North America, from the r5 results in
# results/data_out_2026-emp-c, into results/figures_bovc_us_g1. Loads only the
# function definitions and libraries of the plotting script.
suppressMessages(library(here))
setwd(here("results"))
for (.e in parse("plot_empirical_estimates.R")) {
  is.def <- is.call(.e) && identical(.e[[1]], as.name("<-")) && is.call(.e[[3]]) && identical(.e[[3]][[1]], as.name("function"))
  is.lib <- is.call(.e) && as.character(.e[[1]])[1] %in% c("library", "source")
  if (is.def || is.lib) suppressMessages(eval(.e, globalenv()))
}
for (.t in c("BO", "VC")) {
  export_u_pair_tables(data_dir = "data_out_2026-emp-c", fund_type = .t, region = "North America",
    output_file = file.path("figures_bovc_us_g1", paste0("upair_", .t, "_US")))
}

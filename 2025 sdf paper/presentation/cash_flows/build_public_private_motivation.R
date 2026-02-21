#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

args_all <- commandArgs(trailingOnly = FALSE)
script_flag <- grep("^--file=", args_all, value = TRUE)
if (length(script_flag) == 0) stop("Could not determine script path.")
script_path_raw <- sub("^--file=", "", script_flag[1])
script_path_raw <- gsub("~\\+~", " ", script_path_raw)
script_path <- normalizePath(script_path_raw, mustWork = FALSE)
script_dir <- dirname(script_path)
project_root <- normalizePath(file.path(script_dir, "..", "..", ".."), mustWork = FALSE)

cashflow_path <- file.path(script_dir, "preqin_cashflows_EW_VYP_NAV.csv")
index_path <- file.path(project_root, "2025 sdf paper", "presentation", "nav_returns", "preqin_indices.csv")
out_fig <- file.path(project_root, "2025 sdf paper", "presentation", "figures", "motivating_public_vs_private.pdf")

if (!file.exists(cashflow_path)) stop(sprintf("Missing: %s", cashflow_path))
if (!file.exists(index_path)) stop(sprintf("Missing: %s", index_path))

theme_prof <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "#e6e6e6", linewidth = 0.35),
    panel.grid.major.y = element_line(color = "#e6e6e6", linewidth = 0.35),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, color = "#4d4d4d"),
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.title = element_text(size = 10)
  )

# ------------------------
# Left panel: PE NAV returns vs MKT
# ------------------------
idx <- read.csv(index_path, check.names = FALSE)
idx$date <- as.Date(idx$date)
pe_col <- "IDX00007_PRIVATE_EQUITY_GLOBAL"
if (!pe_col %in% names(idx)) stop("Missing PE index column in preqin_indices.csv")

pe_ret <- idx %>%
  transmute(
    date,
    pe_level = as.numeric(.data[[pe_col]])
  ) %>%
  arrange(date) %>%
  mutate(pe_ret = pe_level / lag(pe_level) - 1) %>%
  filter(!is.na(pe_ret)) %>%
  select(date, pe_ret)

# Q5 factor source: prefer downloaded file from nav_returns, fallback to local cache
q5_candidate_1 <- file.path(project_root, "2025 sdf paper", "presentation", "nav_returns", "q5_factors_monthly_2024_download.csv")
q5_candidate_2 <- file.path(project_root, "2025 sdf paper", "presentation", "nav_returns", "q5_factors_monthly_2022_fallback.csv")
q5_candidate_3 <- file.path(project_root, "r_project", "empirical", "data_prepared_main", "q5_factors_monthly_2022.csv")
q5_path <- NULL
for (p in c(q5_candidate_1, q5_candidate_2, q5_candidate_3)) {
  if (file.exists(p)) {
    q5_path <- p
    break
  }
}
if (is.null(q5_path)) stop("Could not find q5 factor file.")

q5 <- read.csv(q5_path, check.names = FALSE)
if (!all(c("year", "month", "R_MKT", "R_F") %in% names(q5))) stop("q5 file missing year/month/R_MKT/R_F")

mkt_q <- q5 %>%
  transmute(
    month_date = as.Date(sprintf("%04d-%02d-01", as.integer(year), as.integer(month))),
    mkt_excess_m = as.numeric(R_MKT) / 100,
    rf_m = as.numeric(R_F) / 100,
    mkt_total_m = (as.numeric(R_MKT) + as.numeric(R_F)) / 100
  ) %>%
  mutate(date = ceiling_date(month_date, unit = "quarter") - days(1)) %>%
  group_by(date) %>%
  summarise(
    mkt_excess_q = prod(1 + mkt_excess_m, na.rm = TRUE) - 1,
    rf_q = prod(1 + rf_m, na.rm = TRUE) - 1,
    mkt_total_q = prod(1 + mkt_total_m, na.rm = TRUE) - 1,
    .groups = "drop"
  )

ret_comp <- pe_ret %>%
  inner_join(mkt_q, by = "date") %>%
  transmute(
    date,
    `PE NAV return` = pe_ret,
    `q5 Market total return (R_MKT + R_F)` = mkt_total_q
  ) %>%
  pivot_longer(
    cols = c(`PE NAV return`, `q5 Market total return (R_MKT + R_F)`),
    names_to = "series",
    values_to = "ret"
  ) %>%
  group_by(series) %>%
  arrange(date, .by_group = TRUE) %>%
  mutate(cum_return = cumprod(1 + ret) - 1) %>%
  ungroup()

p_left <- ggplot(ret_comp, aes(x = date, y = cum_return, color = series)) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_line(linewidth = 1.2, alpha = 0.95) +
  scale_color_manual(values = c("PE NAV return" = "#1f78b4", "q5 Market total return (R_MKT + R_F)" = "#e31a1c")) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0.02, 0.08))) +
  labs(
    title = "Public Return View",
    subtitle = "Cumulative returns: PE NAV index vs q5 market total return (includes RF)",
    x = NULL,
    y = "Cumulative return"
  ) +
  theme_prof

# ------------------------
# Right panel: cumulative PE contributions/distributions by fund age
# ------------------------
cf <- read.csv(cashflow_path, check.names = FALSE)
cf$date <- as.Date(cf$Date)
cf_pe <- cf %>%
  filter(type == "PE") %>%
  mutate(
    CON = as.numeric(CON),
    DIS = as.numeric(DIS),
    NAV = as.numeric(NAV)
  ) %>%
  group_by(Fund.ID) %>%
  arrange(date, .by_group = TRUE) %>%
  mutate(
    first_date = min(date, na.rm = TRUE),
    age_q = (year(date) - year(first_date)) * 4 + (quarter(date) - quarter(first_date)),
    age_years = age_q / 4
  ) %>%
  ungroup()

cf_cum <- cf_pe %>%
  group_by(Fund.ID) %>%
  arrange(age_q, .by_group = TRUE) %>%
  mutate(
    cum_contrib = cumsum(-CON),      # paid-in capital shown as positive
    cum_distri = cumsum(DIS)
  ) %>%
  ungroup()

cf_path <- cf_cum %>%
  group_by(age_q, age_years) %>%
  summarise(
    `Cumulative contributions` = mean(cum_contrib, na.rm = TRUE),
    `Cumulative distributions` = mean(cum_distri, na.rm = TRUE),
    `Average NAV` = mean(NAV, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(age_years <= 20) %>%
  pivot_longer(
    cols = c(`Cumulative contributions`, `Cumulative distributions`, `Average NAV`),
    names_to = "series",
    values_to = "value"
  )

p_right <- ggplot(cf_path, aes(x = age_years, y = value, color = series, linetype = series)) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_line(linewidth = 1.25, alpha = 0.95) +
  scale_color_manual(values = c(
    "Cumulative contributions" = "#6a3d9a",
    "Cumulative distributions" = "#33a02c",
    "Average NAV" = "#ff7f00"
  )) +
  scale_linetype_manual(values = c(
    "Cumulative contributions" = "solid",
    "Cumulative distributions" = "solid",
    "Average NAV" = "longdash"
  )) +
  scale_x_continuous(limits = c(0, 20), breaks = seq(0, 20, by = 2.5)) +
  scale_y_continuous(labels = number_format(accuracy = 0.1), expand = expansion(mult = c(0.02, 0.08))) +
  labs(
    title = "Private Cash-Flow View",
    subtitle = sprintf("Average cumulative cash flows and NAV over fund age for PE (n=%d)", length(unique(cf_pe$Fund.ID))),
    x = "Fund age (years)",
    y = "Average cumulative amount"
  ) +
  theme_prof

combined <- p_left + p_right + plot_layout(widths = c(1, 1))

ggsave(out_fig, combined, width = 12, height = 5.8)
cat(sprintf("Wrote figure: %s\n", out_fig))

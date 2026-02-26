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
if (length(script_flag) == 0) {
  stop("Could not determine script path.")
}
script_path_raw <- sub("^--file=", "", script_flag[1])
script_path_raw <- gsub("~\\+~", " ", script_path_raw)
script_path <- normalizePath(script_path_raw, mustWork = FALSE)
script_dir <- dirname(script_path)
project_root <- normalizePath(file.path(script_dir, "..", "..", ".."), mustWork = FALSE)

preqin_path <- file.path(script_dir, "preqin_indices.csv")
out_dir <- script_dir
fig_dir <- file.path(script_dir, "..", "figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

series_map <- c(
  "IDX00007_PRIVATE_EQUITY_GLOBAL" = "PE Global",
  "IDX00008_PRIVATE_EQUITY_NORTH_AMERICA" = "PE North America",
  "IDX00009_PRIVATE_EQUITY_EUROPE" = "PE Europe",
  "IDX00012_BUYOUT_GLOBAL" = "Buyout Global",
  "IDX00020_GROWTH_GLOBAL" = "Growth Global",
  "IDX00025_VC_GLOBAL" = "VC Global"
)

if (!file.exists(preqin_path)) {
  stop(sprintf("Missing input file: %s", preqin_path))
}

preqin_raw <- read.csv(preqin_path, check.names = FALSE)
preqin_raw$date <- as.Date(preqin_raw$date)

missing_cols <- setdiff(names(series_map), names(preqin_raw))
if (length(missing_cols) > 0) {
  stop(sprintf("Missing required index columns: %s", paste(missing_cols, collapse = ", ")))
}

nav_long <- preqin_raw %>%
  select(date, all_of(names(series_map))) %>%
  pivot_longer(cols = -date, names_to = "series_id", values_to = "level") %>%
  mutate(
    level = as.numeric(level),
    series = unname(series_map[series_id])
  ) %>%
  arrange(series, date)

nav_returns <- nav_long %>%
  group_by(series, series_id) %>%
  mutate(nav_ret = level / lag(level) - 1) %>%
  filter(!is.na(nav_ret), !is.na(level), !is.na(lag(level))) %>%
  ungroup()

max_lag <- 8

calc_acf_summary <- function(x, lag_max) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  n <- length(x)
  if (n <= lag_max + 1) {
    return(list(
      n_obs = n,
      crit = NA_real_,
      n_sig_lags = NA_integer_,
      max_sig_lag = NA_integer_,
      sig_lags = NA_character_
    ))
  }
  acv <- as.numeric(acf(x, lag.max = lag_max, plot = FALSE)$acf)[-1]
  lags <- seq_len(lag_max)
  crit <- 1.96 / sqrt(n)
  sig <- abs(acv) > crit
  sig_lags <- lags[sig]
  consec_sig <- 0L
  for (k in lags) {
    if (abs(acv[k]) > crit) {
      consec_sig <- consec_sig + 1L
    } else {
      break
    }
  }
  list(
    n_obs = n,
    crit = crit,
    n_sig_lags = sum(sig),
    consec_sig_lags = consec_sig,
    max_sig_lag = ifelse(length(sig_lags) == 0, 0L, max(sig_lags)),
    sig_lags = ifelse(length(sig_lags) == 0, "", paste(sig_lags, collapse = ", "))
  )
}

acf_summary <- nav_returns %>%
  group_by(series, series_id) %>%
  group_modify(~{
    s <- calc_acf_summary(.x$nav_ret, max_lag)
    tibble(
      n_obs = s$n_obs,
      crit = s$crit,
      n_sig_lags = s$n_sig_lags,
      consec_sig_lags = s$consec_sig_lags,
      max_sig_lag = s$max_sig_lag,
      sig_lags = s$sig_lags
    )
  }) %>%
  ungroup() %>%
  arrange(series)

q5_download_path <- file.path(out_dir, "q5_factors_monthly_2024_download.csv")
q5_fallback_copy_path <- file.path(out_dir, "q5_factors_monthly_2022_fallback.csv")
q5_url <- "https://global-q.org/uploads/1/2/2/6/122679606/q5_factors_monthly_2024.csv"
q5_source <- q5_download_path
q5_source_note <- "Downloaded from global-q.org"

ok_download <- FALSE
try({
  utils::download.file(q5_url, q5_download_path, mode = "wb", quiet = TRUE)
  ok_download <- file.exists(q5_download_path) && file.info(q5_download_path)$size > 1000
}, silent = TRUE)

if (!ok_download) {
  fallback_path <- file.path(project_root, "r_project", "empirical", "data_prepared_main", "q5_factors_monthly_2022.csv")
  if (!file.exists(fallback_path)) {
    stop("Could not download q5 factors and local fallback file is missing.")
  }
  q5_source <- fallback_path
  q5_source_note <- "Fallback to local cached global-q file (2022 vintage)"
  file.copy(fallback_path, q5_fallback_copy_path, overwrite = TRUE)
}

q5_raw <- read.csv(q5_source, check.names = FALSE)
required_q5 <- c("year", "month", "R_MKT", "R_F")
if (!all(required_q5 %in% names(q5_raw))) {
  stop("q5 factor file does not contain required columns year, month, R_MKT, R_F.")
}

q5_quarterly <- q5_raw %>%
  transmute(
    month_date = as.Date(sprintf("%04d-%02d-01", as.integer(year), as.integer(month))),
    mkt_ret_m = as.numeric(R_MKT) / 100,
    rf_m = as.numeric(R_F) / 100
  ) %>%
  mutate(qtr_end = ceiling_date(month_date, unit = "quarter") - days(1)) %>%
  group_by(qtr_end) %>%
  summarise(
    mkt_ret_q = prod(1 + mkt_ret_m, na.rm = TRUE) - 1,
    rf_q = prod(1 + rf_m, na.rm = TRUE) - 1,
    .groups = "drop"
  )

build_dimson <- function(df_ret, series_name, lag_count) {
  d0 <- df_ret %>%
    filter(series == series_name) %>%
    transmute(date, y = nav_ret) %>%
    inner_join(q5_quarterly, by = c("date" = "qtr_end")) %>%
    mutate(y = y - rf_q) %>%
    arrange(date)

  if (nrow(d0) < 20) {
    return(tibble(
      series = series_name,
      n_obs = nrow(d0),
      dimson_lags = NA_integer_,
      beta_naive = NA_real_,
      beta_dimson = NA_real_,
      beta_uplift = NA_real_
    ))
  }

  lag_count <- as.integer(max(0, lag_count))

  d <- d0
  for (k in 0:lag_count) {
    d[[sprintf("mkt_lag_%d", k)]] <- dplyr::lag(d$mkt_ret_q, n = k)
  }
  d <- d %>% filter(if_all(starts_with("mkt_lag_"), ~ !is.na(.x)))

  if (nrow(d) < 20) {
    return(tibble(
      series = series_name,
      n_obs = nrow(d),
      dimson_lags = lag_count,
      beta_naive = NA_real_,
      beta_dimson = NA_real_,
      beta_uplift = NA_real_
    ))
  }

  fit_naive <- lm(y ~ mkt_ret_q, data = d0)
  rhs <- paste(sprintf("mkt_lag_%d", 0:lag_count), collapse = " + ")
  fit_dimson <- lm(as.formula(paste("y ~", rhs)), data = d)

  beta_naive <- unname(coef(fit_naive)["mkt_ret_q"])
  beta_dimson <- sum(coef(fit_dimson)[grep("^mkt_lag_", names(coef(fit_dimson)))])

  tibble(
    series = series_name,
    n_obs = nrow(d),
    dimson_lags = lag_count,
    beta_naive = beta_naive,
    beta_dimson = beta_dimson,
    beta_uplift = beta_dimson - beta_naive
  )
}

lag_lookup <- acf_summary %>%
  select(series, consec_sig_lags) %>%
  mutate(consec_sig_lags = as.integer(ifelse(is.na(consec_sig_lags), 0L, pmin(consec_sig_lags, 8L))))

dimson_summary <- lag_lookup %>%
  rowwise() %>%
  do(build_dimson(nav_returns, .$series, .$consec_sig_lags)) %>%
  ungroup() %>%
  arrange(series)

build_dimson_no_intercept <- function(df_ret, series_name, lag_count) {
  d0 <- df_ret %>%
    filter(series == series_name) %>%
    transmute(date, y = nav_ret) %>%
    inner_join(q5_quarterly, by = c("date" = "qtr_end")) %>%
    mutate(y = y - rf_q) %>%
    arrange(date)

  if (nrow(d0) < 20) {
    return(tibble(
      series = series_name,
      n_obs = nrow(d0),
      dimson_lags = NA_integer_,
      beta_naive = NA_real_,
      beta_dimson = NA_real_,
      beta_uplift = NA_real_
    ))
  }

  lag_count <- as.integer(max(0, lag_count))

  d <- d0
  for (k in 0:lag_count) {
    d[[sprintf("mkt_lag_%d", k)]] <- dplyr::lag(d$mkt_ret_q, n = k)
  }
  d <- d %>% filter(if_all(starts_with("mkt_lag_"), ~ !is.na(.x)))

  if (nrow(d) < 20) {
    return(tibble(
      series = series_name,
      n_obs = nrow(d),
      dimson_lags = lag_count,
      beta_naive = NA_real_,
      beta_dimson = NA_real_,
      beta_uplift = NA_real_
    ))
  }

  fit_naive <- lm(y ~ mkt_ret_q - 1, data = d0)
  rhs <- paste(sprintf("mkt_lag_%d", 0:lag_count), collapse = " + ")
  fit_dimson <- lm(as.formula(paste("y ~", rhs, "- 1")), data = d)

  beta_naive <- unname(coef(fit_naive)["mkt_ret_q"])
  beta_dimson <- sum(coef(fit_dimson)[grep("^mkt_lag_", names(coef(fit_dimson)))])

  tibble(
    series = series_name,
    n_obs = nrow(d),
    dimson_lags = lag_count,
    beta_naive = beta_naive,
    beta_dimson = beta_dimson,
    beta_uplift = beta_dimson - beta_naive
  )
}

dimson_summary_no_intercept <- lag_lookup %>%
  rowwise() %>%
  do(build_dimson_no_intercept(nav_returns, .$series, .$consec_sig_lags)) %>%
  ungroup() %>%
  arrange(series)

build_dimson_alpha_mkt <- function(df_ret, series_name, lag_count) {
  d0 <- df_ret %>%
    filter(series == series_name) %>%
    transmute(date, y = nav_ret) %>%
    inner_join(q5_quarterly, by = c("date" = "qtr_end")) %>%
    mutate(y = y - rf_q) %>%
    arrange(date)

  if (nrow(d0) < 20) {
    return(tibble(
      series = series_name,
      n_obs = nrow(d0),
      dimson_lags = as.integer(lag_count),
      alpha_q = NA_real_,
      alpha_q_se = NA_real_,
      alpha_q_ci_low = NA_real_,
      alpha_q_ci_high = NA_real_,
      beta_dimson = NA_real_,
      beta_dimson_se = NA_real_,
      beta_dimson_ci_low = NA_real_,
      beta_dimson_ci_high = NA_real_
    ))
  }

  lag_count <- as.integer(max(0, lag_count))
  d <- d0
  for (k in 0:lag_count) {
    d[[sprintf("mkt_lag_%d", k)]] <- dplyr::lag(d$mkt_ret_q, n = k)
  }
  d <- d %>% filter(if_all(starts_with("mkt_lag_"), ~ !is.na(.x)))

  if (nrow(d) < 20) {
    return(tibble(
      series = series_name,
      n_obs = nrow(d),
      dimson_lags = lag_count,
      alpha_q = NA_real_,
      alpha_q_se = NA_real_,
      alpha_q_ci_low = NA_real_,
      alpha_q_ci_high = NA_real_,
      beta_dimson = NA_real_,
      beta_dimson_se = NA_real_,
      beta_dimson_ci_low = NA_real_,
      beta_dimson_ci_high = NA_real_
    ))
  }

  rhs <- paste(sprintf("mkt_lag_%d", 0:lag_count), collapse = " + ")
  fit <- lm(as.formula(paste("y ~", rhs)), data = d)
  coef_fit <- coef(fit)
  vc <- vcov(fit)
  z <- qnorm(0.975)

  alpha_q <- unname(coef_fit["(Intercept)"])
  alpha_q_se <- sqrt(max(0, vc["(Intercept)", "(Intercept)"]))

  lag_names <- grep("^mkt_lag_", names(coef_fit), value = TRUE)
  beta_dimson <- sum(coef_fit[lag_names])
  beta_dimson_se <- sqrt(max(0, sum(vc[lag_names, lag_names, drop = FALSE])))

  tibble(
    series = series_name,
    n_obs = nrow(d),
    dimson_lags = lag_count,
    alpha_q = alpha_q,
    alpha_q_se = alpha_q_se,
    alpha_q_ci_low = alpha_q - z * alpha_q_se,
    alpha_q_ci_high = alpha_q + z * alpha_q_se,
    beta_dimson = beta_dimson,
    beta_dimson_se = beta_dimson_se,
    beta_dimson_ci_low = beta_dimson - z * beta_dimson_se,
    beta_dimson_ci_high = beta_dimson + z * beta_dimson_se
  )
}

fund_type_series <- c("PE Global", "Buyout Global", "Growth Global", "VC Global")

dimson_alpha_mkt_summary <- lag_lookup %>%
  filter(series %in% fund_type_series) %>%
  rowwise() %>%
  do(build_dimson_alpha_mkt(nav_returns, .$series, .$consec_sig_lags)) %>%
  ungroup() %>%
  mutate(series = factor(series, levels = fund_type_series)) %>%
  arrange(series) %>%
  mutate(
    series_label = paste0(as.character(series), " (L=", dimson_lags, ")"),
    alpha_q_label = ifelse(
      is.finite(alpha_q),
      paste0("alpha_q=", sprintf("%.3f%%", 100 * alpha_q)),
      "alpha_q=NA"
    ),
    series_label_with_alpha = paste0(series_label, "\n", alpha_q_label)
  )

dimson_alpha_all <- lag_lookup %>%
  rowwise() %>%
  do(build_dimson_alpha_mkt(nav_returns, .$series, .$consec_sig_lags)) %>%
  ungroup() %>%
  select(series, alpha_q)

acf_global <- nav_returns %>%
  filter(series == "PE Global") %>%
  arrange(date)

acf_obj <- acf(acf_global$nav_ret, lag.max = max_lag, plot = FALSE)
acf_plot_df <- tibble(
  lag = 1:max_lag,
  acf = as.numeric(acf_obj$acf)[2:(max_lag + 1)]
)
acf_crit <- 1.96 / sqrt(nrow(acf_global))

p_acf <- ggplot(acf_plot_df, aes(x = lag, y = acf)) +
  geom_col(fill = "#1f78b4", width = 0.6) +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.3) +
  geom_hline(yintercept = c(-acf_crit, acf_crit), linetype = "dashed", color = "firebrick") +
  scale_x_continuous(breaks = 1:max_lag) +
  labs(
    title = "ACF of PE Global Quarterly NAV Returns",
    subtitle = sprintf("Dashed lines: +/- 1.96/sqrt(N), N = %d", nrow(acf_global)),
    x = "Lag (quarters)",
    y = "Autocorrelation"
  ) +
  theme_minimal(base_size = 11)

beta_plot_df <- dimson_summary %>%
  left_join(dimson_alpha_all, by = "series") %>%
  mutate(series_label = paste0(series, " (L=", dimson_lags, ")")) %>%
  select(series, series_label, dimson_lags, beta_naive, beta_dimson) %>%
  pivot_longer(cols = c(beta_naive, beta_dimson), names_to = "method", values_to = "beta") %>%
  mutate(method = recode(
    method,
    beta_naive = "Naive OLS with intercept (lag 0)",
    beta_dimson = "Dimson with intercept sum(beta_0...beta_L)"
  ))

beta_label_expr <- dimson_summary %>%
  left_join(dimson_alpha_all, by = "series") %>%
  mutate(
    series_label = paste0(series, " (L=", dimson_lags, ")"),
    alpha_pct = ifelse(is.finite(alpha_q), sprintf("%.3f%%", 100 * alpha_q), "NA"),
    label_expr = paste0("atop('", series_label, "', alpha[q]=='", alpha_pct, "')")
  ) %>%
  select(series_label, label_expr)

p_beta <- ggplot(beta_plot_df, aes(x = series_label, y = beta, fill = method)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  scale_x_discrete(
    labels = function(x) {
      parse(text = beta_label_expr$label_expr[match(x, beta_label_expr$series_label)])
    }
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  coord_flip() +
  scale_fill_manual(values = c(
    "Naive OLS with intercept (lag 0)" = "#b2df8a",
    "Dimson with intercept sum(beta_0...beta_L)" = "#fb9a99"
    )) +
  labs(
    title = "Naive vs Dimson Beta to q5 Market Factor (alpha + MKT)",
    subtitle = "Target is excess PE return: R_PE - R_F; labels include Dimson lag length L and quarterly alpha",
    x = NULL,
    y = "Estimated beta",
    fill = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

combined_plot <- p_acf + p_beta + plot_layout(widths = c(1, 1.2))

beta_plot_df_no_intercept <- dimson_summary_no_intercept %>%
  mutate(series_label = paste0(series, " (L=", dimson_lags, ")")) %>%
  select(series, series_label, dimson_lags, beta_naive, beta_dimson) %>%
  pivot_longer(cols = c(beta_naive, beta_dimson), names_to = "method", values_to = "beta") %>%
  mutate(method = recode(
    method,
    beta_naive = "Naive OLS no intercept (lag 0)",
    beta_dimson = "Dimson no intercept sum(beta_0...beta_L)"
  ))

p_beta_no_intercept <- ggplot(beta_plot_df_no_intercept, aes(x = series_label, y = beta, fill = method)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  coord_flip() +
  scale_fill_manual(values = c(
    "Naive OLS no intercept (lag 0)" = "#33a02c", 
    "Dimson no intercept sum(beta_0...beta_L)" = "#ff7f00"
  )) +
  labs(
    title = "Naive vs Dimson Beta to q5 Market Factor (MKT-only)",
    subtitle = "Target is excess PE return: R_PE - R_F; series labels include Dimson lag length L",
    x = NULL,
    y = "Estimated beta",
    fill = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

combined_plot_no_intercept <- p_acf + p_beta_no_intercept + plot_layout(widths = c(1, 1.2))

twofactor_plot_df <- bind_rows(
  dimson_alpha_mkt_summary %>%
    transmute(
      series_label_with_alpha,
      term = "Alpha (% per quarter)",
      estimate = 100 * alpha_q,
      ci_low = 100 * alpha_q_ci_low,
      ci_high = 100 * alpha_q_ci_high
    ),
  dimson_alpha_mkt_summary %>%
    transmute(
      series_label_with_alpha,
      term = "Dimson beta to MKT",
      estimate = beta_dimson,
      ci_low = beta_dimson_ci_low,
      ci_high = beta_dimson_ci_high
    )
) %>%
  mutate(
    term = factor(term, levels = c("Dimson beta to MKT", "Alpha (% per quarter)")),
    series_label_with_alpha = factor(series_label_with_alpha, levels = rev(unique(dimson_alpha_mkt_summary$series_label_with_alpha)))
  )

p_twofactor <- ggplot(twofactor_plot_df, aes(x = series_label_with_alpha, y = estimate, ymin = ci_low, ymax = ci_high, color = term)) +
  geom_hline(yintercept = 0, color = "grey45", linewidth = 0.3) +
  geom_pointrange(size = 0.45, fatten = 2.3) +
  coord_flip() +
  facet_wrap(~term, ncol = 1, scales = "free_y") +
  scale_color_manual(values = c("Dimson beta to MKT" = "#1f78b4", "Alpha (% per quarter)" = "#e31a1c")) +
  labs(
    title = "Dimson Two-Factor Results by Fund Type",
    subtitle = "Model: (R_PE - R_F)_t = alpha + sum_{k=0}^{L} beta_k MKT_{t-k} + e_t (95% CI)",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9.5),
    axis.text.y = element_text(size = 9)
  )

plot_path <- file.path(fig_dir, "motivating_nav_acf_dimson.pdf")
ggsave(plot_path, combined_plot, width = 12, height = 5.8)
plot_no_intercept_path <- file.path(fig_dir, "motivating_nav_acf_dimson_no_intercept.pdf")
ggsave(plot_no_intercept_path, combined_plot_no_intercept, width = 12, height = 5.8)
plot_twofactor_path <- file.path(fig_dir, "dimson_alpha_mkt_fundtypes_square.pdf")
ggsave(plot_twofactor_path, p_twofactor, width = 7, height = 7)

write.csv(nav_returns, file.path(out_dir, "preqin_nav_returns_quarterly.csv"), row.names = FALSE)
write.csv(acf_summary, file.path(out_dir, "acf_significant_lags.csv"), row.names = FALSE)
write.csv(dimson_summary, file.path(out_dir, "dimson_beta_summary.csv"), row.names = FALSE)
write.csv(dimson_summary_no_intercept, file.path(out_dir, "dimson_beta_summary_no_intercept.csv"), row.names = FALSE)
write.csv(dimson_alpha_mkt_summary, file.path(out_dir, "dimson_alpha_mkt_fundtypes_summary.csv"), row.names = FALSE)

report_lines <- c(
  sprintf("q5_source_note: %s", q5_source_note),
  sprintf("q5_source_path: %s", q5_source),
  sprintf("acf_max_lag_tested: %d", max_lag),
  "",
  "acf_significant_lags_by_series:",
  paste0(
    "  - ", acf_summary$series,
    ": n_sig_lags=", acf_summary$n_sig_lags,
    ", consec_sig_lags=", acf_summary$consec_sig_lags,
    ", sig_lags=[", acf_summary$sig_lags, "]"
  ),
  "",
  "dimson_summary:",
  paste0(
    "  - ", dimson_summary$series,
    ": L=", dimson_summary$dimson_lags,
    ", beta_naive=", sprintf("%.3f", dimson_summary$beta_naive),
    ", beta_dimson=", sprintf("%.3f", dimson_summary$beta_dimson),
    ", uplift=", sprintf("%.3f", dimson_summary$beta_uplift)
  ),
  "",
  "dimson_summary_no_intercept:",
  paste0(
    "  - ", dimson_summary_no_intercept$series,
    ": L=", dimson_summary_no_intercept$dimson_lags,
    ", beta_naive=", sprintf("%.3f", dimson_summary_no_intercept$beta_naive),
    ", beta_dimson=", sprintf("%.3f", dimson_summary_no_intercept$beta_dimson),
    ", uplift=", sprintf("%.3f", dimson_summary_no_intercept$beta_uplift)
  ),
  "",
  "dimson_alpha_mkt_summary:",
  paste0(
    "  - ", dimson_alpha_mkt_summary$series,
    ": L=", dimson_alpha_mkt_summary$dimson_lags,
    ", alpha_q=", sprintf("%.4f", dimson_alpha_mkt_summary$alpha_q),
    ", beta_dimson=", sprintf("%.3f", dimson_alpha_mkt_summary$beta_dimson)
  )
)
writeLines(report_lines, con = file.path(out_dir, "analysis_report.txt"))

cat("Wrote outputs to:\n")
cat(sprintf("- %s\n", file.path(out_dir, "preqin_nav_returns_quarterly.csv")))
cat(sprintf("- %s\n", file.path(out_dir, "acf_significant_lags.csv")))
cat(sprintf("- %s\n", file.path(out_dir, "dimson_beta_summary.csv")))
cat(sprintf("- %s\n", file.path(out_dir, "dimson_beta_summary_no_intercept.csv")))
cat(sprintf("- %s\n", file.path(out_dir, "dimson_alpha_mkt_fundtypes_summary.csv")))
cat(sprintf("- %s\n", file.path(out_dir, "analysis_report.txt")))
cat(sprintf("- %s\n", plot_path))
cat(sprintf("- %s\n", plot_no_intercept_path))
cat(sprintf("- %s\n", plot_twofactor_path))

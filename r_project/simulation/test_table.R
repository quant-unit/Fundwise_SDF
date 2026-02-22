library(dplyr)
library(tidyr)
library(knitr)
options(dplyr.summarise.inform = FALSE)

scenarios <- c(
  "base_case_zero_alpha", "big_n_v_50funds_alpha",
  "big_n_v_50funds_alpha_stdv30", "big_n_v_50funds_alpha_stdv30_shifted"
)

df <- read.csv("/Users/christausch/Library/Mobile Documents/com~apple~CloudDocs/Projects/Fundwise_SDF/r_project/simulation/data_out_2026_new/bias_analysis/2026-02-22_162216_bias_full.csv")

# Filter for the relevant scenarios and horizon 0
df_sub <- df %>%
  filter(max_month == 0 & scenario_id %in% scenarios)

if (nrow(df_sub) == 0) {
  # maybe horizon is 1 for some? Let's check what horizons are available for these
  df_sub <- df %>% filter(scenario_id %in% scenarios)
  print(table(df_sub$scenario_id, df_sub$max_month))
  quit()
}

stats <- df_sub %>%
  group_by(scenario_id) %>%
  summarise(
    n_sim = n(),
    true_beta = mean(true_beta_MKT, na.rm = TRUE),
    true_alpha = mean(true_second, na.rm = TRUE) * 12 * 100, # Assuming alpha in CSV is monthly and we want annualized % for table maybe? Let's print raw first

    mean_beta = mean(est_beta_MKT, na.rm = TRUE),
    sd_beta = sd(est_beta_MKT, na.rm = TRUE),
    p25_beta = quantile(est_beta_MKT, 0.25, na.rm = TRUE),
    p75_beta = quantile(est_beta_MKT, 0.75, na.rm = TRUE),
    mean_alpha = mean(est_second, na.rm = TRUE),
    sd_alpha = sd(est_second, na.rm = TRUE),
    p25_alpha = quantile(est_second, 0.25, na.rm = TRUE),
    p75_alpha = quantile(est_second, 0.75, na.rm = TRUE)
  )

print(stats)

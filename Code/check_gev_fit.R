# Run one section at a time in RStudio, from the WWA_EU_heatwave_2026 project.
# Each call checks Annual and June using bootstrap KS and Ljung-Box.

# 1. Functions and settings ----
source(here::here("Code", "helpers", "fitcheck.R"))
B <- 999L
seed <- 20261007L
lag <- 5L  # Prespecified: first five annual lags.

# 2. ERA5 Tx3x ----
era5_tx_data <- read.csv(here::here("results_era5_tx", "annual_extremes_and_gmst.csv"))
era5_tx_check <- check_gev_fit(era5_tx_data, product = "ERA5", index = "Tx3x",
                              B = B, lag = lag, seed = seed)
era5_tx_check$summary

# 3. CPC Tx3x ----
cpc_tx_data <- read.csv(here::here("results_cpc_tmax", "annual_extremes_and_gmst.csv"))
cpc_tx_check <- check_gev_fit(cpc_tx_data, product = "CPC", index = "Tx3x",
                             B = B, lag = lag, seed = seed + 2L)
cpc_tx_check$summary

# 4. Berkeley Earth Tx3x ----
berkeley_tx_data <- read.csv(here::here("results_berkeley_tx", "annual_extremes_and_gmst.csv"))
berkeley_tx_check <- check_gev_fit(berkeley_tx_data, product = "Berkeley Earth", index = "Tx3x",
                                  B = B, lag = lag, seed = seed + 4L)
berkeley_tx_check$summary

# 5. E-OBS Tx3x ----
eobs_tx_data <- read.csv(here::here("results_eobs_tx", "annual_extremes_and_gmst.csv"))
eobs_tx_check <- check_gev_fit(eobs_tx_data, product = "E-OBS", index = "Tx3x",
                              B = B, lag = lag, seed = seed + 6L)
eobs_tx_check$summary

# 6. CPC Tn3x ----
cpc_tn_data <- read.csv(here::here("results_cpc_tmin", "annual_extremes_and_gmst.csv"))
cpc_tn_check <- check_gev_fit(cpc_tn_data, product = "CPC", index = "Tn3x",
                             B = B, lag = lag, seed = seed + 8L)
cpc_tn_check$summary

# 7. E-OBS Tn3x ----
eobs_tn_data <- read.csv(here::here("results_eobs_tn", "annual_extremes_and_gmst.csv"))
eobs_tn_check <- check_gev_fit(eobs_tn_data, product = "E-OBS", index = "Tn3x",
                              B = B, lag = lag, seed = seed + 10L)
eobs_tn_check$summary

# 8. Combine the 12 results ----
fit_results <- rbind(era5_tx_check$summary, cpc_tx_check$summary,
                    berkeley_tx_check$summary, eobs_tx_check$summary,
                    cpc_tn_check$summary, eobs_tn_check$summary)
rownames(fit_results) <- NULL
fit_results
# KS_p and LB_p are per-fit bootstrap p-values, conditional on successful refits.
# Inspect KS_p_lower/upper and LB_p_lower/upper if refits fail. No multiplicity adjustment.

# Optional examples: inspect any fit's PIT or failure messages.
# cpc_tn_check$details$June$PIT
# table(cpc_tn_check$details$June$failure, useNA = "no")

# 9. Save summary and detailed results (optional) ----
dir.create(here::here("diagnostics"), showWarnings = FALSE)
write.csv(fit_results, here::here("diagnostics", "gev_fitcheck_all.csv"),
          row.names = FALSE)
saveRDS(list(ERA5_Tx3x = era5_tx_check, CPC_Tx3x = cpc_tx_check,
             Berkeley_Tx3x = berkeley_tx_check, EOBS_Tx3x = eobs_tx_check,
             CPC_Tn3x = cpc_tn_check, EOBS_Tn3x = eobs_tn_check),
        here::here("diagnostics", "gev_fitcheck_all.rds"))


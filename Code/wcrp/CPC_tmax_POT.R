# CPC summer Tmax POT attribution: numerical results only.
# Run one numbered section at a time in RStudio and inspect intermediate objects.
# Functions are defined in Code/helpers/pot.R. This script creates no figures.

# ============================================================================
# 1. Paths, helper functions, and settings
# ============================================================================
project_root <- "C:/Users/siyan/WCRP/wcrp_common_attribution_assessment"


daily_file <- file.path(project_root, "outputs/pnw_heatwave",
                        "CPC_reference_daily_tmax.nc")
gmst_file <- file.path(project_root, "data", "processed", "gmst_smoothed.csv")
output_dir <- file.path(project_root, "Results", "CPC_POT")

source(file.path(project_root, "Code", "helpers", "pot.R"))
source(file.path(project_root, "Code", "helpers", "trend.R"))

start_year <- 1979L
end_year <- 2025L
season_start <- "05-15"
season_end <- "09-15"
season_label <- "May15-Sep15"
event_year <- 2021L
factual_year <- 2021L
counterfactual_year <- 1951L
base_temperature <- 30
threshold_probability <- 0.90
quantile_type <- 7L
run_days <- 3L
mk_alpha <- 0.05
# bootstrap_B <- 1000L
# random_seed <- 20261111L
# ci_level <- 0.95

stopifnot(file.exists(daily_file), file.exists(gmst_file),
          requireNamespace("ncdf4", quietly = TRUE), start_year <= end_year,
          threshold_probability > 0, threshold_probability < 1,
          run_days >= 1L, run_days == as.integer(run_days),
          # bootstrap_B >= 0L, bootstrap_B == as.integer(bootstrap_B),
          # ci_level > 0, ci_level < 1,
          length(season_start) == 1L, length(season_end) == 1L)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# 2. Read daily Tmax and select summer BEFORE defining the threshold
# ============================================================================
# CPC_reference_daily_tmax.nc already contains the PNW area-mean daily Tmax.
# The season is May 15 through September 15, inclusive (124 days per year).
# No three-day averaging is applied.
daily <- read_pot_daily(daily_file)
summer <- select_pot_season(daily, start_year, end_year, season_start, season_end)
rownames(summer) <- NULL
stopifnot(nrow(summer) > 0L)
gmst <- read.csv(gmst_file)

# Each complete summer contributes one exposure unit, including zero-event summers.
# Exclude incomplete summers from threshold selection, fitting and bootstrapping.
summer_info <- summarize_pot_summers(summer, gmst, season_start = season_start, season_end = season_end)
years_fit <- summer_info[summer_info$used_in_fit, ]
stopifnot(nrow(years_fit) > 4L, length(unique(years_fit$gmst)) > 1L)
summer$used_in_fit <- summer$year %in% years_fit$year
summer_fit <- summer[summer$used_in_fit, ]
summer_info

# ============================================================================
# 2a. Mann-Kendall diagnostic BEFORE threshold selection and model fitting
# ============================================================================
# Test one seasonal TXx per complete fitting season against calendar year.
# Do not concatenate daily data: that would introduce within-season dependence.
# This diagnostic does not test all forms of nonstationarity or attribute a trend
# to GMST. The PP model remains unchanged regardless of the test's p-value.
mk_input <- data.frame(year = years_fit$year, value = years_fit$txx)
mk_trend <- mann_kendall_trend(
  mk_input$year, mk_input$value, series = "Seasonal TXx",
  value_units = "degC", alpha = mk_alpha)
mk_trend
# Write the diagnostic immediately, even if a later PP fit does not converge.
#write.csv(mk_input, file.path(output_dir, "mann_kendall_input.csv"), row.names = FALSE)
#write.csv(mk_trend, file.path(output_dir, "mann_kendall_results.csv"), row.names = FALSE)

# ============================================================================
# 3. Apply the 30 degC baseline, then calculate the conditional 0.90 quantile
# ============================================================================
# The 30 degC line is a numerical screening level; no plot is generated here.
# First use summer days only, then retain TX > 30 degC, then take quantile 0.90.
# Thus u = Q_0.90(TX | summer, complete summer, TX > 30 degC).
# This is not the 90th percentile of all summer days. Pool the selected summers
# to obtain one fixed threshold, including 2021, as in the original CPC analysis.
above_base <- summer_fit$tx[summer_fit$tx > base_temperature]
stopifnot(length(above_base) > 0L)
u <- unname(quantile(above_base, probs = threshold_probability, type = quantile_type))
summer$exceeds <- is.finite(summer$tx) & summer$tx > u

threshold_summary <- data.frame(
  season = season_label, fitted_summers = nrow(years_fit),
  summer_days_used = nrow(summer_fit), base_C = base_temperature,
  days_above_base = length(above_base), conditional_quantile = threshold_probability,
  pot_threshold_C = u, exceedance_days = sum(summer_fit$tx > u), run_days = run_days)
threshold_summary

# ============================================================================
# 4. Decluster the complete summer daily series
# ============================================================================
# Three consecutive days at or below u separate independent clusters by default.
# Preserve below-threshold days for measuring the separation. Each cluster keeps
# one peak; ties use the earliest date. Separate summers are separate blocks.
peaks <- decluster_pot_summers(summer, years_fit, u, run_days)
years_fit$n_clusters <- tabulate(match(peaks$year, years_fit$year), nrow(years_fit))
threshold_summary$n_clusters <- nrow(peaks)
peaks
years_fit

# ============================================================================
# 5. Define the observed 2021 summer TXx and the two climate conditions
# ============================================================================
# The event is the maximum daily TX within the selected 2021 summer, not the
# annual maximum file. This keeps the event and fitted seasonal exposure aligned.
event_info <- summer_info[summer_info$year == event_year, ]
stopifnot(nrow(event_info) == 1L, event_info$complete)
event_summer <- summer[summer$year == event_year, ]
x_event <- max(event_summer$tx)
event_date <- event_summer$date[which.max(event_summer$tx)]
stopifnot(is.finite(x_event), x_event > u)
g_factual <- gmst$gmst[gmst$year == factual_year]
g_counterfactual <- gmst$gmst[gmst$year == counterfactual_year]
stopifnot(length(g_factual) == 1L, is.finite(g_factual),
          length(g_counterfactual) == 1L, is.finite(g_counterfactual))
reference <- paste0(counterfactual_year, "-like")
x_event
event_date

# ============================================================================
# 6. Fit the nonstationary POT model and calculate attribution point estimates
# ============================================================================
# PP mu changes linearly with GMST; PP sigma and xi remain constant. This also
# allows the cluster rate and implied GPD scale to change with GMST.
# p_factual and p_counterfactual are probabilities of at least one event as hot
# as x_event in ONE SUMMER under each fixed GMST condition. Return periods are
# reciprocals of these probabilities: years, with one summer opportunity per year.
# They describe repeated fixed-climate summers, not a forecast under future warming.
fit_pot <- fit_pp(years_fit, peaks, u)
fit_pot
point_results <- pp_attribution(fit_pot, x_event, g_factual, g_counterfactual, reference)
point_results$index <- "Summer TXx"
point_results$season <- season_label
point_results$probability_basis <- paste("One", season_label, "season per year")
point_results$event_date <- event_date
point_results$n_years <- nrow(years_fit)
point_results$n_clusters <- nrow(peaks)
# Save point estimates directly while bootstrap is disabled.
results <- point_results
point_results
# Keep infinite RR explicitly when the counterfactual event probability is zero.
if (point_results$lambda_counterfactual == 0)
  warning("Counterfactual threshold exceeds the model endpoint; inspect this boundary result.")

# ============================================================================
# 7. Bootstrap confidence intervals for the numerical attribution results
# ============================================================================
# Resample complete summers as blocks, including their GMST and all peaks.
# Bootstrap is disabled below; the saved results contain point estimates only.
# Reference GMST, x_event, u and declustering choices remain fixed in bootstrap.
# bootstrap_output <- bootstrap_pot(
#   years_fit, peaks, u, x_event, g_factual, g_counterfactual,
#   point_results, bootstrap_B, random_seed, ci_level, reference)
# bootstrap_samples <- bootstrap_output$samples
# results <- bootstrap_output$results
# results

# ============================================================================
# 8. Save numerical results and intermediate objects
# ============================================================================
# No plotting functions, graphics devices, or figure files are used.
settings <- list(
  mk_alpha = mk_alpha, trend_test = "Classical Mann-Kendall diagnostic, not a model-selection gate",
  daily_file = daily_file, gmst_file = gmst_file,
  start_year = start_year, end_year = end_year, season_start = season_start, season_end = season_end,
  season = season_label, event_year = event_year, event_date = event_date,
  factual_year = factual_year, counterfactual_year = counterfactual_year,
  base_temperature = base_temperature, threshold_probability = threshold_probability,
  quantile_type = quantile_type, threshold_C = u, run_days = run_days,
  # bootstrap_B = bootstrap_B, random_seed = random_seed, ci_level = ci_level,
  exposure = paste("One complete", season_label, "season per year"),
  model = "PP: mu linear in GMST; constant sigma and xi",
  # bootstrap = "Summer pairs bootstrap; fixed event, threshold, and clusters",
  uncertainty = "Not calculated; bootstrap disabled")
write.csv(results, file.path(output_dir, "cpc_pot_attribution.csv"), row.names = FALSE)
write.csv(threshold_summary, file.path(output_dir, "threshold_summary.csv"), row.names = FALSE)
write.csv(peaks, file.path(output_dir, "cluster_peaks.csv"), row.names = FALSE)
write.csv(summer_info, file.path(output_dir, "summer_completeness_and_gmst.csv"), row.names = FALSE)
write.csv(years_fit, file.path(output_dir, "fitted_summers_and_cluster_counts.csv"), row.names = FALSE)
# write.csv(bootstrap_samples, file.path(output_dir, "bootstrap_samples.csv"), row.names = FALSE)
saveRDS(list(settings = settings, mk_trend = mk_trend, mk_input = mk_input,
             fit = fit_pot, results = results,
             threshold_summary = threshold_summary, summer = summer,
             summer_info = summer_info, years_fit = years_fit,
             # bootstrap = bootstrap_samples,
             peaks = peaks),
        file.path(output_dir, "cpc_pot_fit.rds"))





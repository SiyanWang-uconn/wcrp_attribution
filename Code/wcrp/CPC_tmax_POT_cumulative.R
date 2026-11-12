# CPC summer cumulative heatwave magnitude attribution: exploratory version.
# Run one numbered section at a time in RStudio. No bootstrap or plots are used.
# Every model fit in this file is to event magnitude S, measured in degC day.

# ============================================================================
# 1. Paths, helper functions, and configuration
# ============================================================================
project_root <- "C:/Users/siyan/WCRP/wcrp_common_attribution_assessment"
source(file.path(project_root, "Code", "helpers", "pot.R"))
source(file.path(project_root, "Code", "helpers", "trend.R"))
source(file.path(project_root, "Code", "helpers", "pot_cumulative.R"))

daily_file <- file.path(project_root, "outputs/pnw_heatwave",
                        "CPC_reference_daily_tmax.nc")
gmst_file <- file.path(project_root, "data", "processed", "gmst_smoothed.csv")
output_dir <- file.path(project_root, "Results", "CPC_POT_cumulative")
start_year <- 1979L
end_year <- 2025L
season_start <- "05-15"
season_end <- "09-15"
season_label <- "May15-Sep15"
event_anchor_date <- as.Date("2021-06-29")
factual_year <- 2021L
counterfactual_year <- 1951L
base_temperature <- 30
threshold_probability <- 0.90
quantile_type <- 7L
run_days <- 3L
mk_alpha <- 0.05

# There are TWO thresholds with different units and purposes:
# u_temperature: degC, identifies hot days and determines event boundaries.
# score_threshold: degC day, selects event magnitudes for the PP tail fit.
# Start with 0 to use all positive event scores. This is an exploratory parametric
# model assumption, not a claim that all positive scores satisfy tail asymptotics.
# A positive magnitude threshold can be investigated later with tail diagnostics.
score_threshold <- 0

stopifnot(file.exists(daily_file), file.exists(gmst_file),
          requireNamespace("ncdf4", quietly = TRUE), start_year <= end_year,
          is.finite(score_threshold), score_threshold >= 0,
          threshold_probability > 0, threshold_probability < 1,
          run_days >= 1L, run_days == as.integer(run_days),
          length(event_anchor_date) == 1L, !is.na(event_anchor_date),
          length(season_start) == 1L, length(season_end) == 1L)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# 2. Select May 15-September 15 temperatures and complete seasons
# ============================================================================
daily <- read_pot_daily(daily_file)
summer <- select_pot_season(daily, start_year, end_year, season_start, season_end)
rownames(summer) <- NULL
stopifnot(nrow(summer) > 0L)
gmst <- read.csv(gmst_file)
summer_info <- summarize_pot_summers(summer, gmst, season_start = season_start, season_end = season_end)
years_fit <- summer_info[summer_info$used_in_fit, ]
stopifnot(nrow(years_fit) > 4L, length(unique(years_fit$gmst)) > 1L)
summer$used_in_fit <- summer$year %in% years_fit$year
summer_fit <- summer[summer$used_in_fit, ]
summer_info

# ============================================================================
# 3. Obtain the fixed temperature threshold from summer TX above 30 degC
# ============================================================================
# First select summer, then TX > 30 degC, then calculate conditional quantile 0.90.
# Pool complete fitting summers, including 2021, to define one fixed threshold.
# This screening level is applied numerically; this script does not draw a line.
above_base <- summer_fit$tx[summer_fit$tx > base_temperature]
stopifnot(length(above_base) > 0L)
u_temperature <- unname(quantile(above_base, threshold_probability, type = quantile_type))
summer$exceeds <- is.finite(summer$tx) & summer$tx > u_temperature
u_temperature

# ============================================================================
# 4. Decluster and calculate cumulative magnitude for every event
# ============================================================================
# S = sum(max(TX - u_temperature, 0)) * 1 day over each independent cluster.
# Duration and average excess jointly determine S: S = exceedance_days*mean_excess.
# The default runs rule allows gaps shorter than three below-threshold days.
# Such gap days add zero to S. span_days can therefore exceed exceedance_days.
# Only the May 15-September 15 part of a boundary-crossing event contributes to its score.
events <- cumulative_pot_events(summer, years_fit, u_temperature, run_days)
events

# Select the cluster containing the known PNW event date. Do not automatically
# replace it with the year's largest cumulative score, which could be another event.
event <- events[events$start_date <= event_anchor_date &
                  events$end_date >= event_anchor_date, ]
if (nrow(event) != 1L)
  stop("The anchor date does not identify one retained event; inspect events and settings.")
event_score <- event$cumulative_excess
stopifnot(is.finite(event_score), event_score > score_threshold)
event

# ============================================================================
# 5. Prepare magnitude observations and summer exposure for the new PP fit
# ============================================================================
# This adapter puts S into the existing PP fitter's generic observation slots.
# pp_data$peaks$peak_tx contains cumulative magnitude, NOT temperature, here.
# pp_data$years$txx contains the summer maximum score for initialization only.
# events preserves the original temperatures and all event duration diagnostics.
pp_data <- prepare_cumulative_pp(years_fit, events, score_threshold)
years_magnitude <- pp_data$years
peaks_magnitude <- pp_data$peaks
magnitude_summary <- data.frame(
  season = season_label, n_summers = nrow(years_magnitude),
  temperature_threshold_C = u_temperature, base_temperature_C = base_temperature,
  conditional_quantile = threshold_probability, run_days = run_days,
  n_events = nrow(events), score_threshold_C_days = score_threshold,
  n_events_used = nrow(peaks_magnitude), event_score_C_days = event_score)
magnitude_summary
years_magnitude

# ============================================================================
# 5a. Mann-Kendall diagnostic BEFORE fitting the cumulative-magnitude model
# ============================================================================
# Test one seasonal maximum event score per complete season against calendar year.
# Zero-event seasons have score zero and remain in this series; do not test only
# the selected extreme events, which would omit changes in event occurrence.
# This diagnostic does not test all forms of nonstationarity or attribute a trend
# to GMST. The PP model remains unchanged regardless of the test's p-value.
mk_input <- data.frame(year = years_magnitude$year,
                       value = years_magnitude$txx)
mk_trend <- mann_kendall_trend(
  mk_input$year, mk_input$value, series = "Seasonal maximum cluster cumulative excess",
  value_units = "degC day", alpha = mk_alpha)
mk_trend
# Write the diagnostic immediately, even if a later PP fit does not converge.
#write.csv(mk_input, file.path(output_dir, "mann_kendall_input.csv"), row.names = FALSE)
#write.csv(mk_trend, file.path(output_dir, "mann_kendall_results.csv"), row.names = FALSE)

# ============================================================================
# 6. Fit the nonstationary cumulative-magnitude PP model
# ============================================================================
# Refit all four PP parameters to S, rather than reusing the temperature model.
# mu_S(g) = mu_S0 + alpha_S*(g-g0); sigma_S and xi_S are constant.
# mu_S and sigma_S have units degC day; alpha_S has units degC day per degC GMST.
# The implied conditional distribution above score_threshold is a GPD in S units.
# As in the original fitter, a shape <= -0.5 stops with a nonregular-fit message;
# this exploratory version does not silently change that fitting rule.
fit_magnitude <- fit_pp(years_magnitude, peaks_magnitude, score_threshold)
fit_magnitude
fit_parameters <- data.frame(
  mu_C_days = fit_magnitude$mu, alpha_C_days_per_C_GMST = fit_magnitude$alpha,
  sigma_C_days = fit_magnitude$sigma, shape = fit_magnitude$xi,
  gmst_center = fit_magnitude$g0, nll = fit_magnitude$nll,
  convergence = fit_magnitude$convergence)
fit_parameters

# ============================================================================
# 7. Attribute the 2021 event's cumulative magnitude: point estimates only
# ============================================================================
g_factual <- gmst$gmst[gmst$year == factual_year]
g_counterfactual <- gmst$gmst[gmst$year == counterfactual_year]
stopifnot(length(g_factual) == 1L, is.finite(g_factual),
          length(g_counterfactual) == 1L, is.finite(g_counterfactual))
results <- cumulative_pp_attribution(
  fit_magnitude, event_score, g_factual, g_counterfactual,
  reference = paste0(counterfactual_year, "-like"))
results$index <- "Cluster cumulative positive TX excess"
results$season <- season_label
results$event_anchor_date <- event_anchor_date
results$event_start <- event$start_date
results$event_end <- event$end_date
results$event_peak_tx_C <- event$peak_tx
results$event_exceedance_days <- event$exceedance_days
results$event_mean_excess_C <- event$mean_excess
results$temperature_threshold_C <- u_temperature
results$n_years <- nrow(years_magnitude)
results$n_clusters <- nrow(peaks_magnitude)
results$probability_basis <- paste("At least one event reaching the 2021 cumulative score per",
                                   season_label, "season")
results
# RR/FAR refer to reaching event_score, not reaching 39 degC or summing an entire
# summer's heat. Intensity change is a difference in degC day at equal probability.
# A zero counterfactual probability produces infinite RR; preserve that result.

# ============================================================================
# 8. Save numerical results and intermediate objects
# ============================================================================
settings <- list(
  mk_alpha = mk_alpha, trend_test = "Classical Mann-Kendall diagnostic, not a model-selection gate",
  daily_file = daily_file, gmst_file = gmst_file, start_year = start_year,
  end_year = end_year, season_start = season_start, season_end = season_end, season = season_label,
  event_anchor_date = event_anchor_date, factual_year = factual_year,
  counterfactual_year = counterfactual_year, base_temperature = base_temperature,
  threshold_probability = threshold_probability, quantile_type = quantile_type,
  temperature_threshold_C = u_temperature, run_days = run_days,
  score_threshold_C_days = score_threshold,
  score_definition = "sum(max(TX-u_temperature,0))*1 day per declustered cluster",
  model = "PP fitted to cumulative score; mu linear in GMST; constant sigma and xi",
  exposure = paste("One complete", season_label, "season per year"), uncertainty = "Not calculated")
write.csv(results, file.path(output_dir, "cpc_cumulative_attribution.csv"), row.names = FALSE)
write.csv(fit_parameters, file.path(output_dir, "fit_parameters.csv"), row.names = FALSE)
write.csv(magnitude_summary, file.path(output_dir, "magnitude_summary.csv"), row.names = FALSE)
write.csv(events, file.path(output_dir, "events_cumulative.csv"), row.names = FALSE)
write.csv(summer_info, file.path(output_dir, "summer_completeness_and_gmst.csv"), row.names = FALSE)
write.csv(years_magnitude, file.path(output_dir, "fitted_summers.csv"), row.names = FALSE)
saveRDS(list(settings = settings, mk_trend = mk_trend, mk_input = mk_input,
             fit = fit_magnitude, results = results,
             fit_parameters = fit_parameters, magnitude_summary = magnitude_summary,
             summer = summer, summer_info = summer_info, events = events,
             event = event, years_magnitude = years_magnitude,
             peaks_magnitude = peaks_magnitude),
        file.path(output_dir, "cpc_cumulative_fit.rds"))



# Smith Section 0.6.2 PNW heatwave: MLE method reproduction with existing ERA5 data.
# Run numbered sections in RStudio. This file is not executed by its creation.
# New helper and output names keep all earlier project code and results separate.
# This stage implements MLE, not the paper's later Bayesian Table 0.2 analysis.

# ============================================================================
# 1. Paths, packages, and explicit reconstruction choices
# ============================================================================
project_root <- "C:/Users/siyan/WCRP/wcrp_common_attribution_assessment"
source(file.path(project_root, "Code", "smith_reproduce", "helpers", "smith_pnw_pot.R"))
daily_file <- file.path(project_root, "outputs/pnw_heatwave", "ERA5_reference_daily_tmax.nc")
gmst_file <- file.path(project_root, "data/processed/gmst_smoothed.csv")
output_dir <- file.path(project_root, "Results/ERA5_Smith_PNW_POT_MLE")
fit_start_year <- 1950L
fit_end_year <- 2020L
prediction_year <- 2021L
counterfactual_cooling_C <- 1.2
quantile_tau <- 0.90
threshold_harmonic_order <- 2L


library(texmex)
# The paper's calendar description and numbered-day limits disagree.
# Default: reproduce its explicit PP limits, days 135-255 on a 365-day calendar
# (May 15-September 12). To use May 15-September 15 instead, choose "calendar".
window_mode <- "paper_doy"
season_first_day <- 135L
season_last_day <- 255L
season_first_date <- "05-15"
season_last_date <- "09-15"

# pp.fit parameterizes a full year using npy. Retain its 365-day convention and
# integrate predicted daily rates with 1/npy, including for the 2021 predictions.
# Do not replace the denominator by the number of peaks or multiply by EI again.
npy <- 365
optimizer_maxit <- 50000L
initial_shapes <- c(-0.2, -0.05, 0.1)
count_temperatures_C <- seq(30, 39, by = 0.5)
attribution_temperatures_C <- 34:39

smith_check_packages()
stopifnot(file.exists(daily_file), file.exists(gmst_file),
          fit_start_year < fit_end_year, fit_end_year < prediction_year,
          counterfactual_cooling_C > 0, npy > 0,
          length(initial_shapes) > 0, all(is.finite(initial_shapes)))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# 2. Read existing ERA5/GMST and audit the 1950-2020 fitting seasons
# ============================================================================
# Existing ERA5 and smoothed GMST are used as requested. They are not the exact
# KNMI regional series / historical GMST release used by Smith: numerical equality
# with the paper's endpoint estimates is not assumed or enforced.
daily <- smith_read_daily(daily_file)
gmst <- read.csv(gmst_file)
prepared <- smith_prepare(daily, gmst, fit_start_year, fit_end_year, window_mode,
                           season_first_day, season_last_day,
                           season_first_date, season_last_date)
fit_data <- prepared$data
coverage <- prepared$coverage
gmst_center <- prepared$gmst_center
coverage
#write.csv(coverage, file.path(output_dir, "fitting_season_coverage.csv"), row.names = FALSE)

# ============================================================================
# 3. Estimate the time-varying 0.90 threshold by quantile regression
# ============================================================================
# No 30 degC prefilter is used. The 30 degC levels are for count diagnostics.
# Threshold model: intercept + linear calendar-year trend + annual harmonics;
# the default also includes semiannual harmonics. These details are configurable
# because the paper does not publish the exact quantile-regression formula.
threshold_fit <- smith_threshold(fit_data, quantile_tau, threshold_harmonic_order)
fit_data$threshold_C <- threshold_fit$threshold
fit_data$exceeds <- fit_data$tx > fit_data$threshold_C
threshold_coefficients <- data.frame(
  term = names(stats::coef(threshold_fit$model)),
  estimate = unname(stats::coef(threshold_fit$model)))
threshold_coefficients

# ============================================================================
# 4. Ferro-Segers intervals declustering through texmex
# ============================================================================
# Apply EI to TX-u(t) at zero to accommodate the varying temperature threshold.
# Automatic intervals clustering replaces the earlier fixed three-day runs rule.
# Retain the highest actual TX in each cluster and discard nonpeak exceedances
# from the event term only. All observation days remain in the exposure integral.
declustered <- smith_decluster(fit_data, fit_data$threshold_C)
peaks <- declustered$peaks
extremal_index_summary <- declustered$summary
extremal_index_summary
peaks
write.csv(peaks, file.path(output_dir, "declustered_peaks.csv"), row.names = FALSE)
write.csv(extremal_index_summary, file.path(output_dir, "extremal_index.csv"), row.names = FALSE)

# ============================================================================
# 5. Fit the paper's seasonal PP model using ismev::pp.fit
# ============================================================================
# mu = b0 + bG*(GMST-g0) + b1*sin1 + b2*cos1 + b3*sin2 + b4*cos2
# sigma = s0 + s1*sin1 + s2*cos1, with identity link and positive fitted values.
# xi is constant. All ten parameters are fitted jointly, not copied from old fits.
# The helper independently evaluates the same PP likelihood to check exposure.
model <- smith_fit_pp(fit_data, fit_data$threshold_C, peaks,
                       npy, optimizer_maxit, initial_shapes)
parameter_table <- data.frame(parameter = names(model$coef), estimate = unname(model$coef))
optimizer_audit <- model$audit
parameter_table
optimizer_audit

# ============================================================================
# 6. Compare observed and expected historical declustered exceedance counts
# ============================================================================
count_check <- smith_count_check(model, fit_data, fit_data$threshold_C,
                                 peaks, count_temperatures_C)
count_check
# Thresholds at which any fitted u(t)>x require below-threshold extrapolation
# and are flagged explicitly. Nested cumulative counts are correlated: the
# Pearson components are descriptive and are not assigned chi-square p-values.

# ============================================================================
# 7. Predict 2021 with observed GMST and with GMST reduced by 1.2 degC
# ============================================================================
projection <- smith_season(daily, window_mode, season_first_day, season_last_day,
                            season_first_date, season_last_date)
projection <- projection[projection$year == prediction_year, , drop = FALSE]
expected_days <- unique(coverage$expected_days)
stopifnot(length(expected_days) == 1L, nrow(projection) == expected_days,
          all(is.finite(projection$tx)))
g_factual <- gmst$gmst[gmst$year == prediction_year]
stopifnot(length(g_factual) == 1L, is.finite(g_factual))
factual <- projection
factual$gmst <- g_factual
factual <- smith_features(factual, gmst_center, fit_start_year)
counterfactual <- projection
counterfactual$gmst <- g_factual - counterfactual_cooling_C
counterfactual <- smith_features(counterfactual, gmst_center, fit_start_year)

# The 2021 data define an evaluation target only; they never enter fitting.
event_temperature_C <- max(projection$tx)
event_date <- projection$date[which.max(projection$tx)]
temperatures <- sort(unique(c(attribution_temperatures_C, event_temperature_C)))
results <- smith_attribution(model, factual, counterfactual, temperatures)
results$is_observed_2021_event <- abs(results$temperature_C - event_temperature_C) < 1e-8
results$gmst_factual <- g_factual
results$gmst_counterfactual <- g_factual - counterfactual_cooling_C
results$probability_basis <- "At least one cluster event above x in the selected season"
results

# Report daily and seasonal endpoints, not just an arbitrarily chosen summer day.
daily_parameters <- rbind(smith_daily_parameters(model, factual, "factual"),
                           smith_daily_parameters(model, counterfactual, "counterfactual"))
endpoint_summary <- aggregate(upper_endpoint_C ~ scenario, daily_parameters, max)
endpoint_summary
# Both probabilities may be zero above the fitted endpoints. Keep RR=NA in that
# case: do not replace 0/0 by an infinite RR or force the endpoint above 2021.
# Table 0.2 is Bayesian: these MLE outputs are not its posterior means/intervals.

# ============================================================================
# 8. Save only this stage's numerical results and inspectable intermediate data
# ============================================================================
settings <- list(
  daily_file = daily_file, gmst_file = gmst_file,
  source = "Smith Detection and Attribution, revised April 18 2025, Section 0.6.2",
  data_choice = "Existing regional ERA5 and existing smoothed GMST",
  fit_start_year = fit_start_year, fit_end_year = fit_end_year,
  prediction_year = prediction_year, counterfactual_cooling_C = counterfactual_cooling_C,
  window_mode = window_mode, season_first_day = season_first_day,
  season_last_day = season_last_day, season_first_date = season_first_date,
  season_last_date = season_last_date, quantile_tau = quantile_tau,
  threshold_harmonic_order = threshold_harmonic_order, npy = npy,
  gmst_center = gmst_center, optimizer_maxit = optimizer_maxit,
  initial_shapes = initial_shapes, event_temperature_C = event_temperature_C,
  event_date = event_date, inference = "MLE only; no bootstrap or MCMC")
write.csv(threshold_coefficients, file.path(output_dir, "threshold_coefficients.csv"), row.names = FALSE)
write.csv(fit_data, file.path(output_dir, "daily_thresholds_and_covariates.csv"), row.names = FALSE)
write.csv(parameter_table, file.path(output_dir, "pp_parameters.csv"), row.names = FALSE)
write.csv(optimizer_audit, file.path(output_dir, "optimizer_audit.csv"), row.names = FALSE)
write.csv(count_check, file.path(output_dir, "historical_count_check.csv"), row.names = FALSE)
write.csv(results, file.path(output_dir, "attribution_mle.csv"), row.names = FALSE)
write.csv(daily_parameters, file.path(output_dir, "prediction_daily_parameters.csv"), row.names = FALSE)
write.csv(endpoint_summary, file.path(output_dir, "prediction_endpoints.csv"), row.names = FALSE)
saveRDS(list(settings = settings, model = model, coverage = coverage,
             threshold_fit = threshold_fit, fit_data = fit_data, peaks = peaks,
             declustered = declustered, count_check = count_check,
             factual = factual, counterfactual = counterfactual, results = results,
             endpoint_summary = endpoint_summary),
        file.path(output_dir, "smith_pnw_pot_mle.rds"))


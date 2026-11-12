# Cumulative heatwave magnitude helpers. Definitions only; no analysis runs here.
# Source pot.R before this file to reuse the daily reader and PP likelihood.
# Event magnitude S is sum(max(TX - u, 0)) over one declustered event, in degC day.
# A one-day observation interval is assumed. S is a temperature-based index,
# not physical energy and not a daily maximum temperature.

cumulative_pot_events <- function(summer, years_fit, temperature_threshold, run_days) {
  # Reuse the original runs definition and keep one row for each independent event.
  events <- decluster_pot_summers(summer, years_fit, temperature_threshold, run_days)
  events$span_days <- as.integer(events$end_date - events$start_date) + 1L
  events$cumulative_excess <- NA_real_
  events$mean_excess <- NA_real_
  events$boundary_truncated <- FALSE
  for (i in seq_len(nrow(events))) {
    z <- summer[summer$year == events$year[i] &
                  summer$date >= events$start_date[i] &
                  summer$date <= events$end_date[i], ]
    excess <- pmax(z$tx - temperature_threshold, 0)
    stopifnot(all(is.finite(excess)),
              sum(excess > 0) == events$exceedance_days[i])
    # Short below-threshold gaps may belong to the cluster but contribute zero.
    # Use the sum of daily excesses, not peak excess multiplied by duration.
    events$cumulative_excess[i] <- sum(excess)
    events$mean_excess[i] <- sum(excess) / events$exceedance_days[i]
    season_dates <- summer$date[summer$year == events$year[i]]
    events$boundary_truncated[i] <- events$start_date[i] == min(season_dates) ||
      events$end_date[i] == max(season_dates)
  }
  stopifnot(all(events$cumulative_excess > 0))
  events
}

prepare_cumulative_pp <- function(years_fit, events, score_threshold) {
  # Adapt the existing PP fitter to magnitude units without changing its likelihood.
  # peak_tx becomes S, and txx becomes the largest S in each summer, used ONLY
  # for optimizer initialization. No temperature peak is used in this fit.
  # A summer with no events has maximum S = 0 and still contributes exposure.
  years_data <- years_fit
  years_data$summer_txx_C <- years_data$txx
  years_data$txx <- 0
  years_data$n_clusters <- 0L
  years_data$n_score_exceedances <- 0L
  for (j in seq_len(nrow(years_data))) {
    scores <- events$cumulative_excess[events$year == years_data$year[j]]
    if (length(scores)) years_data$txx[j] <- max(scores)
    years_data$n_clusters[j] <- length(scores)
    years_data$n_score_exceedances[j] <- sum(scores > score_threshold)
  }
  selected <- events[events$cumulative_excess > score_threshold, ]
  if (nrow(selected) < 10L)
    stop("Fewer than 10 events exceed the magnitude threshold; inspect settings.")
  selected$peak_temperature_C <- selected$peak_tx
  selected$peak_tx <- selected$cumulative_excess
  selected$excess <- selected$cumulative_excess - score_threshold
  list(years = years_data, peaks = selected)
}

cumulative_pp_attribution <- function(fit, event_score, gf, gc, reference) {
  # Reuse the same PP probability calculation after refitting in magnitude units.
  # Event score, magnitude threshold and PP scale are all in degC day.
  result <- pp_attribution(fit, event_score, gf, gc, reference)
  names(result)[names(result) == "threshold_C"] <- "event_cumulative_excess_C_days"
  names(result)[names(result) == "pot_threshold_C"] <- "score_threshold_C_days"
  names(result)[names(result) == "gpd_scale_factual"] <- "gpd_scale_factual_C_days"
  names(result)[names(result) == "gpd_scale_counterfactual"] <- "gpd_scale_counterfactual_C_days"
  names(result)[names(result) == "intensity_change"] <- "cumulative_intensity_change_C_days"
  # With a location-only PP change, equal-probability levels shift by alpha*dG.
  # This shift describes event magnitudes, not a temperature change in degC.
  equivalent_cf <- event_score - fit$alpha * (gf - gc)
  result$counterfactual_equivalent_score_C_days <- equivalent_cf
  result$intensity_status <- "Within fitted magnitude tail in both climates"
  if (equivalent_cf <= fit$threshold) {
    # A matching return level below the magnitude threshold is outside the tail
    # model. Do not report its extrapolation as an estimated intensity change.
    result$counterfactual_equivalent_score_C_days <- NA_real_
    result$cumulative_intensity_change_C_days <- NA_real_
    result$intensity_status <- "Matching counterfactual level is outside fitted magnitude tail"
  }
  result
}

# POT attribution helpers. This file defines functions only and runs no analysis.
# One exposure unit is one complete summer, including summers with zero clusters.
# The PP model combines event counts and cluster-peak magnitudes in one likelihood.
# mu(g) = mu + alpha * (g - g0); PP sigma and xi are constant.
# lambda_u(g) = [1 + xi * (u - mu(g)) / sigma]^(-1 / xi).
# The implied GPD scale is sigma_u(g) = sigma + xi * (u - mu(g)).
# A seasonal exceedance probability is 1 - exp(-lambda_x(g)).
# Its reciprocal is a return period in years, with one summer opportunity per year.

read_pot_daily <- function(file) {
  nc <- ncdf4::nc_open(file)
  on.exit(ncdf4::nc_close(nc))
  stopifnot(identical(vapply(nc$var$tmax$dim, function(x) x$name, ""), "time"))
  units <- ncdf4::ncatt_get(nc, "tmax", "units")$value
  stopifnot(units %in% c("degC", "C", "Celsius", "degree_Celsius"))
  tx <- as.numeric(ncdf4::ncvar_get(nc, "tmax"))
  time <- as.numeric(ncdf4::ncvar_get(nc, "time"))
  time_units <- nc$dim$time$units
  stopifnot(grepl("^(hours|days) since ", time_units))
  origin <- as.Date(substr(sub("^(hours|days) since ", "", time_units), 1, 10))
  offset <- time / if (startsWith(time_units, "hours")) 24 else 1
  stopifnot(all(abs(offset - round(offset)) < 1e-7))
  date <- origin + round(offset)
  stopifnot(length(date) == length(tx), !anyDuplicated(date), all(diff(date) == 1))
  data.frame(date = date, tx = tx, year = as.integer(format(date, "%Y")),
             month = as.integer(format(date, "%m")))
}

select_pot_season <- function(daily, start_year, end_year, season_start, season_end) {
  # Boundaries are month-day strings; include both endpoints of a within-year season.
  # May 15-September 15 contains 124 days in both leap and non-leap years.
  stopifnot(length(season_start) == 1L, length(season_end) == 1L,
            !is.na(season_start), !is.na(season_end))
  start <- as.Date(paste0("2000-", season_start), format = "%Y-%m-%d")
  end <- as.Date(paste0("2000-", season_end), format = "%Y-%m-%d")
  stopifnot(!is.na(start), !is.na(end), start <= end,
            format(start, "%m-%d") == season_start,
            format(end, "%m-%d") == season_end)
  month_day <- format(daily$date, "%m-%d")
  daily[daily$year >= start_year & daily$year <= end_year &
          month_day >= season_start & month_day <= season_end, ]
}

summarize_pot_summers <- function(summer, gmst, months = NULL,
                                  season_start = NULL, season_end = NULL) {
  # Completeness concerns the selected summer, not the full calendar year.
  # Missing winter observations therefore do not discard a complete summer.
  stopifnot(all(c("year", "gmst") %in% names(gmst)), !anyDuplicated(gmst$year))
  info <- do.call(rbind, lapply(sort(unique(summer$year)), function(y) {
    z <- summer[summer$year == y, ]
    calendar <- seq(as.Date(sprintf("%d-01-01", y)),
                    as.Date(sprintf("%d-12-31", y)), by = "day")
    # Date limits must match the season selection, including partial months.
    # Retain the optional months argument for existing month-based callers.
    if (!is.null(season_start) && !is.null(season_end)) {
      month_day <- format(calendar, "%m-%d")
      expected <- sum(month_day >= season_start & month_day <= season_end)
    } else {
      stopifnot(!is.null(months), all(months %in% 1:12))
      expected <- sum(as.integer(format(calendar, "%m")) %in% months)
    }
    data.frame(year = y, expected_days = expected, valid_days = sum(is.finite(z$tx)),
               complete = nrow(z) == expected && all(is.finite(z$tx)),
               txx = if (any(is.finite(z$tx))) max(z$tx, na.rm = TRUE) else NA_real_)
  }))
  info <- merge(info, gmst[, c("year", "gmst")], by = "year", all.x = TRUE)
  info$used_in_fit <- info$complete & is.finite(info$gmst)
  info
}

decluster_pot_summers <- function(summer, years_fit, threshold, run_days) {
  # Keep all days within each summer when measuring gaps. Never remove the
  # below-threshold days before declustering. Missing summers are excluded.
  # A new cluster starts after at least r consecutive days at or below u.
  # Clusters at the selected season boundaries use only their observed segment:
  # both counts and peak magnitudes refer to the same fixed summer window.
  clusters <- lapply(years_fit$year, function(y) {
    z <- summer[summer$year == y, ]
    stopifnot(all(is.finite(z$tx)), all(diff(z$date) == 1))
    ii <- which(z$tx > threshold)
    if (!length(ii)) return(NULL)
    group <- cumsum(c(TRUE, diff(ii) > run_days))
    do.call(rbind, lapply(split(ii, group), function(jj) {
      peak <- jj[which.max(z$tx[jj])]
      data.frame(year = y, start_date = min(z$date[jj]), end_date = max(z$date[jj]),
                 peak_date = z$date[peak], peak_tx = z$tx[peak],
                 excess = z$tx[peak] - threshold, exceedance_days = length(jj))
    }))
  })
  peaks <- do.call(rbind, clusters)
  if (is.null(peaks) || nrow(peaks) < 10L)
    stop("Fewer than 10 cluster peaks: inspect threshold and declustering settings.")
  rownames(peaks) <- NULL
  peaks$gmst <- years_fit$gmst[match(peaks$year, years_fit$year)]
  if (nrow(peaks) < 30L) warning("Fewer than 30 cluster peaks: inference is exploratory.")
  peaks
}

# The likelihood integrates the PP rate over one full summer for every row of
# years_data. Zero-event summers must remain in years_data. Seasonal TXx is used
# only to initialize the optimizer; fitting uses all retained cluster peaks.

log_pp_rate <- function(x, mu, sigma, xi) {
  z <- (x - mu) / sigma
  if (abs(xi) < 1e-7) return(-z)
  h <- 1 + xi * z
  result <- rep(if (xi < 0) -Inf else Inf, length(h))
  ok <- h > 0
  result[ok] <- -log(h[ok]) / xi
  result
}

fit_pp <- function(years_data, peak_data, threshold) {
  g0 <- mean(years_data$gmst)
  gy <- years_data$gmst - g0
  gp <- peak_data$gmst - g0
  regression <- lm(txx ~ I(gmst - g0), data = years_data)
  initial_scale <- max(sd(residuals(regression)) * sqrt(6) / pi, 0.1)
  initial_mu <- unname(coef(regression)[1]) - 0.5772157 * initial_scale
  initial_alpha <- unname(coef(regression)[2])
  nll <- function(par) {
    if (any(!is.finite(par))) return(1e100)
    sigma <- exp(par[3])
    xi <- par[4]
    if (!is.finite(sigma) || sigma < 1e-6) return(1e100)
    mu_y <- par[1] + par[2] * gy
    mu_p <- par[1] + par[2] * gp
    log_rate <- log_pp_rate(threshold, mu_y, sigma, xi)
    if (any(is.na(log_rate)) || any(log_rate > 700)) return(1e100)
    z <- (peak_data$peak_tx - mu_p) / sigma
    if (abs(xi) < 1e-7) {
      log_intensity <- -log(sigma) - z
    } else {
      h <- 1 + xi * z
      if (any(!is.finite(h)) || any(h <= 0)) return(1e100)
      log_intensity <- -log(sigma) - (1 + 1 / xi) * log(h)
    }
    value <- sum(exp(log_rate)) - sum(log_intensity)
    if (is.finite(value)) value else 1e100
  }
  candidates <- lapply(c(-0.25, -0.1, 0, 0.1, 0.25), function(shape) {
    tryCatch(optim(c(initial_mu, initial_alpha, log(initial_scale), shape), nll,
                   method = "Nelder-Mead", control = list(maxit = 5000, reltol = 1e-10)),
             error = function(e) NULL)
  })
  good <- vapply(candidates, function(z)
    !is.null(z) && z$convergence == 0L && is.finite(z$value) && z$value < 1e90, TRUE)
  if (!any(good)) stop("PP optimizer failed.")
  best <- candidates[good][[which.min(vapply(candidates[good], function(z) z$value, 0.0))]]
  if (best$par[4] <= -0.5)
    stop("PP shape <= -0.5: nonregular fit; inspect data and model.")
  list(mu = unname(best$par[1]), alpha = unname(best$par[2]),
       sigma = exp(unname(best$par[3])), xi = unname(best$par[4]),
       g0 = g0, threshold = threshold, nll = best$value,
       convergence = best$convergence)
}

log_season_probability <- function(log_rate) {
  # Stable log(1-exp(-exp(log_rate))), including finite-endpoint cases.
  result <- log_rate
  middle <- is.finite(log_rate) & log_rate >= -35 & log_rate <= 35
  result[middle] <- log(-expm1(-exp(log_rate[middle])))
  result[log_rate > 35] <- 0
  result
}

pp_attribution <- function(fit, event, gf, gc, reference) {
  mu_f <- fit$mu + fit$alpha * (gf - fit$g0)
  mu_c <- fit$mu + fit$alpha * (gc - fit$g0)
  rate_f <- log_pp_rate(event, mu_f, fit$sigma, fit$xi)
  rate_c <- log_pp_rate(event, mu_c, fit$sigma, fit$xi)
  log_pf <- log_season_probability(rate_f)
  log_pc <- log_season_probability(rate_c)
  if (!is.finite(log_pf)) stop("Observed event has zero factual probability.")
  log_rr <- log_pf - log_pc
  data.frame(reference = reference,
             delta_G = gf - gc, return_period = exp(-log_pf),
             return_period_counterfactual = exp(-log_pc),
             gmst_factual = gf, gmst_counterfactual = gc, threshold_C = event,
             pot_threshold_C = fit$threshold,
             p_factual = exp(log_pf), p_counterfactual = exp(log_pc),
             log_RR = log_rr, RR = exp(log_rr), FAR = -expm1(-log_rr),
             intensity_change = fit$alpha * (gf - gc),
             lambda_factual = exp(log_pp_rate(fit$threshold, mu_f, fit$sigma, fit$xi)),
             lambda_counterfactual = exp(log_pp_rate(fit$threshold, mu_c, fit$sigma, fit$xi)),
             gpd_scale_factual = if (fit$sigma + fit$xi * (fit$threshold - mu_f) > 0)
               fit$sigma + fit$xi * (fit$threshold - mu_f) else NA_real_,
             gpd_scale_counterfactual = if (fit$sigma + fit$xi * (fit$threshold - mu_c) > 0)
               fit$sigma + fit$xi * (fit$threshold - mu_c) else NA_real_)
}

bootstrap_pot <- function(years_fit, peaks, u, x_event, g_factual,
                          g_counterfactual, point_results, bootstrap_B,
                          random_seed, ci_level, reference) {
  # Resample complete summers with their GMST and all associated cluster peaks.
  # Retain zero-event summers and duplicate peak sets when a summer is sampled
  # more than once. Refit both frequency and magnitude jointly in each replicate.
  # The event magnitude, threshold, declustering and reference GMST are fixed.
  # Intervals exclude uncertainty in threshold, run length and model selection.
  metrics <- c("RR", "FAR", "intensity_change", "return_period",
               "return_period_counterfactual", "p_factual", "p_counterfactual")
  bootstrap_samples <- as.data.frame(matrix(NA_real_, bootstrap_B, length(metrics)))
  names(bootstrap_samples) <- metrics
  bootstrap_samples$success <- rep(FALSE, bootstrap_B)
  bootstrap_samples$status <- rep("not_run", bootstrap_B)
  set.seed(random_seed)
  for (b in seq_len(bootstrap_B)) {
    sampled <- sample.int(nrow(years_fit), nrow(years_fit), replace = TRUE)
    boot_years <- years_fit[sampled, ]
    boot_peaks <- do.call(rbind, lapply(sampled, function(j) {
      peaks[peaks$year == years_fit$year[j], ]
    }))
    result <- tryCatch({
      fitted <- fit_pp(boot_years, boot_peaks, u)
      estimate <- pp_attribution(fitted, x_event, g_factual, g_counterfactual, reference)
      values <- unlist(estimate[metrics], use.names = FALSE)
      if (any(is.na(values))) stop("Undefined bootstrap attribution statistic.")
      values
    }, error = function(e) e)
    if (inherits(result, "error")) {
      bootstrap_samples$status[b] <- conditionMessage(result)
    } else {
      bootstrap_samples[b, metrics] <- as.list(result)
      bootstrap_samples$success[b] <- TRUE
      bootstrap_samples$status[b] <- "ok"
    }
    if (b %% 50L == 0L) message("POT bootstrap ", b, " of ", bootstrap_B)
  }
  bootstrap_success <- sum(bootstrap_samples$success)
  results <- point_results
  results$bootstrap_success <- bootstrap_success
  results$bootstrap_requested <- bootstrap_B
  results$ci_level <- ci_level
  # At least 80% success is required to report intervals; inspect all failures.
  minimum_success <- max(20L, ceiling(0.8 * bootstrap_B))
  for (metric in metrics) {
    interval <- c(NA_real_, NA_real_)
    if (bootstrap_success >= minimum_success) {
      # Type 1 keeps infinite bounds rather than interpolating Inf and finite values.
      interval <- quantile(bootstrap_samples[bootstrap_samples$success, metric],
                           c((1 - ci_level) / 2, (1 + ci_level) / 2),
                           names = FALSE, type = 1)
    }
    results[[paste0(metric, "_low")]] <- interval[1]
    results[[paste0(metric, "_high")]] <- interval[2]
  }
  if (bootstrap_B > 0L && bootstrap_success < minimum_success)
    warning("Too few successful bootstrap fits; confidence intervals are NA.")
  results

  list(results = results, samples = bootstrap_samples)
}




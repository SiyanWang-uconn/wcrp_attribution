# Reusable, inspectable steps for daily gridded temperature attribution.
# Dataset scripts call these functions one step at a time and keep each returned
# object in the Global Environment for inspection.

open_temperature_dataset <- function(file, variable, kelvin_to_celsius = FALSE) {
  nc <- ncdf4::nc_open(file)
  metadata <- nc$var[[variable]]
  dimension_names <- vapply(metadata$dim, function(x) x$name, character(1))
  stopifnot(identical(dimension_names, c("lon", "lat", "time")))

  lon_original <- nc$dim$lon$vals
  lon <- (lon_original + 180) %% 360 - 180
  lat <- nc$dim$lat$vals
  time_values <- nc$dim$time$vals
  time_units <- nc$dim$time$units
  time_origin <- as.Date(substr(sub("^days since ", "", time_units), 1, 10))
  date <- time_origin + time_values

  sample_array <- ncdf4::ncvar_get(
    nc, variable, start = c(1, 1, 1), count = c(-1, -1, 1),
    collapse_degen = FALSE
  )
  sample_grid <- sample_array[, , 1]
  if (kelvin_to_celsius) sample_grid <- sample_grid - 273.15

  list(
    nc = nc,
    metadata = metadata,
    dimension_names = dimension_names,
    lon_original = lon_original,
    lon = lon,
    lat = lat,
    date = date,
    sample_grid = sample_grid
  )
}

find_persistent_footprint <- function(nc, variable, date, nlon, nlat,
                                      minimum_fraction = 0.95) {
  number_of_days <- length(date)
  number_of_cells <- nlon * nlat
  chunk_starts <- seq(1L, number_of_days, by = 366L)
  valid_day_count <- rep(0, number_of_cells)

  for (chunk_start in chunk_starts) {
    chunk_end <- min(chunk_start + 365L, number_of_days)
    chunk_days <- chunk_start:chunk_end
    values <- ncdf4::ncvar_get(
      nc, variable, start = c(1, 1, chunk_start),
      count = c(-1, -1, length(chunk_days)), collapse_degen = FALSE
    )
    value_matrix <- matrix(values, nrow = number_of_cells)
    valid_day_count <- valid_day_count + rowSums(is.finite(value_matrix))
    message("Footprint: ", date[chunk_start], " to ", date[chunk_end])
  }

  valid_fraction <- valid_day_count / number_of_days
  list(
    included = valid_fraction >= minimum_fraction,
    valid_fraction = valid_fraction
  )
}

map_external_mask <- function(lon, lat, mask_file) {
  mask_table <- read.csv(mask_file)
  stopifnot(all(c("lon", "lat", "included") %in% names(mask_table)))
  mask_lon <- unique(mask_table$lon)
  mask_lat <- unique(mask_table$lat)
  mask_grid <- matrix(
    as.logical(mask_table$included),
    nrow = length(mask_lon), ncol = length(mask_lat)
  )
  nearest_lon <- vapply(lon, function(x) which.min(abs(mask_lon - x)), integer(1))
  nearest_lat <- vapply(lat, function(x) which.min(abs(mask_lat - x)), integer(1))
  as.vector(mask_grid[nearest_lon, nearest_lat, drop = FALSE])
}

calculate_daily_regional_temperature <- function(nc, variable, date, lon, lat,
                                                  footprint,
                                                  kelvin_to_celsius = FALSE) {
  number_of_days <- length(date)
  number_of_cells <- length(lon) * length(lat)
  chunk_starts <- seq(1L, number_of_days, by = 366L)
  latitude_weights <- cos(lat * pi / 180)
  grid_weights <- rep(latitude_weights, each = length(lon))
  grid_weights[!footprint] <- 0
  daily_value <- rep(NA_real_, number_of_days)
  coverage <- rep(NA_real_, number_of_days)

  for (chunk_start in chunk_starts) {
    chunk_end <- min(chunk_start + 365L, number_of_days)
    chunk_days <- chunk_start:chunk_end
    values <- ncdf4::ncvar_get(
      nc, variable, start = c(1, 1, chunk_start),
      count = c(-1, -1, length(chunk_days)), collapse_degen = FALSE
    )
    if (kelvin_to_celsius) values <- values - 273.15
    value_matrix <- matrix(values, nrow = number_of_cells)
    valid <- is.finite(value_matrix)
    denominator <- colSums(valid * grid_weights)
    value_matrix[!valid] <- 0
    daily_value[chunk_days] <- colSums(value_matrix * grid_weights) / denominator
    coverage[chunk_days] <- denominator / sum(grid_weights)
    message("Regional mean: ", date[chunk_start], " to ", date[chunk_end])
  }

  daily_value[coverage < 0.95] <- NA_real_
  data.frame(date = date, value = daily_value, spatial_coverage = coverage)
}

calculate_three_day_extremes <- function(daily_data) {
  daily_data$year <- as.integer(format(daily_data$date, "%Y"))
  daily_data$month <- as.integer(format(daily_data$date, "%m"))
  daily_data$value3 <- as.numeric(
    stats::filter(daily_data$value, rep(1 / 3, 3), sides = 1)
  )

  years <- sort(unique(daily_data$year))
  annual <- data.frame()
  june <- data.frame()

  for (selected_year in years) {
    annual_start <- as.Date(sprintf("%d-01-01", selected_year))
    annual_end <- as.Date(sprintf("%d-12-31", selected_year))
    annual_block <- daily_data[
      daily_data$date >= annual_start & daily_data$date <= annual_end,
    ]
    expected <- as.integer(annual_end - annual_start) + 1L
    complete <- nrow(annual_block) == expected &&
      mean(is.finite(annual_block$value)) >= 0.95
    candidates <- annual_block$value3[annual_block$date >= annual_start + 2]
    extreme <- if (complete) max(candidates, na.rm = TRUE) else NA_real_
    annual <- rbind(annual, data.frame(
      year = selected_year, season = "Annual", complete = complete, X = extreme
    ))

    june_start <- as.Date(sprintf("%d-06-01", selected_year))
    june_end <- as.Date(sprintf("%d-06-30", selected_year))
    june_block <- daily_data[
      daily_data$date >= june_start & daily_data$date <= june_end,
    ]
    complete <- nrow(june_block) == 30L &&
      mean(is.finite(june_block$value)) >= 0.95
    candidates <- june_block$value3[june_block$date >= june_start + 2]
    extreme <- if (complete) max(candidates, na.rm = TRUE) else NA_real_
    june <- rbind(june, data.frame(
      year = selected_year, season = "June", complete = complete, X = extreme
    ))
  }

  list(daily = daily_data, annual = annual, june = june,
       combined = rbind(annual, june))
}

prepare_smoothed_gmst <- function(gmst_file) {
  raw <- read.csv(gmst_file, skip = 1, check.names = FALSE,
                  na.strings = c("***", "****", "*****"))
  months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
  monthly <- as.matrix(raw[, months])
  storage.mode(monthly) <- "double"
  annual_value <- rowMeans(monthly)
  annual_value[rowSums(is.finite(monthly)) != 12] <- NA_real_
  observed <- data.frame(year = as.integer(raw$Year), annual = annual_value)
  observed <- observed[is.finite(observed$annual), ]
  endpoint_trend <- lm(annual ~ year, data = tail(observed, 15))
  complete <- merge(
    data.frame(year = seq(min(observed$year), 2027L)), observed,
    by = "year", all.x = TRUE
  )
  complete$extrapolated <- !is.finite(complete$annual)
  complete$annual[complete$extrapolated] <- predict(
    endpoint_trend, newdata = complete[complete$extrapolated, ]
  )
  complete$gmst <- NA_real_
  for (row_number in seq_len(nrow(complete))) {
    y <- complete$year[row_number]
    values <- complete$annual[match((y - 2):(y + 1), complete$year)]
    if (!anyNA(values)) complete$gmst[row_number] <- mean(values)
  }
  list(raw = raw, observed = observed, endpoint_trend = endpoint_trend,
       complete = complete,
       gmst_2026 = complete$gmst[complete$year == 2026])
}

join_extremes_and_gmst <- function(extremes, gmst_complete) {
  joined <- merge(extremes, gmst_complete[, c("year", "gmst")],
                  by = "year", all.x = TRUE)
  joined$used_in_fit <- joined$complete & joined$year <= 2025 &
    is.finite(joined$gmst)
  list(
    all = joined,
    annual = joined[joined$season == "Annual" & joined$used_in_fit, ],
    june = joined[joined$season == "June" & joined$used_in_fit, ]
  )
}

bootstrap_attribution <- function(fit, data, gmst_factual, gmst_counterfactual, B,
                                  seed = 2026) {
  stopifnot(
    all(c("year", "tmax", "gmst") %in% names(data)),
    nrow(data) > 0L, all(is.finite(data$tmax)),
    all(is.finite(data$gmst)), sum(data$year == 2021, na.rm = TRUE) == 1L,
    length(gmst_factual) == 1L, is.finite(gmst_factual),
    length(gmst_counterfactual) == 1L, is.finite(gmst_counterfactual),
    length(B) == 1L, is.finite(B), B >= 1, B == floor(B),
    length(seed) == 1L, is.finite(seed)
  )
  set.seed(seed)
  # Store three statistics for the single counterfactual climate.
  simulations <- matrix(
    NA_real_, nrow = 3, ncol = B,
    dimnames = list(c("RR", "FAR", "intensity"), NULL)
  )
  location <- fit$mu + fit$alpha * (data$gmst - fit$g0)
  for (b in seq_len(B)) {
    simulated <- qgev_local(runif(nrow(data)), location, fit$sigma, fit$xi)
    new_fit <- tryCatch(fit_gev(simulated, data$gmst), error = function(e) NULL)
    if (!is.null(new_fit)) {
      # Keep the observed 2021 event threshold fixed across bootstrap fits.
      result <- attribution(new_fit, gmst_factual, gmst_counterfactual, data)
      simulations[, b] <- c(result$RR, result$FAR, result$intensity_change)
    }
    if (b %% 20 == 0) message("Bootstrap ", b, " of ", B)
  }
  simulations
}

add_bootstrap_intervals <- function(point_results, annual_bootstrap,
                                     level = 0.95) {
  stopifnot(is.matrix(annual_bootstrap),
            all(c("RR", "FAR", "intensity") %in% rownames(annual_bootstrap)),
            length(level) == 1L, is.finite(level), level > 0, level < 1)
  alpha <- (1 - level) / 2
  limits <- c(alpha, 1 - alpha)
  result <- point_results
  result$RR_low <- NA_real_
  result$RR_high <- NA_real_
  result$FAR_low <- NA_real_
  result$FAR_high <- NA_real_
  result$intensity_change_low <- NA_real_
  result$intensity_change_high <- NA_real_
  result$bootstrap_success <- NA_integer_

  for (i in seq_len(nrow(result))) {
    boot <- annual_bootstrap
    complete <- colSums(is.na(boot)) == 0
    result$bootstrap_success[i] <- sum(complete)
    if (!any(complete)) next
    for (statistic in c("RR", "FAR", "intensity")) {
      values <- boot[statistic, complete]
      # Empirical order statistics preserve valid infinite bounds.
      interval <- quantile(values, limits, na.rm = TRUE, names = FALSE, type = 1)
      column_stem <- if (statistic == "intensity") "intensity_change" else statistic
      result[i, paste0(column_stem, "_low")] <- interval[1]
      result[i, paste0(column_stem, "_high")] <- interval[2]
    }
  }
  result
}

write_bootstrap_samples <- function(annual_bootstrap, june_bootstrap, file) {
  to_data_frame <- function(x, season) {
    data.frame(
      iteration = seq_len(ncol(x)),
      season = season,
      as.data.frame(t(x)),
      check.names = FALSE
    )
  }
  samples <- rbind(
    to_data_frame(annual_bootstrap, "Annual"),
    to_data_frame(june_bootstrap, "June")
  )
  write.csv(samples, file, row.names = FALSE)
  invisible(samples)
}

# RR inference uses the same likelihood and parameterization as fit_gev().
# The threshold and both GMST values are held fixed during inference.
# No uncertainty in the observed threshold or GMST is propagated here.
rr_ci_setup <- function(fit, data, gmst_factual, gmst_past,
                        threshold, level) {
  stopifnot(all(c("tmax", "gmst") %in% names(data)),
            nrow(data) > 4L, all(is.finite(data$tmax)),
            all(is.finite(data$gmst)),
            length(gmst_factual) == 1L, is.finite(gmst_factual),
            length(gmst_past) == 1L, is.finite(gmst_past),
            length(threshold) == 1L, is.finite(threshold),
            length(level) == 1L, level > 0, level < 1)
  theta <- c(mu = unname(fit$mu), alpha = unname(fit$alpha),
             log_sigma = log(unname(fit$sigma)), xi = unname(fit$xi))
  stopifnot(all(is.finite(theta)), is.finite(fit$g0))
  centered_gmst <- data$gmst - fit$g0
  nll <- function(par) {
    if (any(!is.finite(par))) return(1e100)
    sigma <- exp(par[3])
    if (!is.finite(sigma) || sigma < 1e-6) return(1e100)
    z <- (data$tmax - par[1] - par[2] * centered_gmst) / sigma
    xi <- par[4]
    if (any(!is.finite(z))) return(1e100)
    if (abs(xi) < 1e-7) {
      value <- sum(log(sigma) + z + exp(-z))
    } else {
      h <- 1 + xi * z
      if (any(h <= 0)) return(1e100)
      value <- sum(log(sigma) + (1 + 1 / xi) * log(h) +
                     exp(-log(h) / xi))
    }
    if (is.finite(value)) value else 1e100
  }
  log_rr <- function(par) {
    mu_f <- par[1] + par[2] * (gmst_factual - fit$g0)
    mu_p <- par[1] + par[2] * (gmst_past - fit$g0)
    unname(log_survival(threshold, mu_f, exp(par[3]), par[4]) -
      log_survival(threshold, mu_p, exp(par[3]), par[4]))
  }
  list(theta = theta, nll = nll, log_rr = log_rr)
}

# Delta intervals use the inverse observed Hessian and log(RR) gradient.
# Pass the same point_results used by add_bootstrap_intervals() to obtain
# identical table columns. Diagnostics are stored in the delta_diagnostics
# attribute. bootstrap_success is NA because no bootstrap is performed.
rr_ci_delta <- function(fit, data, gmst_factual, gmst_past,
                        threshold = data$tmax[data$year==2021], level = 0.95,
                        point_results = NULL) {
  setup <- rr_ci_setup(fit, data, gmst_factual, gmst_past, threshold, level)
  theta <- setup$theta
  estimate <- setup$log_rr(theta)
  delta_g <- gmst_factual - gmst_past
  if (is.null(point_results)) {
    mu_f <- fit$mu + fit$alpha * (gmst_factual - fit$g0)
    mu_p <- fit$mu + fit$alpha * (gmst_past - fit$g0)
    log_pf <- log_survival(threshold, mu_f, fit$sigma, fit$xi)
    log_pc <- log_survival(threshold, mu_p, fit$sigma, fit$xi)
    point_results <- data.frame(
      reference = "1951-like", delta_G = delta_g,
      return_period = exp(-log_pf), gmst_factual = gmst_factual,
      threshold_C = threshold, p_factual = exp(log_pf),
      p_counterfactual = exp(log_pc), log_RR = estimate,
      RR = exp(estimate), FAR = -expm1(-estimate),
      intensity_change = unname(fit$alpha * delta_g)
    )
  }
  stopifnot(is.data.frame(point_results), nrow(point_results) == 1L,
            all(c("RR", "FAR", "intensity_change") %in% names(point_results)))
  # Reject a table referring to a different event or climate comparison.
  for (column in c("threshold_C", "gmst_factual", "delta_G")) {
    expected <- switch(column, threshold_C = threshold,
                       gmst_factual = gmst_factual, delta_G = delta_g)
    if (column %in% names(point_results))
      stopifnot(isTRUE(all.equal(as.numeric(point_results[[column]]),
                                as.numeric(expected))))
  }
  result <- point_results
  result$RR_low <- NA_real_
  result$RR_high <- NA_real_
  result$FAR_low <- NA_real_
  result$FAR_high <- NA_real_
  result$intensity_change_low <- NA_real_
  result$intensity_change_high <- NA_real_
  result$bootstrap_success <- NA_integer_
  diagnostics <- list(method = "delta", ci_level = level,
                      log_RR_se = NA_real_, intensity_change_se = NA_real_,
                      status = "not_computed")
  finish <- function() {
    attr(result, "delta_diagnostics") <- diagnostics
    result
  }
  if (gmst_factual == gmst_past) {
    result$RR_low <- result$RR_high <- 1
    result$FAR_low <- result$FAR_high <- 0
    result$intensity_change_low <- result$intensity_change_high <- 0
    diagnostics$log_RR_se <- diagnostics$intensity_change_se <- 0
    diagnostics$status <- "identical_climates"
    return(finish())
  }
  if (theta[4] <= -0.5) {
    diagnostics$status <- "nonregular_shape"
    return(finish())
  }
  calculation <- tryCatch({
    hessian <- optimHess(theta, setup$nll)
    if (any(!is.finite(hessian)) ||
        min(eigen(hessian, symmetric = TRUE, only.values = TRUE)$values) <= 0)
      stop("Observed Hessian is not positive definite")
    covariance <- chol2inv(chol(hessian))
    intensity_se <- abs(delta_g) * sqrt(covariance[2, 2])
    if (!is.finite(estimate)) {
      return_value <- list(log_rr_se = NA_real_, intensity_se = intensity_se)
    } else {
      step <- .Machine$double.eps^(1 / 3) * pmax(abs(theta), 1)
      gradient <- vapply(seq_along(theta), function(j) {
        plus <- minus <- theta
        plus[j] <- plus[j] + step[j]
        minus[j] <- minus[j] - step[j]
        (setup$log_rr(plus) - setup$log_rr(minus)) / (2 * step[j])
      }, 0.0)
      variance <- as.numeric(t(gradient) %*% covariance %*% gradient)
      if (!is.finite(variance) || variance < 0) stop("Invalid delta variance")
      return_value <- list(log_rr_se = sqrt(variance),
                           intensity_se = intensity_se)
    }
    return_value
  }, error = function(e) e)
  if (inherits(calculation, "error")) {
    diagnostics$status <- conditionMessage(calculation)
    return(finish())
  }
  diagnostics$log_RR_se <- calculation$log_rr_se
  diagnostics$intensity_change_se <- calculation$intensity_se
  z <- qnorm((1 + level) / 2)
  intensity_limits <- unname(fit$alpha * delta_g) +
    c(-1, 1) * z * calculation$intensity_se
  result$intensity_change_low <- intensity_limits[1]
  result$intensity_change_high <- intensity_limits[2]
  if (!is.finite(estimate)) {
    diagnostics$status <- "nonfinite_log_RR"
    return(finish())
  }
  limits <- estimate + c(-1, 1) * z * calculation$log_rr_se
  result$RR_low <- exp(limits[1])
  result$RR_high <- exp(limits[2])
  # FAR = 1 - 1/RR is monotone; transform the log(RR) interval.
  result$FAR_low <- -expm1(-limits[1])
  result$FAR_high <- -expm1(-limits[2])
  diagnostics$status <- "ok"
  finish()
}

# Profile likelihood interval, using the asymptotic chi-square(1) LR cutoff.
# At each fixed log(RR), optimize three nuisance parameters. Parameterizing
# factual/counterfactual exceedance probabilities enforces RR exactly.
# search_distance limits the search on the log(RR) scale on either side.
# An unreached endpoint is NA with status search_limit. An infinite upper
# endpoint requires a feasible zero-counterfactual-probability boundary fit
# inside the LR confidence set; optimization failure never implies infinity.
# The returned profile table allows inspection of LR values and convergence.
rr_ci_lrt <- function(fit, data, gmst_factual, gmst_past,
                      threshold = data$tmax[data$year==2021], level = 0.95,
                      search_distance = 40, maxit = 3000L,
                      root_tolerance = 1e-4) {
  stopifnot(is.finite(search_distance), search_distance > 0,
            maxit > 0, root_tolerance > 0)
  setup <- rr_ci_setup(fit, data, gmst_factual, gmst_past, threshold, level)
  theta <- setup$theta
  estimate <- setup$log_rr(theta)
  baseline <- setup$nll(theta)
  cutoff <- qchisq(level, df = 1)
  result <- data.frame(
    RR = exp(estimate), RR_lrt_low = NA_real_, RR_lrt_high = NA_real_,
    log_RR = estimate, ci_level = level,
    lrt_low_status = "not_computed", lrt_high_status = "not_computed"
  )
  trace <- data.frame(log_RR = numeric(), nll = numeric(),
                      LR = numeric(), convergence = integer())
  finish <- function() list(interval = result, profile = trace, cutoff = cutoff)
  if (!is.finite(estimate) || theta[4] <= -0.5) {
    result$lrt_low_status <- result$lrt_high_status <-
      "nonfinite_log_RR_or_nonregular_shape"
    return(finish())
  }
  if (gmst_factual == gmst_past) {
    result$RR_lrt_low <- result$RR_lrt_high <- 1
    result$lrt_low_status <- result$lrt_high_status <- "identical_climates"
    return(finish())
  }
  log_pf <- log_survival(
    threshold, fit$mu + fit$alpha * (gmst_factual - fit$g0),
    fit$sigma, fit$xi
  )
  # Stable standardized quantile from a log exceedance probability.
  standard_quantile <- function(log_p, xi) {
    log_hazard <- if (log_p < -35) log_p else log(-log1p(-exp(log_p)))
    if (abs(xi) < 1e-7) -log_hazard else expm1(-xi * log_hazard) / xi
  }
  to_theta <- function(nuisance, target) {
    # softplus keeps both exceedance probabilities strictly below one.
    eta <- nuisance[1]
    softplus <- max(eta, 0) + log1p(exp(-abs(eta)))
    log_f <- -softplus - max(0, -target)
    log_p <- log_f - target
    sigma <- exp(nuisance[2])
    xi <- nuisance[3]
    mu_f <- threshold - sigma * standard_quantile(log_f, xi)
    mu_p <- threshold - sigma * standard_quantile(log_p, xi)
    slope <- (mu_f - mu_p) / (gmst_factual - gmst_past)
    c(mu_f - slope * (gmst_factual - fit$g0), slope, nuisance[2], xi)
  }
  profile_at <- function(target) {
    cached <- which(abs(trace$log_RR - target) < 1e-10)
    if (length(cached)) return(trace$LR[cached[1]] - cutoff)
    objective <- function(nuisance) {
      if (any(!is.finite(nuisance))) return(1e100)
      setup$nll(to_theta(nuisance, target))
    }
    softplus_start <- max(-log_pf - max(0, -target), 0.01)
    eta_start <- if (softplus_start > 35) softplus_start else
      log(expm1(softplus_start))
    starts <- lapply(unique(c(theta[4], -0.1, 0.1)), function(shape)
      c(eta_start, theta[3], shape))
    fits <- lapply(starts, function(start) tryCatch(
      optim(start, objective, method = "Nelder-Mead",
            control = list(maxit = maxit, reltol = 1e-10)),
      error = function(e) NULL
    ))
    good <- vapply(fits, function(x)
      !is.null(x) && x$convergence == 0 && is.finite(x$value) &&
        x$value < 1e90, TRUE)
    if (!any(good)) stop("Profile optimization failed at log(RR) = ", target)
    best <- fits[good][[which.min(vapply(fits[good], function(x) x$value, 0.0))]]
    # A substantially better constrained fit means the supplied MLE needs review.
    if (best$value < baseline - 1e-4)
      stop("Profile improves on supplied MLE; refit the unconstrained model")
    lr <- 2 * (best$value - baseline)
    trace <<- rbind(trace, data.frame(log_RR = target, nll = best$value,
                                     LR = lr, convergence = best$convergence))
    lr - cutoff
  }
  infinite_upper_boundary <- function() {
    # For xi < 0, endpoint_past = mu_past - sigma / xi.
    # Set endpoint_past strictly below the threshold, giving p_past = 0.
    boundary_theta <- function(par) {
      sigma <- exp(par[2])
      xi <- -exp(par[3])
      endpoint_past <- threshold - exp(par[4])
      mu_past <- endpoint_past + sigma / xi
      c(mu_past - par[1] * (gmst_past - fit$g0), par[1], par[2], xi)
    }
    objective <- function(par) {
      if (any(!is.finite(par))) return(1e100)
      candidate <- boundary_theta(par)
      value <- setup$nll(candidate)
      mu_f <- candidate[1] + candidate[2] * (gmst_factual - fit$g0)
      log_f <- log_survival(threshold, mu_f, exp(candidate[3]), candidate[4])
      if (!is.finite(log_f)) return(1e100)
      value
    }
    starts <- list()
    for (multiplier in c(0.5, 1, 1.5, 2)) {
      slope <- theta[2] * multiplier
      shifted_max <- max(data$tmax - slope * (data$gmst - gmst_past))
      gap <- (threshold - shifted_max) / 2
      if (gap > 0) {
        for (shape in c(-0.1, -0.25, -0.4)) {
          # Choose sigma so that the boundary location fits the data scale.
          starts[[length(starts) + 1L]] <-
            c(slope, theta[3], log(-shape), log(gap))
        }
      }
    }
    fits <- lapply(starts, function(start) tryCatch(
      optim(start, objective, method = "Nelder-Mead",
            control = list(maxit = maxit, reltol = 1e-10)),
      error = function(e) NULL
    ))
    good <- vapply(fits, function(x)
      !is.null(x) && x$convergence == 0 && is.finite(x$value) &&
        x$value < 1e90, TRUE)
    if (!any(good)) return(FALSE)
    best <- fits[good][[which.min(vapply(fits[good], function(x) x$value, 0.0))]]
    if (best$value < baseline - 1e-4) return(FALSE)
    candidate <- boundary_theta(best$par)
    if (!is.infinite(setup$log_rr(candidate)) || setup$log_rr(candidate) < 0)
      return(FALSE)
    lr <- 2 * (best$value - baseline)
    trace <<- rbind(trace, data.frame(log_RR = Inf, nll = best$value,
                                     LR = lr, convergence = best$convergence))
    lr <= cutoff
  }
  endpoint <- function(direction) {
    tryCatch({
      inside <- estimate
      if (profile_at(inside) >= 0) stop("Profile does not reproduce supplied MLE")
      distance <- min(0.5, search_distance)
      repeat {
        outside <- estimate + direction * distance
        if (profile_at(outside) >= 0) {
          bracket <- sort(c(inside, outside))
          root <- uniroot(profile_at, bracket, tol = root_tolerance)$root
          return(list(value = exp(root), status = "ok"))
        }
        if (distance >= search_distance)
          return(list(value = NA_real_, status = "search_limit"))
        inside <- outside
        distance <- min(2 * distance, search_distance)
      }
    }, error = function(e) list(value = NA_real_, status = conditionMessage(e)))
  }
  lower <- endpoint(-1)
  upper <- endpoint(1)
  if (upper$status == "search_limit" &&
      isTRUE(tryCatch(infinite_upper_boundary(), error = function(e) FALSE))) {
    upper <- list(value = Inf, status = "unbounded_likelihood_set")
  }
  result$RR_lrt_low <- lower$value
  result$RR_lrt_high <- upper$value
  result$lrt_low_status <- lower$status
  result$lrt_high_status <- upper$status
  trace <- trace[order(trace$log_RR), ]
  finish()
}

# Parametric bootstrap used for Figure 7 return-level confidence bands.
bootstrap_return_level_curves <- function(fit, data, gmst_factual, gmst_past,
                                          probabilities, B,
                                          seed = 20260922) {
  set.seed(seed)
  current_curves <- matrix(NA_real_, B, length(probabilities))
  past_curves <- matrix(NA_real_, B, length(probabilities))
  success <- logical(B)
  fitted_location <- fit$mu + fit$alpha * (data$gmst - fit$g0)

  for (b in seq_len(B)) {
    simulated <- qgev_local(
      runif(nrow(data)), fitted_location, fit$sigma, fit$xi
    )
    new_fit <- tryCatch(fit_gev(simulated, data$gmst), error = function(e) NULL)
    if (!is.null(new_fit)) {
      current_location <- new_fit$mu +
        new_fit$alpha * (gmst_factual - new_fit$g0)
      past_location <- new_fit$mu +
        new_fit$alpha * (gmst_past - new_fit$g0)
      current_curves[b, ] <- qgev_local(
        probabilities, current_location, new_fit$sigma, new_fit$xi
      )
      past_curves[b, ] <- qgev_local(
        probabilities, past_location, new_fit$sigma, new_fit$xi
      )
      success[b] <- TRUE
    }
    if (b %% 20 == 0) message("Figure 7 bootstrap ", b, " of ", B)
  }

  confidence_band <- function(curves) {
    apply(curves[success, , drop = FALSE], 2, quantile,
          probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE, type = 1)
  }

  list(
    current_curves = current_curves,
    past_curves = past_curves,
    current_ci = confidence_band(current_curves),
    past_ci = confidence_band(past_curves),
    success = success
  )
}

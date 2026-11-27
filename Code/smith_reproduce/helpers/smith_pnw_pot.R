# Helpers for Smith (revised April 2025), Section 0.6.2: PNW POT MLE.
# Definitions only. No existing project helper is sourced or modified.
# Package references:
# https://search.r-project.org/CRAN/refmans/ismev/html/pp.fit.html
# https://search.r-project.org/CRAN/refmans/texmex/html/extremalIndex.html

smith_check_packages <- function() {
  required <- c("ncdf4", "quantreg", "texmex", "ismev")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "),
                            ". Install these separately before running the analysis.")
}

smith_read_daily <- function(file) {
  nc <- ncdf4::nc_open(file)
  on.exit(ncdf4::nc_close(nc))
  stopifnot(identical(vapply(nc$var$tmax$dim, function(z) z$name, ""), "time"))
  stopifnot(ncdf4::ncatt_get(nc, "tmax", "units")$value %in%
              c("degC", "C", "Celsius", "degree_Celsius"))
  units <- nc$dim$time$units
  stopifnot(grepl("^(hours|days) since ", units))
  origin <- as.Date(substr(sub("^(hours|days) since ", "", units), 1, 10))
  offset <- as.numeric(ncdf4::ncvar_get(nc, "time")) /
    if (startsWith(units, "hours")) 24 else 1
  stopifnot(all(is.finite(offset)), all(abs(offset - round(offset)) < 1e-7))
  d <- data.frame(date = origin + round(offset),
                   tx = as.numeric(ncdf4::ncvar_get(nc, "tmax")))
  stopifnot(!anyDuplicated(d$date), all(diff(d$date) == 1))
  d$year <- as.integer(format(d$date, "%Y"))
  d$month_day <- format(d$date, "%m-%d")
  leap <- d$year %% 4 == 0 & (d$year %% 100 != 0 | d$year %% 400 == 0)
  # Use a reproducible 365-day calendar: dates after February lose the leap day.
  d$season_day <- as.integer(format(d$date, "%j")) -
    as.integer(leap & as.integer(format(d$date, "%m")) > 2)
  d
}

smith_season <- function(d, mode, first_day, last_day, first_date, last_date) {
  stopifnot(mode %in% c("paper_doy", "calendar"))
  if (mode == "paper_doy") {
    stopifnot(first_day >= 1, last_day <= 365, first_day <= last_day)
    keep <- d$season_day >= first_day & d$season_day <= last_day
  } else {
    a <- as.Date(paste0("2001-", first_date), format = "%Y-%m-%d")
    b <- as.Date(paste0("2001-", last_date), format = "%Y-%m-%d")
    stopifnot(!is.na(a), !is.na(b), a <= b,
              format(a, "%m-%d") == first_date, format(b, "%m-%d") == last_date)
    keep <- d$month_day >= first_date & d$month_day <= last_date
  }
  d[keep & d$month_day != "02-29", , drop = FALSE]
}

smith_prepare <- function(daily, gmst, start_year, end_year, mode,
                           first_day, last_day, first_date, last_date) {
  stopifnot(all(c("year", "gmst") %in% names(gmst)), !anyDuplicated(gmst$year))
  d <- smith_season(daily, mode, first_day, last_day, first_date, last_date)
  d <- d[d$year >= start_year & d$year <= end_year, , drop = FALSE]
  expected <- if (mode == "paper_doy") last_day - first_day + 1L else
    as.integer(as.Date(paste0("2001-", last_date)) -
                 as.Date(paste0("2001-", first_date))) + 1L
  coverage <- do.call(rbind, lapply(seq.int(start_year, end_year), function(y) {
    z <- d[d$year == y, , drop = FALSE]
    g <- gmst$gmst[match(y, gmst$year)]
    data.frame(year = y, expected_days = expected, available_days = nrow(z),
               valid_days = sum(is.finite(z$tx)), gmst = g,
               used_in_fit = nrow(z) == expected && all(is.finite(z$tx)) && is.finite(g))
  }))
  d <- d[d$year %in% coverage$year[coverage$used_in_fit], , drop = FALSE]
  stopifnot(length(unique(d$year)) >= 10L)
  d$gmst <- gmst$gmst[match(d$year, gmst$year)]
  center <- mean(coverage$gmst[coverage$used_in_fit])
  d <- smith_features(d, center, start_year)
  rownames(d) <- NULL
  list(data = d, coverage = coverage, gmst_center = center)
}

smith_features <- function(d, gmst_center, trend_origin) {
  phase <- 2 * pi * (d$season_day - 1) / 365
  d$sin1 <- sin(phase)
  d$cos1 <- cos(phase)
  d$sin2 <- sin(2 * phase)
  d$cos2 <- cos(2 * phase)
  d$gmst_centered <- d$gmst - gmst_center
  d$trend_decades <- (d$year - trend_origin) / 10
  d
}

smith_threshold <- function(d, tau, harmonic_order) {
  stopifnot(tau > 0, tau < 1, harmonic_order %in% c(1L, 2L))
  # The paper does not supply its exact threshold design matrix. Use calendar
  # year as the long-term trend, with explicit selectable seasonal harmonics.
  terms <- c("trend_decades", "sin1", "cos1",
              if (harmonic_order == 2L) c("sin2", "cos2"))
  model <- quantreg::rq(stats::reformulate(terms, response = "tx"),
                        data = d, tau = tau, method = "fn")
  u <- as.numeric(stats::predict(model, newdata = d))
  stopifnot(length(u) == nrow(d), all(is.finite(u)))
  list(model = model, threshold = u)
}

smith_decluster <- function(d, threshold) {
  # Standardize the exceedance decision, not the physical temperature scale.
  # Concatenate complete seasonal windows for the intervals estimate; this
  # convention is explicit because the paper does not specify winter handling.
  # Cluster boundaries are always forced between years after automatic clustering.
  residual <- d$tx - threshold
  ei <- texmex::extremalIndex(residual, threshold = 0)
  dc <- texmex::declust(ei)
  index <- as.integer(dc$exceedanceTimes)
  stopifnot(identical(index, which(residual > 0)), length(index) == length(dc$clusters))
  groups <- interaction(d$year[index], dc$clusters, drop = TRUE)
  peaks <- do.call(rbind, lapply(split(index, groups), function(ii) {
    # Select the highest actual temperature, not the largest threshold residual.
    p <- ii[which.max(d$tx[ii])]
    data.frame(row = p, year = d$year[p], date = d$date[p],
               tx = d$tx[p], threshold = threshold[p], excess = residual[p],
               start_date = min(d$date[ii]), end_date = max(d$date[ii]),
               exceedance_days = length(ii))
  }))
  peaks <- peaks[order(peaks$date), , drop = FALSE]
  rownames(peaks) <- NULL
  stopifnot(nrow(peaks) > 12L)
  summary <- data.frame(extremal_index = ei$EIintervals,
                        automatic_gap_days = dc$r, exceedance_days = length(index),
                        clusters_before_year_split = dc$nClusters,
                        retained_clusters = nrow(peaks),
                        time_axis = "Concatenated seasons; clusters split at year boundaries")
  list(peaks = peaks, summary = summary, extremal_index = ei, declustering = dc)
}

smith_design <- function(d) {
  list(mu = cbind(intercept = 1, gmst_centered = d$gmst_centered,
                   sin1 = d$sin1, cos1 = d$cos1, sin2 = d$sin2, cos2 = d$cos2),
       sigma = cbind(intercept = 1, sin1 = d$sin1, cos1 = d$cos1))
}

smith_parameters <- function(coef, d) {
  design <- smith_design(d)
  list(mu = drop(design$mu %*% coef[1:6]),
       sigma = drop(design$sigma %*% coef[7:9]), xi = unname(coef[10]))
}

smith_log_tail <- function(x, mu, sigma, xi) {
  stopifnot(all(is.finite(sigma)), all(sigma > 0), is.finite(xi))
  z <- (x - mu) / sigma
  if (abs(xi) < 1e-7) return(-z)
  h <- 1 + xi * z
  out <- rep(if (xi < 0) -Inf else Inf, length(h))
  ok <- h > 0
  out[ok] <- -log(h[ok]) / xi
  out
}

smith_log_sum <- function(x) {
  if (anyNA(x)) return(NA_real_)
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}

smith_nll <- function(coef, d, threshold, peak_rows, npy) {
  # Same Poisson likelihood as pp.fit, with daily exposure 1/npy and a stable
  # xi=0 limit. No extremal-index multiplier is applied after declustering.
  if (any(!is.finite(coef))) return(1e100)
  p <- smith_parameters(coef, d)
  if (any(!is.finite(p$sigma)) || any(p$sigma <= 0)) return(1e100)
  hu <- 1 + p$xi * (threshold - p$mu) / p$sigma
  if (any(!is.finite(hu)) || any(hu <= 0)) return(1e100)
  log_a <- smith_log_tail(threshold, p$mu, p$sigma, p$xi)
  z <- (d$tx[peak_rows] - p$mu[peak_rows]) / p$sigma[peak_rows]
  if (abs(p$xi) < 1e-7) {
    log_intensity <- -log(p$sigma[peak_rows]) - z
  } else {
    h <- 1 + p$xi * z
    if (any(!is.finite(h)) || any(h <= 0)) return(1e100)
    log_intensity <- -log(p$sigma[peak_rows]) - (1 + 1 / p$xi) * log(h)
  }
  value <- sum(exp(log_a)) / npy - sum(log_intensity)
  if (is.finite(value)) value else 1e100
}

smith_fit_pp <- function(d, threshold, peaks, npy, maxit, shapes) {
  covariates <- as.matrix(d[, c("gmst_centered", "sin1", "cos1", "sin2", "cos2")])
  # Keep every exposure day. Replacing nonpeak observations by their threshold
  # censors only the event term. Passing peaks alone would discard exposure and
  # corrupt event rates; keeping raw nonpeaks would incorrectly include them.
  xdat <- threshold
  xdat[peaks$row] <- d$tx[peaks$row]
  selected <- which(xdat > threshold)
  stopifnot(identical(selected, sort(as.integer(peaks$row))))
  design <- smith_design(d)
  u_coef <- as.numeric(qr.solve(design$mu, threshold))
  scale0 <- max(stats::sd(peaks$tx), 0.5)
  annualized_rate <- nrow(peaks) / (nrow(d) / npy)
  starts <- lapply(shapes, function(shape) {
    shift <- if (abs(shape) < 1e-7) scale0 * log(annualized_rate) else
      -scale0 * expm1(-shape * log(annualized_rate)) / shape
    mu0 <- u_coef
    mu0[1] <- mu0[1] + shift
    tryCatch(ismev::pp.fit(
      xdat = xdat, threshold = threshold, npy = npy, ydat = covariates,
      mul = 1:5, sigl = 2:3, shl = NULL,
      mulink = identity, siglink = identity, shlink = identity,
      muinit = mu0, siginit = c(scale0, 0, 0), shinit = shape,
      show = FALSE, method = "Nelder-Mead", maxit = maxit, reltol = 1e-10),
      error = function(e) list(conv = NA_integer_, nllh = Inf, error = conditionMessage(e)))
  })
  audit <- do.call(rbind, lapply(seq_along(starts), function(i) {
    z <- starts[[i]]
    data.frame(start_shape = shapes[i], convergence = z$conv, nll = z$nllh,
               error = if (is.null(z$error)) "" else z$error)
  }))
  good <- which(vapply(starts, function(z)
    !is.na(z$conv) && z$conv == 0L && is.finite(z$nllh) &&
      z$nllh < 1e6 && !is.null(z$mle) && all(is.finite(z$mle)), logical(1)))
  if (!length(good)) {
    print(audit)
    stop("No pp.fit start converged. Inspect the audit; no alternative model was substituted.")
  }
  best_index <- good[which.min(vapply(starts[good], function(z) z$nllh, 0.0))]
  best <- starts[[best_index]]
  checked <- smith_nll(best$mle, d, threshold, peaks$row, npy)
  stopifnot(is.finite(checked), abs(checked - best$nllh) < 1e-4 * (1 + abs(checked)))
  names(best$mle) <- c(paste0("mu_", colnames(design$mu)),
                       paste0("sigma_", colnames(design$sigma)), "xi")
  if (best$mle[10] <= -0.5)
    warning("Shape <= -0.5: the MLE is retained, but ordinary Hessian intervals are not reported.")
  list(fit = best, audit = audit, coef = best$mle, npy = npy, checked_nll = checked)
}

smith_expected_count <- function(coef, d, temperature, npy) {
  p <- smith_parameters(coef, d)
  log_rate <- smith_log_sum(smith_log_tail(temperature, p$mu, p$sigma, p$xi)) - log(npy)
  c(log_rate = log_rate, expected_count = exp(log_rate))
}

smith_count_check <- function(model, d, threshold, peaks, temperatures) {
  do.call(rbind, lapply(temperatures, function(x) {
    expected <- smith_expected_count(model$coef, d, x, model$npy)["expected_count"]
    observed <- sum(peaks$tx > x)
    below <- sum(threshold > x)
    data.frame(temperature_C = x, observed_cluster_count = observed,
               expected_cluster_count = unname(expected),
               below_fitted_threshold_days = below,
               threshold_coverage = if (below == 0) "Within fitted tail" else "Includes below-threshold extrapolation",
               pearson_component = if (expected > 0 && is.finite(expected))
                 (observed - expected)^2 / expected else NA_real_)
  }))
}

smith_log_probability <- function(log_rate) {
  if (is.na(log_rate)) return(NA_real_)
  if (log_rate < -35) return(log_rate)
  if (log_rate > 35) return(0)
  log(-expm1(-exp(log_rate)))
}

smith_attribution <- function(model, factual, counterfactual, temperatures) {
  do.call(rbind, lapply(temperatures, function(x) {
    lf <- smith_expected_count(model$coef, factual, x, model$npy)["log_rate"]
    lc <- smith_expected_count(model$coef, counterfactual, x, model$npy)["log_rate"]
    pf <- smith_log_probability(lf)
    pc <- smith_log_probability(lc)
    both_zero <- is.infinite(pf) && pf < 0 && is.infinite(pc) && pc < 0
    ratio <- if (both_zero) NA_real_ else pf - pc
    data.frame(temperature_C = x, expected_events_factual = exp(lf),
               expected_events_counterfactual = exp(lc),
               p_factual = exp(pf), p_counterfactual = exp(pc),
               return_period_factual_years = exp(-pf),
               return_period_counterfactual_years = exp(-pc),
               log_RR = ratio, RR = exp(ratio), FAR = -expm1(-ratio),
               inverse_RR = exp(-ratio),
               ratio_status = if (both_zero) "Undefined: both probabilities are zero" else
                 if (pc == -Inf) "Infinite RR: counterfactual probability is zero" else
                   if (pf == -Inf) "Zero RR: factual probability is zero" else "Defined")
  }))
}

smith_daily_parameters <- function(model, d, scenario) {
  p <- smith_parameters(model$coef, d)
  data.frame(date = d$date, season_day = d$season_day, scenario = scenario,
             gmst = d$gmst, mu_C = p$mu, sigma_C = p$sigma, xi = p$xi,
             upper_endpoint_C = if (p$xi < 0) p$mu - p$sigma / p$xi else Inf)
}



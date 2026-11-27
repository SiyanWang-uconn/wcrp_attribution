# Smith-inspired seasonal location margins and a logistic extreme-value copula.
# Definitions only. Source bivariate_pot_gmst.R first for numerical utilities.
# Censored window working likelihood, not a declustered Poisson likelihood.

bs_design <- function(date, gmst, g0, location_order, scale_order) {
  date <- as.Date(date)
  year <- as.integer(format(date, "%Y"))
  leap <- year %% 4 == 0 & (year %% 100 != 0 | year %% 400 == 0)
  day <- as.integer(format(date, "%j")) -
    as.integer(leap & as.integer(format(date, "%m")) > 2)
  phase <- 2 * pi * (day - 1) / 365
  harmonics <- function(order) {
    z <- matrix(numeric(), nrow = length(date), ncol = 0)
    if (order > 0) for (k in seq_len(order)) {
      z <- cbind(z, sin(k * phase), cos(k * phase))
      colnames(z)[(ncol(z) - 1):ncol(z)] <- paste0(c("sin", "cos"), k)
    }
    z
  }
  list(mu = cbind(intercept = 1, gmst = gmst - g0,
                  harmonics(location_order)), season_scale = harmonics(scale_order))
}

bs_parameters <- function(par, design, station) {
  mu_keys <- paste0(station, "_mu_", colnames(design$mu))
  mu <- drop(design$mu %*% par[mu_keys])
  sigma <- rep(exp(par[paste0(station, "_logsigma")]), length(mu))
  if (ncol(design$season_scale)) sigma <- sigma + drop(
    design$season_scale %*%
      par[paste0(station, "_sigma_", colnames(design$season_scale))])
  list(mu = mu, sigma = sigma, xi = unname(par[paste0(station, "_xi")]))
}

bs_margin <- function(x, p, npy = 365) {
  # F_window(x) = exp(-t(x)/npy); t is the Smith PP tail measure.
  # npy is a parameterization convention, not an independent-window count.
  z <- (x - p$mu) / p$sigma
  h <- 1 + p$xi * z
  valid <- is.finite(h) & h > 0 & is.finite(p$sigma) & p$sigma > 0
  loga <- rep(NA_real_, length(z))
  if (abs(p$xi) < 1e-7) loga[valid] <- -z[valid] - log(npy)
  else loga[valid] <- -log1p(p$xi * z[valid]) / p$xi - log(npy)
  # Above a finite upper endpoint the CDF is one and density zero.
  upper <- !valid & p$xi < 0 & is.finite(z) & h <= 0 & p$sigma > 0
  loga[upper] <- -Inf
  logF <- -exp(loga)
  logsurvival <- rep(NA_real_, length(logF))
  defined <- !is.na(logF)
  logsurvival[defined] <- bp_log1mexp(logF[defined])
  logdensity <- rep(-Inf, length(z))
  logdensity[valid] <- logF[valid] + loga[valid] - log(p$sigma[valid]) -
    if (abs(p$xi) < 1e-7) 0 else log(h[valid])
  list(loga = loga, logF = logF, logsurvival = logsurvival,
       logdensity = logdensity, valid = valid)
}

bs_loglik <- function(par, data, threshold, design, npy) {
  if (any(!is.finite(par))) return(-Inf)
  alpha <- plogis(unname(par["alpha_logit"]))
  if (alpha < 1e-5 || alpha >= 1) return(-Inf)
  margins <- lapply(1:2, function(j) {
    p <- bs_parameters(par, design, c("A", "B")[j])
    if (any(!is.finite(p$sigma)) || any(p$sigma <= 0)) return(NULL)
    bs_margin(pmax(data[[c("X", "Y")[j]]], threshold[j]), p, npy)
  })
  if (any(vapply(margins, is.null, logical(1)))) return(-Inf)
  a <- margins[[1]]; b <- margins[[2]]
  if (any(!a$valid) || any(!b$valid) ||
      any(!is.finite(c(a$loga, b$loga)))) return(-Inf)
  cop <- bp_logistic(a$loga, b$loga, alpha)
  ia <- data$X > threshold[1]; ib <- data$Y > threshold[2]
  ll <- -cop$V
  ll[ia] <- ll[ia] + (cop$logVa + a$logdensity - a$logF)[ia]
  ll[ib] <- ll[ib] + (cop$logVb + b$logdensity - b$logF)[ib]
  both <- ia & ib
  ll[both] <- ll[both] + bp_softplus(log(1 / alpha - 1) - cop$logV[both])
  if (any(!is.finite(ll))) -Inf else sum(ll)
}

fit_bivariate_smith <- function(data, threshold, location_order = 2L,
                                scale_order = 1L, npy = 365,
                                maxit = 3000L) {
  stopifnot(all(c("X", "Y", "gmst", "date", "year") %in% names(data)),
    nrow(data) >= 50, length(threshold) == 2,
    all(is.finite(threshold)), all(is.finite(as.matrix(data[c("X", "Y", "gmst")]))),
    !anyNA(as.Date(data$date)), !anyDuplicated(as.Date(data$date)),
    sd(data$gmst) > 0, location_order %in% 0:2, scale_order %in% 0:1,
    length(npy) == 1, is.finite(npy), npy > 0)
  g0 <- mean(tapply(data$gmst, data$year, mean))
  design <- bs_design(data$date, data$gmst, g0, location_order, scale_order)
  par <- numeric()
  for (j in 1:2) {
    station <- c("A", "B")[j]
    x <- data[[c("X", "Y")[j]]]
    excess <- x[x > threshold[j]] - threshold[j]
    if (length(excess) < 20 || length(excess) == length(x))
      stop("Each margin needs at least 20 exceedances and some nonexceedances.")
    sigma0 <- max(mean(excess), 0.2)
    zeta <- length(excess) / length(x)
    mu <- setNames(rep(0, ncol(design$mu)),
                   paste0(station, "_mu_", colnames(design$mu)))
    mu[1] <- threshold[j] + sigma0 * log(npy * -log1p(-zeta))
    scale_season <- setNames(rep(0, ncol(design$season_scale)),
      paste0(station, "_sigma_", colnames(design$season_scale)))
    par <- c(par, mu, setNames(log(sigma0), paste0(station, "_logsigma")),
      scale_season, setNames(0, paste0(station, "_xi")))
  }
  par <- c(par, alpha_logit = qlogis(0.65))
  objective <- function(theta) {
    ll <- bs_loglik(theta, data, threshold, design, npy)
    if (is.finite(ll)) -ll else 1e100
  }
  starts <- c(0.35, 0.65, 0.9)
  candidates <- lapply(starts, function(alpha) {
    initial <- par; initial["alpha_logit"] <- qlogis(alpha)
    tryCatch(optim(initial, objective, method = "BFGS",
      control = list(maxit = maxit, reltol = 1e-9)), error = function(e)
        list(convergence = 999L, value = Inf, message = conditionMessage(e)))
  })
  ok <- vapply(candidates, function(z)
    z$convergence == 0 && is.finite(z$value) && z$value < 1e90, logical(1))
  if (!any(ok)) stop("No start converged. Inspect design, support and thresholds.")
  best <- candidates[ok][[which.min(vapply(candidates[ok], function(z) z$value, 0.0))]]
  list(coefficients = best$par, logLik = -best$value, threshold = threshold,
    g0 = g0, location_order = location_order, scale_order = scale_order,
    npy = npy, gmst_range = range(data$gmst), candidates = candidates,
    n = nrow(data), inference = "Censored window working likelihood; no uncertainty fitted")
}

predict_bivariate_smith <- function(fit, gmst, date, target) {
  stopifnot(length(target) == 2, all(is.finite(target)),
    all(target >= fit$threshold), length(gmst) > 0, all(is.finite(gmst)),
    length(date) %in% c(1L, length(gmst)), !anyNA(as.Date(date)))
  date <- rep(as.Date(date), length.out = length(gmst))
  design <- bs_design(date, gmst, fit$g0, fit$location_order, fit$scale_order)
  margins <- lapply(1:2, function(j) {
    p <- bs_parameters(fit$coefficients, design, c("A", "B")[j])
    if (any(!is.finite(p$sigma)) || any(p$sigma <= 0))
      stop("Nonpositive predicted scale; inspect seasonal extrapolation.")
    bs_margin(rep(target[j], length(gmst)), p, fit$npy)
  })
  a <- margins[[1]]; b <- margins[[2]]
  if (anyNA(c(a$loga, b$loga))) stop("Target below modeled lower endpoint.")
  alpha <- plogis(unname(fit$coefficients["alpha_logit"]))
  logjoint <- rep(-Inf, length(gmst))
  ok <- is.finite(a$loga) & is.finite(b$loga)
  if (any(ok)) {
    cop <- bp_logistic(a$loga[ok], b$loga[ok], alpha)
    logs <- bp_logadd(a$loga[ok], b$loga[ok])
    logdelta <- logs + bp_log1mexp(pmin(cop$logV - logs, 0))
    logextra <- -exp(logs) + bp_log_expm1_from_log(logdelta)
    logjoint[ok] <- bp_logadd(a$logsurvival[ok] + b$logsurvival[ok], logextra)
    logjoint[ok] <- pmin(logjoint[ok], a$logsurvival[ok], b$logsurvival[ok])
  }
  logba <- logab <- rep(NA_real_, length(gmst))
  ka <- is.finite(a$logsurvival); kb <- is.finite(b$logsurvival)
  logba[ka] <- logjoint[ka] - a$logsurvival[ka]
  logab[kb] <- logjoint[kb] - b$logsurvival[kb]
  data.frame(date = date, gmst = gmst, target_A_C = target[1], target_B_C = target[2],
    p_A = exp(a$logsurvival), p_B = exp(b$logsurvival), p_joint = exp(logjoint),
    p_B_given_A = exp(logba), p_A_given_B = exp(logab),
    PMF = exp(logjoint - a$logsurvival - b$logsurvival),
    log_p_A = a$logsurvival, log_p_B = b$logsurvival, log_p_joint = logjoint,
    log_p_B_given_A = logba, log_p_A_given_B = logab,
    extrapolated_gmst = gmst < fit$gmst_range[1] | gmst > fit$gmst_range[2])
}

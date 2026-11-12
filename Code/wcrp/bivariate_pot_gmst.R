# Nonstationary bivariate POT with GPD margins and a logistic extreme-value copula.
# This is a censored upper-tail model, not a conditional multivariate GP density.
# Source this file in RStudio; sourcing creates no files and fits no station data.
# With overlapping windows, the objective is a working/composite likelihood.
# Its optimizer Hessian must not be used as an independence-based covariance.

# ---- Interactive configuration and example (not run automatically) ----
PROJECT_DIR <- "C:/Users/siyan/WCRP/wcrp_common_attribution_assessment"
PAIRED_FILE <- file.path(PROJECT_DIR, "Results/Seattle_Princeton_bivariate",
                         "paired_3day_JJA_windows.csv")
GMST_FILE <- file.path(PROJECT_DIR, "data/processed/gmst_smoothed.csv")
POT_THRESHOLDS <- c(A = 29.80, B = 33.3683333333333)
TARGET_TEMPERATURES <- c(A = 40.3666666666667, B = 41.00)
RATE_GMST <- c(TRUE, TRUE)
SCALE_GMST <- c(FALSE, FALSE)
# Example after sourcing this file:
# paired <- read.csv(PAIRED_FILE)
# gmst <- read.csv(GMST_FILE)
# d <- prepare_bivariate_pot(paired, gmst)
# fit <- fit_bivariate_pot(d, POT_THRESHOLDS, RATE_GMST, SCALE_GMST)
# fit$coefficients
# predict_bivariate_pot(fit, gmst = 1.0, target = TARGET_TEMPERATURES)
# compare_bivariate_pot(fit, TARGET_TEMPERATURES, factual_gmst, counterfactual_gmst)
# Reference GMST values must be selected explicitly; no historical year is
# automatically treated as a world without anthropogenic forcing.

# ---- Numerical helpers ----
bp_softplus <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))
bp_logadd <- function(a, b) {
  m <- pmax(a, b)
  out <- m + log1p(exp(-abs(a - b)))
  out[is.infinite(m) & m < 0] <- -Inf
  out
}
bp_log1mexp <- function(x) {
  if (any(x > 1e-12, na.rm = TRUE)) stop("log1mexp requires x <= 0.")
  x <- pmin(x, 0)
  out <- x
  k <- x < -log(2)
  out[k] <- log1p(-exp(x[k]))
  out[!k] <- log(-expm1(x[!k]))
  out
}
bp_log_expm1_from_log <- function(logx) {
  out <- logx
  k <- is.finite(logx) & logx >= -20
  x <- exp(logx[k])
  out[k] <- ifelse(x > 50, x + log1p(-exp(-x)), log(expm1(pmin(x, 50))))
  out
}
bp_switch <- function(x, label) {
  if (!is.logical(x) || !length(x) %in% c(1L, 2L) || anyNA(x))
    stop(label, " must be one or two logical values.")
  rep(x, length.out = 2L)
}

# ---- Explicit pairing and GMST preparation ----
prepare_bivariate_pot <- function(paired, gmst) {
  if (!all(c("year", "X_C", "Y_C") %in% names(paired)) ||
      !all(c("year", "gmst") %in% names(gmst))) stop("Required columns are missing.")
  if (anyNA(gmst$year) || anyDuplicated(gmst$year)) stop("GMST years must be unique.")
  d <- data.frame(year = paired$year, X = paired$X_C, Y = paired$Y_C,
                  gmst = gmst$gmst[match(paired$year, gmst$year)])
  if ("date" %in% names(paired)) d$date <- as.Date(paired$date)
  valid <- is.finite(d$X) & is.finite(d$Y)
  if (any(!is.finite(d$gmst[valid]))) stop("GMST is missing for valid station pairs.")
  message("Excluded ", sum(!valid), " incomplete station windows; no infilling.")
  d[valid, , drop = FALSE]
}

# ---- Parameter layout; slopes are per degree of centered GMST ----
bp_spec <- function(rate_gmst, scale_gmst, alpha) {
  rate_gmst <- bp_switch(rate_gmst, "rate_gmst")
  scale_gmst <- bp_switch(scale_gmst, "scale_gmst")
  if (!is.null(alpha) && (length(alpha) != 1L || !is.finite(alpha) ||
                          alpha <= 0 || alpha > 1)) stop("Fixed alpha must be in (0, 1].")
  keys <- character()
  for (j in c("A", "B")) {
    i <- match(j, c("A", "B"))
    keys <- c(keys, paste0("rate_", j, "_intercept"),
      if (rate_gmst[i]) paste0("rate_", j, "_gmst"),
      paste0("logscale_", j, "_intercept"),
      if (scale_gmst[i]) paste0("logscale_", j, "_gmst"), paste0("shape_", j))
  }
  if (is.null(alpha)) keys <- c(keys, "alpha_logit")
  list(keys = keys, rate_gmst = rate_gmst, scale_gmst = scale_gmst, alpha = alpha)
}
bp_parameters <- function(par, gmst, g0, spec) {
  g <- gmst - g0
  z <- lapply(c("A", "B"), function(j) {
    i <- match(j, c("A", "B"))
    eta <- rep(unname(par[paste0("rate_", j, "_intercept")]), length(g))
    ls <- rep(unname(par[paste0("logscale_", j, "_intercept")]), length(g))
    if (spec$rate_gmst[i]) eta <- eta + par[paste0("rate_", j, "_gmst")] * g
    if (spec$scale_gmst[i]) ls <- ls + par[paste0("logscale_", j, "_gmst")] * g
    list(logzeta = -bp_softplus(-eta), logscale = ls,
         shape = unname(par[paste0("shape_", j)]))
  })
  list(A = z[[1]], B = z[[2]],
       alpha = if (is.null(spec$alpha)) plogis(unname(par["alpha_logit"])) else spec$alpha)
}

# ---- Tail margins; defined at or above the fixed modeling threshold ----
bp_margin <- function(x, u, p) {
  q <- (x - u) * exp(-p$logscale)
  support <- 1 + p$shape * q
  valid <- is.finite(q) & is.finite(support) & support > 0
  loggp <- rep(-Inf, length(q))
  if (p$shape == 0) loggp[valid] <- -q[valid]
  else loggp[valid] <- -log1p(p$shape * q[valid]) / p$shape
  ls <- p$logzeta + loggp
  lf <- p$logzeta - p$logscale + (1 + p$shape) * loggp
  lf[!valid] <- -Inf
  logF <- bp_log1mexp(ls)
  # a = -log(F); use its asymptotic value when survival is extremely small.
  loga <- log(-logF)
  tiny <- ls < -36
  loga[tiny] <- ls[tiny]
  list(logsurvival = ls, logdensity = lf, logF = logF,
       loga = loga, support = valid)
}
bp_logistic <- function(loga, logb, alpha) {
  logV <- alpha * bp_logadd(loga / alpha, logb / alpha)
  rminus1 <- 1 / alpha - 1
  logVa <- rminus1 * (loga - logV)
  logVb <- rminus1 * (logb - logV)
  if (alpha == 1) {
    logVa[] <- 0
    logVb[] <- 0
    correction <- rep(0, length(logV))
  } else correction <- bp_softplus(log(rminus1) - logV)
  list(logV = logV, V = exp(logV), logVa = logVa, logVb = logVb,
       log_density_factor = logVa + logVb + correction)
}

# ---- Four-state censored log likelihood ----
bivariate_pot_loglik <- function(par, data, threshold, g0, spec, individual = FALSE) {
  if (length(par) != length(spec$keys) || any(!is.finite(par))) return(-Inf)
  names(par) <- spec$keys
  p <- bp_parameters(par, data$gmst, g0, spec)
  if (!is.finite(p$alpha) || p$alpha <= 1e-6 || p$alpha > 1) return(-Inf)
  ia <- data$X > threshold[1]
  ib <- data$Y > threshold[2]
  a <- bp_margin(pmax(data$X, threshold[1]), threshold[1], p$A)
  b <- bp_margin(pmax(data$Y, threshold[2]), threshold[2], p$B)
  if (any(!a$support[ia]) || any(!b$support[ib]) ||
      any(!is.finite(a$loga)) || any(!is.finite(b$loga))) return(-Inf)
  c <- bp_logistic(a$loga, b$loga, p$alpha)
  ll <- -c$V
  ll[ia] <- ll[ia] + (c$logVa + a$logdensity - a$logF)[ia]
  ll[ib] <- ll[ib] + (c$logVb + b$logdensity - b$logF)[ib]
  both <- ia & ib
  if (p$alpha < 1) ll[both] <- ll[both] +
    bp_softplus(log(1 / p$alpha - 1) - c$logV[both])
  if (any(!is.finite(ll))) return(-Inf)
  if (individual) ll else sum(ll)
}

# ---- Maximum working-likelihood fit ----
fit_bivariate_pot <- function(data, threshold, rate_gmst = TRUE, scale_gmst = FALSE,
                              alpha = NULL, g0 = NULL, start = NULL,
                              control = list(maxit = 1500, reltol = 1e-9)) {
  if (!all(c("X", "Y", "gmst") %in% names(data))) stop("Use columns X, Y and gmst.")
  if (!all(vapply(data[, c("X", "Y", "gmst")], is.numeric, logical(1))) ||
      any(!is.finite(as.matrix(data[, c("X", "Y", "gmst")]))))
    stop("Supply complete finite numeric station pairs and GMST; NA is not censoring.")
  if (nrow(data) < 20L || length(threshold) != 2L || any(!is.finite(threshold)))
    stop("Need at least 20 pairs and two finite thresholds.")
  spec <- bp_spec(rate_gmst, scale_gmst, alpha)
  if (any(c(spec$rate_gmst, spec$scale_gmst)) && sd(data$gmst) <= 0)
    stop("GMST must vary when slopes are estimated.")
  if (is.null(g0)) g0 <- mean(data$gmst)
  if (length(g0) != 1L || !is.finite(g0)) stop("g0 must be finite and scalar.")
  par <- setNames(rep(0, length(spec$keys)), spec$keys)
  counts <- integer(2)
  for (i in 1:2) {
    j <- c("A", "B")[i]
    x <- data[[c("X", "Y")[i]]]
    excess <- x[x > threshold[i]] - threshold[i]
    counts[i] <- length(excess)
    if (counts[i] < 5L || counts[i] == nrow(data))
      stop("Each station needs at least five exceedances and some nonexceedances.")
    par[paste0("rate_", j, "_intercept")] <- qlogis(counts[i] / nrow(data))
    par[paste0("logscale_", j, "_intercept")] <- log(mean(excess))
  }
  if (!is.null(start)) {
    if (is.null(names(start)) || !setequal(names(start), spec$keys))
      stop("start must have exactly the named parameters in spec$keys.")
    par <- start[spec$keys]
  }
  nll <- function(theta) {
    ll <- bivariate_pot_loglik(theta, data, threshold, g0, spec)
    if (is.finite(ll)) -ll else 1e100
  }
  starts <- if (is.null(alpha) && is.null(start)) c(0.35, 0.65, 0.90) else NA_real_
  candidates <- lapply(starts, function(value) {
    initial <- par
    if (is.finite(value)) initial["alpha_logit"] <- qlogis(value)
    tryCatch(optim(initial, nll, method = "BFGS", control = control),
             error = function(e) NULL)
  })
  good <- vapply(candidates, function(z) !is.null(z) && z$convergence == 0 &&
    is.finite(z$value) && z$value < 1e90, logical(1))
  if (!any(good)) stop("No optimizer start converged; inspect thresholds, support and starts.")
  eligible <- candidates[good]
  best <- eligible[[which.min(vapply(eligible, function(z) z$value, numeric(1)))]]
  names(best$par) <- spec$keys
  structure(list(coefficients = best$par, threshold = unname(threshold), g0 = g0,
    spec = spec, logLik = -best$value, convergence = best$convergence,
    n = nrow(data), exceedances = setNames(counts, c("A", "B")),
    gmst_range = range(data$gmst), optimizer_starts = candidates,
    inference = "Working likelihood; use whole-year refitting bootstrap for dependent windows."),
    class = "bivariate_pot_gmst")
}

# ---- Stable joint and conditional window probabilities ----
predict_bivariate_pot <- function(fit, gmst, target) {
  if (!inherits(fit, "bivariate_pot_gmst")) stop("Supply a fitted bivariate POT model.")
  if (!length(gmst) || any(!is.finite(gmst)) || length(target) != 2L ||
      any(!is.finite(target)) || any(target < fit$threshold))
    stop("Supply finite GMST and two target values at or above modeling thresholds.")
  p <- bp_parameters(fit$coefficients, gmst, fit$g0, fit$spec)
  a <- bp_margin(rep(target[1], length(gmst)), fit$threshold[1], p$A)
  b <- bp_margin(rep(target[2], length(gmst)), fit$threshold[2], p$B)
  logjoint <- rep(-Inf, length(gmst))
  ok <- is.finite(a$logsurvival) & is.finite(b$logsurvival)
  if (any(ok)) {
    la <- a$loga[ok]
    lb <- b$loga[ok]
    cop <- bp_logistic(la, lb, p$alpha)
    logs <- bp_logadd(la, lb)
    # joint = survival_A * survival_B + H - F_A * F_B.
    # The second term is nonnegative for the logistic model.
    if (p$alpha == 1) logextra <- rep(-Inf, sum(ok)) else {
      logdelta <- logs + bp_log1mexp(pmin(cop$logV - logs, 0))
      logextra <- -exp(logs) + bp_log_expm1_from_log(logdelta)
    }
    logjoint[ok] <- bp_logadd(a$logsurvival[ok] + b$logsurvival[ok], logextra)
    logjoint[ok] <- pmin(logjoint[ok], a$logsurvival[ok], b$logsurvival[ok])
  }
  cb <- ca <- rep(NA_real_, length(gmst))
  ka <- is.finite(a$logsurvival)
  kb <- is.finite(b$logsurvival)
  cb[ka] <- exp(logjoint[ka] - a$logsurvival[ka])
  ca[kb] <- exp(logjoint[kb] - b$logsurvival[kb])
  data.frame(gmst = gmst, target_A = target[1], target_B = target[2],
    p_A = exp(a$logsurvival), p_B = exp(b$logsurvival), p_joint = exp(logjoint),
    p_B_given_A = cb, p_A_given_B = ca,
    log_p_A = a$logsurvival, log_p_B = b$logsurvival, log_p_joint = logjoint,
    log_p_B_given_A = ifelse(ka, logjoint - a$logsurvival, NA_real_),
    log_p_A_given_B = ifelse(kb, logjoint - b$logsurvival, NA_real_),
    extrapolated_gmst = gmst < fit$gmst_range[1] | gmst > fit$gmst_range[2])
}

# ---- Probability ratios for explicitly chosen climate states ----
compare_bivariate_pot <- function(fit, target, factual_gmst, counterfactual_gmst) {
  if (length(factual_gmst) != 1L || length(counterfactual_gmst) != 1L)
    stop("Choose one GMST value for each climate state.")
  z <- predict_bivariate_pot(fit, c(factual_gmst, counterfactual_gmst), target)
  keys <- c("p_A", "p_B", "p_joint", "p_B_given_A", "p_A_given_B")
  lf <- unlist(z[1, paste0("log_", keys)], use.names = FALSE)
  lc <- unlist(z[2, paste0("log_", keys)], use.names = FALSE)
  data.frame(quantity = keys, factual_gmst = factual_gmst,
    counterfactual_gmst = counterfactual_gmst,
    p_factual = unlist(z[1, keys], use.names = FALSE),
    p_counterfactual = unlist(z[2, keys], use.names = FALSE),
    log_PR = lf - lc, PR = exp(lf - lc),
    extrapolated_gmst = any(z$extrapolated_gmst))
}

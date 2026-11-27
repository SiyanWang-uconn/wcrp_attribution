# Smith Section 0.6.1 GEV helpers. Definitions only; no analysis runs on sourcing.

smith_gev_annual <- function(daily, gmst, first_year, last_year) {
  stopifnot(!anyDuplicated(gmst$year))
  do.call(rbind, lapply(seq.int(first_year, last_year), function(y) {
    z <- daily[daily$year == y, ]
    expected <- as.integer(as.Date(sprintf("%d-01-01", y + 1L)) - as.Date(sprintf("%d-01-01", y)))
    g <- gmst$gmst[match(y, gmst$year)]
    valid <- is.finite(z$tx)
    data.frame(year = y, expected_days = expected, valid_days = sum(valid),
      txx = if (any(valid)) max(z$tx[valid]) else NA_real_,
      peak_date = if (any(valid)) z$date[which.max(replace(z$tx, !valid, -Inf))] else as.Date(NA),
      gmst = g, used_in_fit = nrow(z) == expected && all(valid) && is.finite(g))
  }))
}

smith_gev_nll <- function(b, data, g0) {
  if (length(b) != 4L || any(!is.finite(b)) || b[3] <= 0) return(Inf)
  z <- (data$txx - b[1] - b[2] * (data$gmst - g0)) / b[3]
  xi <- b[4]
  if (any(!is.finite(z))) return(Inf)
  if (abs(xi) < 1e-7) value <- sum(log(b[3]) + z + exp(-z)) else {
    h <- 1 + xi * z
    if (any(h <= 0)) return(Inf)
    value <- sum(log(b[3]) + (1 + 1 / xi) * log(h) + exp(-log(h) / xi))
  }
  if (is.finite(value)) value else Inf
}

smith_gev_fit <- function(data, shapes, maxit) {
  stopifnot(nrow(data) > 5L, all(is.finite(data$txx)), all(is.finite(data$gmst)))
  g0 <- mean(data$gmst)
  gc <- data$gmst - g0
  regression <- lm(data$txx ~ gc)
  scale0 <- max(stats::sd(stats::residuals(regression)) * sqrt(6) / pi, 0.1)
  mu0 <- c(unname(stats::coef(regression)[1]) - 0.5772157 * scale0,
            unname(stats::coef(regression)[2]))
  starts <- lapply(shapes, function(shape) {
    tryCatch(ismev::gev.fit(data$txx, ydat = matrix(gc, ncol = 1), mul = 1,
      sigl = NULL, shl = NULL, muinit = mu0, siginit = scale0, shinit = shape,
      show = FALSE, maxit = maxit, method = "Nelder-Mead", reltol = 1e-10),
      error = function(e) list(conv = NA_integer_, nllh = Inf, error = conditionMessage(e)))
  })
  audit <- do.call(rbind, lapply(seq_along(starts), function(i) {
    z <- starts[[i]]
    data.frame(initial_shape = shapes[i], convergence = z$conv, nll = z$nllh,
               error = if (is.null(z$error)) "" else z$error)
  }))
  good <- which(vapply(starts, function(z)
    !is.na(z$conv) && z$conv == 0 && is.finite(z$nllh) && z$nllh < 1e6 &&
      !is.null(z$mle) && all(is.finite(z$mle)), logical(1)))
  if (!length(good)) { print(audit); stop("No GEV MLE start converged.") }
  best_index <- good[which.min(vapply(starts[good], function(z) z$nllh, 0.0))]
  best <- starts[[best_index]]
  checked <- smith_gev_nll(best$mle, data, g0)
  stopifnot(is.finite(checked), abs(checked - best$nllh) < 1e-5 * (1 + abs(checked)))
  names(best$mle) <- c("mu_centered", "beta_gmst", "sigma", "xi")
  if (best$mle[4] <= -0.5) warning("Nonregular shape: interpret Hessian inference cautiously.")
  list(coef = best$mle, g0 = g0, data = data, fit = best, audit = audit, checked_nll = checked)
}

smith_gev_table_parameters <- function(model, case) {
  # Undo GMST centering for the paper's beta0 and its standard error.
  transform <- diag(4)
  transform[1, 2] <- -model$g0
  estimate <- drop(transform %*% model$coef)
  variance <- transform %*% model$fit$cov %*% t(transform)
  se <- ifelse(diag(variance) >= 0, sqrt(pmax(diag(variance), 0)), NA_real_)
  data.frame(case = case, parameter = c("beta0", "beta_gmst", "sigma", "xi", "NLLH"),
             estimate = c(estimate, model$checked_nll), standard_error = c(se, NA_real_),
             n_years = nrow(model$data))
}

smith_gev_paper_table <- function() {
  rbind(data.frame(case = "exclude_2021", parameter = c("beta0", "beta_gmst", "sigma", "xi", "NLLH"),
    paper_estimate = c(29.720, 1.758, 1.729, -0.469, 132.2),
    paper_standard_error = c(0.286, 0.548, 0.174, 0.085, NA_real_)),
    data.frame(case = "include_2021", parameter = c("beta0", "beta_gmst", "sigma", "xi", "NLLH"),
    paper_estimate = c(29.242, 2.595, 1.694, -0.128, 145.1),
    paper_standard_error = c(0.284, 0.641, 0.148, 0.055, NA_real_)))
}

smith_gev_log_hazard <- function(x, mu, sigma, xi) {
  n <- max(length(x), length(mu), length(sigma), length(xi))
  x <- rep(x, length.out = n)
  mu <- rep(mu, length.out = n)
  sigma <- rep(sigma, length.out = n)
  xi <- rep(xi, length.out = n)
  stopifnot(all(is.finite(sigma)), all(sigma > 0))
  z <- (x - mu) / sigma
  out <- -z
  nonzero <- abs(xi) >= 1e-7
  h <- 1 + xi * z
  good <- nonzero & h > 0
  out[good] <- -log(h[good]) / xi[good]
  out[nonzero & h <= 0 & xi < 0] <- -Inf
  out[nonzero & h <= 0 & xi > 0] <- Inf
  out
}

smith_gev_log_probability <- function(x, mu, sigma, xi) {
  log_a <- smith_gev_log_hazard(x, mu, sigma, xi)
  out <- log_a
  mid <- is.finite(log_a) & log_a >= -35 & log_a <= 35
  out[mid] <- log(-expm1(-exp(log_a[mid])))
  out[log_a > 35] <- 0
  out
}

smith_gev_endpoints <- function(model, gf, gc, case) {
  climates <- c(factual = gf, counterfactual = gc)
  do.call(rbind, lapply(names(climates), function(scenario) {
    g <- climates[[scenario]]
    b <- model$coef
    endpoint <- if (b[4] < 0) b[1] + b[2] * (g - model$g0) - b[3] / b[4] else Inf
    gradient <- if (b[4] < 0) c(1, g - model$g0, -1 / b[4], b[3] / b[4]^2) else rep(NA_real_, 4)
    variance <- if (b[4] < 0) drop(t(gradient) %*% model$fit$cov %*% gradient) else NA_real_
    data.frame(case = case, scenario = scenario, gmst = g, upper_endpoint_C = endpoint,
               delta_standard_error_C = if (is.finite(variance) && variance >= 0) sqrt(variance) else NA_real_)
  }))
}

smith_gev_attribution <- function(model, gf, gc, temperatures, case) {
  b <- model$coef
  pf <- smith_gev_log_probability(temperatures, b[1] + b[2] * (gf - model$g0), b[3], b[4])
  pc <- smith_gev_log_probability(temperatures, b[1] + b[2] * (gc - model$g0), b[3], b[4])
  log_rr <- pf - pc
  both <- pf == -Inf & pc == -Inf
  log_rr[both] <- NA_real_
  data.frame(case = case, temperature_C = temperatures,
    p_factual = exp(pf), p_counterfactual = exp(pc), RR = exp(log_rr),
    inverse_RR = exp(-log_rr), FAR = -expm1(-log_rr),
    return_period_factual_years = exp(-pf), return_period_counterfactual_years = exp(-pc),
    intensity_change_C = b[2] * (gf - gc), both_probabilities_zero = both)
}

smith_gev_gof <- function(model, case) {
  b <- model$coef
  mu <- b[1] + b[2] * (model$data$gmst - model$g0)
  u <- exp(-exp(smith_gev_log_hazard(model$data$txx, mu, b[3], b[4])))
  n <- length(u)
  sorted <- sort(pmin(pmax(u, .Machine$double.eps), 1 - .Machine$double.eps))
  i <- seq_len(n)
  statistics <- data.frame(case = case,
    KS = max(c(i / n - sorted, sorted - (i - 1) / n)),
    Cramer_von_Mises = 1 / (12 * n) + sum((sorted - (2 * i - 1) / (2 * n))^2),
    Anderson_Darling = -n - sum((2 * i - 1) * (log(sorted) + log1p(-rev(sorted)))) / n,
    p_value_status = "Not calculated: fitted-parameter simulation calibration is required")
  qq <- data.frame(case = case, theoretical_uniform = stats::ppoints(n),
    observed_uniform = sorted, theoretical_gumbel = -log(-log(stats::ppoints(n))),
    observed_gumbel = -log(-log(sorted)))
  list(statistics = statistics, qq = qq)
}

smith_gev_coordinates <- function(model) {
  eta <- model$coef
  eta[3] <- log(eta[3])
  names(eta)[3] <- "log_sigma"
  jac <- diag(4)
  jac[3, 3] <- 1 / model$coef[3]
  cov <- jac %*% model$fit$cov %*% jac
  cov <- (cov + t(cov)) / 2
  e <- eigen(cov, symmetric = TRUE)
  stopifnot(all(is.finite(e$values)), max(e$values) > 0)
  repaired <- e$vectors %*% diag(pmax(e$values, max(e$values) * 1e-10)) %*% t(e$vectors)
  list(center = eta, transform = t(chol(repaired)))
}

smith_gev_target <- function(model) {
  # Stated GEV prior: flat beta0, betaG, log(sigma), xi; |xi|<1.
  # Centering the intercept has unit Jacobian. No positive-betaG restriction.
  compiler::cmpfun(function(eta) {
    if (length(eta) != 4L || any(!is.finite(eta)) || abs(eta[4]) >= 1) return(-Inf)
    raw <- eta
    raw[3] <- exp(raw[3])
    if (!is.finite(raw[3]) || raw[3] <= 0) return(-Inf)
    value <- smith_gev_nll(raw, model$data, model$g0)
    if (is.finite(value)) -value else -Inf
  })
}

smith_gev_diagnostics <- function(chains) {
  do.call(rbind, lapply(seq_len(4), function(j) {
    x <- do.call(cbind, lapply(chains, function(z) {
      v <- z$eta_draws[, j]
      if (j == 3L) exp(v) else v
    }))
    smith_am_diag(x, c("mu_centered", "beta_gmst", "sigma", "xi")[j])
  }))
}

smith_gev_posterior_predict <- function(chains, model, gf, gc, temperatures) {
  n <- nrow(chains[[1]]$eta_draws)
  nc <- length(chains)
  nt <- length(temperatures)
  lf <- lc <- array(NA_real_, c(n, nc, nt))
  endpoint_f <- endpoint_c <- slope <- matrix(NA_real_, n, nc)
  for (j in seq_len(nc)) {
    eta <- chains[[j]]$eta_draws
    sigma <- exp(eta[, 3])
    xi <- eta[, 4]
    mf <- eta[, 1] + eta[, 2] * (gf - model$g0)
    mc <- eta[, 1] + eta[, 2] * (gc - model$g0)
    slope[, j] <- eta[, 2]
    for (k in seq_len(nt)) {
      lf[, j, k] <- smith_gev_log_probability(temperatures[k], mf, sigma, xi)
      lc[, j, k] <- smith_gev_log_probability(temperatures[k], mc, sigma, xi)
    }
    endpoint_f[, j] <- endpoint_c[, j] <- Inf
    bounded <- xi < 0
    endpoint_f[bounded, j] <- mf[bounded] - sigma[bounded] / xi[bounded]
    endpoint_c[bounded, j] <- mc[bounded] - sigma[bounded] / xi[bounded]
  }
  list(temperatures = temperatures, log_p_factual = lf, log_p_counterfactual = lc,
       endpoint_factual = endpoint_f, endpoint_counterfactual = endpoint_c, gmst_slope = slope)
}

smith_gev_posterior_summary <- function(prediction) {
  rows <- diags <- list()
  nr <- nd <- 0L
  positive <- prediction$gmst_slope > 0
  for (k in seq_along(prediction$temperatures)) {
    lf <- prediction$log_p_factual[, , k, drop = FALSE]
    lc <- prediction$log_p_counterfactual[, , k, drop = FALSE]
    dim(lf) <- dim(lc) <- dim(positive)
    ratio <- smith_am_inverse_ratio(lf, lc)
    quantities <- list(p_factual = exp(lf), p_counterfactual = exp(lc),
      inverse_ratio_author_convention = ratio$value,
      supports_factual = 1 * (lf > -Inf), supports_counterfactual = 1 * (lc > -Inf))
    for (name in names(quantities)) {
      nd <- nd + 1L
      diags[[nd]] <- smith_am_diag(quantities[[name]], paste(prediction$temperatures[k], name, sep = ":"))
    }
    for (subset in c("all_draws", "positive_gmst")) {
      selected <- if (subset == "positive_gmst") positive else matrix(TRUE, nrow(lf), ncol(lf))
      r <- ratio$value[selected]
      q <- if (length(r)) stats::quantile(r, c(0.025, 0.975), names = FALSE, type = 1) else c(NA_real_, NA_real_)
      nr <- nr + 1L
      rows[[nr]] <- data.frame(temperature_C = prediction$temperatures[k], ratio_subset = subset,
        E_pF = mean(exp(lf)), E_pC = mean(exp(lc)),
        E_inverse_RR = if (length(r)) mean(r) else NA_real_, q0_025 = q[1], q0_975 = q[2],
        Pr_inverse_RR_lt_0_01 = if (length(r)) mean(r < 0.01) else NA_real_,
        posterior_Pr_positive_gmst = mean(positive), fraction_both_zero = mean(ratio$both_zero),
        n_support_factual = sum(lf > -Inf), n_support_counterfactual = sum(lc > -Inf),
        n_draws = length(lf), n_ratio_draws = length(r))
    }
  }
  list(summary = do.call(rbind, rows), diagnostics = do.call(rbind, diags))
}

smith_gev_endpoint_summary <- function(prediction, temperatures) {
  do.call(rbind, lapply(c("factual", "counterfactual"), function(scenario) {
    endpoint <- prediction[[paste0("endpoint_", scenario)]]
    q <- stats::quantile(endpoint, c(0.025, 0.5, 0.975), names = FALSE, type = 1)
    do.call(rbind, lapply(temperatures, function(x) {
      data.frame(scenario = scenario, temperature_C = x, q0_025 = q[1], median = q[2], q0_975 = q[3],
        posterior_Pr_endpoint_above_x = mean(endpoint > x),
        posterior_Pr_unbounded_endpoint = mean(is.infinite(endpoint)))
    }))
  }))
}

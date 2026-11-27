# Adaptive-Metropolis helpers for the existing Smith seasonal PP MLE.
# Sourcing defines functions only. No packages are attached, no data are fitted.

smith_am_coordinates <- function(mle) {
  raw <- mle$model$coef
  stopifnot(length(raw) == 10L, raw[7] > 0)
  eta <- raw
  eta[7] <- log(raw[7])
  names(eta)[7] <- "log_sigma_intercept"
  # MLE Hessian covariance initializes proposals only; it is NOT the posterior.
  jac <- diag(10)
  jac[7, 7] <- 1 / raw[7]
  cov <- jac %*% mle$model$fit$cov %*% jac
  cov <- (cov + t(cov)) / 2
  stopifnot(all(is.finite(cov)))
  eig <- eigen(cov, symmetric = TRUE)
  stopifnot(max(eig$values) > 0)
  floor <- max(eig$values) * 1e-10
  repaired <- eig$vectors %*% diag(pmax(eig$values, floor)) %*% t(eig$vectors)
  list(center = eta, transform = t(chol(repaired)),
       original_covariance_eigenvalues = eig$values)
}

smith_am_target <- function(mle, prior_measure) {
  stopifnot(prior_measure %in% c("flat_log_scale_intercept", "flat_scale_coefficients"))
  d <- mle$fit_data
  days <- sort(unique(d$season_day))
  years <- sort(unique(d$year))
  di <- match(d$season_day, days)
  yi <- match(d$year, years)
  one_day <- d[match(days, d$season_day), ]
  g <- d$gmst_centered[match(years, d$year)]
  threshold <- matrix(NA_real_, length(days), length(years))
  threshold[cbind(di, yi)] <- d$threshold_C
  # All fitted seasons must be complete; no interpolation or exposure loss.
  stopifnot(all(is.finite(threshold)))
  peak <- mle$peaks$row
  peak_day <- di[peak]
  peak_year <- yi[peak]
  peak_tx <- d$tx[peak]
  sin1 <- one_day$sin1
  cos1 <- one_day$cos1
  sin2 <- one_day$sin2
  cos2 <- one_day$cos2
  npy <- mle$model$npy

  # Prior is flat in (six mu coefficients, log(sigma0), sigma_sin, sigma_cos, xi),
  # restricted to |xi|<1 and positive scales/support in the fitted season.
  # This is an explicit PP reconstruction, NOT a fully specified prior in Smith.
  # Alternative flat raw scale coefficients requires the exp(eta7) Jacobian.
  target <- function(eta) {
    if (length(eta) != 10L || any(!is.finite(eta)) || abs(eta[10]) >= 1) return(-Inf)
    scale0 <- exp(eta[7])
    if (!is.finite(scale0) || scale0 <= 0) return(-Inf)
    sc <- scale0 + eta[8] * sin1 + eta[9] * cos1
    if (any(!is.finite(sc)) || any(sc <= 0)) return(-Inf)
    mu_day <- eta[1] + eta[3] * sin1 + eta[4] * cos1 +
      eta[5] * sin2 + eta[6] * cos2
    mu <- outer(mu_day, eta[2] * g, "+")
    xi <- eta[10]
    z <- (threshold - mu) / sc
    zp <- (peak_tx - mu_day[peak_day] - eta[2] * g[peak_year]) / sc[peak_day]
    if (any(!is.finite(z)) || any(!is.finite(zp))) return(-Inf)
    if (abs(xi) < 1e-7) {
      integral <- sum(exp(-z)) / npy
      intensity <- -log(sc[peak_day]) - zp
    } else {
      h <- 1 + xi * z
      hp <- 1 + xi * zp
      if (any(h <= 0) || any(hp <= 0)) return(-Inf)
      integral <- sum(exp(-log(h) / xi)) / npy
      intensity <- -log(sc[peak_day]) - (1 + 1 / xi) * log(hp)
    }
    lp <- -integral + sum(intensity)
    if (prior_measure == "flat_scale_coefficients") lp <- lp + eta[7]
    if (is.finite(lp)) lp else -Inf
  }
  compiler::cmpfun(target)
}

smith_am_run <- function(log_target, coordinates, n_iter, burn_in, thin, seed,
                          adapt_start, adapt_every, freeze_after_burn,
                          mixture_probability, covariance_epsilon,
                          checkpoint_every, progress_every,
                          checkpoint_file = NULL, resume = FALSE, input_key = NULL) {
  # Haario-style Gaussian random-walk adaptation in fixed whitened coordinates.
  # A small fixed Gaussian component remains available throughout the chain.
  # Rejected states are included in both adaptation and posterior storage.
  p <- length(coordinates$center)
  stopifnot(n_iter > burn_in, burn_in >= 0, thin >= 1,
            all(c(n_iter, burn_in, thin) == as.integer(c(n_iter, burn_in, thin))),
            adapt_start >= 2L, adapt_every >= 1L, checkpoint_every >= 1L,
            progress_every >= 1L, mixture_probability > 0, mixture_probability < 1,
            covariance_epsilon > 0)
  n_keep <- floor((n_iter - burn_in) / thin)
  stopifnot(n_keep >= 2L)
  signature <- list(version = "smith-am-1", input = input_key,
    center = coordinates$center, transform = coordinates$transform,
    n_iter = n_iter, burn_in = burn_in, thin = thin, seed = seed,
    adapt_start = adapt_start, adapt_every = adapt_every,
    freeze_after_burn = freeze_after_burn, mixture_probability = mixture_probability,
    covariance_epsilon = covariance_epsilon)
  has_checkpoint <- !is.null(checkpoint_file) && file.exists(checkpoint_file)
  if (has_checkpoint && !resume) stop("Checkpoint exists. Enable resume or choose a new output directory.")
  if (has_checkpoint) {
    s <- readRDS(checkpoint_file)
    if (!identical(s$signature, signature)) stop("Checkpoint configuration/input mismatch.")
    assign(".Random.seed", s$rng_state, envir = .GlobalEnv)
  } else {
    set.seed(seed)
    z <- rep(0, p)
    lp <- log_target(coordinates$center)
    # Dispersed feasible starts across chains; no positive-GMST prior constraint.
    found <- FALSE
    for (attempt in seq_len(100L)) {
      z0 <- rnorm(p)
      eta0 <- coordinates$center + drop(coordinates$transform %*% z0)
      lp0 <- log_target(eta0)
      if (is.finite(lp0)) { z <- z0; lp <- lp0; found <- TRUE; break }
    }
    if (!found) stop("Could not initialize a dispersed feasible chain.")
    s <- list(signature = signature, iteration = 0L, z = z, log_posterior = lp,
      n_history = 1L, running_mean = z, m2 = matrix(0, p, p),
      proposal_chol = diag(2.38 / sqrt(p), p), accepted = 0L,
      accepted_sampling = 0L, kept = 0L,
      eta_draws = matrix(NA_real_, n_keep, p,
                         dimnames = list(NULL, names(coordinates$center))),
      log_posterior_draws = rep(NA_real_, n_keep))
  }
  if (s$iteration >= n_iter) return(s)
  for (iteration in seq.int(s$iteration + 1L, n_iter)) {
    if (runif(1) < mixture_probability) {
      step <- rnorm(p, sd = 0.1 / sqrt(p))
    } else {
      step <- drop(rnorm(p) %*% s$proposal_chol)
    }
    proposed_z <- s$z + step
    proposed_eta <- coordinates$center + drop(coordinates$transform %*% proposed_z)
    proposed_lp <- log_target(proposed_eta)
    if (is.finite(proposed_lp) && log(runif(1)) < proposed_lp - s$log_posterior) {
      s$z <- proposed_z
      s$log_posterior <- proposed_lp
      s$accepted <- s$accepted + 1L
      if (iteration > burn_in) s$accepted_sampling <- s$accepted_sampling + 1L
    }
    # Numerically stable online covariance using every visited state.
    s$n_history <- s$n_history + 1L
    delta <- s$z - s$running_mean
    s$running_mean <- s$running_mean + delta / s$n_history
    s$m2 <- s$m2 + tcrossprod(delta, s$z - s$running_mean)
    if (iteration >= adapt_start && iteration %% adapt_every == 0L &&
        (!freeze_after_burn || iteration <= burn_in)) {
      empirical <- s$m2 / (s$n_history - 1L)
      proposal <- (2.38^2 / p) *
        ((empirical + t(empirical)) / 2 + diag(covariance_epsilon, p))
      s$proposal_chol <- chol(proposal)
    }
    if (iteration > burn_in && (iteration - burn_in) %% thin == 0L) {
      s$kept <- s$kept + 1L
      s$eta_draws[s$kept, ] <- coordinates$center +
        drop(coordinates$transform %*% s$z)
      s$log_posterior_draws[s$kept] <- s$log_posterior
    }
    s$iteration <- iteration
    if (iteration %% progress_every == 0L)
      message("Chain seed ", seed, ": ", iteration, "/", n_iter,
              "; acceptance ", round(s$accepted / iteration, 3))
    if (!is.null(checkpoint_file) &&
        (iteration %% checkpoint_every == 0L || iteration == n_iter)) {
      s$rng_state <- get(".Random.seed", envir = .GlobalEnv)
      saveRDS(s, checkpoint_file)
    }
  }
  s
}

smith_am_raw <- function(eta) {
  raw <- eta
  raw[, 7] <- exp(eta[, 7])
  colnames(raw)[7] <- "sigma_intercept"
  raw
}

smith_am_diag <- function(x, variable) {
  stopifnot(is.matrix(x))
  finite <- all(is.finite(x))
  constant <- finite && length(unique(as.vector(x))) == 1L
  assess <- finite && !constant && nrow(x) >= 20L && ncol(x) >= 2L
  safe <- function(fun) if (assess)
    tryCatch(suppressWarnings(fun(x)), error = function(e) NA_real_) else NA_real_
  data.frame(variable = variable, mean = mean(x), sd = stats::sd(as.vector(x)),
    rhat = safe(posterior::rhat), ess_bulk = safe(posterior::ess_bulk),
    ess_tail = safe(posterior::ess_tail), mcse_mean = safe(posterior::mcse_mean),
    diagnostic_status = if (!finite) "Nonfinite values; diagnostics undefined" else
      if (constant) "Constant across retained draws; mixing cannot be established" else
        if (!assess) "Too few draws or chains" else "Calculated")
}

smith_am_parameter_diagnostics <- function(chains) {
  n_keep <- unique(vapply(chains, function(z) nrow(z$eta_draws), integer(1)))
  stopifnot(length(n_keep) == 1L)
  raw <- lapply(chains, function(z) smith_am_raw(z$eta_draws))
  do.call(rbind, lapply(seq_len(ncol(raw[[1]])), function(j) {
    x <- do.call(cbind, lapply(raw, function(z) z[, j]))
    smith_am_diag(x, colnames(raw[[1]])[j])
  }))
}

smith_am_integrals <- function(mu, sigma, xi, temperatures, npy) {
  vapply(temperatures, function(x) {
    smith_log_sum(smith_log_tail(x, mu, sigma, xi)) - log(npy)
  }, 0.0)
}

smith_am_predict <- function(chains, mle, temperatures) {
  n <- nrow(chains[[1]]$eta_draws)
  nc <- length(chains)
  nt <- length(temperatures)
  log_f <- log_c <- array(NA_real_, c(n, nc, nt))
  endpoint_f <- endpoint_c <- slope <- matrix(NA_real_, n, nc)
  f <- smith_design(mle$factual)
  c <- smith_design(mle$counterfactual)
  for (j in seq_len(nc)) {
    raw <- smith_am_raw(chains[[j]]$eta_draws)
    slope[, j] <- raw[, 2]
    for (i in seq_len(n)) {
      b <- raw[i, ]
      sf <- drop(f$sigma %*% b[7:9])
      sc <- drop(c$sigma %*% b[7:9])
      stopifnot(all(sf > 0), all(sc > 0))
      mf <- drop(f$mu %*% b[1:6])
      mc <- drop(c$mu %*% b[1:6])
      log_f[i, j, ] <- smith_am_integrals(mf, sf, b[10], temperatures, mle$model$npy)
      log_c[i, j, ] <- smith_am_integrals(mc, sc, b[10], temperatures, mle$model$npy)
      endpoint_f[i, j] <- if (b[10] < 0) max(mf - sf / b[10]) else Inf
      endpoint_c[i, j] <- if (b[10] < 0) max(mc - sc / b[10]) else Inf
    }
    message("Posterior prediction: chain ", j, "/", nc)
  }
  list(log_count_factual = log_f, log_count_counterfactual = log_c,
       endpoint_factual = endpoint_f, endpoint_counterfactual = endpoint_c,
       gmst_slope = slope, temperatures = temperatures)
}

smith_am_log_probabilities <- function(x) {
  out <- x
  mid <- is.finite(x) & x >= -35 & x <= 35
  out[mid] <- log(-expm1(-exp(x[mid])))
  out[x > 35] <- 0
  out
}

smith_am_inverse_ratio <- function(log_f, log_c) {
  # Follow the author's explicit posterior convention: 0/0 maps to zero.
  # This is not a mathematical definition of 0/0. Save the affected fraction.
  # Negative-GMST-slope draws are retained and may yield infinite inverse ratios.
  both <- log_f == -Inf & log_c == -Inf
  log_ratio <- log_c - log_f
  log_ratio[both] <- -Inf
  list(value = exp(log_ratio), both_zero = both)
}

smith_am_summarize <- function(prediction, min_positive_draws) {
  summaries <- diagnostics <- list()
  counter <- dc <- 0L
  positive <- prediction$gmst_slope > 0
  for (k in seq_along(prediction$temperatures)) {
    temperature <- prediction$temperatures[k]
    count_f <- prediction$log_count_factual[, , k, drop = FALSE]
    count_c <- prediction$log_count_counterfactual[, , k, drop = FALSE]
    dim(count_f) <- dim(count_c) <- dim(prediction$gmst_slope)
    probability_f <- smith_am_log_probabilities(count_f)
    probability_c <- smith_am_log_probabilities(count_c)
    for (basis in c("expected_events", "at_least_one_probability")) {
      lf <- if (basis == "expected_events") count_f else probability_f
      lc <- if (basis == "expected_events") count_c else probability_c
      vf <- exp(lf)
      vc <- exp(lc)
      ratio <- smith_am_inverse_ratio(lf, lc)
      supports_f <- lf > -Inf
      supports_c <- lc > -Inf
      quantities <- list(factual = vf, counterfactual = vc,
                         inverse_ratio_author_convention = ratio$value,
                         supports_factual = 1 * supports_f,
                         supports_counterfactual = 1 * supports_c)
      for (name in names(quantities)) {
        dc <- dc + 1L
        diagnostic <- smith_am_diag(quantities[[name]], paste(basis, temperature, name, sep = ":"))
        diagnostic$temperature_C <- temperature
        diagnostic$basis <- basis
        diagnostics[[dc]] <- diagnostic
      }
      for (subset in c("all_draws", "positive_gmst")) {
        selected <- if (subset == "positive_gmst") positive else
          matrix(TRUE, nrow(positive), ncol(positive))
        r <- ratio$value[selected]
        qs <- if (length(r)) stats::quantile(r, c(0.025, 0.975), names = FALSE, type = 1) else
          c(NA_real_, NA_real_)
        counter <- counter + 1L
        summaries[[counter]] <- data.frame(
          temperature_C = temperature, basis = basis, ratio_subset = subset,
          E_pF = mean(vf), E_pC = mean(vc),
          E_inverse_RR = if (length(r)) mean(r) else NA_real_,
          q0_025 = qs[1], q0_975 = qs[2],
          Pr_inverse_RR_lt_0_01 = if (length(r)) mean(r < 0.01) else NA_real_,
          ratio_of_posterior_means = mean(vc) / mean(vf),
          posterior_Pr_positive_gmst = mean(positive), n_draws = length(vf),
          n_ratio_draws = length(r), fraction_both_zero = mean(ratio$both_zero),
          n_support_factual = sum(supports_f), n_support_counterfactual = sum(supports_c),
          rare_tail_status = if (sum(supports_f) < min_positive_draws ||
                                sum(supports_c) < min_positive_draws)
            "Too few supporting draws; long-chain assessment required" else
              "Inspect ESS, Rhat and MCSE; support count alone is insufficient")
      }
    }
  }
  list(summary = do.call(rbind, summaries), diagnostics = do.call(rbind, diagnostics))
}

smith_am_paper_table <- function() {
  # Transcribed from Table 0.2, not targets imposed on the sampler.
  data.frame(temperature_C = 34:39,
    paper_E_pF = c(0.122, 0.031, 4.8e-3, 5.0e-4, 4.9e-5, 5.3e-6),
    paper_E_pC = c(0.0022, 2.1e-4, 2.0e-5, 2.7e-6, 4.8e-7, 1.0e-7),
    paper_E_inverse_RR = c(0.020, 4.8e-3, 1.4e-3, 6.4e-4, 3.6e-4, 2.1e-5),
    paper_q0_025 = 0, paper_q0_975 = c(0.132, 0.053, 0.0013, 0.00083, 0, 0),
    paper_Pr_inverse_RR_lt_0_01 = c(0.65, 0.901, 0.971, 0.988, 0.994, 0.996))
}


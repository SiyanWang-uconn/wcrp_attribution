# Function-only helper: check_gev_fit(data) checks Annual and June.
# KS assesses PIT distribution; Ljung-Box assesses serial correlation of normal PIT scores.
# Both p-values use independent simulation followed by GEV refitting.

gev_pit <- function(x, g, fit) {
  z <- (x - fit$mu - fit$alpha * (g - fit$g0)) / fit$sigma
  if (abs(fit$xi) < 1e-7) return(exp(-exp(-z)))
  h <- 1 + fit$xi * z
  u <- rep(if (fit$xi < 0) 1 else 0, length(x))
  inside <- h > 0
  u[inside] <- exp(-exp(-log(h[inside]) / fit$xi))
  u
}

ks_distance <- function(u) {
  stopifnot(length(u) > 0, all(is.finite(u)), all(u >= 0 & u <= 1))
  u <- sort(u)
  n <- length(u)
  max(seq_len(n) / n - u, u - (seq_len(n) - 1) / n)
}

pit_normal_scores <- function(u) {
  # Bound numerical 0/1 values to avoid infinite normal scores.
  qnorm(pmin(pmax(u, .Machine$double.eps), 1 - .Machine$double.eps))
}

ljung_box_statistic <- function(u, lag = 5L) {
  z <- pit_normal_scores(u)
  stopifnot(all(is.finite(z)), sd(z) > 0, length(z) > lag)
  unname(stats::Box.test(z, lag = lag, type = "Ljung-Box", fitdf = 0)$statistic)
  # Only Q is used. Parameter-estimation effects are calibrated by refitting,
  # not by subtracting the four GEV parameters from Ljung-Box degrees of freedom.
}

check_gev_fit <- function(data, product = "", index = "", B = 999L,
                          seed = 20261007L, lag = 5L,
                          gev_file = here::here("Code", "helpers", "gev.R")) {
  gev <- new.env(parent = environment())
  sys.source(gev_file, envir = gev)
  if (is.character(data) && length(data) == 1L) data <- read.csv(data)
  stopifnot(is.data.frame(data),
    all(c("year", "season", "X", "gmst", "used_in_fit") %in% names(data)),
    length(B) == 1, is.finite(B), B >= 1, B == as.integer(B),
    length(seed) == 1, is.finite(seed), seed >= 0,
    seed < .Machine$integer.max, seed == as.integer(seed),
    length(lag) == 1, is.finite(lag), lag >= 1, lag == as.integer(lag))
  seasons <- c("Annual", "June")
  selected <- lapply(seasons, function(s) {
    d <- data[data$season %in% s & data$used_in_fit %in% TRUE, ]
    d <- d[order(d$year), ]
    if (nrow(d) < max(10, lag + 1)) stop(s, ": insufficient fitting years")
    stopifnot(!anyDuplicated(d$year), all(is.finite(d$year)),
      all(d$year == as.integer(d$year)), all(is.finite(d$X)),
      all(is.finite(d$gmst)), sd(d$gmst) > 0)
    d
  })
  details <- setNames(vector("list", 2L), seasons)
  for (j in seq_along(seasons)) {
    d <- selected[[j]]
    x <- d$X
    g <- d$gmst
    continuous <- all(diff(d$year) == 1)
    if (!continuous) warning(seasons[j], ": missing years; Ljung-Box is skipped, KS still runs.")
    fit <- gev$fit_gev(x, g)
    u <- gev_pit(x, g, fit)
    observed <- c(KS = ks_distance(u),
                  LB = if (continuous) ljung_box_statistic(u, lag) else NA_real_)
    mu <- fit$mu + fit$alpha * (g - fit$g0)
    boot <- matrix(NA_real_, B, 2L, dimnames = list(NULL, c("KS", "LB")))
    failure <- rep(NA_character_, B)
    set.seed(seed + j - 1L)
    message(product, " ", index, " ", seasons[j], ": ", nrow(d), " years")
    for (b in seq_len(B)) {
      simulated <- tryCatch({
        xb <- gev$qgev_local(runif(length(x)), mu, fit$sigma, fit$xi)
        fb <- gev$fit_gev(xb, g)
        ub <- gev_pit(xb, g, fb)
        list(stats = c(ks_distance(ub), if (continuous) ljung_box_statistic(ub, lag) else NA_real_),
             error = NA_character_)
      }, error = function(e) list(stats = c(NA_real_, NA_real_), error = conditionMessage(e)))
      boot[b, ] <- simulated$stats
      failure[b] <- simulated$error
      if (b %% 100L == 0L) message("Fit-check bootstrap: ", b, "/", B)
    }
    summary <- data.frame(product, index, season = seasons[j],
      first_year = min(d$year), last_year = max(d$year), n = length(x),
      years_continuous = continuous, KS_D = observed[["KS"]], LB_Q = observed[["LB"]], LB_lag = lag,
      B_requested = B, B_success = sum(is.na(failure)), B_failed = sum(!is.na(failure)),
      seed = seed + j - 1L)
    for (test in c("KS", "LB")) {
      ok <- is.finite(boot[, test])
      m <- sum(ok)
      k <- if (is.finite(observed[test])) sum(boot[ok, test] >= observed[test]) else NA_real_
      summary[[paste0(test, "_p")]] <- if (m > 0) (1 + k) / (1 + m) else NA_real_
      summary[[paste0(test, "_p_lower")]] <- if (is.finite(k)) (1 + k) / (1 + B) else NA_real_
      summary[[paste0(test, "_p_upper")]] <- if (is.finite(k)) (1 + k + B - m) / (1 + B) else NA_real_
    }
    details[[j]] <- list(summary = summary, fit = fit, year = d$year,
      PIT = u, normal_scores = pit_normal_scores(u), bootstrap = boot, failure = failure)
    if (any(!is.na(failure))) warning(seasons[j], ": ", sum(!is.na(failure)),
      " bootstrap failures; inspect p-value bounds and failure messages.")
  }
  summary <- do.call(rbind, lapply(details, function(z) z$summary))
  rownames(summary) <- NULL
  list(summary = summary, details = details)
}
# Per-fit p-values are conditional on successful replicates and not multiplicity-adjusted.
# Large p-values prove neither model correctness nor independence of every kind.

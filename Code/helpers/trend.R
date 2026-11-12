# Trend diagnostics for one seasonal observation per year.
# Definitions only: sourcing this file does not run tests or change any data.
# Mann-Kendall tests a monotonic temporal trend, not every form of nonstationarity.
# A significant result is not proof of a causal GMST effect; a nonsignificant
# result does not establish stationarity. Do not select the PP model by p alone.

mann_kendall_trend <- function(year, value, series, value_units, alpha = 0.05) {
  stopifnot(length(year) == length(value), is.numeric(year), is.numeric(value),
            length(alpha) == 1L, is.finite(alpha), alpha > 0, alpha < 1)
  keep <- is.finite(year) & is.finite(value)
  omitted <- sum(!keep)
  year <- year[keep]
  value <- value[keep]
  order_year <- order(year)
  year <- year[order_year]
  value <- value[order_year]
  n <- length(year)
  stopifnot(n >= 4L, !anyDuplicated(year))

  # With unique ordered years, Kendall's association with time is the MK trend
  # test. Use a two-sided asymptotic test with ties and continuity correction.
  # No additional R packages are required. A constant series has no trend.
  if (length(unique(value)) == 1L) {
    tau <- 0
    z <- 0
    p <- 1
  } else {
    test <- stats::cor.test(year, value, method = "kendall", alternative = "two.sided",
                            exact = FALSE, continuity = TRUE)
    tau <- unname(test$estimate)
    z <- unname(test$statistic)
    p <- test$p.value
  }

  # Sen's slope is the median of all pairwise slopes. Actual year differences
  # are used, so excluded/incomplete years do not distort the slope units.
  pairs <- utils::combn(seq_len(n), 2L)
  slopes <- (value[pairs[2L, ]] - value[pairs[1L, ]]) /
    (year[pairs[2L, ]] - year[pairs[1L, ]])
  sen <- stats::median(slopes)

  # This is a descriptive dependence diagnostic, not a second significance test.
  # Remove the Sen trend only for this diagnostic; MK uses the original values.
  # Match actual adjacent years rather than treating gaps as one-year lags.
  residual <- value - sen * (year - min(year))
  next_index <- match(year + 1, year)
  adjacent <- which(!is.na(next_index))
  lag1 <- NA_real_
  if (length(adjacent) >= 3L) {
    left <- residual[adjacent]
    right <- residual[next_index[adjacent]]
    if (stats::sd(left) > 0 && stats::sd(right) > 0)
      lag1 <- stats::cor(left, right)
  }
  significant <- p < alpha
  direction <- if (!significant) "No significant monotonic trend" else
    if (tau > 0) "Increasing" else "Decreasing"

  # Classical MK assumes independence between years. The p-value below is NOT
  # corrected for serial correlation; inspect dependence before interpretation.
  data.frame(
    series = series, n = n, n_omitted = omitted,
    start_year = min(year), end_year = max(year),
    kendall_tau = tau, z = z, p_value = p, alpha = alpha,
    significant = significant, trend = direction,
    sen_slope_per_year = sen, sen_slope_per_decade = 10 * sen,
    value_units = value_units, detrended_lag1_correlation = lag1,
    adjacent_year_pairs = length(adjacent),
    method = "Two-sided asymptotic Mann-Kendall; ties and continuity correction",
    dependence_adjustment = "None; inspect detrended lag-1 correlation")
}

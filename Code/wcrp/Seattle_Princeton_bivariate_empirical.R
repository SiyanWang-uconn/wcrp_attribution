# Empirical joint heat exceedances at Seattle and Princeton.
# Run sections line by line in RStudio. Uses base R only.

# ---- Configuration ----
PROJECT_DIR <- "C:/Users/siyan/WCRP/wcrp_common_attribution_assessment"
INPUT_DIR <- file.path(PROJECT_DIR, "outputs/pnw_heatwave/ghcn_stations/seattle_princeton")
OUTPUT_DIR <- file.path(PROJECT_DIR, "Results/Seattle_Princeton_bivariate")
STATION_IDS <- c("USW00024233", "CAN01126510")
START_YEAR <- 1948L
END_YEAR <- 2025L
BASELINE_YEARS <- c(1961L, 1990L)
SUMMER_MONTHS <- 6:8
WINDOW_DAYS <- 3L
THRESHOLD_PROB <- 0.95
MIN_YEAR_COVERAGE <- 0.90
EVENT_START <- as.Date("2021-06-26")
BOOTSTRAP_REPS <- 2000L
RANDOM_SEED <- 20261111L

# ---- Read QC-screened observations; do not remove missing dates ----
read_station <- function(id) {
  f <- file.path(INPUT_DIR, paste0(id, "_daily_tmax_1948_2025.csv"))
  d <- read.csv(f, na.strings = c("", "NA"), stringsAsFactors = FALSE)
  d$date <- as.Date(d$date)
  stopifnot(!anyNA(d$date), !anyDuplicated(d$date), all(d$station_id == id))
  flagged <- !is.na(d$quality_flag) & nzchar(d$quality_flag)
  stopifnot(!any(flagged & !is.na(d$tmax_C)))
  d[, c("date", "tmax_C")]
}
a <- read_station(STATION_IDS[1])
b <- read_station(STATION_IDS[2])
names(a)[2] <- "tmax_A_C"
names(b)[2] <- "tmax_B_C"
calendar <- data.frame(date = seq(as.Date(paste0(START_YEAR, "-01-01")),
  as.Date(paste0(END_YEAR, "-12-31")), by = "day"))
daily <- merge(merge(calendar, a, by = "date", all.x = TRUE),
  b, by = "date", all.x = TRUE)
stopifnot(all(diff(daily$date) == 1))

# ---- Forward rolling means, labeled by their first day ----
forward_mean <- function(x, k) {
  ans <- rep(NA_real_, length(x))
  for (i in seq_len(length(x) - k + 1L)) {
    z <- x[i:(i + k - 1L)]
    if (all(is.finite(z))) ans[i] <- mean(z)
  }
  ans
}
daily$X_C <- forward_mean(daily$tmax_A_C, WINDOW_DAYS)
daily$Y_C <- forward_mean(daily$tmax_B_C, WINDOW_DAYS)
daily$window_end <- daily$date + WINDOW_DAYS - 1L
daily$year <- as.integer(format(daily$date, "%Y"))
keep <- vapply(daily$date, function(t) {
  dates <- as.Date(t, origin = "1970-01-01") + seq.int(0L, WINDOW_DAYS - 1L)
  all(as.integer(format(dates, "%m")) %in% SUMMER_MONTHS) &&
    length(unique(format(dates, "%Y"))) == 1L
}, logical(1))
w <- daily[keep, ]
w$valid_pair <- is.finite(w$X_C) & is.finite(w$Y_C)
baseline <- w$valid_pair & w$year >= BASELINE_YEARS[1] &
  w$year <= BASELINE_YEARS[2]
stopifnot(sum(baseline) > 0)
u_A <- unname(quantile(w$X_C[baseline], THRESHOLD_PROB, type = 7))
u_B <- unname(quantile(w$Y_C[baseline], THRESHOLD_PROB, type = 7))
w$exceed_A <- ifelse(w$valid_pair, w$X_C > u_A, NA)
w$exceed_B <- ifelse(w$valid_pair, w$Y_C > u_B, NA)
w$state <- ifelse(w$valid_pair,
  paste0(as.integer(w$exceed_A), as.integer(w$exceed_B)), NA)

# ---- Annual four-state counts and completeness ----
annual <- do.call(rbind, lapply(split(w, w$year), function(z) {
  s <- z$state[z$valid_pair]
  data.frame(year = z$year[1], expected_windows = nrow(z),
    valid_windows = length(s), coverage = length(s) / nrow(z),
    n00 = sum(s == "00"), n10 = sum(s == "10"),
    n01 = sum(s == "01"), n11 = sum(s == "11"))
}))
rownames(annual) <- NULL
annual$joint_probability <- annual$n11 / annual$valid_windows
annual$adequate_coverage <- annual$coverage >= MIN_YEAR_COVERAGE
stopifnot(all(rowSums(annual[, c("n00", "n10", "n01", "n11")]) ==
  annual$valid_windows))

# ---- Probabilities and whole-year bootstrap intervals ----
metrics <- function(z) {
  n <- sum(z$valid_windows)
  counts <- colSums(z[, c("n00", "n10", "n01", "n11"), drop = FALSE])
  pA <- (counts["n10"] + counts["n11"]) / n
  pB <- (counts["n01"] + counts["n11"]) / n
  pAB <- counts["n11"] / n
  unname(c(pA, pB, pAB, pA * pB,
    if (pA * pB > 0) pAB / (pA * pB) else NA_real_,
    if (pA > 0) pAB / pA else NA_real_,
    if (pB > 0) pAB / pB else NA_real_))
}
metric_names <- c("p_A", "p_B", "p_AB", "p_independent", "PMF",
  "p_B_given_A", "p_A_given_B")
set.seed(RANDOM_SEED)
summarize_years <- function(z, label) {
  boot <- replicate(BOOTSTRAP_REPS,
    metrics(z[sample.int(nrow(z), nrow(z), replace = TRUE), , drop = FALSE]))
  ci <- t(apply(boot, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE))
  data.frame(analysis = label, metric = metric_names, estimate = metrics(z),
    ci_lower = ci[, 1], ci_upper = ci[, 2], years = nrow(z),
    valid_windows = sum(z$valid_windows))
}
summary_table <- rbind(
  summarize_years(annual, "1948-2025_all_valid_windows"),
  summarize_years(annual[annual$adequate_coverage, ], "coverage_at_least_90pct"),
  summarize_years(annual[annual$year >= 1961 & annual$year <= 1990, ], "1961-1990"),
  summarize_years(annual[annual$year >= 1991, ], "1991-2025"))
thresholds <- data.frame(station_id = STATION_IDS,
  station = c("Seattle-Tacoma Airport", "Princeton A, BC"),
  baseline_start = BASELINE_YEARS[1], baseline_end = BASELINE_YEARS[2],
  percentile = THRESHOLD_PROB, threshold_C = c(u_A, u_B),
  paired_baseline_windows = sum(baseline), window_days = WINDOW_DAYS)
event <- w[w$date >= as.Date("2021-06-23") & w$date <= as.Date("2021-07-02"), ]
target <- w[w$date == EVENT_START, ]
stopifnot(nrow(target) == 1L, target$valid_pair)

# ---- Save empirical results ----
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
write.csv(w, file.path(OUTPUT_DIR, "paired_3day_JJA_windows.csv"), row.names = FALSE)
write.csv(annual, file.path(OUTPUT_DIR, "annual_counts.csv"), row.names = FALSE)
write.csv(thresholds, file.path(OUTPUT_DIR, "thresholds.csv"), row.names = FALSE)
write.csv(summary_table, file.path(OUTPUT_DIR, "empirical_summary.csv"), row.names = FALSE)
write.csv(event, file.path(OUTPUT_DIR, "event_2021_windows.csv"), row.names = FALSE)
writeLines(c("# Seattle-Princeton empirical joint heat analysis", "",
  "QC-screened daily Tmax; no infilling, homogenization or spatial averaging here.",
  "Forward 3-day means are labeled by their first day. All dates must be in JJA.",
  "Thresholds use paired baseline windows, 1961-1990, quantile type 7, strict >.",
  "Thresholds are pooled summer percentiles, not calendar-day percentiles.",
  "A 3-day mean exceedance does not require all three daily values to exceed a threshold.",
  "PMF = p_AB / (p_A * p_B), using observed marginals in each analysis period.",
  "Window counts are not counts of independent heatwave events.",
  "95% intervals: 2000 whole-year bootstrap replicates; baseline thresholds stay fixed.",
  "Intervals omit threshold uncertainty and do not model long-term nonstationarity.",
  "Full-period probabilities describe a historical mixture, not the current climate.",
  "Coverage sensitivity excludes years below 90% paired rolling-window coverage.",
  "The target June 26-28, 2021 window is a candidate, not an optimized joint peak.",
  "No EVT model, annual return period, significance test or attribution PR is fitted."),
  file.path(OUTPUT_DIR, "README.md"))

print(thresholds, row.names = FALSE)
print(summary_table, row.names = FALSE)
print(annual[!annual$adequate_coverage, ], row.names = FALSE)
print(event[, c("date", "window_end", "X_C", "Y_C", "state")], row.names = FALSE)
cat("Results:", OUTPUT_DIR, "\n")

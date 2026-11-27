# Station-based Smith-inspired bivariate tail model.
# Run numbered sections interactively in RStudio; no automatic package installs.
# Existing scripts/results are preserved. This is not Smith's exact PP likelihood.

# ---- 1. Paths and configuration ----
project_root <- "C:/Users/siyan/WCRP/wcrp_common_attribution_assessment"
input_dir <- file.path(project_root, "Results/Seattle_Princeton_bivariate")
paired_file <- file.path(input_dir, "paired_3day_JJA_windows.csv")
threshold_file <- file.path(input_dir, "thresholds.csv")
gmst_file <- file.path(project_root, "data/processed/gmst_smoothed.csv")
output_dir <- file.path(project_root, "Results/Seattle_Princeton_bivariate_Smith_location")
station_ids <- c("USW00024233", "CAN01126510")
fit_years <- c(1950L, 2020L)
minimum_paired_coverage <- 0.90
expected_JJA_windows <- 90L
event_start <- as.Date("2021-06-26")
factual_year <- 2021L
reference_year <- 1951L
location_harmonic_order <- 2L
scale_harmonic_order <- 1L
npy <- 365
optimizer_maxit <- 3000L
# "event" uses June 26-28, 2021 temperatures; "baseline" uses baseline thresholds.
target_mode <- "event"

# ---- 2. Load definitions and paired station observations ----
# The older helper supplies stable arithmetic and data preparation only.
# Its example configuration uses uppercase names and does not fit or save data.
source(file.path(project_root, "Code/wcrp/bivariate_pot_gmst.R"))
source(file.path(project_root, "Code/wcrp/bivariate_smith_location.R"))
paired <- read.csv(paired_file, stringsAsFactors = FALSE)
paired$date <- as.Date(paired$date)
gmst <- read.csv(gmst_file)
stopifnot(!anyNA(paired$date), !anyDuplicated(paired$date),
  all(paired$year == as.integer(format(paired$date, "%Y"))),
  all(as.integer(format(paired$date, "%m")) %in% 6:8),
  !anyDuplicated(gmst$year), target_mode %in% c("event", "baseline"))
threshold_table <- read.csv(threshold_file, stringsAsFactors = FALSE)
stopifnot(!anyDuplicated(threshold_table$station_id),
  all(station_ids %in% threshold_table$station_id))
threshold <- threshold_table$threshold_C[match(station_ids, threshold_table$station_id)]
complete <- is.finite(paired$X_C) & is.finite(paired$Y_C)
coverage <- aggregate(as.integer(complete), list(year = paired$year), sum)
names(coverage)[2] <- "paired_windows"
coverage$coverage <- coverage$paired_windows / expected_JJA_windows
coverage$used <- coverage$year >= fit_years[1] & coverage$year <= fit_years[2] &
  coverage$coverage >= minimum_paired_coverage
selected <- paired[paired$year %in% coverage$year[coverage$used], ]
fit_data <- prepare_bivariate_pot(selected, gmst)
stopifnot(length(unique(fit_data$year)) >= 20)
event_row <- paired[paired$date == event_start, ]
stopifnot(nrow(event_row) == 1, all(is.finite(unlist(event_row[c("X_C", "Y_C")]))))
target <- if (target_mode == "event") unname(c(event_row$X_C, event_row$Y_C)) else threshold
coverage
data.frame(station = station_ids, threshold_C = threshold, target_C = target)

# ---- 3. Joint fit with seasonal Smith-style marginal parameters ----
# mu_j = intercept + beta_G * centered_GMST + seasonal harmonics.
# sigma_j = positive baseline + seasonal harmonics (identity scale).
# Shape and logistic alpha are constant. No rate-GMST term is estimated.
# All four censoring states are retained, including neither station exceeding.
fit <- fit_bivariate_smith(fit_data, threshold,
  location_order = location_harmonic_order, scale_order = scale_harmonic_order,
  npy = npy, maxit = optimizer_maxit)
parameter_table <- data.frame(parameter = names(fit$coefficients),
  estimate = unname(fit$coefficients))
optimizer_table <- do.call(rbind, lapply(seq_along(fit$candidates), function(i) {
  z <- fit$candidates[[i]]
  data.frame(start = i, convergence = z$convergence, nll = z$value)
}))
parameter_table
optimizer_table

# ---- 4. Conditional and joint probabilities at the same calendar phase ----
# Compare climate states at June 26, not different seasonal dates.
# These are probabilities per paired window conditional on calendar phase.
prediction_years <- c(factual_year, reference_year)
prediction_gmst <- gmst$gmst[match(prediction_years, gmst$year)]
stopifnot(all(is.finite(prediction_gmst)))
probabilities <- cbind(climate_year = prediction_years,
  predict_bivariate_smith(fit, prediction_gmst, event_start, target))
stopifnot(all(probabilities$p_joint <= pmin(probabilities$p_A, probabilities$p_B)),
  all(abs(probabilities$p_B_given_A * probabilities$p_A - probabilities$p_joint) < 1e-12),
  all(abs(probabilities$p_A_given_B * probabilities$p_B - probabilities$p_joint) < 1e-12))
quantities <- c("p_A", "p_B", "p_joint", "p_B_given_A", "p_A_given_B")
lf <- unlist(probabilities[1, paste0("log_", quantities)], use.names = FALSE)
lr <- unlist(probabilities[2, paste0("log_", quantities)], use.names = FALSE)
comparison <- data.frame(quantity = quantities,
  p_factual = unlist(probabilities[1, quantities], use.names = FALSE),
  p_reference = unlist(probabilities[2, quantities], use.names = FALSE),
  log_PR = lf - lr, PR = exp(lf - lr))
# Zero conditioning-event probabilities give undefined conditional quantities.
# 0/0 PRs remain NaN; do not replace these with finite numbers.
probabilities
comparison

# ---- 5. Inspect seasonal scale and endpoint diagnostics before saving ----
diagnostic_dates <- seq(as.Date("2021-06-01"), as.Date("2021-08-29"), by = "day")
grid <- expand.grid(date = as.character(diagnostic_dates), gmst = prediction_gmst)
design <- bs_design(as.Date(grid$date), grid$gmst, fit$g0,
  fit$location_order, fit$scale_order)
diagnostics <- do.call(rbind, lapply(1:2, function(j) {
  p <- bs_parameters(fit$coefficients, design, c("A", "B")[j])
  data.frame(station = station_ids[j], grid,
    mu = p$mu, sigma = p$sigma, xi = p$xi,
    upper_endpoint_C = if (p$xi < 0) p$mu - p$sigma / p$xi else Inf)
}))
stopifnot(all(is.finite(diagnostics$sigma)), all(diagnostics$sigma > 0))
summary(diagnostics)
# Review threshold adequacy, seasonal identifiability, support and dependence.
# No independent-window Hessian intervals or annual return periods are reported.
# Whole-year refitting bootstrap is a separate stage, not launched here.

# ---- 6. Save to a new result directory ----
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(parameter_table, file.path(output_dir, "parameters.csv"), row.names = FALSE)
write.csv(optimizer_table, file.path(output_dir, "optimizer_starts.csv"), row.names = FALSE)
write.csv(coverage, file.path(output_dir, "coverage.csv"), row.names = FALSE)
write.csv(probabilities, file.path(output_dir, "probabilities.csv"), row.names = FALSE)
write.csv(comparison, file.path(output_dir, "climate_comparison.csv"), row.names = FALSE)
write.csv(diagnostics, file.path(output_dir, "seasonal_diagnostics.csv"), row.names = FALSE)
saveRDS(list(fit = fit, fit_data = fit_data, station_ids = station_ids,
  coverage = coverage, target = target, target_mode = target_mode,
  event_start = event_start, prediction_years = prediction_years,
  probabilities = probabilities, comparison = comparison,
  interpretation = "Station 3-day window working likelihood; point estimates only; historical GMST contrast"),
  file.path(output_dir, "smith_location_bivariate.rds"))

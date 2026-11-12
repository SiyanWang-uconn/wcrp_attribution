# Seattle-Princeton nonstationary bivariate POT: default model.
# Run sections interactively in RStudio. No figures or output files are created.
# Station rates vary with GMST; GPD scales, shapes and dependence are constant.

# ---- 1. Configuration ----
project_root <- "C:/Users/siyan/WCRP/wcrp_common_attribution_assessment"
source(file.path(project_root, "Code/wcrp/bivariate_pot_gmst.R"))
paired_file <- file.path(project_root, "Results/Seattle_Princeton_bivariate",
                         "paired_3day_JJA_windows.csv")
threshold_file <- file.path(project_root, "Results/Seattle_Princeton_bivariate",
                            "thresholds.csv")
gmst_file <- file.path(project_root, "data/processed/gmst_smoothed.csv")
station_ids <- c("USW00024233", "CAN01126510")
event_start <- as.Date("2021-06-26")
factual_year <- 2021L
reference_year <- 1951L
prediction_years <- c(factual_year, reference_year)
rate_gmst <- c(TRUE, TRUE)
scale_gmst <- c(FALSE, FALSE)

# ---- 2. Paired windows, exact thresholds and annual GMST ----
paired <- read.csv(paired_file, stringsAsFactors = FALSE)
paired$date <- as.Date(paired$date)
stopifnot(!anyDuplicated(paired$date))
gmst <- read.csv(gmst_file)
threshold_table <- read.csv(threshold_file, stringsAsFactors = FALSE)
stopifnot(!anyDuplicated(threshold_table$station_id),
          all(station_ids %in% threshold_table$station_id))
threshold <- threshold_table$threshold_C[
  match(station_ids, threshold_table$station_id)]
# Keep all complete pairs, including both-below and singly-exceeding windows.
# Incomplete pairs are explicitly excluded; GMST gaps are not silently dropped.
fit_data <- prepare_bivariate_pot(paired, gmst)
target_row <- paired[paired$date == event_start, , drop = FALSE]
stopifnot(nrow(target_row) == 1L,
          all(is.finite(unlist(target_row[, c("X_C", "Y_C")]))))
target <- unname(c(target_row$X_C, target_row$Y_C))
data_summary <- data.frame(
  first_year = min(fit_data$year), last_year = max(fit_data$year),
  years = length(unique(fit_data$year)), paired_windows = nrow(fit_data),
  threshold_A_C = threshold[1], threshold_B_C = threshold[2],
  target_A_C = target[1], target_B_C = target[2])
data_summary

# ---- 3. Fit the default working-likelihood model ----
fit <- fit_bivariate_pot(fit_data, threshold,
  rate_gmst = rate_gmst, scale_gmst = scale_gmst)
parameter_table <- data.frame(parameter = names(fit$coefficients),
                               estimate = unname(fit$coefficients))
optimizer_table <- do.call(rbind, lapply(seq_along(fit$optimizer_starts), function(i) {
  z <- fit$optimizer_starts[[i]]
  if (is.null(z)) return(data.frame(start = i, convergence = NA_integer_, nll = NA_real_))
  data.frame(start = i, convergence = z$convergence, nll = z$value)
}))
stopifnot(fit$convergence == 0L, is.finite(fit$logLik))
parameter_table
optimizer_table
# GPD scale is exp(logscale_*_intercept); alpha is plogis(alpha_logit).
# Rate slopes apply to logit(exceedance probability), not directly to temperature.
p0 <- bp_parameters(fit$coefficients, fit$g0, fit$g0, fit$spec)
model_summary <- data.frame(gmst_center = fit$g0, log_likelihood = fit$logLik,
  scale_A_C = exp(p0$A$logscale), scale_B_C = exp(p0$B$logscale),
  shape_A = p0$A$shape, shape_B = p0$B$shape, alpha = p0$alpha,
  upper_endpoint_A_C = if (p0$A$shape < 0)
    threshold[1] - exp(p0$A$logscale) / p0$A$shape else Inf,
  upper_endpoint_B_C = if (p0$B$shape < 0)
    threshold[2] - exp(p0$B$logscale) / p0$B$shape else Inf)
model_summary

# ---- 4. Predict at explicitly identified observed climate states ----
reference_gmst <- gmst$gmst[match(prediction_years, gmst$year)]
stopifnot(all(is.finite(reference_gmst)))
probabilities <- cbind(year = prediction_years,
  predict_bivariate_pot(fit, gmst = reference_gmst, target = target))
stopifnot(all(probabilities$p_joint >= 0),
  all(probabilities$p_joint <= pmin(probabilities$p_A, probabilities$p_B)),
  all(abs(probabilities$p_B_given_A * probabilities$p_A -
            probabilities$p_joint) < 1e-12),
  all(abs(probabilities$p_A_given_B * probabilities$p_B -
            probabilities$p_joint) < 1e-12))
probabilities

# These are window probabilities, not annual probabilities or return periods.
# The reference below is a 1951-like GMST climate, not an anthropogenic-free world.
# Predictions are preliminary: threshold sensitivity, marginal/dependence model
# adequacy, coverage sensitivity and refitting bootstrap uncertainty remain to do.
# Whole-year bootstrap is needed for uncertainty with overlapping windows;
# do not use an independent-observation Hessian as the final covariance estimate.

# ---- 5. Probability attribution: 2021 versus a 1951-like GMST climate ----
# Keep the modeling thresholds, target temperatures and dependence model fixed.
# Only GMST changes. This is an observed-GMST climate contrast, not a complete
# separation of anthropogenic forcing from natural influences.
factual_gmst <- gmst$gmst[match(factual_year, gmst$year)]
reference_gmst <- gmst$gmst[match(reference_year, gmst$year)]
stopifnot(length(factual_gmst) == 1L, is.finite(factual_gmst),
          length(reference_gmst) == 1L, is.finite(reference_gmst))
attribution <- compare_bivariate_pot(fit, target, factual_gmst, reference_gmst)
names(attribution)[names(attribution) == "counterfactual_gmst"] <- "reference_gmst"
names(attribution)[names(attribution) == "p_counterfactual"] <- "p_reference"
attribution$factual_year <- factual_year
attribution$reference_year <- reference_year
attribution$comparison <- paste0(factual_year, " vs ", reference_year, "-like")
# Conditional PRs must equal the joint PR divided by the corresponding marginal PR.
rr <- setNames(attribution$PR, attribution$quantity)
stopifnot(all(is.finite(attribution$PR)),
  abs(rr["p_B_given_A"] - rr["p_joint"] / rr["p_A"]) < 1e-10,
  abs(rr["p_A_given_B"] - rr["p_joint"] / rr["p_B"]) < 1e-10)
attribution
# Point estimates only. No confidence intervals have been calculated.
# Constant dependence does not imply constant conditional probabilities.
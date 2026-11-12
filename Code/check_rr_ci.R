# Run section by section in RStudio from the Code directory.
project_root <- "C:/Users/siyan/WCRP/wcrp_common_attribution_assessment"
code_dir <- file.path(project_root, "Code")
confidence_level <- 0.95
datasets_to_check <- c("CPC", "ERA5")

source(file.path(code_dir, "helpers", "gev.R"))
source(file.path(code_dir, "helpers", "workflow.R"))
gmst <- read.csv(file.path(project_root, "data/processed/gmst_smoothed.csv"))
gmst_factual <- gmst$gmst[gmst$year == 2021]
gmst_past <- gmst$gmst[gmst$year == 1951]
rr_ci_results <- list()

for (dataset in datasets_to_check) {
  data <- read.csv(file.path(project_root, "Results", dataset,
                             "annual_extremes_and_gmst.csv"))
  data <- data[data$used_in_fit, ]
  fit <- fit_gev(data$tmax, data$gmst)
  print(dataset)
  print(attribution(fit, gmst_factual, gmst_past, data$tmax))
  delta <- rr_ci_delta(fit, data, gmst_factual, gmst_past,
                       level = confidence_level)
  lrt <- rr_ci_lrt(fit, data, gmst_factual, gmst_past,
                   level = confidence_level)
  print(delta)
  print(lrt$interval)
  rr_ci_results[[dataset]] <- cbind(
    dataset = dataset, delta,
    lrt$interval[, c("RR_lrt_low", "RR_lrt_high",
                     "lrt_low_status", "lrt_high_status")]
  )
  stopifnot(delta$delta_status == "ok",
            delta$RR_delta_low < delta$RR,
            delta$RR_delta_high > delta$RR,
            lrt$interval$lrt_low_status == "ok",
            lrt$interval$RR_lrt_low < delta$RR)
  # Successful endpoints must actually satisfy the LR cutoff.
  for (bound in c("RR_lrt_low", "RR_lrt_high")) {
    value <- lrt$interval[[bound]]
    if (is.finite(value)) {
      nearest <- which.min(abs(lrt$profile$log_RR - log(value)))
      stopifnot(abs(lrt$profile$LR[nearest] - lrt$cutoff) < 0.01)
    }
  }
  same <- rr_ci_lrt(fit, data, gmst_factual, gmst_factual)$interval
  stopifnot(same$RR_lrt_low == 1, same$RR_lrt_high == 1)
  if (is.infinite(lrt$interval$RR_lrt_high)) {
    boundary <- lrt$profile[is.infinite(lrt$profile$log_RR), ]
    stopifnot(nrow(boundary) == 1L, boundary$LR <= lrt$cutoff)
  }
}

# Inspect this table or export it with write.csv() to your chosen location.
rr_ci_summary <- do.call(rbind, rr_ci_results)
rownames(rr_ci_summary) <- NULL
print(rr_ci_summary)

# Check a well-identified simulated case with two finite LR endpoints.
set.seed(20261029)
simulated_data <- data.frame(gmst = seq(-1, 1, length.out = 600))
simulated_data$tmax <- qgev_local(
  runif(nrow(simulated_data)), 25 + simulated_data$gmst, 2, 0.1
)
simulated_fit <- fit_gev(simulated_data$tmax, simulated_data$gmst)
simulated_delta <- rr_ci_delta(simulated_fit, simulated_data, 1, -1,
                                threshold = 30)
simulated_lrt <- rr_ci_lrt(simulated_fit, simulated_data, 1, -1,
                          threshold = 30)
stopifnot(simulated_delta$delta_status == "ok",
          simulated_lrt$interval$lrt_low_status == "ok",
          simulated_lrt$interval$lrt_high_status == "ok",
          simulated_lrt$interval$RR_lrt_low < simulated_delta$RR,
          simulated_lrt$interval$RR_lrt_high > simulated_delta$RR)
for (bound in c("RR_lrt_low", "RR_lrt_high")) {
  target <- log(simulated_lrt$interval[[bound]])
  nearest <- which.min(abs(simulated_lrt$profile$log_RR - target))
  stopifnot(abs(simulated_lrt$profile$LR[nearest] - simulated_lrt$cutoff) < 0.01)
}
message("RR CI checks passed.")

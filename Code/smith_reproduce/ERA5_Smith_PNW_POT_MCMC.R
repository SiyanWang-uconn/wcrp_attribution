# Smith Table 0.2: Bayesian seasonal PP inference from the existing MLE snapshot.
# Run section by section in RStudio. Creating this file does not run any chain.
# No existing script, helper, or MLE output is modified. No bootstrap is used.

# ============================================================================
# 1. Paths, prior reconstruction, and sampler configuration
# ============================================================================
project_root <- "C:/Users/siyan/WCRP/wcrp_common_attribution_assessment"
source(file.path(project_root, "Code/smith_reproduce/helpers/smith_pnw_pot.R"))
source(file.path(project_root, "Code/smith_reproduce/helpers/smith_pnw_mcmc.R"))
mle_file <- file.path(project_root, "Results/ERA5_Smith_PNW_POT_MLE/smith_pnw_pot_mle.rds")
output_root <- file.path(project_root, "Results/ERA5_Smith_PNW_POT_MCMC")

# Start with a pilot to assess mixing and runtime. Do not treat a pilot table as
# a reproduction of rare probabilities near 39 degrees C.
# The paper reports 20 million iterations, half discarded, thin=10 for its GEV
# analysis; PP is described as similar without a full sampler specification.
run_profile <- "pilot"
n_chains <- 4L
if (run_profile == "pilot") {
  n_iter <- 50000L
  burn_in <- 25000L
  thin <- 10L
} else if (run_profile == "paper_length") {
  n_iter <- 20000000L
  burn_in <- 10000000L
  thin <- 10L
} else {
  stop("Choose pilot or paper_length.")
}
base_seed <- 20261117L
adapt_start <- 1000L
adapt_every <- 50L
freeze_after_burn <- FALSE
mixture_probability <- 0.05
covariance_epsilon <- 1e-6
checkpoint_every <- if (run_profile == "pilot") 10000L else 500000L
progress_every <- if (run_profile == "pilot") 5000L else 100000L
resume_chains <- TRUE

# The seasonal PP prior is not fully stated in the paper. Explicit default:
# uniform in six location coefficients, log(scale intercept), two raw seasonal
# scale coefficients and xi; |xi|<1; positive scale and valid fitted-tail support.
# This extends the stated GEV prior but is not uniquely implied by the paper.
# Alternative: flat_scale_coefficients, including its transformation Jacobian.
prior_measure <- "flat_log_scale_intercept"
# The article's caption calls pF/pC expected exceedances; its prose calls them
# probabilities. BOTH conventions are calculated and saved, with explicit labels.
table_basis <- "expected_events"
ratio_subset_for_comparison <- "positive_gmst"
min_positive_draws <- 100L
rhat_limit <- 1.01
ess_minimum <- 400

stopifnot(file.exists(mle_file), requireNamespace("posterior", quietly = TRUE),
          n_chains >= 2L, n_iter > burn_in, thin >= 1L,
          table_basis %in% c("expected_events", "at_least_one_probability"),
          ratio_subset_for_comparison %in% c("all_draws", "positive_gmst"))
output_dir <- file.path(output_root, paste(run_profile, prior_measure, sep = "_"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# 2. Load the saved MLE, thresholds, peaks and fixed prediction climates
# ============================================================================
# Do not source the MLE main script or rerun threshold selection/declustering.
# 2021 remains excluded from fitting. Factual/counterfactual states are exactly
# those in the existing MLE snapshot, including its 1.2-degree GMST difference.
mle <- readRDS(mle_file)
stopifnot(all(c("model", "fit_data", "peaks", "factual", "counterfactual") %in% names(mle)),
          max(mle$fit_data$year) <= 2020,
          abs(unique(mle$factual$gmst - mle$counterfactual$gmst) - 1.2) < 1e-8)
coordinates <- smith_am_coordinates(mle)
log_target <- smith_am_target(mle, prior_measure)
# Verify the accelerated daily-grid likelihood at the original MLE before MCMC.
prior_at_mle <- if (prior_measure == "flat_scale_coefficients") coordinates$center[7] else 0
stopifnot(abs(log_target(coordinates$center) - prior_at_mle + mle$model$checked_nll) < 1e-7)
input_key <- list(mle_md5 = unname(tools::md5sum(mle_file)), prior = prior_measure)

# ============================================================================
# 3. Run or resume independent adaptive-Metropolis chains
# ============================================================================
# This is the computationally expensive section. Run deliberately in RStudio.
# Checkpoints include RNG state, adaptation state, retained draws and input hash.
# Matching completed checkpoints are reused; incompatible ones are rejected.
chains <- vector("list", n_chains)
for (chain in seq_len(n_chains)) {
  message("Starting/resuming chain ", chain, "/", n_chains)
  chains[[chain]] <- smith_am_run(
    log_target, coordinates, n_iter, burn_in, thin, base_seed + chain,
    adapt_start, adapt_every, freeze_after_burn, mixture_probability,
    covariance_epsilon, checkpoint_every, progress_every,
    checkpoint_file = file.path(output_dir, sprintf("chain_%02d.rds", chain)),
    resume = resume_chains, input_key = input_key)
}

# ============================================================================
# 4. Assess parameter mixing before interpreting posterior tail predictions
# ============================================================================
parameter_diagnostics <- smith_am_parameter_diagnostics(chains)
parameter_diagnostics
acceptance <- data.frame(
  chain = seq_len(n_chains), seed = base_seed + seq_len(n_chains),
  total_acceptance = vapply(chains, function(z) z$accepted / n_iter, 0.0),
  sampling_acceptance = vapply(chains, function(z) z$accepted_sampling / (n_iter - burn_in), 0.0),
  retained_draws = vapply(chains, function(z) z$kept, integer(1)))
acceptance
parameter_mixing_ok <- all(is.finite(parameter_diagnostics$rhat)) &&
  all(parameter_diagnostics$rhat < rhat_limit) &&
  all(parameter_diagnostics$ess_bulk >= ess_minimum) &&
  all(parameter_diagnostics$ess_tail >= ess_minimum)
if (!parameter_mixing_ok)
  warning("Parameter diagnostics do not meet the configured criteria; all tables are preliminary.")

# ============================================================================
# 5. Predict each retained posterior draw at 34-39 C and the observed event
# ============================================================================
# Store integrated event rates separately from 1-exp(-rate) probabilities.
# Structural zeros are determined in log space rather than by numerical underflow.
# The GMST coefficient is not constrained positive during sampling.
temperatures_C <- sort(unique(c(34:39, mle$settings$event_temperature_C)))
prediction <- smith_am_predict(chains, mle, temperatures_C)
posterior_tables <- smith_am_summarize(prediction, min_positive_draws)
posterior_summary <- posterior_tables$summary
generated_diagnostics <- posterior_tables$diagnostics
posterior_summary$parameter_mixing_ok <- parameter_mixing_ok
posterior_summary$run_profile <- run_profile
posterior_summary$prior_measure <- prior_measure

# ============================================================================
# 6. Build a clearly labeled comparison with the published Table 0.2
# ============================================================================
# E(pF) and E(pC) always use all draws. Ratio summaries are saved BOTH over all
# draws and conditional on positive GMST slope. Conditioning is not a prior.
# Author convention: inverse RR=0 when both probabilities are structurally zero.
# Its affected fraction is reported; an unconditional mean may legitimately be Inf.
# Average draw-wise ratios, not the ratio of posterior mean probabilities.
comparison_rows <- posterior_summary[
  posterior_summary$basis == table_basis &
    posterior_summary$ratio_subset == ratio_subset_for_comparison &
    posterior_summary$temperature_C %in% 34:39, ]
table_0_2_comparison <- merge(comparison_rows, smith_am_paper_table(),
                              by = "temperature_C", all.x = TRUE, sort = TRUE)
table_0_2_comparison

# Inspect generated_diagnostics for rare tail probabilities AND zero/nonzero
# support indicators. Parameter convergence alone does not establish accurate
# tail means. A constant-zero sample gives undefined convergence diagnostics,
# not proof that the posterior probability is exactly zero.
generated_diagnostics

# ============================================================================
# 7. Save new numerical summaries without changing MLE files
# ============================================================================
settings <- list(
  mle_file = mle_file, input_key = input_key, run_profile = run_profile,
  n_chains = n_chains, n_iter = n_iter, burn_in = burn_in, thin = thin,
  base_seed = base_seed, adapt_start = adapt_start, adapt_every = adapt_every,
  freeze_after_burn = freeze_after_burn, mixture_probability = mixture_probability,
  covariance_epsilon = covariance_epsilon, prior_measure = prior_measure,
  table_basis = table_basis, ratio_subset_for_comparison = ratio_subset_for_comparison,
  ratio_zero_convention = "pC/pF=0 for structural 0/0, following the author",
  rhat_limit = rhat_limit, ess_minimum = ess_minimum,
  min_positive_draws = min_positive_draws,
  prior_limitation = "Seasonal PP prior reconstructed; not fully specified in the paper")
write.csv(parameter_diagnostics, file.path(output_dir, "parameter_diagnostics.csv"), row.names = FALSE)
write.csv(acceptance, file.path(output_dir, "acceptance_rates.csv"), row.names = FALSE)
write.csv(generated_diagnostics, file.path(output_dir, "generated_quantity_diagnostics.csv"), row.names = FALSE)
write.csv(posterior_summary, file.path(output_dir, "posterior_summaries_all_conventions.csv"), row.names = FALSE)
write.csv(table_0_2_comparison, file.path(output_dir, "table_0_2_comparison.csv"), row.names = FALSE)
saveRDS(list(settings = settings, mle_settings = mle$settings,
             parameter_diagnostics = parameter_diagnostics, acceptance = acceptance,
             generated_diagnostics = generated_diagnostics,
             posterior_summary = posterior_summary,
             comparison = table_0_2_comparison, prediction = prediction),
        file.path(output_dir, "posterior_predictions_and_summary.rds"))
# Full parameter draws remain in the per-chain checkpoint files.


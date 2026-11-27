# Smith Section 0.6.1 GEV posterior analysis, using the stated uniform GEV prior.
# Run sections deliberately in RStudio; creating the file starts no chains.

# ---- 1. Paths, chain settings and interpretation settings ----
project_root <- "C:/Users/siyan/WCRP/wcrp_common_attribution_assessment"
smith_dir <- file.path(project_root, "Code/smith_reproduce")
source(file.path(smith_dir, "helpers/smith_pnw_mcmc.R"))
source(file.path(smith_dir, "helpers/smith_pnw_gev.R"))
mle_file <- file.path(project_root, "Results/ERA5_Smith_PNW_GEV_MLE/smith_pnw_gev_mle.rds")
output_root <- file.path(project_root, "Results/ERA5_Smith_PNW_GEV_MCMC")
fit_case <- "exclude_2021"
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
} else stop("Choose pilot or paper_length.")
base_seed <- 20261118L
adapt_start <- 1000L
adapt_every <- 50L
freeze_after_burn <- FALSE
mixture_probability <- 0.05
covariance_epsilon <- 1e-6
checkpoint_every <- if (run_profile == "pilot") 10000L else 500000L
progress_every <- if (run_profile == "pilot") 5000L else 100000L
resume_chains <- TRUE
rhat_limit <- 1.01
ess_minimum <- 400

stopifnot(file.exists(mle_file), requireNamespace("posterior", quietly = TRUE),
          fit_case %in% c("exclude_2021", "include_2021"), n_chains >= 2L)
output_dir <- file.path(output_root, paste(fit_case, run_profile, sep = "_"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 2. Stated GEV prior and saved annual-data likelihood ----
# Uniform beta0, betaG, log(sigma), xi; |xi|<1. No positive-betaG prior.
# An affine centered intercept has unit Jacobian, so centering preserves the prior.
snapshot <- readRDS(mle_file)
model <- snapshot$fits[[fit_case]]
coordinates <- smith_gev_coordinates(model)
log_target <- smith_gev_target(model)
stopifnot(abs(log_target(coordinates$center) + model$checked_nll) < 1e-7)
input_key <- list(mle_md5 = unname(tools::md5sum(mle_file)), fit_case = fit_case,
                  prior = "Uniform beta0,betaG,log_sigma,xi; abs(xi)<1")

# ---- 3. Adaptive-Metropolis chains using the same generic sampler as PP ----
# MLE covariance initializes proposals; annual observations determine the posterior.
# Rejections are retained. Checkpoint input/configuration must match for resuming.
# The paper_length option is per chain, not the sum across all chains.
chains <- vector("list", n_chains)
for (chain in seq_len(n_chains)) {
  chains[[chain]] <- smith_am_run(log_target, coordinates, n_iter, burn_in, thin,
    base_seed + chain, adapt_start, adapt_every, freeze_after_burn,
    mixture_probability, covariance_epsilon, checkpoint_every, progress_every,
    checkpoint_file = file.path(output_dir, sprintf("chain_%02d.rds", chain)),
    resume = resume_chains, input_key = input_key)
}

# ---- 4. Parameter convergence and Monte Carlo uncertainty ----
parameter_diagnostics <- smith_gev_diagnostics(chains)
acceptance <- data.frame(chain = seq_len(n_chains),
  total_acceptance = vapply(chains, function(z) z$accepted / n_iter, 0.0),
  sampling_acceptance = vapply(chains, function(z) z$accepted_sampling / (n_iter - burn_in), 0.0),
  retained_draws = vapply(chains, function(z) z$kept, integer(1)))
parameter_diagnostics
acceptance
check_values <- as.matrix(parameter_diagnostics[, c("rhat", "ess_bulk", "ess_tail")])
parameter_mixing_ok <- all(is.finite(check_values)) &&
  all(parameter_diagnostics$rhat < rhat_limit) &&
  all(parameter_diagnostics$ess_bulk >= ess_minimum) &&
  all(parameter_diagnostics$ess_tail >= ess_minimum)
if (!parameter_mixing_ok) warning("GEV parameter diagnostics do not pass; outputs are preliminary.")

# ---- 5. Posterior annual event probabilities and endpoint probabilities ----
# Evaluate the paper's temperatures AND this dataset's observed 2021 maximum.
temperatures <- sort(unique(c(snapshot$settings$paper_previous_record_C,
  snapshot$settings$paper_event_temperature_C, snapshot$settings$observed_event_temperature_C)))
prediction <- smith_gev_posterior_predict(chains, model,
  snapshot$settings$g_factual, snapshot$settings$g_counterfactual, temperatures)
tables <- smith_gev_posterior_summary(prediction)
posterior_summary <- tables$summary
generated_diagnostics <- tables$diagnostics
endpoint_summary <- smith_gev_endpoint_summary(prediction, temperatures)
posterior_summary$parameter_mixing_ok <- parameter_mixing_ok
posterior_summary$fit_case <- fit_case
posterior_summary$run_profile <- run_profile
posterior_summary
endpoint_summary
# Means use all draws. Inverse ratios are summarized with all draws and conditional
# on positive betaG; structural 0/0 is set to zero under the author's convention.
# An all-zero tail sample is not proof of zero posterior event probability.

# ---- 6. Compare with Section 0.6.1 text, not POT Table 0.2 ----
# References below concern the excluding-2021 model; never impose them on MCMC.
paper_reference <- data.frame(temperature_C = c(33.98, 39.6),
  paper_E_pF = c(0.095, 3.1e-6), paper_E_pC = c(0.0043, 1.8e-7),
  paper_E_inverse_RR_given_positive_gmst = c(0.0055, 7.7e-5),
  paper_q0_025 = c(0, NA_real_), paper_q0_975 = c(0.06, NA_real_),
  paper_Pr_inverse_RR_lt_0_01 = c(0.93, 0.9991))
comparison <- merge(posterior_summary[posterior_summary$ratio_subset == "positive_gmst", ],
                     paper_reference, by = "temperature_C", all.x = TRUE)
comparison$paper_reference_applicable <- fit_case == "exclude_2021"
comparison
# At 39.6 C the paper reports Pr(endpoint>x)=0.0088 factual and 0.0025
# counterfactual. Compare endpoint_summary; these are different from E(p).

# ---- 7. Save numerical summaries; chain files retain all parameter draws ----
settings <- list(mle_file = mle_file, input_key = input_key, fit_case = fit_case,
  run_profile = run_profile, n_chains = n_chains, n_iter = n_iter,
  burn_in = burn_in, thin = thin, base_seed = base_seed,
  adapt_start = adapt_start, adapt_every = adapt_every, freeze_after_burn = freeze_after_burn,
  mixture_probability = mixture_probability, covariance_epsilon = covariance_epsilon,
  prior = "Uniform beta0,betaG,log_sigma,xi; abs(xi)<1",
  ratio_zero_convention = "Author convention: structural 0/0 inverse ratio set to zero",
  probability_basis = "Annual maximum exceedance probability")
write.csv(parameter_diagnostics, file.path(output_dir, "parameter_diagnostics.csv"), row.names = FALSE)
write.csv(acceptance, file.path(output_dir, "acceptance_rates.csv"), row.names = FALSE)
write.csv(generated_diagnostics, file.path(output_dir, "generated_quantity_diagnostics.csv"), row.names = FALSE)
write.csv(posterior_summary, file.path(output_dir, "posterior_summary.csv"), row.names = FALSE)
write.csv(endpoint_summary, file.path(output_dir, "posterior_endpoint_summary.csv"), row.names = FALSE)
write.csv(comparison, file.path(output_dir, "gev_paper_comparison.csv"), row.names = FALSE)
saveRDS(list(settings = settings, mle_settings = snapshot$settings,
  prediction = prediction, posterior_summary = posterior_summary,
  endpoint_summary = endpoint_summary, comparison = comparison,
  parameter_diagnostics = parameter_diagnostics, generated_diagnostics = generated_diagnostics),
  file.path(output_dir, "smith_gev_posterior_summary.rds"))

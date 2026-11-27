# Smith Table 0.2: MCMC reconstruction

Use the new `ERA5_Smith_PNW_POT_MCMC.R` main script section by section in RStudio.
It reads the existing `smith_pnw_pot_mle.rds`; it does not rerun or edit MLE code.
It sources the existing PP definitions and the new `helpers/smith_pnw_mcmc.R`.
Only the posterior R package is additionally required; it is already installed.

## What the paper states and what remains unspecified

Section 0.6.1 states a uniform prior in the GEV location coefficients, log scale
and shape, restricted to absolute shape below one; it describes adaptive Metropolis
with 20 million iterations, half discarded and every tenth retained. Section
0.6.2 says the seasonal PP Bayesian analysis follows similar lines, without a
complete PP prior or MCMC implementation. Those GEV settings are therefore a
transparent reconstruction for PP, not a verified exact PP configuration.

The default PP prior is uniform in six location coefficients, log(scale
intercept), the two raw seasonal scale coefficients, and shape, with |shape|<1.
The scale intercept is positive; the seasonal scale and fitted-threshold support
must be positive on all fitted days. This is an improper flat prior; posterior
propriety is not established here. Exact agreement with the published table
cannot be claimed from the paper alone. Check for parameter drift and perform
prior sensitivity before interpreting very small tail means. The alternative
flat raw scale-coefficient measure is available, with its required Jacobian.
The MLE covariance only whitens/initializes proposals; it is not a prior or a
substitute for posterior sampling. No prior forces the GMST coefficient positive.

## Algorithm and workflow

Independent chains use dispersed feasible starts, online covariance adaptation,
a Gaussian random walk scaled by 2.38^2/d, and a 5% fixed small-step component.
Adaptation continues by default; it can be frozen after burn-in. The daily-grid
likelihood matches the existing PP likelihood and is verified at the MLE.
All rejections remain as repeated states. Checkpoints preserve RNG and adaptation
state; resuming requires identical input hash and sampler configuration.

The default pilot is four chains of 50,000 iterations, with 25,000 discarded and
thin=10. It is for runtime/mixing checks, not publication-quality rare-tail means.
`paper_length` selects 20 million iterations, half discarded, thin=10 PER CHAIN;
with four chains this is 80 million iterations and can be very expensive.
Nothing runs automatically until you execute the relevant main-script section.

## Table conventions and diagnostics

Table 0.2's caption says expected exceedances, while its text uses probabilities.
Both expected event counts and at-least-one-event probabilities are saved.
The default comparison uses expected counts; switch `table_basis` explicitly.
Means of factual/counterfactual quantities use all draws. Inverse-ratio statistics
are supplied both unconditionally and conditional on positive GMST slope. The
latter follows the conditioning explicitly described for the preceding GEV
example, but its exact reuse for every PP table column is not fully specified.
The author maps structural 0/0 inverse ratios to zero; affected fractions are
reported. Negative-slope draws are not dropped from chains. Infinite unconditional
ratios are not truncated. Draw-wise ratio means are never replaced by ratios of
means. Quantiles use type 1 to preserve atoms at zero and infinite values.

Rank-normalized Rhat, bulk/tail ESS and MCSE are reported using posterior for both
parameters and generated quantities. Constant-zero tail samples cannot establish
mixing or prove zero posterior risk. Check nonzero-support indicators, supporting
draw counts, chain-to-chain agreement and sensitivity to runtime/prior choices.
No result is declared reproduced based only on proximity to the published table.
Existing ERA5/GMST, threshold harmonics, seasonal window and declustering choices
retain the differences documented for the MLE stage. No old files are modified.

References: Smith, revised April 18 2025, Sections 0.6.1-0.6.2, Table 0.2;
Haario, Saksman and Tamminen (2001), Bernoulli 7, 223-242;
https://mc-stan.org/posterior/reference/diagnostics.html.


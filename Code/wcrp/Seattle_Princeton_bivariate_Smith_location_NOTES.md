# Station bivariate extension with Smith-inspired location parameters

## Scope and execution

Open `Seattle_Princeton_bivariate_Smith_location.R` in RStudio and run numbered
sections interactively. Configuration is at the top. The script uses existing
paired Seattle-Princeton JJA three-day means, not regional ERA5. Existing code
and outputs are preserved. New outputs go to
`Results/Seattle_Princeton_bivariate_Smith_location/` when section 6 is executed.
Implementation and numerical checks do not establish a completed station fit.

## Model and relationship to Smith

For each station j, mu_j is linear in centered annual GMST with annual and
semiannual seasonal harmonics. Sigma_j has an identity link, a positive baseline
and annual harmonics; xi_j is constant. These covariate structures follow the
local Smith seasonal PP reconstruction. Harmonic orders are configurable.

Define t_j(x) = [1 + xi_j (x - mu_j) / sigma_j]^(-1/xi_j), with its Gumbel limit.
The window margin is F_j(x) = exp[-t_j(x)/365]. The joint CDF is
H(x,y) = exp[-{a^(1/alpha) + b^(1/alpha)}^alpha], where
a = -log F_A(x), b = -log F_B(y), and 0 < alpha < 1.
Alpha is constant across dates and GMST. Independence is a limiting case;
the current optimizer fits the interior and does not select the boundary model.
The logistic family assumes asymptotic dependence for interior alpha; its tail
adequacy must be assessed rather than assumed from optimizer convergence.

The factor 365 fixes the PP-style parameterization. It does not turn overlapping
three-day windows into independent daily events, and it does not justify annual
return periods. This is a Smith-inspired censored window model, not a reproduction
of Smith's declustered Poisson likelihood or an annual-maxima GEV fit.

The four likelihood contributions are H(u_A,u_B), H_x(x,u_B), H_y(u_A,y), and
H_xy(x,y). Subthreshold observations are censored. All complete paired windows
are retained; separate marginal declustering would destroy synchronous pairing.
Time dependence makes the product a working likelihood. No Hessian standard
errors, bootstrap intervals, Bayesian inference or annual risks are implemented.

## Event, training and climate comparison

Default training is 1950-2020, excluding the target year. Years require at least
90% complete paired windows out of 90 possible JJA forward three-day windows.
The modeling thresholds remain the earlier fixed 1961-1990 station percentiles;
this version does not copy Smith's changing quantile-regression thresholds.
Threshold adequacy across season and GMST needs review. Missing pairs are dropped
without infilling; missing GMST for retained observations stops execution.

The default target uses observed June 26-28, 2021 station three-day temperatures.
Set target_mode to "baseline" to evaluate conditional exceedances at the baseline
thresholds instead. Predictions hold the calendar phase at June 26 and compare
2021 GMST with 1951 GMST. They are date-specific paired-window probabilities,
not probabilities averaged over summer. The reference is a historical-GMST
contrast, not an anthropogenic-free counterfactual.

Outputs retain both marginal risks, joint risk, both conditional risks, PMF and
climate-state probability ratios. Constant alpha does not imply constant
conditional probabilities. Zero-probability conditioning events yield undefined
conditional probabilities; zero/zero ratios are not assigned an artificial value.

## Review before scientific use

This 21-parameter nonstationary combination is a proposed implementation, not a
model reproduced from a single reference. Algebra checks do not validate the
mapping from Smith's PP event process to paired-window probabilities. Establish
that mapping and justify the window likelihood before treating this version as
the selected Smith-based bivariate approach.

Inspect optimizer starts, support and seasonal scale diagnostics. Seasonal
harmonics on JJA-only data may be weakly identified; the default complexity is
not a validated model choice. Simplified harmonics are available for sensitivity
analysis. Joint optimization re-estimates both margins rather than copying ERA5
coefficients. Whole-year refitting bootstrap, threshold sensitivity and tail
dependence diagnostics belong to a later authorized stage.

## Reference provenance

The earlier chat "研究双变量分析" identified the classical bivariate threshold
route through evd::fbvpot, with Joe, Smith and Weissman (1992), *Bivariate
Threshold Methods for Extremes*, JRSS B 54, 171-183, as a point-process reference.
The package documentation specifically attributes its censored likelihood to
Coles (2001), Section 8.3.1, and its Poisson alternative to Coles and Tawn (1991)
and Joe et al. (1992). These are related likelihoods, not interchangeable recipes.
https://search.r-project.org/CRAN/refmans/evd/html/fbvpot.html

The same chat separately discussed strict multivariate GP theory: Rootzen and
Tajvidi (2006); Rootzen, Segers and Wadsworth (2018); and Kiriliouk, Rootzen,
Segers and Wadsworth (2019), *Peaks Over Thresholds Modeling With Multivariate
Generalized Pareto Distributions*, Technometrics 61, 123-135.
https://doi.org/10.1080/00401706.2018.1462738
That conditional multivariate GP likelihood is not the joint window CDF
likelihood implemented here or in the earlier custom helper.

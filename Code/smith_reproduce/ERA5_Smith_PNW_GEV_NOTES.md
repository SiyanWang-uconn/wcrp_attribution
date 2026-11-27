# Smith PNW GEV reconstruction

Run ERA5_Smith_PNW_GEV_MLE.R, then ERA5_Smith_PNW_GEV_MCMC.R section by section
in RStudio. Implementation and the authorized 2026-11-18 MLE/MCMC pilot are now complete.
Existing regional ERA5 daily Tmax and prepared GMST are used, not the exact KNMI
releases in Smith, revised April 18 2025, Section 0.6.1.

## MLE: Table 0.1 and Figure 1 numerical diagnostics

Use full-calendar annual TXx. Fit 1950-2020 and 1950-2021 separately. Location is
linear in GMST; scale and shape are constant. ismev::gev.fit estimates both models.
GMST is centered internally, then beta0 and its covariance are transformed back
for comparison with the paper. Multi-start audits and independent likelihood
checks are retained. Outputs include endpoints with delta-method errors, annual
risks and intensity shifts. The endpoint is not constrained above the 2021 event.

Figure 1 QQ coordinates and KS/Cramer-von Mises/Anderson-Darling statistics are
saved numerically, without drawing figures. Simulated GOF p-values are not claimed
reproduced: the author says the calibration will be described elsewhere.
Parameter-estimation effects preclude naive fixed-parameter p-values. Hessian
inference is approximate, especially when the shape is in a nonregular region.

## Bayesian: Figure 2 and Section 0.6.1 text

Default input is the excluding-2021 GEV MLE. The prior follows the stated flat
measure in beta0, betaG, log(sigma), xi, restricted to |xi|<1. No Gaussian prior or
positive-GMST constraint is added. This unbounded flat prior is improper; following
the article does not itself prove posterior propriety. MLE covariance only
whitens proposals; the original annual observations determine the likelihood.
The generic adaptive-Metropolis/checkpoint/diagnostics implementation is reused
from the PP reconstruction without changing that helper.

Default pilot: four chains of 50,000 iterations, half discarded and thin=10.
paper_length: 20 million iterations PER CHAIN, half discarded, thin=10. The exact
proposal/adaptation settings are not supplied by Smith; sampler identity is not
claimed. Run-budget and convergence assessment are separate execution decisions.

Outputs include annual exceedance probabilities, inverse ratios both unconditional
and conditional on positive betaG, endpoints, endpoint-above-temperature
probabilities, Rhat/ESS/MCSE and numerical paper references. Endpoint exceedance
probability differs from posterior mean event probability. The author's 0/0
inverse-ratio convention is explicit and its affected fraction is reported.
GEV results are compared with Section 0.6.1, not the PP Bayesian Table 0.2.
Parameter convergence in a pilot does not establish accurate rare-event means.

## Inputs, outputs, and validation

Inputs remain under outputs/pnw_heatwave and data/processed. New output locations:
Results/ERA5_Smith_PNW_GEV_MLE and Results/ERA5_Smith_PNW_GEV_MCMC.
Existing PP result/checkpoint paths are unchanged. Required packages: ncdf4,
ismev and posterior. No package is installed automatically.

Validation on 2026-11-18: all Smith R files parse; source targets exist; analytical
GEV likelihood matches evd; centering, prior/support, endpoint, and a small
artificial GEV fit were checked. Existing PP checkpoint input hashes and proposal
coordinate signatures still match. No formal ERA5 GEV fit or posterior run occurred.

## Authorized saved run: 2026-11-18

Both MLE cases and the default four-chain pilot have now run. See WORKFLOW.md
and Results/ERA5_Smith_PNW_GEV_MCMC/exclude_2021_pilot/reproduction_review.md.
Parameter mixing meets the configured criteria, but rare-tail sampling is limited
and the posterior values do not match Smith exactly. No paper-length chain or
figure generation was performed. The earlier validation paragraph records the
implementation-stage checks, not the later run status.


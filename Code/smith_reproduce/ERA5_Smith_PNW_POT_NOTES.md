# Smith PNW POT: MLE reproduction stage

Run `ERA5_Smith_PNW_POT_MLE.R` section by section in RStudio. It uses only the new
`Code/smith_reproduce/helpers/smith_pnw_pot.R` helper. No earlier scripts/helpers are modified.
The script creates numerical outputs only when you execute it.

## Source and data

Richard L. Smith, *Detection and Attribution of Extreme Weather Events*, revised
April 18, 2025, Section 0.6.2 (PDF pages 26-28; printed pages 24-26).
Use the existing regional ERA5 daily Tmax NetCDF and smoothed GMST CSV.
These differ from the historical KNMI series and spatial bounds cited by Smith;
matching the paper's numerical estimates exactly is not guaranteed.
Fit 1950-2020; use 2021 only for prediction and the observed evaluation target.

## Method implemented

- Quantile regression at tau=0.90, using seasonal harmonics and a linear year trend.
- Ferro-Segers extremal index and automatic intervals clustering with texmex.
- Actual temperature maximum retained within each cluster.
- ismev PP fit with GMST, annual and semiannual harmonics in location; annual
  harmonics in scale with identity link; constant shape.
- Full daily exposure retained while nonpeak exceedances are censored.
- Historical observed/predicted cluster counts at and above 30 degrees C.
- Predictions under GMST(2021) and GMST(2021)-1.2 degrees C.
- Daily/seasonal upper endpoints and MLE probabilities, RR and FAR.

## Explicit reconstruction choices

The paper does not publish its exact threshold formula, harmonic order,
declustering treatment of winter gaps, complete code, or daily integration rule.
The default threshold has annual and semiannual harmonics and a linear calendar
year trend. EI is estimated on concatenated complete seasons; clusters cannot
cross years. Cluster maxima use actual temperature rather than TX-u residual.
The 365-day phase removes leap-day displacement. pp.fit uses npy=365, and expected
counts use the matching sum of daily tail measures divided by 365. No additional
extremal-index multiplier is applied after declustering.

The paper calls the season May 15-September 15 but specifies days 135-255 for PP
fitting. Default `window_mode="paper_doy"` follows the latter, May 15-September 12
on the 365-day calendar. Choose `window_mode="calendar"` for May 15-September 15.
This discrepancy is intentionally not hidden.

Required packages: ncdf4, quantreg, texmex, ismev. Missing packages are reported;
none are installed automatically. Optimization uses three explicit initial shapes.
No fallback to an earlier fixed-threshold model is performed if optimization fails.

This stage has no MCMC, posterior inference, bootstrap, or graphs. The paper's
Table 0.2 requires a separate Bayesian stage. MLE 0/0 ratios are left undefined;
the paper's posterior zero-ratio convention is not silently applied to MLE.
Nested count-check Pearson components are descriptive, not independent tests.


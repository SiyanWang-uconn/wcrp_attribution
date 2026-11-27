# WCRP | 2021 Pacific Northwest heatwave

This directory contains the R code for the regional heatwave attribution and
Seattle–Princeton joint heat analyses. The [WCRP research workflow](https://github.com/SiyanWang-uconn/wcrp_attribution/blob/main/WORKFLOW.md)
records the latest research progress, open questions, and planned work.

## Research status

Snapshot from the research workflow updated on 2026-11-27.

| Workstream | Status |
| --- | --- |
| Regional and station inputs | Regional temperature, station, and smoothed GMST inputs are prepared. |
| Regional and paired-station analyses | Initial GEV, POT, cumulative heat, and empirical Seattle–Princeton results are available. The Smith-inspired station location model is implemented and awaits its first fit. |
| Smith GEV and seasonal point-process reproduction | MLE and Bayesian pilot outputs are saved. Numerical reproduction and rare-tail precision remain under review. |
| Next checks | Review MCMC chains individually and compare the GEV and POT upper bounds. GMST effects in scale and shape and goodness-of-fit calibration remain open for discussion. |

## Directory guide

| Directory | Contents |
| --- | --- |
| [`wcrp/`](wcrp/) | Regional temperature extremes and Seattle–Princeton analyses. |
| [`smith_reproduce/`](smith_reproduce/README.md) | Smith reproduction scripts, method notes, and [POT MLE report](smith_reproduce/documentation/Smith_POT_MLE_Reproduction_Theory_and_Results.md). |
| [`helpers/`](helpers/) | Shared fitting and diagnostics functions. |
| [`config/`](config/) | Paths and run settings for the regional scripts. |

Run analysis scripts section by section in RStudio from the WCRP checkout.
Set the paths at the top of the relevant script or in `config/paths.R`.
The scripts use inputs under `outputs/pnw_heatwave/` and `data/processed/` in
the WCRP repository. Saved fits and diagnostics are kept under `Results/`.

## Updating both copies

Maintain the code and this README in [WCRP `Code/`](https://github.com/SiyanWang-uconn/wcrp_attribution/tree/main/Code).
When a research milestone changes, update the WCRP research workflow and this
summary in the same change. Commit the intended files by path and push to
WCRP `main`. The [Sync WCRP PNW code workflow](https://github.com/SiyanWang-uconn/wcrp_attribution/actions/workflows/sync-pnw-code.yml)
copies tracked code and documentation changes, including deletions, to
[`extreme-event-attribution/WCRP_PNW_heatwave_2021/`](https://github.com/SiyanWang-uconn/extreme-event-attribution/tree/main/WCRP_PNW_heatwave_2021).
The mirror excludes `.bak` files. After a successful run, use
`git pull --ff-only` in a local checkout of `extreme-event-attribution`
to refresh its copy.

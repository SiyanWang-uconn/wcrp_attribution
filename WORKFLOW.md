# WCRP research workflow and discussion notes

Updated: 2026-11-27

## 1. Research focus and current progress

The current case study is the 2021 Pacific Northwest heatwave, with regional
analyses and Seattle-Princeton joint heat as complementary branches.

Status vocabulary:
- **Implemented**: code exists; a complete run is not yet confirmed.
- **Outputs present**: saved results exist; scientific validation may remain.
- **Pilot**: exploratory analysis; further assessment is needed.
- **Review needed**: a specific methodological or numerical issue is unresolved.
- **Proposed**: an idea retained for discussion, not an agreed task.
- **Planned**: a selected next direction; work has not started.
- **Completed**: the stated task is finished, within its stated scope.

| Workstream | Status | Brief note |
| --- | --- | --- |
| Temperature data and GMST | Outputs present | Regional and station inputs prepared. |
| Initial GEV, POT and cumulative heat analyses | Outputs present; Pilot | Exploratory fits and summaries available. |
| Empirical bivariate analysis | Outputs present | Paired-window summaries available. |
| Nonstationary bivariate POT | Implemented; Pilot; Review needed | New station Smith-inspired location model implemented; earlier pilot remains under review. |
| Smith GEV/PP reproduction | Outputs present; Pilot; Review needed | MLE fits and Bayesian pilots saved; numerical reproduction remains incomplete. |
| Smith POT MLE theory and results report | Completed | English PDF and editable Markdown document the saved reproduction; no refitting. |
| Causal framing and upper-bound exploration | Planned | Research directions selected; advisor requests a GEV/POT upper-bound comparison. |
| Smith model extensions and assessment | Planned; Proposed | Advisor recommends GMST parameter extensions, chain-wise diagnostics and model selection; GoF remains a discussion item. |

## 2. Code and output map

| Location | Contents |
| --- | --- |
| `outputs/pnw_heatwave/` | Regional temperature and station inputs |
| `data/processed/` | Prepared GMST |
| `Code/wcrp/` | Initial regional and bivariate analyses |
| `Code/helpers/` | Shared routines for initial analyses |
| `Code/smith_reproduce/` | Smith GEV/PP scripts, helpers and method notes |
| `Results/` | Saved fits, diagnostics and comparisons |

## 3. Analysis branches and discussion points

The branches below include completed exploration and ideas retained for meeting
discussion. Open questions are not automatically commitments to execute work.

### A. Inputs and GMST

**Status:** Outputs present. Circulation extension: Proposed.

Regional ERA5, CPC and Berkeley inputs, station data, and smoothed GMST are
available. Check input consistency when comparing methods. Circulation analysis
is a possible extension; pressure data availability alone is not an analysis.

### B. Annual extremes and GEV

**Status:** Outputs present; Review needed.

Initial ERA5/CPC annual-extreme fits and bootstrap outputs are available.
Compare event definitions and climate states before interpreting differences.
The existing Figure 7 summary has a probability/return-period labeling issue.

### C. Regional peak POT / PP

**Status:** Outputs present; Review needed.

Initial analyses use declustered daily temperature peaks and GMST-dependent PP
models. Their event definition differs from three-day means. The earlier POT
approach needs reassessment alongside the Smith reconstruction.

### D. Cumulative heat magnitude

**Status:** Outputs present; Pilot.

Exploratory fits use heat accumulated above a temperature threshold within
events. Discussion: does this capture a useful aspect of event severity, and
is the magnitude-tail threshold appropriate?

### E. Seattle-Princeton empirical joint heat

**Status:** Outputs present.

Paired summer three-day means, exceedance frequencies and joint-dependence
summaries are available. Discussion: what joint event best serves the research
question? Overlapping windows must be accounted for in uncertainty estimates.

### F. Nonstationary bivariate POT

**Status:** Implemented; Pilot; Review needed.

A preliminary GMST-dependent model and attribution outputs exist. Its validity
is unresolved. Reconsider the formulation from the Smith reproduction before
further extension; retain both GEV and POT/PP as options.
The development target retains conditional co-exceedance probabilities in both
directions, their climate-state ratios, and joint and marginal probabilities.
Extend Smith's marginal climate response with a separately specified
tail-dependence model. Distinguish concurrent paired-window heat from same-season
componentwise maxima; exact event thresholds and model choice remain open.
The new station script `Code/wcrp/Seattle_Princeton_bivariate_Smith_location.R`
implements seasonal location-GMST margins with a censored window likelihood;
no station fit has been run for this version. See its method notes for the
distinction from Smith's declustered PP likelihood and remaining model checks.
The PP-to-window mapping and reference justification remain unresolved; the
21-parameter combination is a candidate, not an accepted theoretical extension.

### G. Smith seasonal PP: MLE

**Status:** Outputs present; Review needed.

Reconstruction fits and diagnostics are saved. The approach uses a changing
seasonal threshold, declustering and a nonstationary PP model. Agreement with
the source depends on unresolved input and implementation conventions.
Advisor recommendation (recorded 2026-11-27): consider GMST effects in scale ($\sigma$) and shape
($\xi$), in addition to location ($\mu$); assess the alternatives through model
selection, with `zhang2024leveraging` as the suggested method reference.

Documentation completed 2026-11-19; rendering updated 2026-11-27. The theory and
saved-result report, Markdown, native LaTeX source and builder are colocated in
`Code/smith_reproduce/documentation/`. The PDF uses native mathematical fonts;
the report covers threshold estimation, declustering, PP likelihood and saved
attribution results. No statistical models were rerun; scientific validation
remains under review.

### H. Smith seasonal PP: Bayesian pilot

**Status:** Pilot; Review needed.

Four-chain pilot outputs are available. Parameter diagnostics are satisfactory,
but rare-tail precision is insufficient. Discussion: which discrepancies need
explanation, and would prior checks or longer runs be informative?
Advisor-requested next check (recorded 2026-11-27): inspect trace plots and convergence separately
for each chain in both the PP and GEV pilots; this review is Planned.

### I. Smith annual GEV

**Status:** Outputs present; Pilot; Review needed.

MLE fits excluding and including 2021 and an excluding-2021 Bayesian pilot are
saved. MLE results are broadly close to Smith; Bayesian numerical discrepancies
and rare-tail uncertainty remain unresolved.
The matched excluding-2021 MLE comparison shows higher seasonal PP endpoints
than annual GEV endpoints, but both remain below the observed 2021 extreme.
This is a descriptive comparison, not evidence of superior GOF; see
`Results/ERA5_Smith_POT_GEV_comparison/endpoint_comparison_review.md`.
Advisor-requested update (recorded 2026-11-27): present the GEV/POT upper-bound comparison explicitly,
building on this saved comparison; preparation of the update is Planned.

### J. Smith input comparability

**Status:** Review needed.

The local 2021 temperature differs slightly from the paper's evaluation value.
Exact source-series and GMST conventions remain to be compared. A spatial-domain
difference has not been established from the quoted bounds alone.

### K. Daily temperature aggregation

**Status:** Proposed (unverified explanation).

Instantaneous versus interval-maximum temperature and daily window conventions
are possible sources of input differences. This is an unverified explanation,
not a demonstrated cause of the reproduction discrepancy.

### L. Goodness-of-fit calibration

**Status:** Completed (reference review); Proposed (calibration extension).

The KS reference highlights parameter-estimation and dependence effects on test
calibration. Current GEV outputs include statistics without naive KS p-values.
Advisor question (recorded 2026-11-27): "Do GoF?" Retain as Proposed; the assessment scope and method
remain to be decided. Refitting and dependence-aware calibration are relevant
considerations if GoF is pursued.

### M. Research directions for discussion

**Status:** Planned.

- **Causal framing:** "Given the weather situation like this one, how much hotter did human influence make this？"
  Develop the causal-inference framing around this question.
- **Bivariate development:** build on the Smith reproduction and reassess the
  previous POT approach, retaining conditional co-exceedance as the focus and
  keeping GEV and POT/PP options open.
- **Upper bounds:** read the literature on finite endpoints and uncertainty;
  explore whether the issue offers a methodological contribution. Distinguish
  statistical endpoints from physical limits. Add a GEV/POT upper-bound
  comparison to the next research update (advisor recommendation recorded
  2026-11-27).
- **Model development and selection** (advisor recommendation recorded
  2026-11-27): consider GMST in $\sigma$ and $\xi$,
  beyond the current $\mu$ response, and assess whether more covariates are
  warranted using the method in `zhang2024leveraging` (reference to be identified).

## 4. Interpretation reminders

Saved outputs establish that an analysis was run, not that its conclusions are
validated. Keep event metrics, fit periods and climate comparisons explicit.
A 1951-like climate state is not an anthropogenic-free counterfactual. Parameter
convergence alone does not establish rare-tail accuracy.

## 5. Current plans

| Direction | Status | Recorded / updated |
| --- | --- | --- |
| Weather-conditioned causal intensity attribution | Planned; framing to be developed | 2026-11-18 (plan recorded) |
| Smith-based bivariate GEV or POT/PP | Implemented; Review needed; station location-model extension not yet fitted | 2026-11-19 (implementation update) |
| Upper-bound literature exploration | Planned; review not started | 2026-11-18 (plan recorded) |
| GEV/POT upper-bound comparison for the research update | Planned; existing comparison available | 2026-11-27 (advisor recommendation recorded) |
| GMST effects in scale and shape | Proposed; advisor recommendation to consider | 2026-11-27 (advisor recommendation recorded) |
| GoF assessment | Proposed; scope and method to discuss | 2026-11-27 (advisor recommendation recorded) |
| MCMC trace plots and convergence for each chain | Planned; chain-wise review pending | 2026-11-27 (advisor recommendation recorded) |
| Model selection and additional covariates | Planned; identify and review `zhang2024leveraging` method | 2026-11-27 (advisor recommendation recorded) |


## 6. Document scope

Keep a brief record of research progress, open ideas and questions useful for
supervisor discussions. Retain the status vocabulary and brief status labels
in the progress table and analysis branches. Tentative ideas may be retained
without becoming tasks.
Detailed settings and diagnostics belong in method notes and result reviews;
routine operations and conversation-by-conversation logs are omitted.
Record dates beside substantive recommendations, decisions and task updates;
label recording dates separately from meeting or completion dates. Do not infer
missing dates or add a separate chronological log for routine edits.

# CHANGELOG


## v0.1.0 (2026-11-27)

### Build System

- Add repository Make targets and usage notes
  ([`f42e676`](https://github.com/SiyanWang-uconn/wcrp_attribution/commit/f42e676aa6896291e412f5fbf1d24584b8196706))

Add a default help target, locked Python setup, lint and formatting checks, R syntax validation,
  packaging, and explicit data and report entry points. Document the RStudio workflow and local PDF
  builder requirements.

### Chores

- Archive legacy EU heatwave scripts locally
  ([`9dffc37`](https://github.com/SiyanWang-uconn/wcrp_attribution/commit/9dffc37d6a31e37285f398a9fd581aa52398a2d5))

Remove 11 legacy temperature, GMST, and GEV scripts from the shared Code directory. Identical copies
  are retained in the ignored Code/previous_EU_heatwave directory of the current workspace.

- Keep local research artifacts out of Git
  ([`a8537ee`](https://github.com/SiyanWang-uconn/wcrp_attribution/commit/a8537ee9fa9c2029137fba54894517430504c4a5))

Ignore R and RStudio session state, the archived EU heatwave scripts, and new generated results.
  Share only the Smith Markdown report from its documentation directory. Existing tracked results
  remain in version control.

### Documentation

- Record attribution workflow and research priorities
  ([`4cce5b8`](https://github.com/SiyanWang-uconn/wcrp_attribution/commit/4cce5b8911c123facb7a7b7998c8afd4d6d01927))

Capture the status of the PNW research branches and distinguish implementations, saved runs, and
  reviewed findings. Record the causal framing question, endpoint comparison, and planned diagnostic
  and model-selection discussions.

- Share the Smith seasonal POT MLE report
  ([`f2f6f74`](https://github.com/SiyanWang-uconn/wcrp_attribution/commit/f2f6f7478e74e41e81948b52ded99364e7f9f695))

Document the seasonal threshold, declustering, point-process likelihood, and attribution estimates
  from saved runs. Include the editable Markdown report in the repository.

### Features

- Add Smith PNW GEV and point-process reproduction
  ([`7e29a6b`](https://github.com/SiyanWang-uconn/wcrp_attribution/commit/7e29a6b2ee7e14009b2c3a33027e9ce9011ad7b8))

Add annual GEV and seasonal POT point-process scripts for maximum likelihood fits and Bayesian
  pilots. Include shared likelihood and sampling helpers, RStudio run guidance, and notes on input
  conventions, checkpoint dependencies, and limitations of the saved pilots.

- Add Smith-inspired bivariate station model
  ([`1658986`](https://github.com/SiyanWang-uconn/wcrp_attribution/commit/16589869dfce7c477b61ebf9e977c024d2093239))

Implement seasonal GMST-dependent location margins and a censored paired-window likelihood for
  Seattle and Princeton. Document the distinction from a declustered point-process likelihood and
  the likelihood and model checks still open. The station fit has not been run.


## v0.0.0 (2026-11-12)

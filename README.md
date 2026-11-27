# wcrp-common-attribution-assessment

Python project scaffold managed with `uv`.

## Getting started

```bash
make setup
make hooks
```

Run `make` from the repository root to list the GNU Make targets. Use
`make check` for Python lint and formatting checks and R syntax validation.
If `Rscript` is outside `PATH`, set `RSCRIPT` to its executable path, for example
`make check-r "RSCRIPT=C:/Program Files/R/R-4.5.1/bin/Rscript.exe"` in PowerShell.

Run the analysis scripts section by section in RStudio. The data targets are
opt-in and may download inputs or replace generated files. The ERA5 targets
require access to the remote data source. `make report-smith` requires the local
Python report builder in `Code/smith_reproduce/documentation/` and `pdflatex`.
The Markdown report is available in that directory after cloning.

## Updating the PNW code mirror

Edit `Code/` here, stage the intended files by path, commit, and push to `main`.
The **Sync WCRP PNW code** workflow mirrors the tracked files to
[`WCRP_PNW_heatwave_2021/`](https://github.com/SiyanWang-uconn/extreme-event-attribution/tree/main/WCRP_PNW_heatwave_2021),
excluding `.bak` files. Make code changes in this repository. The local
checkout of `extreme-event-attribution` can then be refreshed with `git pull --ff-only`.
Check the workflow run in the **Actions** tab. Use **Run workflow** on `main`
when you need to sync the current snapshot again.

## Development

- Runtime and development dependencies are managed in `pyproject.toml` with `uv`.
- Code quality checks run through `pre-commit` using `ruff`.
- Versioning follows semantic versioning and conventional commits.

### Conventional commit examples

- `feat: add attribution data loader`
- `fix: correct baseline anomaly calculation`
- `feat!: remove deprecated processing path`
- `chore: update contributor guidance`

Changes merged into `main` trigger the GitHub release workflow, which calculates the next semantic version from the commit history and writes it back to `pyproject.toml`.

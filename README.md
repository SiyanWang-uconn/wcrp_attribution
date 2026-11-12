# wcrp-common-attribution-assessment

Python project scaffold managed with `uv`.

## Getting started

```bash
uv sync --group dev
uv run pre-commit install
```

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

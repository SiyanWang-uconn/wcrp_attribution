.DEFAULT_GOAL := help

UV ?= uv
RSCRIPT ?= Rscript
PYTHON_PACKAGE := wcrp_common_attribution_assessment
SMITH_REPORT_BUILDER := Code/smith_reproduce/documentation/build_pot_mle_report.py

.PHONY: help setup hooks check check-python check-r format build \
        data-era5 data-cpc data-berkeley data-ghcn data-seattle-princeton \
        data-pakistan pnw-event-definitions report-smith

help:
	@echo "setup                   Install the locked Python environment"
	@echo "hooks                   Install the Git pre-commit hook"
	@echo "check                   Check Python lint and formatting, and parse R scripts"
	@echo "check-python            Check Python lint and formatting"
	@echo "check-r                 Parse R scripts under Code without running them"
	@echo "format                  Apply Ruff formatting to src"
	@echo "build                   Build the Python source distribution and wheel"
	@echo "data-era5               Produce ERA5 PNW reference data"
	@echo "data-cpc                Produce CPC PNW reference data"
	@echo "data-berkeley           Produce Berkeley Earth PNW reference data"
	@echo "data-ghcn               Prepare the PNW GHCN station dataset"
	@echo "data-seattle-princeton  Extract Seattle and Princeton station series"
	@echo "data-pakistan           Run the Pakistan rainfall script"
	@echo "pnw-event-definitions   Calculate ERA5 PNW heatwave event definitions"
	@echo "report-smith            Render the Smith PDF using the local builder and pdflatex"

setup:
	"$(UV)" sync --locked --group dev

hooks:
	"$(UV)" run --no-sync pre-commit install

check: check-python check-r

check-python:
	"$(UV)" run --no-sync ruff check src
	"$(UV)" run --no-sync ruff format --check src

check-r:
	"$(RSCRIPT)" --vanilla -e "invisible(lapply(list.files('Code', pattern = '[.]R$$', recursive = TRUE, full.names = TRUE), function(path) parse(file = path)))"

format:
	"$(UV)" run --no-sync ruff format src

build:
	"$(UV)" build

data-era5:
	"$(UV)" run --no-sync python -m $(PYTHON_PACKAGE).pnw_heatwave_era5_reference_data

data-cpc:
	"$(UV)" run --no-sync python -m $(PYTHON_PACKAGE).pnw_heatwave_cpc_reference_data

data-berkeley:
	"$(UV)" run --no-sync python -m $(PYTHON_PACKAGE).pnw_heatwave_berkeley_reference_data

data-ghcn:
	"$(UV)" run --no-sync python -m $(PYTHON_PACKAGE).pnw_heatwave_ghcn_station_data

data-seattle-princeton:
	"$(UV)" run --no-sync python -m $(PYTHON_PACKAGE).pnw_heatwave_select_ghcn_stations

data-pakistan:
	"$(UV)" run --no-sync python -m $(PYTHON_PACKAGE).pakistan_rainfall

pnw-event-definitions:
	"$(UV)" run --no-sync python -m $(PYTHON_PACKAGE).pnw_heatwave

report-smith: $(SMITH_REPORT_BUILDER)
	"$(UV)" run --no-sync python $(SMITH_REPORT_BUILDER)

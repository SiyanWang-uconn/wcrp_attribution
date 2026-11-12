"""Retrieve CPC daily maximum temperatures and build PNW reference series.

Run this file directly, or run its sections interactively. Downloads use NOAA's
NetCDF Subset Service, so only the requested region is transferred. No account
or additional dependencies beyond this project's environment are required.
"""

from datetime import date
from io import BytesIO
from pathlib import Path

import numpy as np
import requests
import xarray as xr

# %% Configuration
START_YEAR = 1979
END_YEAR = 2025
PNW_LATITUDE_BOUNDS = (52.0, 45.0)  # North, south in degrees north.
PNW_LONGITUDE_BOUNDS = (237.0, 241.0)  # West, east in degrees east (0-360).
OUTPUT_DIR = Path(__file__).resolve().parents[2] / "outputs" / "pnw_heatwave"
CACHE_DIR = OUTPUT_DIR / "cpc_cache"
USE_CACHE = True
REQUEST_TIMEOUT = (30, 180)  # Connection and read timeouts, in seconds.
CPC_SUBSET_URL = "https://psl.noaa.gov/thredds/ncss/grid/Datasets/cpc_global_temp"
CPC_SOURCE_URL = "https://psl.noaa.gov/data/gridded/data.cpc.globaltemp.html"


# %% Retrieval and processing helpers
def select_domain(dataset, year):
    """Trim the subset service's rounded grid and time bounds precisely."""
    north, south = PNW_LATITUDE_BOUNDS
    west, east = PNW_LONGITUDE_BOUNDS
    tx = dataset["tmax"].sel(time=slice(f"{year}-01-01", f"{year}-12-31"))
    tx = tx.where(
        (tx.lat >= south) & (tx.lat <= north) & (tx.lon >= west) & (tx.lon <= east),
        drop=True,
    )
    if any(tx.sizes[dim] == 0 for dim in ("time", "lat", "lon")):
        raise ValueError(f"CPC subset for {year} contains no data in the requested domain.")
    expected_days = int(
        (np.datetime64(f"{year + 1}-01-01") - np.datetime64(f"{year}-01-01"))
        / np.timedelta64(1, "D")
    )
    expected_times = np.arange(np.datetime64(f"{year}-01-01"), np.datetime64(f"{year + 1}-01-01"))
    if tx.sizes["time"] != expected_days or not np.array_equal(tx.time.values, expected_times):
        raise ValueError(f"CPC data for {year} do not contain a complete daily calendar.")
    if tx.attrs.get("units") != "degC":
        raise ValueError(f"Unexpected CPC temperature units: {tx.attrs.get('units')!r}")
    missing_days = tx.isnull().all(dim=["lat", "lon"])
    if bool(missing_days.any()):
        dates = np.datetime_as_string(tx.time.values[missing_days.values], unit="D")
        print(
            f"WARNING: CPC {year} has {len(dates)} day(s) with no valid temperatures: "
            f"{', '.join(dates)}. Retaining these days as NaN; annual maxima use "
            "available days and may underestimate the full-year maximum.",
            flush=True,
        )
    return tx


def fetch_year(session, year):
    """Load one complete year, using a region-specific local cache when available."""
    north, south = PNW_LATITUDE_BOUNDS
    west, east = PNW_LONGITUDE_BOUNDS
    cache_path = CACHE_DIR / (f"tmax_{year}_lat{south:g}_{north:g}_lon{west:g}_{east:g}.nc")
    if USE_CACHE and cache_path.exists():
        print(f"Loading cached CPC data for {year}...", flush=True)
        with xr.open_dataset(cache_path, engine="scipy") as dataset:
            return select_domain(dataset, year).load()

    url = f"{CPC_SUBSET_URL}/tmax.{year}.nc"
    params = {
        "var": "tmax",
        "north": north,
        "south": south,
        "west": (west + 180) % 360 - 180,
        "east": (east + 180) % 360 - 180,
        "time_start": f"{year}-01-01T00:00:00Z",
        "time_end": f"{year}-12-31T00:00:00Z",
        "accept": "netcdf",
    }
    print(f"Downloading CPC PNW daily maximum temperatures for {year}...", flush=True)
    try:
        with session.get(url, params=params, timeout=REQUEST_TIMEOUT) as response:
            response.raise_for_status()
            with xr.open_dataset(BytesIO(response.content), engine="scipy") as dataset:
                tx = select_domain(dataset, year).load()
    except requests.RequestException as exc:
        raise RuntimeError(
            f"Could not retrieve CPC data for {year} from {url}. "
            "Check the network and data availability, then rerun; completed years are cached."
        ) from exc

    if USE_CACHE:
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        temporary_path = cache_path.with_suffix(".nc.part")
        tx.to_netcdf(temporary_path, engine="scipy")
        temporary_path.replace(cache_path)
    return tx


def domain_mean(tx):
    """Average daily maxima over valid grid cells using cosine-latitude weights."""
    weights = np.cos(np.deg2rad(tx.lat))
    return tx.weighted(weights).mean(dim=["lat", "lon"], skipna=True, keep_attrs=True)


# %% Build and save reference series
def main():
    """Retrieve each year and save daily and annual PNW temperature series."""
    if not 1979 <= START_YEAR <= END_YEAR < date.today().year:
        raise ValueError("Choose complete calendar years from 1979 through last year.")
    north, south = PNW_LATITUDE_BOUNDS
    west, east = PNW_LONGITUDE_BOUNDS
    if not -90 <= south < north <= 90 or not 0 <= west < east < 360:
        raise ValueError("Use north > south and 0 <= west < east < 360 for the PNW domain.")
    if west < 180 <= east:
        raise ValueError("The configured domain must not cross the 180-degree meridian.")

    daily_series = []
    with requests.Session() as session:
        for year in range(START_YEAR, END_YEAR + 1):
            tx_grid = fetch_year(session, year)
            daily_series.append(domain_mean(tx_grid))

    TX = xr.concat(daily_series, dim="time").sortby("time").rename("tmax")
    TXx = TX.resample(time="1YE").max(skipna=True, keep_attrs=True)
    TXx = TXx.assign_coords(
        valid_days=TX.resample(time="1YE").count(),
        expected_days=xr.ones_like(TX, dtype="int32").resample(time="1YE").sum(),
    )
    TXx.valid_days.attrs = {
        "long_name": "Number of days with a valid area-mean temperature in the year",
        "units": "days",
    }
    TXx.expected_days.attrs = {
        "long_name": "Number of calendar days expected in the year",
        "units": "days",
    }
    shared_attrs = {
        "source_dataset": "CPC Global Unified Temperature V1.0",
        "source_variable": "tmax (daily maximum temperature)",
        "source_url": CPC_SOURCE_URL,
        "source_subset_url": CPC_SUBSET_URL,
        "retrieval_date": date.today().isoformat(),
        "retrieval_date_note": "Execution date; cached data may have been retrieved earlier.",
        "spatial_domain": f"{south:g}-{north:g} N, {west:g}-{east:g} E",
        "spatial_weighting": "Cosine of latitude over valid grid cells within the bounds",
        "missing_data_handling": (
            "Missing grid cells excluded from spatial means; days with no valid cells "
            "retained as NaN. Annual maxima use available days without interpolation. "
            "Incomplete years may underestimate the full-year maximum; compare "
            "valid_days with expected_days in the annual output."
        ),
        "reference_period": f"{START_YEAR}-{END_YEAR}",
        "processing_notes": (
            "CPC tmax is already a daily maximum in degrees Celsius; no hourly "
            "resampling or Kelvin conversion was applied. Grid cell centres must lie "
            "within the configured bounds. CPC reports a 06 UTC to 06 UTC maximum period."
        ),
    }
    TX.attrs.update(
        {
            **shared_attrs,
            "units": "degC",
            "long_name": "PNW area-mean daily maximum temperature",
            "temporal_aggregation": "Daily maximum at each grid cell, then spatial mean",
            "cell_methods": "time: maximum within days area: mean",
        }
    )
    TXx.attrs.update(
        {
            **shared_attrs,
            "units": "degC",
            "long_name": "Annual maximum of PNW area-mean daily maximum temperature",
            "temporal_aggregation": "Annual maximum of the area-mean daily maximum temperature",
            "cell_methods": "time: maximum within days area: mean time: maximum over days",
        }
    )
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    daily_path = OUTPUT_DIR / "CPC_reference_daily_tmax.nc"
    daily_duplicate_path = OUTPUT_DIR / "CPC_reference_daily_daily_tmax.nc"
    annual_path = OUTPUT_DIR / "CPC_reference_tmax.nc"
    daily_output = TX.to_dataset().assign_coords(
        valid_days=TX.notnull().astype("int32"),
        expected_days=xr.ones_like(TX, dtype="int32"),
    )
    # Daily records use the same fields as the annual output, with daily counts.
    daily_output.valid_days.attrs = {
        "long_name": "Number of days with a valid area-mean temperature in the daily record",
        "units": "days",
    }
    daily_output.expected_days.attrs = {
        "long_name": "Number of calendar days expected in the daily record",
        "units": "days",
    }
    daily_output["tmax_date"] = TX.time.where(TX.notnull())
    daily_output.tmax_date.attrs = {
        "long_name": "Date of regional daily maximum temperature",
        "missing_data_handling": "Missing when the daily regional temperature is missing",
    }
    daily_output.tmax_date.encoding = {
        "units": "days since 1970-01-01",
        "dtype": "float64",
        "_FillValue": np.nan,
    }
    daily_output.to_netcdf(daily_path, engine="scipy")
    daily_output.to_netcdf(daily_duplicate_path, engine="scipy")
    annual_output = TXx.to_dataset()
    annual_output["tmax_date"] = TX.resample(time="1YE").map(
        lambda year: year.time.where(year == year.max(skipna=True)).min(skipna=True)
    )
    annual_output.tmax_date.attrs = {
        "long_name": "Date of annual maximum regional daily maximum temperature",
        "tie_handling": "Earliest date when multiple days share the annual maximum",
        "missing_data_handling": "Missing when the entire year has no valid temperatures",
    }
    annual_output.tmax_date.encoding = {
        "units": "days since 1970-01-01",
        "dtype": "float64",
        "_FillValue": np.nan,
    }
    annual_output.to_netcdf(annual_path, engine="scipy")
    print(f"Saved daily reference series: {daily_path}", flush=True)
    print(f"Saved identical daily reference series: {daily_duplicate_path}", flush=True)
    print(f"Saved annual maximum reference series: {annual_path}", flush=True)
    return TX, TXx


if __name__ == "__main__":
    TX, TXx = main()

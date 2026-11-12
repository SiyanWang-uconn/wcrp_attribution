"""Download Berkeley Earth daily TMAX and build PNW reference series.

Run directly or execute sections interactively. The source is an experimental
1-degree land dataset. Global files are grouped by decade (hundreds of MB per
file); completed downloads are cached. Requires h5netcdf and h5py in addition to
the project's dependencies. Source: https://berkeleyearth.org/data/
"""

from datetime import date
from pathlib import Path

import numpy as np
import requests
import xarray as xr

# %% Configuration
START_YEAR = 1979  # Set to 1880 to include the full historical record.
END_YEAR = 2023  # Last complete year in the verified November 2024 release.
PNW_LATITUDE_BOUNDS = (52.0, 45.0)  # North, south in degrees north.
PNW_LONGITUDE_BOUNDS = (-123.0, -119.0)  # West, east in degrees east (-180 to 180).
OUTPUT_DIR = Path(__file__).resolve().parents[2] / "outputs" / "pnw_heatwave"
CACHE_DIR = OUTPUT_DIR / "berkeley_cache"
REQUEST_TIMEOUT = (30, 180)  # Connection and read timeouts, in seconds.
SOURCE_URL = "https://berkeleyearth.org/data/"
DATA_BASE_URL = "https://berkeley-earth-temperature.s3.us-west-1.amazonaws.com/Global/Gridded"


# %% Download helpers
def download_decade(session, decade):
    """Download one global source file, reusing completed local downloads."""
    filename = f"Complete_TMAX_Daily_LatLong1_{decade}.nc"
    path = CACHE_DIR / filename
    if path.exists():
        print(f"Using cached Berkeley Earth file: {filename}", flush=True)
        return path

    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    temporary_path = path.with_suffix(".nc.part")
    print(f"Downloading {filename} (global decade file)...", flush=True)
    try:
        with session.get(
            f"{DATA_BASE_URL}/{filename}", stream=True, timeout=REQUEST_TIMEOUT
        ) as response:
            response.raise_for_status()
            total = int(response.headers.get("Content-Length", 0))
            if total:
                print(f"Download size: {total / 1024**2:.1f} MiB", flush=True)
            downloaded = 0
            next_progress = 64 * 1024**2
            with temporary_path.open("wb") as output:
                for chunk in response.iter_content(chunk_size=1024**2):
                    output.write(chunk)
                    downloaded += len(chunk)
                    if downloaded >= next_progress:
                        print(f"  Downloaded {downloaded / 1024**2:.0f} MiB", flush=True)
                        next_progress += 64 * 1024**2
            if total and downloaded != total:
                raise OSError(
                    f"Incomplete download: expected {total}, received {downloaded} bytes."
                )
    except requests.RequestException as exc:
        raise RuntimeError(
            f"Berkeley Earth download failed for {decade}. "
            "Rerun to reuse completed decade files; access errors may require checking "
            f"the current official download links at {SOURCE_URL}."
        ) from exc
    temporary_path.replace(path)
    return path


# %% Temperature reconstruction and area averaging
def process_decade(dataset):
    """Reconstruct absolute daily TMAX before calculating the regional mean."""
    required = {"temperature", "climatology", "year", "month", "day", "day_of_year"}
    missing = required.difference(dataset.variables)
    if missing:
        raise ValueError(f"Berkeley Earth file is missing variables: {sorted(missing)}")
    north, south = PNW_LATITUDE_BOUNDS
    west, east = PNW_LONGITUDE_BOUNDS
    region = dataset.sel(latitude=slice(south, north), longitude=slice(west, east))
    if not region.sizes["latitude"] or not region.sizes["longitude"]:
        raise ValueError("No Berkeley Earth grid cell centres fall within the configured bounds.")

    years = region.year.values.astype(int)
    indices = np.flatnonzero((years >= START_YEAR) & (years <= END_YEAR))
    if not indices.size:
        raise ValueError("The requested years are not present in the source file.")
    selected = region.isel(time=indices)
    years = selected.year.values.astype(int)
    months = selected.month.values.astype(int)
    days = selected.day.values.astype(int)
    times = np.array(
        [f"{y:04d}-{m:02d}-{d:02d}" for y, m, d in zip(years, months, days, strict=True)],
        dtype="datetime64[ns]",
    )
    selected = selected.assign_coords(time=times)

    # Use Berkeley's supplied climate-day mapping, including its leap-year handling.
    climate_days = selected.day_of_year.values
    if not np.all(np.isfinite(climate_days)) or not np.all(
        climate_days == climate_days.astype(int)
    ):
        raise ValueError("Berkeley Earth day_of_year contains invalid climatology indices.")
    indices = climate_days.astype(int) - 1
    if np.any(indices < 0) or np.any(indices >= region.sizes["day_number"]):
        raise ValueError("Berkeley Earth day_of_year falls outside the available climatology.")
    for name in ("temperature", "climatology"):
        if selected[name].attrs.get("units") not in ("degree C", "degrees C", "degC"):
            raise ValueError(f"Unexpected units for {name}: {selected[name].attrs.get('units')!r}")

    anomaly = selected.temperature.load()
    climate = region.climatology.load().isel(
        day_number=xr.DataArray(indices, dims="time", coords={"time": times})
    )
    absolute_tmax = anomaly + climate
    weights = np.cos(np.deg2rad(absolute_tmax.latitude))
    daily = (
        absolute_tmax.weighted(weights)
        .mean(dim=["latitude", "longitude"], skipna=True)
        .rename("tmax")
    )
    if bool(daily.isnull().any()):
        print(
            f"WARNING: {int(daily.isnull().sum())} day(s) in this source subset "
            "have no valid regional temperatures; retaining NaN.",
            flush=True,
        )
    return daily


# %% Build and save reference series
def main():
    """Retrieve source files, reconstruct daily temperatures, and save results."""
    try:
        import h5netcdf  # noqa: F401
        import h5py  # noqa: F401
    except ImportError as exc:
        raise ImportError(
            "Berkeley Earth requires h5netcdf and h5py. Install them into the "
            "Python environment used to run this script."
        ) from exc
    if not 1880 <= START_YEAR <= END_YEAR <= 2023:
        raise ValueError("Choose complete years between 1880 and 2023 for the verified release.")
    north, south = PNW_LATITUDE_BOUNDS
    west, east = PNW_LONGITUDE_BOUNDS
    if not -90 <= south < north <= 90 or not -180 <= west < east <= 180:
        raise ValueError("Use north > south and -180 <= west < east <= 180.")

    daily_series = []
    source_histories = []
    with requests.Session() as session:
        for decade in range(START_YEAR // 10 * 10, END_YEAR // 10 * 10 + 1, 10):
            path = download_decade(session, decade)
            print(f"Processing PNW temperatures from {path.name}...", flush=True)
            with xr.open_dataset(path, engine="h5netcdf", decode_times=False) as dataset:
                daily_series.append(process_decade(dataset))
                source_histories.append(f"{path.name}: {dataset.attrs.get('history', 'unknown')}")

    TX = xr.concat(daily_series, dim="time").sortby("time")
    expected_times = np.arange(
        np.datetime64(f"{START_YEAR}-01-01"), np.datetime64(f"{END_YEAR + 1}-01-01")
    )
    if not np.array_equal(TX.time.values, expected_times):
        raise ValueError("Source files do not contain a unique, complete requested daily calendar.")
    TXx = TX.resample(time="1YE").max(skipna=True)
    TXx = TXx.assign_coords(
        valid_days=TX.resample(time="1YE").count(),
        expected_days=xr.ones_like(TX, dtype="int32").resample(time="1YE").sum(),
    )
    TXx.valid_days.attrs = {"long_name": "Days with a valid regional mean", "units": "days"}
    TXx.expected_days.attrs = {"long_name": "Calendar days expected in the year", "units": "days"}
    shared_attrs = {
        "source_dataset": "Berkeley Earth experimental global daily land TMAX, 1 degree",
        "source_url": SOURCE_URL,
        "source_file_url_template": f"{DATA_BASE_URL}/Complete_TMAX_Daily_LatLong1_<decade>.nc",
        "source_histories": "\n".join(source_histories),
        "execution_date": date.today().isoformat(),
        "units": "degC",
        "spatial_domain": f"{south:g}-{north:g} N, {west:g}-{east:g} E",
        "spatial_weighting": "Cosine of latitude over valid cells with centres within bounds",
        "reference_period": f"{START_YEAR}-{END_YEAR}",
        "temperature_reconstruction": (
            "Absolute daily TMAX = temperature anomaly + 1951-1980 daily climatology; "
            "climatology indexed using the source day_of_year variable minus one."
        ),
        "missing_data_handling": (
            "Missing cells excluded from spatial means; all-missing days retained as NaN. "
            "Annual maxima use available days and may be underestimated in incomplete years; "
            "compare valid_days with expected_days. No interpolation applied."
        ),
    }
    TX.attrs = {
        **shared_attrs,
        "long_name": "PNW area-mean daily maximum temperature",
        "temporal_aggregation": "Daily grid-cell maximum, then spatial mean",
    }
    TXx.attrs = {
        **shared_attrs,
        "long_name": "Annual maximum of PNW area-mean daily maximum temperature",
        "temporal_aggregation": "Annual maximum of regional daily TMAX",
    }
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
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    daily_output = TX.to_dataset().assign_coords(
        valid_days=TX.notnull().astype("int32"),
        expected_days=xr.ones_like(TX, dtype="int32"),
    )
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
    for series, filename in (
        (daily_output, "Berkeley_reference_daily_tmax.nc"),
        (annual_output, "Berkeley_reference_tmax.nc"),
    ):
        path = OUTPUT_DIR / filename
        series.to_netcdf(path, engine="scipy")
        print(f"Saved {path}", flush=True)
    return TX, TXx


if __name__ == "__main__":
    TX, TXx = main()

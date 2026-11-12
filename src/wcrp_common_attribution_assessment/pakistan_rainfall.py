import hashlib
import netrc
import tempfile
from datetime import date
from pathlib import Path

import aiohttp
import numpy as np
import regionmask
import xarray as xr

from wcrp_common_attribution_assessment import constants

MM_PER_M = 1000.0

# Source selection and local storage. Edit these settings before running.
ERA5_URL = "https://data.earthdatahub.destine.eu/era5/reanalysis-era5-single-levels-v0.zarr"
START_YEAR = 1990
END_YEAR = 2022
LATITUDE_BOUNDS = (60, -10)
LONGITUDE_BOUNDS = (50, 130)
RAW_DATA_DIR = constants.OUTPUTS_DIR / "pakistan_rainfall"
NETRC_FILE = constants.PROJECT_ROOT / "_netrc"


def open_raw_precipitation():
    """Save the unaggregated source subset once, then open the local copy."""
    selection = {
        "latitude": slice(*LATITUDE_BOUNDS),
        "longitude": slice(*LONGITUDE_BOUNDS),
        "valid_time": slice(str(START_YEAR), str(END_YEAR)),
    }
    cache_key = hashlib.sha256(repr((ERA5_URL, selection)).encode("utf-8")).hexdigest()[:16]
    cache_path = RAW_DATA_DIR / f"era5_tp_hourly_{START_YEAR}_{END_YEAR}_{cache_key}.zarr"
    if not cache_path.exists():
        RAW_DATA_DIR.mkdir(parents=True, exist_ok=True)
        client_kwargs = {"trust_env": True}
        if NETRC_FILE.is_file():
            credentials = netrc.netrc(str(NETRC_FILE)).authenticators(
                "data.earthdatahub.destine.eu"
            )
            if credentials is not None:
                login, account, password = credentials
                client_kwargs["auth"] = aiohttp.BasicAuth(login or account, password)
        print(f"Saving original hourly ERA5 precipitation to {cache_path}...", flush=True)
        with xr.open_dataset(
            ERA5_URL,
            storage_options={"client_kwargs": client_kwargs},
            chunks={},
            engine="zarr",
        ) as remote:
            raw = remote[["tp"]].sel(**selection)
            raw.attrs.update(
                source_url=ERA5_URL,
                retrieval_date=date.today().isoformat(),
                source_selection=repr(selection),
                processing="Subset only; original units, grid and hourly values retained.",
            )
            for variable in raw.variables.values():
                variable.encoding = {}
            raw = raw.chunk({"valid_time": 744, "latitude": 81, "longitude": 121})
            # Publish only a completed download; interrupted stores are never reused.
            partial_path = Path(tempfile.mkdtemp(
                prefix="era5_tp_", suffix=".partial", dir=RAW_DATA_DIR
            ))
            raw.to_zarr(partial_path, mode="w", zarr_format=2, consolidated=True)
            partial_path.replace(cache_path)
    print(f"Reading original precipitation from {cache_path}...", flush=True)
    return xr.open_dataset(cache_path, engine="zarr", chunks={})


if __name__ == "__main__":
    # To ensure a common dataset source is used, we use the DestinE Earth Data Hub ARCO zarr.
    # Note that to be able to access this dataset, you will need to set up an Earthdata login
    # and configure your environment. See the DestinE Earth Data Hub documentation for details.
    era5 = open_raw_precipitation()

    print(
        "Selecting ERA5 total precipitation over the Pakistan rainfall analysis region...",
        flush=True,
    )
    era5_tp = (
        era5.tp
        .coarsen(valid_time=24)  # use coarsen for efficiency
        .sum()
    )

    print("Building Pakistan land mask...", flush=True)
    pk_regions = regionmask.defined_regions.natural_earth_v5_1_2.countries_10
    pk_key = pk_regions.map_keys("PK")
    pk_mask = pk_regions.mask(era5_tp.longitude, era5_tp.latitude, flag=None)
    pk_mask = pk_mask == pk_key

    print("Computing Pakistan area-mean daily total precipitation...", flush=True)
    # Convert hourly accumulations to daily totals before deriving the seasonal mean
    # and Rx5d series.
    era5_tp_pk = (
        era5_tp.where(pk_mask)
        .weighted(np.cos(np.deg2rad(era5_tp.latitude)))
        .mean(dim=["latitude", "longitude"])
        .compute()
    )
    era5.close()

    print("Computing JJA accumulated precipitation and annual Rx5d precipitation...", flush=True)
    pr_jja = (
        era5_tp_pk.sel(valid_time=era5_tp_pk.valid_time.dt.month.isin([6, 7, 8]))
        .groupby("valid_time.year")
        .sum()
    )
    rx5d = (
        era5_tp_pk.rolling(valid_time=5, center=True, min_periods=5)
        .sum()
        .groupby("valid_time.year")
        .max()
    )

    print(
        "Calculating event definitions for the 2022 Pakistan rainfall JJA accumulated "
        "precipitation...",
        flush=True,
    )
    # 1. Absolute definition
    pr_jja_absolute = pr_jja.sel(year=2022).squeeze()

    # 2. Seasonal mean anomaly definition
    pr_jja_anom = pr_jja_absolute - pr_jja.sel(year=slice(1991, 2020)).mean()

    # 3. Seasonal mean anomaly as percentage of mean
    pr_jja_anom_percent = pr_jja_absolute / pr_jja.sel(year=slice(1991, 2020)).mean()

    jja_output_lines = [
        (
            "Absolute event definition "
            "(2022 June-August accumulation of the Pakistan area-mean daily total precipitation, "
            "in mm): "
            f"{pr_jja_absolute.item() * MM_PER_M:.2f} mm"
        ),
        (
            "Seasonal-accumulation anomaly event definition "
            "(2022 June-August accumulation minus the 1991-2020 mean June-August Pakistan "
            "area-mean precipitation accumulation, in mm): "
            f"{pr_jja_anom.item() * MM_PER_M:.2f} mm"
        ),
        (
            "Seasonal-accumulation standardized event definition "
            "(2022 June-August accumulation divided by the 1991-2020 mean June-August "
            "Pakistan area-mean precipitation accumulation): "
            f"{pr_jja_anom_percent.item():.2%}"
        ),
    ]

    for line in jja_output_lines:
        print(line, flush=True)

    print(
        "Calculating event definitions for the 2022 Pakistan rainfall Rx5d precipitation...",
        flush=True,
    )
    # 1. Absolute definition
    rx5d_absolute = rx5d.sel(year=2022).squeeze()

    # 2. Annual maximum anomaly definition
    rx5d_anom = rx5d_absolute - rx5d.sel(year=slice(1991, 2020)).mean()

    # 3. Annual maximum anomaly as percentage of mean
    rx5d_anom_percent = rx5d_absolute / rx5d.sel(year=slice(1991, 2020)).mean()

    rx5d_output_lines = [
        (
            "Absolute event definition "
            "(2022 annual maximum of the Pakistan area-mean 5-day precipitation sum, in mm): "
            f"{rx5d_absolute.item() * MM_PER_M:.2f} mm"
        ),
        (
            "Annual-maximum anomaly event definition "
            "(2022 annual maximum minus the 1991-2020 mean annual maximum of the same "
            "Pakistan area-mean 5-day precipitation sum, in mm): "
            f"{rx5d_anom.item() * MM_PER_M:.2f} mm"
        ),
        (
            "Annual-maximum standardized event definition "
            "(2022 annual maximum divided by the 1991-2020 mean annual maximum of the same "
            "Pakistan area-mean 5-day precipitation sum): "
            f"{rx5d_anom_percent.item():.2%}"
        ),
    ]

    for line in rx5d_output_lines:
        print(line, flush=True)

    output_lines = jja_output_lines + [""] + rx5d_output_lines

    output_path = constants.OUTPUTS_DIR / "pakistan_rainfall.event_definitions"
    print(f"Writing event definitions to {output_path}...", flush=True)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(output_lines) + "\n")
    print("Done.", flush=True)

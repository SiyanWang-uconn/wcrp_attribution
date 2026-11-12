import hashlib
import netrc
import tempfile
from datetime import date
from pathlib import Path

import aiohttp
import numpy as np
import xarray as xr

# %% Configuration
START_YEAR = 1940
END_YEAR = 2025
ERA5_SINGLE_LEVELS_URL = (
    "https://data.earthdatahub.destine.eu/era5/reanalysis-era5-single-levels-v0.zarr"
)
ERA5_PRESSURE_LEVELS_URL = "https://data.earthdatahub.destine.eu/era5/era5-pressure-levels-v0.zarr"
OUTPUT_DIR = Path(__file__).resolve().parents[2] / "outputs" / "pnw_heatwave"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
CACHE_DIR = OUTPUT_DIR / "era5_cache"
NETRC_FILE = Path(__file__).resolve().parents[2] / "_netrc"

# Bounds use descending latitude and the ERA5 0-360 longitude convention.
PNW_LATITUDE_BOUNDS = (52, 45)
PNW_LONGITUDE_BOUNDS = (237, 241)
ANCILLARY_LATITUDE_BOUNDS = (70, 30)
ANCILLARY_LONGITUDE_BOUNDS = (180, 300)

CLIMATOLOGY_MONTHS = [6, 7, 8]
REFERENCE_DATE_RANGE = (str(START_YEAR), str(END_YEAR))
CLIMATOLOGY_DATE_RANGE = ("1991", "2020")
EVENT_DATE_RANGE = ("2021-06-26", "2021-06-30")
RETRIEVAL_DATE = date.today().isoformat()
KELVIN_OFFSET = 273.15

# %% Retrieval helper
def open_cached_subset(name, url, variables, selection, periods):
    """Cache unaggregated ERA5 values, then read from the local store."""
    # Different source selections use separate caches.
    cache_key = hashlib.sha256(
        repr((url, variables, selection, periods)).encode("utf-8")
    ).hexdigest()[:16]
    cache_path = CACHE_DIR / f"{name}_{cache_key}.zarr"
    if not cache_path.exists():
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        print(f"Downloading raw ERA5 subset to {cache_path}...", flush=True)
        # Earth Data Hub credentials are required only when creating a new cache.
        client_kwargs = {"trust_env": True}
        if NETRC_FILE.is_file():
            credentials = netrc.netrc(str(NETRC_FILE)).authenticators(
                "data.earthdatahub.destine.eu"
            )
            if credentials is not None:
                login, account, password = credentials
                client_kwargs["auth"] = aiohttp.BasicAuth(login or account, password)
        with xr.open_dataset(
            url,
            storage_options={"client_kwargs": client_kwargs},
            chunks={},
            engine="zarr",
        ) as remote:
            subset = remote[variables].sel(**selection)
            pieces = []
            for start, end, months in periods:
                piece = subset.sel(valid_time=slice(start, end))
                if months is not None:
                    piece = piece.sel(valid_time=piece.valid_time.dt.month.isin(months))
                pieces.append(piece)
            raw = xr.concat(pieces, dim="valid_time").sortby("valid_time")
            raw.attrs.update(
                source_url=url,
                retrieval_date=RETRIEVAL_DATE,
                cache_selection=repr((variables, selection, periods)),
                processing="Subset only; original units, grid and time resolution retained.",
            )
            # Bound memory use and discard source-specific storage encodings.
            for variable in raw.variables.values():
                variable.encoding = {}
            raw = raw.chunk({"valid_time": 744, "latitude": 81, "longitude": 121})
            # Publish only complete caches. Interrupted writes leave unused .partial folders.
            partial_path = Path(tempfile.mkdtemp(
                prefix=f"{name}_", suffix=".partial", dir=CACHE_DIR
            ))
            raw.to_zarr(partial_path, mode="w", zarr_format=2, consolidated=True)
            partial_path.replace(cache_path)
    print(f"Reading local ERA5 cache: {cache_path}", flush=True)
    return xr.open_dataset(cache_path, engine="zarr", chunks={})


# %% Build and save reference series
era5 = open_cached_subset(
    "pnw_mx2t",
    ERA5_SINGLE_LEVELS_URL,
    ["mx2t"],
    {"latitude": slice(*PNW_LATITUDE_BOUNDS), "longitude": slice(*PNW_LONGITUDE_BOUNDS)},
    [(*REFERENCE_DATE_RANGE, None)],
)

## REFERENCE SERIES ##

print("Computing domain-mean daily maximum 2 m temperatures for the PNW region...", flush=True)
TX = (
    era5.mx2t.sel(  # select maximum temperature variable
        latitude=slice(*PNW_LATITUDE_BOUNDS),
        longitude=slice(*PNW_LONGITUDE_BOUNDS),
    )  # select chosen PNW domain
    .resample(valid_time="1D")
    .max()  # resample to daily maximum at each gridpoint
    .weighted(np.cos(np.deg2rad(era5.latitude.sel(latitude=slice(*PNW_LATITUDE_BOUNDS)))))
    .mean(dim=["latitude", "longitude"])  # average over the domain
    .compute()  # trigger computation and load into memory
)

TX = (TX - KELVIN_OFFSET).rename(valid_time="time").rename("tmax")
TX = TX.drop_vars(["number", "surface"], errors="ignore")
TX.attrs["units"] = "degC"

print("Computing annual maxima and their dates from the daily time series...", flush=True)
TXx = TX.resample(time="1YE").max()
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
TXx_date = TX.resample(time="1YE").map(
    lambda daily_values: daily_values.idxmax(dim="time", skipna=True)
)

TXx.attrs.update(
    {
        "units": "degC",
        "long_name": "Annual maximum of PNW area-mean daily maximum temperature",
        "reference_period": f"{START_YEAR}-{END_YEAR}",
        "source_dataset": "ERA5 reanalysis, single levels",
        "source_variable": "mx2t (maximum temperature at 2 metres)",
        "source_url": ERA5_SINGLE_LEVELS_URL,
        "retrieval_date": era5.attrs["retrieval_date"],
        "processing": (
            "Selected the PNW domain (45-52 N, 119-123 W), resampled to daily maxima, "
            "computed a cosine-latitude-weighted spatial mean, converted kelvin to degrees "
            "Celsius, and resampled to annual maxima."
        ),
        "spatial_domain": "Pacific Northwest: 45-52 N, 119-123 W",
        "spatial_weighting": "Cosine of latitude",
        "temporal_aggregation": "Annual maximum of daily maximum 2 m temperature",
        "processing_notes": (
            "Data were loaded from the DestinE Earth Data Hub ARCO Zarr store using xarray."
        ),
    }
)

TX.attrs.update(TXx.attrs)
TX.attrs.update(
    long_name="PNW area-mean daily maximum temperature",
    temporal_aggregation="Daily maximum at each grid cell, then spatial mean",
    processing=(
        "Selected the PNW domain, resampled to daily maxima, computed a "
        "cosine-latitude-weighted spatial mean, and converted kelvin to degrees Celsius."
    ),
    cell_methods="time: maximum within days area: mean",
)
TXx.attrs["cell_methods"] = "time: maximum within days area: mean time: maximum over days"
daily_path = OUTPUT_DIR / "ERA5_reference_daily_tmax.nc"
annual_path = OUTPUT_DIR / "ERA5_reference_mx2t.nc"
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
daily_output["year"] = TX.time.dt.year
daily_output.year.attrs = {"long_name": "Calendar year of daily temperature"}
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
reference = TXx.to_dataset()
reference["year"] = TXx.time.dt.year
reference["year"].attrs["long_name"] = "Calendar year of annual maximum"
reference["tmax_date"] = TXx_date
reference["tmax_date"].attrs = (
    {
        "long_name": "Date of annual maximum regional mean daily maximum temperature",
        "tie_handling": "Earliest date when multiple days share the annual maximum",
        "missing_data_handling": "Missing when the entire year has no valid temperatures",
        "description": (
            "First date in UTC if the annual maximum occurs on multiple days; "
            "missing if all daily values in the year are missing."
        ),
    }
)
reference["tmax_date"].encoding = {
    "units": "days since 1970-01-01",
    "dtype": "float64",
    "_FillValue": np.nan,
}
reference.sel(time=slice(*REFERENCE_DATE_RANGE)).to_netcdf(annual_path, engine="scipy")
print(f"Saved daily reference series: {daily_path}", flush=True)
print(f"Saved annual maximum reference series: {annual_path}", flush=True)

# %% Ancillary data

# SURFACE VARIABLES #

era5_t3_domain = open_cached_subset(
    "surface",
    ERA5_SINGLE_LEVELS_URL,
    ["t2m", "msl", "tcc", "tcwv", "swvl1"],
    {
        "latitude": slice(*ANCILLARY_LATITUDE_BOUNDS),
        "longitude": slice(*ANCILLARY_LONGITUDE_BOUNDS),
    },
    [(*CLIMATOLOGY_DATE_RANGE, CLIMATOLOGY_MONTHS), (*EVENT_DATE_RANGE, None)],
)

era5_t3_domain_clim = (
    era5_t3_domain.sel(valid_time=era5_t3_domain.valid_time.dt.month.isin(CLIMATOLOGY_MONTHS))
    .sel(valid_time=slice(*CLIMATOLOGY_DATE_RANGE))
    .mean(dim="valid_time")
)  # compute climatology for June-August 1991-2020

era5_t3_domain_event = (
    era5_t3_domain.sel(valid_time=slice(*EVENT_DATE_RANGE))
    .resample(valid_time="1D")
    .mean(dim="valid_time")
)  # compute mean for the event period

era5_t3_domain_clim = era5_t3_domain_clim.get(["t2m", "msl", "tcc", "tcwv", "swvl1"]).compute()

era5_t3_domain_event = era5_t3_domain_event.get(["t2m", "msl", "tcc", "tcwv", "swvl1"]).compute()

# Convert temperature only; retain all other surface variables in their native units.
for surface_data in (era5_t3_domain_clim, era5_t3_domain_event):
    temperature_attrs = dict(surface_data["t2m"].attrs)
    surface_data["t2m"] = surface_data["t2m"] - KELVIN_OFFSET
    surface_data["t2m"].attrs = {
        **temperature_attrs,
        "units": "degC",
        "unit_conversion": "Original ERA5 temperature in kelvin minus 273.15",
    }

era5_t3_domain_common_attrs = {
    "source_dataset": "ERA5 reanalysis, single levels",
    "source_variables": "t2m, msl, tcc, tcwv, swvl1",
    "source_url": ERA5_SINGLE_LEVELS_URL,
    "retrieval_date": era5_t3_domain.attrs["retrieval_date"],
    "spatial_domain": "30-70 N, 180-300 E (180 W-60 W)",
    "spatial_weighting": "None; variables retained on the native ERA5 grid",
    "processing_notes": (
        "Data were loaded from the DestinE Earth Data Hub ARCO Zarr store using xarray."
    ),
}

era5_t3_domain_clim.attrs.update(
    {
        **era5_t3_domain_common_attrs,
        "temporal_aggregation": "Mean over June-August for 1991-2020",
        "processing": (
            "Selected the 30-70 N, 180-300 E domain, retained June-August 1991-2020, "
            "computed the mean over valid_time, selected t2m, msl, tcc, tcwv, "
            "and swvl1, "
            "and loaded the result into memory."
        ),
    }
)

era5_t3_domain_event.attrs.update(
    {
        **era5_t3_domain_common_attrs,
        "temporal_aggregation": "Daily mean over 2021-06-26 to 2021-06-30",
        "processing": (
            "Selected the 30-70 N, 180-300 E domain and the 2021-06-26 to "
            "2021-06-30 event period, resampled to daily means, selected t2m, "
            "msl, tcc, tcwv, and swvl1, and loaded the result "
            "into memory."
        ),
    }
)

era5_t3_domain_clim.to_netcdf(OUTPUT_DIR / "ERA5_reference_T3_sfc_climate.nc")
era5_t3_domain_event.rename(valid_time="time").to_netcdf(OUTPUT_DIR / "ERA5_reference_T3_sfc.nc")

# PRESSURE-LEVEL VARIABLES #

print("Opening ERA5 pressure-level dataset...", flush=True)
era5_pl = open_cached_subset(
    "z500",
    ERA5_PRESSURE_LEVELS_URL,
    ["z"],
    {
        "isobaricInhPa": 500,
        "latitude": slice(*ANCILLARY_LATITUDE_BOUNDS),
        "longitude": slice(*ANCILLARY_LONGITUDE_BOUNDS),
    },
    [(*CLIMATOLOGY_DATE_RANGE, CLIMATOLOGY_MONTHS), (*EVENT_DATE_RANGE, None)],
)

# The cached pressure-level subset already contains only 500 hPa over the chosen domain.
era5_z500_domain = era5_pl.z

print("Computing 1991-2020 climate and 5-day event daily averages...", flush=True)
era5_z500_domain_clim = (
    era5_z500_domain.sel(valid_time=era5_z500_domain.valid_time.dt.month.isin(CLIMATOLOGY_MONTHS))
    .sel(valid_time=slice(*CLIMATOLOGY_DATE_RANGE))
    .mean(dim="valid_time")
)

era5_z500_domain_event = (
    era5_z500_domain.sel(valid_time=slice(*EVENT_DATE_RANGE))
    .resample(valid_time="1D")
    .mean(dim="valid_time")
)

era5_z500_domain_clim = era5_z500_domain_clim.compute()
era5_z500_domain_event = era5_z500_domain_event.compute()

era5_z500_domain_common_attrs = {
    "source_dataset": "ERA5 reanalysis, pressure levels",
    "source_variable": "z at 500 hPa (geopotential)",
    "source_url": ERA5_PRESSURE_LEVELS_URL,
    "pressure_level": "500 hPa",
    "retrieval_date": era5_pl.attrs["retrieval_date"],
    "spatial_domain": "30-70 N, 180-300 E (180 W-60 W)",
    "spatial_weighting": "None; geopotential retained on the native ERA5 grid",
    "processing_notes": (
        "Data were loaded from the DestinE Earth Data Hub pressure-level ARCO Zarr "
        "store using xarray. The z variable is retained in native ERA5 geopotential "
        "units; no conversion to geopotential height "
        "in metres was applied."
    ),
}

era5_z500_domain_clim.attrs.update(
    {
        **era5_z500_domain_common_attrs,
        "temporal_aggregation": "Mean over June-August for 1991-2020",
        "processing": (
            "Selected z at 500 hPa over the 30-70 N, 180-300 E domain, "
            "retained June-August 1991-2020, "
            "and computed the mean over valid_time."
        ),
    }
)

era5_z500_domain_event.attrs.update(
    {
        **era5_z500_domain_common_attrs,
        "temporal_aggregation": "Daily mean over 2021-06-26 to 2021-06-30",
        "processing": (
            "Selected z at 500 hPa over the 30-70 N, 180-300 E domain and the "
            "2021-06-26 to 2021-06-30 "
            "event period, then resampled to daily means."
        ),
    }
)

era5_z500_domain_clim.to_netcdf(OUTPUT_DIR / "ERA5_reference_T3_pl_climate.nc")
era5_z500_domain_event.rename(valid_time="time").to_netcdf(OUTPUT_DIR / "ERA5_reference_T3_pl.nc")

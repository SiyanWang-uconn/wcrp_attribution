import numpy as np
import xarray as xr

from wcrp_common_attribution_assessment import constants

if __name__ == "__main__":
    print("Opening ERA5 dataset...", flush=True)
    # To ensure a common dataset source is used, we use the DestinE Earth Data Hub ARCO zarr.
    # Note that to be able to access this dataset, you will need to set up an Earthdata login
    # and configure your environment. See the DestinE Earth Data Hub documentation for details.
    era5 = xr.open_dataset(
        "https://data.earthdatahub.destine.eu/era5/reanalysis-era5-single-levels-v0.zarr",
        storage_options={"client_kwargs": {"trust_env": True}},
        chunks={},
        engine="zarr",
    )

    print("Computing domain-mean daily maximum 2 m temperatures for the PNW region...", flush=True)
    TX = (
        era5.mx2t.sel(  # select maximum temperature variable
            latitude=slice(52, 45),
            longitude=slice((-123 - 180) % 360 + 180, (-119 - 180) % 360 + 180),
        )  # select chosen PNW domain
        .resample(valid_time="1D")
        .max()  # resample to daily maximum at each gridpoint
        .weighted(np.cos(np.deg2rad(era5.latitude.sel(latitude=slice(52, 45)))))
        .mean(dim=["latitude", "longitude"])  # average over the domain
        .compute()  # trigger computation and load into memory
    )

    print("Computing annual maxima from the daily time series...", flush=True)
    TXx = TX.resample(valid_time="1YE").max()

    print("Calculating event definitions for the 2021 PNW heatwave...", flush=True)
    # 1. Absolute definition
    absolute = TXx.sel(valid_time="2021").squeeze().values

    # 2. Annual maximum anomaly definition
    TXx_anom = (absolute - TXx.sel(valid_time=slice("1991", "2020")).mean()).squeeze().values

    # 3. Annual maximum standard deviation definition
    TXx_stdev = (TXx_anom / TXx.sel(valid_time=slice("1991", "2020")).std()).squeeze().values

    # 4. Climatological anomaly definition
    TX_dayofyear = TX.sel(valid_time="2021").idxmax().dt.dayofyear.squeeze().values
    era5_climatology_TX = TX.sel(
        valid_time=TX.valid_time.dt.dayofyear.isin(np.arange(TX_dayofyear - 15, TX_dayofyear + 15))
    ).sel(valid_time=slice("1991", "2020"))
    TX_anom = (absolute - era5_climatology_TX.mean()).squeeze().values

    # 5. Climatological standard deviation definition
    TX_stdev = (TX_anom / era5_climatology_TX.std()).squeeze().values

    output_lines = [
        (
            "Absolute event definition "
            "(2021 annual maximum of the domain-mean daily maximum 2 m temperature, in degC): "
            f"{absolute - constants.ZERO_CELSIUS_IN_KELVIN:.2f}°C"
        ),
        (
            "Annual-maximum anomaly event definition "
            "(2021 annual maximum minus the 1991-2020 mean annual maximum of the same "
            "domain-mean daily maximum 2 m temperature, in degC): "
            f"{TXx_anom:.2f}°C"
        ),
        (
            "Annual-maximum standardized event definition "
            "(annual-maximum anomaly divided by the 1991-2020 standard deviation of annual "
            "maxima for the same domain-mean daily maximum 2 m temperature): "
            f"{TXx_stdev:.2f}"
        ),
        (
            "Climatological anomaly event definition "
            "(2021 annual maximum minus the centered 30-day 1991-2020 climatological "
            "mean of the same domain-mean daily maximum 2 m temperature, "
            f"in degC): {TX_anom:.2f}°C"
        ),
        (
            "Climatological standardized event definition "
            "(climatological anomaly divided by the centered 30-day 1991-2020 "
            "climatological standard deviation of the same domain-mean daily maximum "
            "2 m temperature): "
            f"{TX_stdev:.2f}"
        ),
    ]

    for line in output_lines:
        print(line, flush=True)

    output_path = constants.OUTPUTS_DIR / "pnw_heatwave_local.event_definitions"
    print(f"Writing event definitions to {output_path}...", flush=True)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(output_lines) + "\n")
    print("Done.", flush=True)

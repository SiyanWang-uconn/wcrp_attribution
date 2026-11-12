"""Retrieve NOAA GHCN-Daily station Tmax for the PNW, preserving source files.

Run sections interactively or execute this module. Set DOWNLOAD_OBSERVATIONS
to False to inspect the station inventory without downloading daily records.
Documentation: https://www.ncei.noaa.gov/pub/data/ghcn/daily/readme.txt
"""

import calendar
import csv
import gzip
import json
from datetime import date
from pathlib import Path

import requests

# %% Configuration
START_YEAR = 1979
END_YEAR = 2025
LATITUDE_BOUNDS = (45.0, 52.0)  # South, north.
LONGITUDE_BOUNDS = (-123.0, -119.0)  # West, east.
OUTPUT_DIR = Path(__file__).resolve().parents[2] / "outputs" / "pnw_heatwave" / "ghcn_stations"
SOURCE_DIR = OUTPUT_DIR / "source"
BASE_URL = "https://www.ncei.noaa.gov/pub/data/ghcn/daily"
DOWNLOAD_OBSERVATIONS = True
REQUEST_TIMEOUT = (30, 180)


# %% Source retrieval and station selection
def download(session, relative_path):
    """Reuse complete local files; publish new downloads only after success."""
    path = SOURCE_DIR / relative_path
    if path.exists():
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".part")
    url = f"{BASE_URL}/{relative_path}"
    print(f"Downloading {url}", flush=True)
    with session.get(url, timeout=REQUEST_TIMEOUT, stream=True) as response:
        response.raise_for_status()
        with temporary.open("wb") as handle:
            for block in response.iter_content(chunk_size=1024 * 1024):
                handle.write(block)
    temporary.replace(path)
    return path


def select_stations(stations_path, inventory_path):
    """Select stations with Tmax inventory overlapping the region and years."""
    south, north = LATITUDE_BOUNDS
    west, east = LONGITUDE_BOUNDS
    candidates = {}
    with stations_path.open(encoding="utf-8") as handle:
        for line in handle:
            latitude, longitude = float(line[12:20]), float(line[21:30])
            if south <= latitude <= north and west <= longitude <= east:
                elevation = float(line[31:37])
                candidates[line[:11]] = {
                    "station_id": line[:11], "name": line[41:71].strip(),
                    "latitude": latitude, "longitude": longitude,
                    "elevation_m": "" if elevation == -999.9 else elevation,
                    "state": line[38:40].strip(), "country_code": line[:2],
                }
    selected = []
    with inventory_path.open(encoding="utf-8") as handle:
        for line in handle:
            station_id = line[:11]
            if station_id not in candidates or line[31:35] != "TMAX":
                continue
            first, last = int(line[36:40]), int(line[41:45])
            if first <= END_YEAR and last >= START_YEAR:
                selected.append({**candidates[station_id],
                                 "tmax_first_year": first, "tmax_last_year": last})
    return sorted(selected, key=lambda row: row["station_id"])


def write_table(path, rows):
    if not rows:
        raise ValueError(f"No rows to write to {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


# %% Daily Tmax extraction
def tmax_rows(path, start_year=START_YEAR, end_year=END_YEAR):
    """Retain reported values and all flags; blank QC flags define usable Tmax.

    Missing source values (-9999) and nonblank quality flags yield blank tmax_C.
    tmax_reported_C retains flagged measurements for inspection. No infilling,
    homogenization, station merging, or spatial averaging is applied.
    """
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            year, month = int(line[11:15]), int(line[15:17])
            if line[17:21] != "TMAX" or not start_year <= year <= end_year:
                continue
            for day in range(1, calendar.monthrange(year, month)[1] + 1):
                offset = 21 + (day - 1) * 8
                value = int(line[offset:offset + 5])
                measurement, quality, source = (
                    line[offset + 5:offset + 8].ljust(3)
                )
                reported = "" if value == -9999 else value / 10.0
                usable = value != -9999 and quality == " "
                yield {
                    "station_id": line[:11], "date": date(year, month, day).isoformat(),
                    "tmax_C": reported if quality == " " else "",
                    "tmax_date": date(year, month, day).isoformat() if usable else "",
                    "valid_days": int(usable), "expected_days": 1,
                    "tmax_reported_C": reported,
                    "measurement_flag": measurement.strip(), "quality_flag": quality.strip(),
                    "source_flag": source.strip(),
                }


def gather_observations(session, stations):
    """Save individual source files and a combined compressed daily table."""
    destination = OUTPUT_DIR / "station_daily_tmax.csv.gz"
    temporary = destination.with_suffix(destination.suffix + ".part")
    fields = ["station_id", "date", "tmax_C", "tmax_date", "valid_days", "expected_days",
              "tmax_reported_C",
              "measurement_flag", "quality_flag", "source_flag"]
    counts = []
    with gzip.open(temporary, "wt", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for number, station in enumerate(stations, 1):
            station_id = station["station_id"]
            print(f"Station {number}/{len(stations)}: {station_id}", flush=True)
            source_path = download(session, f"all/{station_id}.dly")
            rows = list(tmax_rows(source_path))
            writer.writerows(rows)
            counts.append({
                "station_id": station_id, "reported_calendar_rows": len(rows),
                "valid_tmax_days": sum(row["tmax_C"] != "" for row in rows),
                "quality_flagged_days": sum(bool(row["quality_flag"]) for row in rows),
            })
    temporary.replace(destination)
    write_table(OUTPUT_DIR / "station_coverage.csv", counts)
    print(f"Saved {destination}", flush=True)


# %% Run the selected stage
def main():
    if not START_YEAR <= END_YEAR < date.today().year:
        raise ValueError("Choose a start year and complete end year before the current year.")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    with requests.Session() as session:
        station_path = download(session, "ghcnd-stations.txt")
        inventory_path = download(session, "ghcnd-inventory.txt")
        version_path = download(session, "ghcnd-version.txt")
        download(session, "readme.txt")
        stations = select_stations(station_path, inventory_path)
        write_table(OUTPUT_DIR / "selected_stations.csv", stations)
        manifest = {
            "source_url": BASE_URL, "execution_date": date.today().isoformat(),
            "ghcnd_version": version_path.read_text(encoding="utf-8").strip(),
            "latitude_bounds": LATITUDE_BOUNDS, "longitude_bounds": LONGITUDE_BOUNDS,
            "start_year": START_YEAR, "end_year": END_YEAR, "station_count": len(stations),
            "selection": "TMAX inventory overlaps requested years; no completeness filter",
            "cached_source_policy": "Existing files reused; may predate execution date",
            "daily_table": "Only days in source TMAX monthly records; absent months omitted",
            "quality_policy": "Nonblank quality flags excluded from tmax_C; reported values kept",
        }
        (OUTPUT_DIR / "manifest.json").write_text(
            json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
        )
        print(f"Selected {len(stations)} stations; saved {OUTPUT_DIR / 'selected_stations.csv'}")
        if DOWNLOAD_OBSERVATIONS:
            gather_observations(session, stations)
        else:
            print("Inventory stage complete. Observation downloads are disabled.")
    return stations


if __name__ == "__main__":
    stations = main()

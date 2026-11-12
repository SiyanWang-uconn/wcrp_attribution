"""Export the selected Seattle-Tacoma and Princeton Airport daily Tmax records."""

import csv
import json
from datetime import date, timedelta
from pathlib import Path

from wcrp_common_attribution_assessment.pnw_heatwave_ghcn_station_data import tmax_rows

# %% Configuration
DATA_DIR = Path(__file__).resolve().parents[2] / "outputs" / "pnw_heatwave" / "ghcn_stations"
OUTPUT_DIR = DATA_DIR / "seattle_princeton"
COMMON_START_YEAR = 1948
END_YEAR = 2025
STATION_IDS = ["USW00024233", "CAN01126510"]


def read_table(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_table(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    stations = {r["station_id"]: r for r in read_table(DATA_DIR / "selected_stations.csv")}
    selected = [dict(stations[sid]) for sid in STATION_IDS]
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    yearly = []
    for station in selected:
        station_id = station["station_id"]
        source_rows = list(tmax_rows(
            DATA_DIR / "source" / "all" / f"{station_id}.dly",
            start_year=int(station["tmax_first_year"]), end_year=END_YEAR,
        ))
        first_valid = min(r["date"] for r in source_rows if r["valid_days"])
        start_date = date.fromisoformat(first_valid)
        expected = (date(END_YEAR + 1, 1, 1) - start_date).days
        by_date = {r["date"]: r for r in source_rows}
        assert len(by_date) == len(source_rows), "Duplicate source dates"
        rows = []
        for offset in range(expected):
            day = (start_date + timedelta(days=offset)).isoformat()
            row = by_date.get(day, {
                "station_id": station_id, "date": day, "tmax_C": "", "tmax_date": "",
                "valid_days": 0, "expected_days": 1, "tmax_reported_C": "",
                "measurement_flag": "", "quality_flag": "", "source_flag": "",
            })
            rows.append(row)
        station["export_start_date"] = first_valid
        station["export_end_date"] = f"{END_YEAR}-12-31"
        station["valid_days"] = sum(r["valid_days"] for r in rows)
        station["expected_days"] = expected
        station["coverage_percent"] = 100 * station["valid_days"] / expected
        common_rows = [r for r in rows if r["date"] >= f"{COMMON_START_YEAR}-01-01"]
        station["common_period_valid_days"] = sum(r["valid_days"] for r in common_rows)
        station["common_period_expected_days"] = len(common_rows)
        station["common_period_coverage_percent"] = (
            100 * station["common_period_valid_days"] / len(common_rows)
        )
        write_table(OUTPUT_DIR / f"{station_id}_daily_tmax.csv", rows)
        write_table(OUTPUT_DIR / f"{station_id}_daily_tmax_{COMMON_START_YEAR}_{END_YEAR}.csv",
                    common_rows)
        for year in range(start_date.year, END_YEAR + 1):
            year_rows = [r for r in rows if r["date"].startswith(str(year))]
            valid = sum(r["valid_days"] for r in year_rows)
            yearly.append({"station_id": station_id, "year": year, "valid_days": valid,
                           "expected_days": len(year_rows),
                           "coverage_percent": 100 * valid / len(year_rows)})
        event = [r for r in rows if "2021-06-26" <= r["date"] <= "2021-06-30"]
        station["heatwave_valid_days_June26_30_2021"] = sum(r["valid_days"] for r in event)
        print(station)
        print("June 26-30, 2021:", [(r["date"], r["tmax_C"]) for r in event])
    write_table(OUTPUT_DIR / "selected_stations.csv", selected)
    write_table(OUTPUT_DIR / "annual_coverage.csv", yearly)
    (OUTPUT_DIR / "selection_notes.json").write_text(json.dumps({
        "station_ids": STATION_IDS, "end_year": END_YEAR,
        "common_start_year": COMMON_START_YEAR,
        "selection": "Seattle-Tacoma and Princeton Airport, explicitly selected by user",
        "export_period": "Earliest QC-valid observation through end year for each station",
        "restriction": "Only previously downloaded stations within 45-52 N, 123-119 W",
        "daily_calendar": "Full calendar; absent source dates added as missing, no interpolation",
    }, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

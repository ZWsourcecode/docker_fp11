#!/usr/bin/env python3
"""
Script to calculate daily fossil-fuel CO2 (ffco2) for stations over a date range,
and/or merge monthly CSVs. Uses `shared.cal_ffco2` and `shared.merge_months`.

Modes:
  - `daily`: generate daily `ffco2` CSVs only
  - `merge`: run monthly merges only
  - `both` : run daily calculation then merge months (default)

Usage:
    python gen_ffco2.py --start YYYY-MM-DD --end YYYY-MM-DD --domain global|eu|both --mode daily|merge|both [--station STATION]

Examples:
    python gen_ffco2.py --start 2025-08-01 --end 2025-08-31 --domain global --mode daily
    python gen_ffco2.py --start 2025-08-01 --end 2025-09-30 --domain both --station hfd100 --mode both
"""
from shared import cal_ffco2, merge_months
import argparse
import datetime as dt
import pandas as pd
import os
from pathlib import Path

# Configuration
IN_PATH = "/flexpart/postprocess/"
MONTH_PATH = IN_PATH + "flexpartweb/"
PATH_FF = IN_PATH + "flux/"  # fossil fuel emission files

DEFAULT_START = None
DEFAULT_END = None

# Station configuration files
STATION_CONFIGS = {
    "global": "stations.conf",
    "eu": "stations_eu.conf"
}


def load_station_config(config_file):
    """Load station config similar to `gen_month14C.py`.
    Returns a DataFrame with columns [station, lon, lat, z, cpuid]
    """
    if not os.path.exists(config_file):
        raise FileNotFoundError(f"Station config file not found: {config_file}")

    with open(config_file, "r") as f:
        lines = [line.strip() for line in f if line.strip() and not line.strip().startswith('#')]

    data = []
    for line in lines:
        parts = line.split()
        if len(parts) >= 5:
            data.append(parts)
        else:
            print(f"Warning: Skipping malformed line: {line}", flush=True)

    if not data:
        raise ValueError(f"No valid station data found in {config_file}")

    return pd.DataFrame(data, columns=["station", "lon", "lat", "z", "cpuid"])


def iterate_days(start_date, end_date):
    """Yield datetime.date objects from start_date to end_date inclusive."""
    d = start_date
    while d <= end_date:
        yield d
        d = d + dt.timedelta(days=1)


def process_station_for_period(station, start_date, end_date, domain, mode='both'):
    """Calculate ffco2 for each day in the date range for one station.
    After computing days (if mode includes `daily`), call merge_months for each year encountered
    when mode includes `merge`.
    """
    station_upper = station.upper()
    print(f"\n=== Processing station {station_upper} domain={domain} from {start_date} to {end_date} ===", flush=True)

    years = set(d.year for d in iterate_days(start_date, end_date))
    day_count = 0

    if mode in ('daily', 'both'):
        for day in iterate_days(start_date, end_date):
            Year = day.year
            Month = day.month
            Day = day.day
            # ATT_PATH as used by postprocess.py
            # Use different station folders for global vs eu
            base_folder = "stations" if domain == "global" else "stations_eu"
            ATT_PATH = MONTH_PATH + base_folder + "/" + station_upper + "/" + str(Year) + "/" + str(Month).zfill(2)
            Path(ATT_PATH).mkdir(parents=True, exist_ok=True)
            try:
                cal_ffco2(PATH_FF, ATT_PATH, station_upper, Year, Month, Day, domain=domain)
                print(f"  ✓ {station_upper} {day} done", flush=True)
            except Exception as e:
                print(f"  ✗ Error calculating ffco2 for {station_upper} {day}: {e}", flush=True)
            day_count += 1

    print(f"Processed {day_count} day(s) for {station_upper}", flush=True)

    # Merge months for each year encountered if requested
    if mode in ('merge', 'both'):
        for yr in sorted(years):
            # create a date object in that year to pass to merge_months
            date_for_year = dt.date(yr, 1, 1)
            try:
                merge_months(MONTH_PATH, station_upper, date_for_year, domain=domain, prefix="ffco2")
                print(f"Merged months for {station_upper} year {yr}", flush=True)
            except Exception as e:
                print(f"Error merging months for {station_upper} year {yr}: {e}", flush=True)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Calculate daily fossil-fuel CO2 (ffco2) for stations over a period and merge months",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --start 2025-08-01 --end 2025-08-31 --domain global
  %(prog)s --start 2025-08-01 --end 2025-09-30 --domain both --station hfd100
"""
    )

    parser.add_argument('--start', required=True, help='Start date YYYY-MM-DD')
    parser.add_argument('--end', required=True, help='End date YYYY-MM-DD')
    parser.add_argument('--domain', choices=['global', 'eu', 'both'], default='global', help='Domain to process')
    parser.add_argument('--mode', choices=['daily', 'merge', 'both'], default='both', help='Operation mode: daily (calc only), merge (merge only), both (default)')
    parser.add_argument('--station', '-s', type=str, help='Optional: process only a specific station code (case-insensitive)')
    parser.add_argument('--list-stations', action='store_true', help='List available stations and exit')
    return parser.parse_args()


def list_stations():
    print("Available stations:", flush=True)
    for domain, config_file in STATION_CONFIGS.items():
        try:
            station_conf = load_station_config(config_file)
            stations = station_conf["station"].tolist()
            print(f"\n{domain.upper()} domain ({len(stations)} stations):", flush=True)
            for i, station in enumerate(stations, 1):
                print(f"  {i:2d}. {station}", flush=True)
        except Exception as e:
            print(f"\n{domain.upper()} domain: Error loading {config_file} - {e}", flush=True)


def main():
    args = parse_arguments()

    if args.list_stations:
        list_stations()
        return

    # parse dates
    try:
        start_date = dt.datetime.strptime(args.start, "%Y-%m-%d").date()
        end_date = dt.datetime.strptime(args.end, "%Y-%m-%d").date()
    except ValueError:
        print("Error: invalid date format. Use YYYY-MM-DD", flush=True)
        return

    if start_date > end_date:
        print("Error: start date must be <= end date", flush=True)
        return

    domains = [args.domain] if args.domain != 'both' else ['global', 'eu']
    mode = args.mode

    # If station specified, only process that station
    if args.station:
        stations_to_run = [args.station]
        # ensure case-insensitive match could be added, but we pass as-is
    else:
        # Load stations from selected domains
        stations_to_run = []
        for d in domains:
            try:
                conf = load_station_config(STATION_CONFIGS[d])
                stations_to_run.extend(conf['station'].tolist())
            except Exception as e:
                print(f"Could not load stations for domain {d}: {e}", flush=True)

    if not stations_to_run:
        print("No stations to process.", flush=True)
        return

    for domain in domains:
        # filter stations per domain
        domain_conf = load_station_config(STATION_CONFIGS[domain])
        domain_station_list = [s for s in stations_to_run if s in domain_conf['station'].tolist() or s.lower() in [x.lower() for x in domain_conf['station'].tolist()]]
        if not domain_station_list:
            print(f"No stations found for domain {domain}", flush=True)
            continue

        for station in domain_station_list:
            process_station_for_period(station, start_date, end_date, domain, mode=mode)

    print("All done.", flush=True)


if __name__ == '__main__':
    main()

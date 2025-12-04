import glob
import os
import pandas as pd
import datetime as dt
from pathlib import Path
import numpy as np
import xarray as xr
from os.path import isfile, join, exists

def merge_month(in_path, station, date, domain="global", prefix='delta_14C'):
    """Merge daily CSV files for a station/month into a single monthly CSV.

    This function is generic and can merge different kinds of CSV outputs produced
    by the processing pipeline. Use the `prefix` argument to specify the file
    name prefix to look for (e.g. 'delta_14C' or 'ffco2').

    Args:
        in_path: Path or string to data root (e.g. app_dir / 'data/v11')
        station: station code (case-insensitive)
        date: a datetime.date or datetime-like with .year and .month
        domain: 'global' or 'eu' to select file name pattern
        prefix: filename prefix to match (default: 'delta_14C')

    Returns:
        pd.DataFrame: concatenated and sorted DataFrame (empty DataFrame if no files).
    """

    in_path = Path(in_path)
    station_str = str(station)
    station_lower = station_str.lower()
    Year = str(date.year)
    Month = str(date.month).zfill(2)

    # Build pattern and outfile depending on domain
    if domain == "global":
        datapath = in_path / "stations" / station_str / Year / Month
        outpath = in_path / "stations" / station_str / Year
        pattern = f"{datapath}/**/{prefix}_{station_lower}_[0-9]*.csv"
        outfile = outpath / f"{prefix}_{station_lower}_month_{Year}{Month}.csv"
    elif domain == "eu":
        datapath = in_path / "stations_eu" / station_str / Year / Month
        outpath = in_path / "stations_eu" / station_str / Year
        pattern = f"{datapath}/**/{prefix}_{station_lower}_[0-9]*eu.csv"
        outfile = outpath / f"{prefix}_{station_lower}_month_{Year}{Month}eu.csv"
    else:
        raise ValueError("domain must be 'global' or 'eu'")
    outpath.mkdir(parents=True, exist_ok=True)

    # Find all matching CSV files recursively
    csv_files = glob.glob(pattern, recursive=True)

    if not csv_files:
        print(f"No files found for station={station_str}, path={datapath}")
        return pd.DataFrame()

    all_dataframes = []
    for csv_file in sorted(csv_files):
        try:
            df = pd.read_csv(csv_file)
            # Add source file info if needed
            df['source_file'] = os.path.basename(csv_file)
            all_dataframes.append(df)
        except Exception as e:
            print(f"Error reading {csv_file}: {e}")
            continue

    if not all_dataframes:
        print("No readable CSV files found.")
        return pd.DataFrame()

    # Concatenate all DataFrames
    combined_df = pd.concat(all_dataframes, ignore_index=True)

    # Normalize and sort by datetime if present
    col = 'datetime'
    if col in combined_df.columns:
        try:
            combined_df[col] = pd.to_datetime(combined_df[col], errors='coerce')
            combined_df = combined_df.sort_values(col).reset_index(drop=True)
        except Exception as e:
            print(f"Error processing {col}: {e}")

    # Save the combined DataFrame to a new CSV file
    try:
        combined_df.to_csv(outfile, index=False)
    except Exception as e:
        print(f"Could not write output file {outfile}: {e}")
    print(f"Merged {len(csv_files)} files into {os.path.basename(outfile)}")
    # return combined_df

def merge_months(in_path, station, date, domain="global", prefix='delta_14C'):
    """For each month in the year, merge daily CSV files for a station/month into a single monthly CSV.

    This function forwards the `prefix` argument to `merge_month` so callers can
    choose which file prefix to merge (e.g. 'delta_14C' or 'ffco2').

    Args:
        in_path: Path or string to data root (e.g. app_dir / 'data/v11')
        station: station code (case-insensitive)
        date: a datetime.date or datetime-like with .year used to select the year
        domain: 'global' or 'eu' to select which station subfolder to check ('stations' vs 'stations_eu')
        prefix: filename prefix to match (default: 'delta_14C')
    """

    in_path = Path(in_path)
    station_str = str(station)
    Year = str(date.year)

    # Choose base folder depending on domain to match merge_month behavior
    base_folder = "stations" if domain == "global" else "stations_eu"
    for month in range(1, 13):
        month_z = str(month).zfill(2)
        monthpath = in_path / base_folder / station_str / Year / month_z
        # print(f"{monthpath}, monthpath.exists()={monthpath.exists()}")
        if not monthpath.exists():
            continue  # Skip months with no data directory
        month_date = date.replace(month=month, day=1)
        merge_month(in_path, station, month_date, domain=domain, prefix=prefix)

def merge_year(in_path, station, date, domain="global", overwrite=False, prefix='delta_14C'):
    """
    Merge the daily CSVs (or create them) in a year into a single yearly CSV.

    This works for different file-prefixes; use `prefix` to select which files
    to include in the yearly merge (default: 'delta_14C').

    Args:
        in_path: Path or string to data root (e.g. app_dir / 'data/v11')
        station: station code (case-insensitive)
        date: a datetime.date or datetime-like with .year and .month
        domain: 'global' or 'eu' to select file name pattern
        overwrite: if False and yearly outfile exists, the function returns the existing file
        prefix: filename prefix to match (default: 'delta_14C')

    Returns:
        pd.DataFrame: concatenated and sorted DataFrame (empty DataFrame if no files).
    """

    in_path = Path(in_path)
    station_str = str(station)
    station_lower = station_str.lower()
    Year = str(date.year)

    # Determine yearly outfile name based on domain
    if domain == "global":
        datapath = in_path / "stations" / station_str / Year
        outpath = in_path / "stations" / station_str / Year
        pattern = f"{datapath}/**/{prefix}_{station_lower}_[0-9]*.csv"
        outfile = outpath / f"{prefix}_{station_lower}_year_{Year}.csv"
    elif domain == "eu":
        datapath = in_path / "stations_eu" / station_str / Year 
        outpath = in_path / "stations_eu" / station_str / Year
        pattern = f"{datapath}/**/{prefix}_{station_lower}_[0-9]*eu.csv"
        outfile = outpath / f"{prefix}_{station_lower}_year_{Year}eu.csv"
    else:
        raise ValueError("domain must be 'global' or 'eu'")
    outpath.mkdir(parents=True, exist_ok=True)

    if outfile.exists() and not overwrite:
        # Load and return existing file to avoid recomputing
        try:
            return pd.read_csv(outfile)
        except Exception as e:
            print(f"Could not read existing yearly file {outfile}: {e}")
            # continue and attempt to rebuild
    
    # Find all matching CSV files recursively
    csv_files = glob.glob(pattern, recursive=True)

    if not csv_files:
        print(f"No files found for station={station_str}, path={datapath}")
        return pd.DataFrame()

    all_dataframes = []
    for csv_file in sorted(csv_files):
        try:
            df = pd.read_csv(csv_file)
            # Add source file info if needed
            df['source_file'] = os.path.basename(csv_file)
            all_dataframes.append(df)
        except Exception as e:
            print(f"Error reading {csv_file}: {e}")
            continue

    if not all_dataframes:
        print("No readable CSV files found.")
        return pd.DataFrame()

    # Concatenate all DataFrames
    combined_df = pd.concat(all_dataframes, ignore_index=True)

    # Sort by datetime if present
    col = 'datetime'
    if col in combined_df.columns:
        try:
            combined_df[col] = pd.to_datetime(combined_df[col], errors='coerce')
            combined_df = combined_df.sort_values(col).reset_index(drop=True)
        except Exception as e:
            print(f"Error processing {col}: {e}")

    # Save the combined DataFrame to a new CSV file
    try:
        combined_df.to_csv(outfile, index=False)
    except Exception as e:
        print(f"Could not write output file {outfile}: {e}")
    print(f"Merged {len(csv_files)} files into {os.path.basename(outfile)}")
    # return combined_df

def cal_ffco2(PATH_FF, ATT_PATH, Station, Year, Month, Day, domain="global"):
    """Calculate fossil-fuel CO2 at a station for one day using
    GCP emissions and a merged FLEXPART footprint file."""
    # choose ff emission file 
    if domain == "global":
        ds = xr.open_dataset(PATH_FF + "GCP_GridFED_ff_1deg_2024.nc")
    elif domain == "eu":
        ds = xr.open_dataset(PATH_FF + "GCP_GridFED_ff_0.1deg_2024.nc")
    darray_FF = ds["emission"]      
    darray_FF = darray_FF.astype(np.float64)

    df_ffco2 = pd.DataFrame(columns=['datetime','ffco2'])
    Hour_columns = [str(i).zfill(2) for i in range(24)]
    for strHour in Hour_columns:
        folder_Flexpart = f"{Year}x{str(Month).zfill(2)}x{str(Day).zfill(2)}x{strHour}"
        fp = join(ATT_PATH, folder_Flexpart, 'foot')
        try:
            darray_flexpart = xr.open_dataarray(fp, engine='netcdf4')
        except Exception as e:
            print(f'Missing footprint file {fp}: {e}')
            df_ffco2.loc[strHour] = [pd.NaT, np.nan]
            continue
        # ensure same lon convention: convert >180 to -180..
        try:
            lon = darray_flexpart['lon']
            darray_flexpart = darray_flexpart.assign_coords(lon = xr.where(lon>180, lon-360, lon))
            darray_flexpart = darray_flexpart.sortby('lon')
        except Exception:
            pass
        # get 2D footprint slice (time dimension may be present)
        if 'time' in darray_flexpart.dims:
            fp2d = darray_flexpart.isel(time=0)
        else:
            fp2d = darray_flexpart
        # align FF grid to footprint lon/lat if necessary by reindexing
        try:
            ff_on_fp = darray_FF.interp(lon=fp2d['lon'], lat=fp2d['lat'], method='nearest')
        except Exception:
            # fall back to simple broadcasting if coords already align
            ff_on_fp = darray_FF
        # multiply: ff (kg CO2 cell-1 day-1 or micromol m-2 s-1 depending) * footprint (ppm per (micromol m-2 s-1))
        try:
            arr = ff_on_fp.values * fp2d.values
            total = float(np.nansum(arr))
        except Exception as e:
            print('Multiplication error:', e)
            total = np.nan
        dt = np.datetime_as_string(darray_flexpart['time'].values[0], unit='s')
        df_ffco2.loc[strHour] = [dt, round(total, 3)]
    # write csv if requested
    simulate_date = f"{Year}{str(Month).zfill(2)}{str(Day).zfill(2)}"
    Filename_ff = f"ffco2_{Station}_{simulate_date}.csv"
    outpath = join(ATT_PATH, Filename_ff)
    df_ffco2.to_csv(outpath, index=False)

    ds.close()
    # return df_ffco2
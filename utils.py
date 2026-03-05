"""
Utilities for the Eswatini HIV model
"""

import os
import re
import numpy as np
import pandas as pd
import sciris as sc


def set_font(size=None):
    sc.options(fontsize=size)


def parse_unaids_value(s):
    """
    Parse a UNAIDS cell like '1.2 [1.0 - 1.5]' or '2900 [2400 - 3700]'.
    Returns the central estimate as a float, or NaN if unparseable.
    Handles '<200' style entries and space-separated thousands (e.g., '11 000').
    """
    if pd.isna(s):
        return np.nan
    s = str(s).strip()
    if not s or s == '...':
        return np.nan

    # Extract the part before the bracket
    bracket = s.find('[')
    central = s[:bracket].strip() if bracket >= 0 else s.strip()

    # Handle '<200' → treat as the number itself
    central = central.lstrip('<')

    # Remove spaces used as thousand separators (e.g., '11 000' → '11000')
    central = central.replace(' ', '')

    try:
        return float(central)
    except ValueError:
        return np.nan


def extract_unaids(unaids_dir='raw_data/UNAIDS', country='Eswatini'):
    """
    Read all UNAIDS Excel files and extract the row for the given country.

    Each file is named like 'hiv_new_infections.xlsx' where the stem maps
    to a sim result name (with '.' replacing the first '_'): 'hiv.new_infections'.

    Returns a dict of {result_name: pd.Series indexed by year}.
    """
    results = {}
    for fname in sorted(os.listdir(unaids_dir)):
        if not fname.endswith('.xlsx'):
            continue

        stem = fname.replace('.xlsx', '')
        # Convert filename to sim result name: first underscore becomes a dot
        col_name = stem.replace('_', '.', 1)

        df = pd.read_excel(os.path.join(unaids_dir, fname))
        row = df[df['Country'].astype(str).str.contains(country, case=False, na=False)]
        if len(row) == 0:
            print(f'  Warning: {country} not found in {fname}, skipping')
            continue

        row = row.iloc[0]
        year_cols = [c for c in df.columns if c != 'Country']
        series = pd.Series(
            {int(yr): parse_unaids_value(row[yr]) for yr in year_cols},
            dtype=float,
        )
        series.index.name = 'time'
        series.name = col_name
        results[col_name] = series

    return results


def format_calibration_data(prev_file='raw_data/SWAZILAND_nationalprevalence_all_updatedPHIA3.csv',
                            inc_file='raw_data/Swaziland_Incidence_Data2.csv',
                            unaids_dir='raw_data/UNAIDS',
                            outfile='data/eswatini_hiv_calib.csv'):
    """
    Read raw prevalence, incidence, and UNAIDS data and format for STIsim calibration.

    Column names must match what sim.to_df(sep='.') produces:
        - Native HIV results: hiv.prevalence_15_49, hiv.n_infected, etc.
        - hiv_epi analyzer results: hiv_epi.prevalence_15_24, hiv_epi.n_infected_15_49, etc.
        - PHIA age/sex prevalence:
            - Native bins (15-35): hiv.prevalence_f_15_20, etc.
            - Analyzer bins (35+): hiv_epi.prevalence_f_35_40, etc.
        - PHIA incidence: hiv_epi.incidence_f_15_49, hiv_epi.incidence_15_100, etc.
    """
    sex_map = {0: 'm', 1: 'f'}

    # Native HIV result columns (from the disease module's default age bins)
    native_prevalence_bins = {(15, 20), (20, 25), (25, 30), (30, 35)}

    # UNAIDS columns that map to native sim results (no age suffix or overall)
    native_unaids = {'hiv.n_infected', 'hiv.n_infected_f', 'hiv.new_infections',
                     'hiv.new_deaths', 'hiv.prevalence_15_49'}

    # --- UNAIDS data ---
    print('Reading UNAIDS data...')
    unaids = extract_unaids(unaids_dir)
    unaids_df = pd.DataFrame(unaids)
    unaids_df.index.name = 'time'

    # Prevalence is in percent in UNAIDS, convert to proportion
    prev_cols = [c for c in unaids_df.columns if 'prevalence' in c]
    for c in prev_cols:
        unaids_df[c] = unaids_df[c] / 100

    # Rename non-native UNAIDS columns to use hiv_epi prefix
    rename = {}
    for c in unaids_df.columns:
        if c not in native_unaids:
            new_name = c.replace('hiv.', 'hiv_epi.', 1)
            rename[c] = new_name
    unaids_df = unaids_df.rename(columns=rename)

    # --- PHIA survey prevalence: age/sex-specific ---
    print('Reading PHIA prevalence data...')
    prev = pd.read_csv(prev_file)

    prev_rows = {}
    for _, row in prev.iterrows():
        year = int(row['Year'])
        sex = sex_map[row['Gender']]
        lo = int(row['start age'])
        hi = lo + 5
        # Use native prefix for bins the disease module tracks, analyzer prefix for others
        if (lo, hi) in native_prevalence_bins:
            col = f'hiv.prevalence_{sex}_{lo}_{hi}'
        else:
            col = f'hiv_epi.prevalence_{sex}_{lo}_{hi}'
        prev_rows.setdefault(year, {})[col] = row['NationalPrevalence']

    prev_df = pd.DataFrame.from_dict(prev_rows, orient='index')
    prev_df.index.name = 'time'

    # --- PHIA survey incidence: sex-specific ---
    print('Reading PHIA incidence data...')
    inc = pd.read_csv(inc_file)

    inc_rows = {}
    for _, row in inc.iterrows():
        year = int(row['Year'])
        gender = row['Gender']
        lo = int(row['Startage'])
        hi = int(row['Endage'])
        val = row['Incidence'] / 100  # percentage to proportion

        if gender in sex_map:
            col = f'hiv_epi.incidence_{sex_map[gender]}_{lo}_{hi}'
        else:
            col = f'hiv_epi.incidence_{lo}_{hi}'
        inc_rows.setdefault(year, {})[col] = val

    inc_df = pd.DataFrame.from_dict(inc_rows, orient='index')
    inc_df.index.name = 'time'

    # --- Combine all sources ---
    calib = unaids_df.join(prev_df, how='outer').join(inc_df, how='outer')
    calib = calib.sort_index().reset_index()

    # Reorder columns: time, native hiv results, then analyzer results
    native_cols = sorted([c for c in calib.columns if c.startswith('hiv.') or c == 'n_alive'])
    analyzer_cols = sorted([c for c in calib.columns if c.startswith('hiv_epi.')])
    col_order = ['time'] + native_cols + analyzer_cols
    calib = calib[col_order]

    # Save
    calib.to_csv(outfile, index=False)
    print(f'\nSaved calibration data to {outfile} ({len(calib)} rows, {len(calib.columns)} columns)')
    print(f'  Years: {int(calib.time.min())}-{int(calib.time.max())}')
    print(f'  Columns:')
    for c in calib.columns[1:]:
        n_vals = calib[c].notna().sum()
        print(f'    {c} ({n_vals} data points)')

    return calib


if __name__ == '__main__':
    df = format_calibration_data()

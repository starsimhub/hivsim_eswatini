"""Build HIV-deleted background mortality rates for experiment 016.

`data/eswatini_deaths.csv` holds all-cause mortality rates, which for Eswatini
include AIDS deaths — adult mortality rises ~5.5x to a 2005 peak and returns to
baseline by 2025. stisim feeds these to `ss.Deaths`, which kills agents
regardless of HIV status, while the HIV module separately kills agents via
`p_hiv_death` and the `ti_zero` AIDS pathway. That double-counts HIV mortality.

This module constructs a non-AIDS counterfactual by log-linear interpolation
between the pre-epidemic (1985) and post-epidemic (2025) rates, per age and sex.
The interpolated line carries the genuine secular trend — falling child
mortality and so on — and only the observed excess above it is treated as AIDS.

Stated assumptions, which the SUMMARY must repeat:
  - 1985 is AIDS-free (very nearly true) and 2025 is AIDS-free (not quite —
    residual AIDS mortality persists, so this slightly under-deletes).
  - Non-AIDS mortality moved smoothly between the endpoints. Where it did not,
    residual gets misattributed to AIDS.
  - Only years strictly between the endpoints are modified. Rates at and beyond
    2025 are projections with no AIDS hump and are left untouched.

Sanity check available without running the model: the implied AIDS share of
all-cause mortality peaks at ~84% in the mid-30s and falls to 0% at 80+, which
is the right shape for AIDS. See README.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pandas as pd

BASE_YEAR = 1985   # pre-epidemic anchor
END_YEAR = 2025    # post-epidemic anchor


def build_hiv_deleted(deaths_df: pd.DataFrame,
                      base_year: int = BASE_YEAR,
                      end_year: int = END_YEAR) -> pd.DataFrame:
    """Return a copy of `deaths_df` with the AIDS hump removed.

    Rates strictly between `base_year` and `end_year` are replaced by a
    log-linear interpolation of the two anchors, per (Sex, AgeStart). Rates are
    only ever lowered — if the observed rate is already at or below the
    counterfactual, it is kept, so bins AIDS never touched are unchanged.
    """
    df = deaths_df.copy()
    anchors = df[df.Time.isin([base_year, end_year])]
    lo = anchors[anchors.Time == base_year].set_index(['Sex', 'AgeStart'])['Value']
    hi = anchors[anchors.Time == end_year].set_index(['Sex', 'AgeStart'])['Value']

    mid = df.Time.between(base_year, end_year, inclusive='neither')
    idx = pd.MultiIndex.from_frame(df.loc[mid, ['Sex', 'AgeStart']])
    w = ((df.loc[mid, 'Time'] - base_year) / (end_year - base_year)).values

    # Log-linear: constant proportional change between anchors. Mortality trends
    # are multiplicative, so this is the right interpolation; linear would
    # overstate the counterfactual mid-period and under-attribute to AIDS.
    a = lo.reindex(idx).values
    b = hi.reindex(idx).values
    with np.errstate(divide='ignore', invalid='ignore'):
        cf = np.exp((1 - w) * np.log(a) + w * np.log(b))
    cf = np.where(np.isfinite(cf), cf, df.loc[mid, 'Value'].values)

    # Never raise a rate: only the excess above the counterfactual is AIDS.
    df.loc[mid, 'Value'] = np.minimum(df.loc[mid, 'Value'].values, cf)
    return df


def deleted_fraction(deaths_df: pd.DataFrame, hiv_deleted: pd.DataFrame) -> pd.DataFrame:
    """Per row, how much mortality was removed — the audit trail for the SUMMARY.

    This is where the construction's assumptions are visible: bins with a large
    deleted fraction are the ones the method is doing the most work on.
    """
    m = deaths_df.merge(hiv_deleted, on=['Time', 'Sex', 'AgeStart'],
                        suffixes=('_all_cause', '_non_aids'))
    m['deleted_rate'] = m['Value_all_cause'] - m['Value_non_aids']
    m['aids_share'] = np.where(m['Value_all_cause'] > 0,
                               m['deleted_rate'] / m['Value_all_cause'], 0.0)
    return m


def make_datafolder(src: Path, dest: Path, hiv_deleted: pd.DataFrame,
                    deaths_filename: str) -> Path:
    """Create an alternate datafolder with HIV-deleted mortality.

    Copies every CSV from `src` so stisim finds whatever demographic input it
    asks for, then overwrites the deaths file. Only demographics are read from
    the datafolder; the model's other inputs stay pointed at data/.
    """
    dest.mkdir(parents=True, exist_ok=True)
    for f in src.glob('*.csv'):
        shutil.copy2(f, dest / f.name)
    hiv_deleted.to_csv(dest / deaths_filename, index=False)
    return dest

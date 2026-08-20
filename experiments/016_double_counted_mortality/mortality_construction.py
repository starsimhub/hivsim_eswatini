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


def _warp(u: np.ndarray, method: str, par: float) -> np.ndarray:
    """Map elapsed fraction u in [0,1] to interpolation weight in [0,1].

    The weight is what decides how much of the decline has happened by a given
    year, and therefore how much of the observed rate is called AIDS. It is the
    single most consequential assumption in the construction, and it is *not*
    identifiable from this data — there are only two AIDS-free anchors, with
    every year between them contaminated.

      loglinear/linear : w = u, constant (proportional / absolute) change
      power            : w = u**par. par > 1 delays the decline (counterfactual
                         stays high mid-period, so LESS is attributed to AIDS);
                         par < 1 front-loads it.
      sigmoid          : symmetric S-curve, slow-fast-slow. Note this is a no-op
                         at the midpoint year — a symmetric sigmoid passes
                         through the midpoint at u = 0.5, exactly where linear
                         does. It changes the shoulders, not the peak.
    """
    if method in ('loglinear', 'linear'):
        return u
    if method == 'power':
        return u ** par
    if method == 'sigmoid':
        f = lambda x: 1.0 / (1.0 + np.exp(-par * (x - 0.5)))
        return (f(u) - f(0.0)) / (f(1.0) - f(0.0))
    raise ValueError(f'unknown method: {method}')


def build_hiv_deleted(deaths_df: pd.DataFrame,
                      base_year: int = BASE_YEAR,
                      end_year: int = END_YEAR,
                      method: str = 'loglinear',
                      par: float = 1.0) -> pd.DataFrame:
    """Return a copy of `deaths_df` with the AIDS hump removed.

    Rates strictly between `base_year` and `end_year` are replaced by an
    interpolation of the two anchors, per (Sex, AgeStart). Rates are only ever
    lowered — if the observed rate is already at or below the counterfactual it
    is kept, so bins AIDS never touched are unchanged.

    `method='loglinear'` (the default, and what exp 016 ran) interpolates
    geometrically: constant *proportional* change, i.e. exponential decay.
    `method='linear'` interpolates arithmetically, which sits above the
    geometric curve mid-period and so attributes less to AIDS.
    """
    df = deaths_df.copy()
    anchors = df[df.Time.isin([base_year, end_year])]
    lo = anchors[anchors.Time == base_year].set_index(['Sex', 'AgeStart'])['Value']
    hi = anchors[anchors.Time == end_year].set_index(['Sex', 'AgeStart'])['Value']

    mid = df.Time.between(base_year, end_year, inclusive='neither')
    idx = pd.MultiIndex.from_frame(df.loc[mid, ['Sex', 'AgeStart']])
    u = ((df.loc[mid, 'Time'] - base_year) / (end_year - base_year)).values
    w = _warp(u, method, par)

    a = lo.reindex(idx).values
    b = hi.reindex(idx).values
    with np.errstate(divide='ignore', invalid='ignore'):
        if method == 'linear':
            cf = (1 - w) * a + w * b
        else:
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

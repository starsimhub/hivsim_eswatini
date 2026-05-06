"""
Experiment dashboards: two summary figures per run.

  Dashboard 1 — Fit vs targets (2x3, 6 panels):
    A: Sex-stratified HIV incidence (15-49) vs PHIA targets w/ 95% CI bars
    B: Prevalence by age at 4 survey years — Female
    C: Prevalence by age at 4 survey years — Male
    D: ART coverage by age group — Female
    E: ART coverage by age group — Male
    F: VMMC coverage by age (sim lines + data targets at 2007/2016/2021)

  Dashboard 2 — Sexual network attributes (2x2):
    A: Age at sexual debut — CDF by sex
    B: Risk-group composition by sex (stacked bars)
    C: Age-pairing heatmap — continuous (female age, male age) density
    D: Lifetime partnership distribution by sex (capped at 30)

Usage:
    python plot_dashboard.py                    # 10 seeds, label=current
    python plot_dashboard.py --n_seeds 3        # quick
    python plot_dashboard.py --label 004_beta   # used in title + filename
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LogNorm
import starsim as ss
from run_sims import make_sim

FIGURES_DIR = 'figures'
DATA_DIR = 'data'
CALIB_DIR = 'calibration_data'

AGE_BINS_YOUNG = [(15, 20), (20, 25), (25, 30), (30, 35)]
AGE_BINS_OLD = [(35, 40), (40, 45), (45, 50), (50, 55), (55, 60), (60, 65)]
ALL_BINS = AGE_BINS_YOUNG + AGE_BINS_OLD

ART_BINS = [(15, 25), (25, 35), (35, 45)]
VMMC_BINS = [(15, 25), (25, 35), (35, 45), (45, 65)]  # broader bins for plotting; data CSV uses 5-yr bins
VMMC_DATA_BINS = [(10,15),(15,20),(20,25),(25,30),(30,35),(35,40),(40,45),(45,50),(50,55),(55,60),(60,65)]
PREV_YEARS = [2007, 2011, 2016, 2021]
VMMC_SURVEY_YEARS = [2007, 2016, 2021]
NETWORK_SNAPSHOT_YEAR = 2020
LIFETIME_PARTNERS_CAP = 30
PARTNERS_YEAR_CAP = 10  # cap for "partners in last year" histogram
RISK_SNAPSHOT_YEAR = 2020
AGE_GAP_BINS = [(15, 24), (25, 34), (35, 49)]  # female age groups for age-gap labels

F_COLOR = '#d62728'
M_COLOR = '#1f77b4'
F_ART_SHADES = ['#fbb4b9', '#f768a1', '#ae017e']
M_ART_SHADES = ['#9ecae1', '#4292c6', '#08519c']
M_VMMC_SHADES = ['#bdd7e7', '#6baed6', '#2171b5', '#08306b']  # 4 shades for 4 VMMC bins
# Shades by survey year (lightest = earliest)
F_YEAR_SHADES = {2007: '#fcae91', 2011: '#fb6a4a', 2016: '#cb181d', 2021: '#67000d'}
M_YEAR_SHADES = {2007: '#bdd7e7', 2011: '#6baed6', 2016: '#2171b5', 2021: '#08306b'}

BAND_ALPHA = 0.18


# ── Analyzers ────────────────────────────────────────────────────────────────

class ARTbyAgeSex(ss.Analyzer):
    def init_results(self):
        super().init_results()
        results = []
        for sex in ['f', 'm']:
            for lo, hi in ART_BINS:
                results.append(ss.Result(f'p_art_{sex}_{lo}_{hi}', dtype=float, scale=False))
        self.define_results(*results)

    def step(self):
        sim = self.sim
        ti = self.ti
        hiv = sim.diseases.hiv
        ppl = sim.people
        alive = ppl.alive
        for sex_key, sex_bool in [('f', ppl.female), ('m', ppl.male)]:
            for lo, hi in ART_BINS:
                in_bin = alive & sex_bool & (ppl.age >= lo) & (ppl.age < hi)
                n_inf = (in_bin & hiv.infected).count()
                if n_inf > 0:
                    n_art = (in_bin & hiv.on_art).count()
                    self.results[f'p_art_{sex_key}_{lo}_{hi}'][ti] = n_art / n_inf
                else:
                    self.results[f'p_art_{sex_key}_{lo}_{hi}'][ti] = np.nan


class VMMCPrevByAge(ss.Analyzer):
    """Track male circumcision prevalence by age bin over time."""
    def init_results(self):
        super().init_results()
        results = []
        for lo, hi in VMMC_BINS:
            results.append(ss.Result(f'p_circ_{lo}_{hi}', dtype=float, scale=False))
        self.define_results(*results)

    def step(self):
        sim = self.sim
        ti = self.ti
        ppl = sim.people
        if 'vmmc' not in sim.interventions:
            for lo, hi in VMMC_BINS:
                self.results[f'p_circ_{lo}_{hi}'][ti] = np.nan
            return
        circ = sim.interventions['vmmc'].circumcised
        male = ppl.alive & ppl.male
        for lo, hi in VMMC_BINS:
            in_bin = male & (ppl.age >= lo) & (ppl.age < hi)
            n_total = in_bin.count()
            if n_total > 0:
                self.results[f'p_circ_{lo}_{hi}'][ti] = (in_bin & circ).count() / n_total
            else:
                self.results[f'p_circ_{lo}_{hi}'][ti] = np.nan


class PartnershipSnapshot(ss.Analyzer):
    """Capture (female age, male age) pairs + new partners in last year at snapshot year."""
    def __init__(self, year=NETWORK_SNAPSHOT_YEAR, network='structuredsexual', **kwargs):
        super().__init__(**kwargs)
        self.year = year
        self.network = network
        self.f_ages = None
        self.m_ages = None
        self.prev_year_lp = None  # lifetime_partners array at year - 1, UID-indexed
        self.partners_last_year_f = None
        self.partners_last_year_m = None

    def init_results(self):
        super().init_results()
        self.define_results(ss.Result('_dummy', dtype=float, scale=False))

    def step(self):
        sim = self.sim
        year = float(sim.t.yearvec[sim.ti])
        net = sim.networks[self.network]
        ppl = sim.people

        # Capture lifetime_partners baseline at year - 1
        if self.prev_year_lp is None and (self.year - 1 - 0.05 <= year < self.year - 1 + 0.1):
            self.prev_year_lp = np.array(net.lifetime_partners)
            return

        # Capture partnership pairs + partners-in-last-year at snapshot year
        if self.f_ages is None and (self.year - 0.05 <= year < self.year + 0.1):
            active = net.edges.dur > 1
            p1 = net.p1[active]
            p2 = net.p2[active]
            self.m_ages = np.array(ppl.age[p1])
            self.f_ages = np.array(ppl.age[p2])

            if self.prev_year_lp is not None:
                current_lp = np.array(net.lifetime_partners)
                # Align lengths if pop changed — pad baseline with zeros for new UIDs
                n = len(current_lp)
                if len(self.prev_year_lp) < n:
                    padded = np.zeros(n)
                    padded[:len(self.prev_year_lp)] = self.prev_year_lp
                    baseline = padded
                else:
                    baseline = self.prev_year_lp[:n]
                diff = np.clip(current_lp - baseline, 0, None)
                participant = np.array((net.participant & ppl.alive))
                female = np.array(ppl.female)
                male = np.array(ppl.male)
                self.partners_last_year_f = diff[participant & female]
                self.partners_last_year_m = diff[participant & male]


class RiskComposition(ss.Analyzer):
    """Snapshot risk-group composition (incl. FSW/client) at a single year."""
    def __init__(self, year=RISK_SNAPSHOT_YEAR, network='structuredsexual', **kwargs):
        super().__init__(**kwargs)
        self.year = int(year)
        self.network = network
        self.snapshot = None

    def init_results(self):
        super().init_results()
        self.define_results(ss.Result('_dummy', dtype=float, scale=False))

    def step(self):
        if self.snapshot is not None:
            return
        sim = self.sim
        year = float(sim.t.yearvec[sim.ti])
        if not (self.year - 0.05 <= year < self.year + 0.1):
            return

        net = sim.networks[self.network]
        ppl = sim.people
        participant = net.participant & ppl.alive

        snap = {}
        for sex_label, sex_bool, role_label, role_attr in [
            ('Female', ppl.female, 'FSW', net.fsw),
            ('Male', ppl.male, 'Client', net.client),
        ]:
            base = sex_bool & participant
            counts = {}
            counts[role_label] = int((base & role_attr).count())
            for rg_i, rg_label in [(0, 'Low'), (1, 'Medium'), (2, 'High')]:
                counts[rg_label] = int((base & (net.risk_group == rg_i) & ~role_attr).count())
            counts['total'] = int(base.count())
            snap[sex_label] = counts
        self.snapshot = snap


# ── Helpers ──────────────────────────────────────────────────────────────────

def annual_mean(yearvec, vals):
    yi = np.floor(yearvec).astype(int)
    uy = np.unique(yi)
    return uy, np.array([np.nanmean(vals[yi == yr]) for yr in uy])


def annual_last(yearvec, vals):
    yi = np.floor(yearvec).astype(int)
    uy = np.unique(yi)
    return uy, np.array([vals[yi == yr][-1] for yr in uy])


def run_sims(n_seeds):
    sims = []
    for seed in range(1, n_seeds + 1):
        print(f'  seed {seed}/{n_seeds}...')
        analyzers = [ARTbyAgeSex(),
                     VMMCPrevByAge(),
                     PartnershipSnapshot(year=NETWORK_SNAPSHOT_YEAR),
                     RiskComposition(year=RISK_SNAPSHOT_YEAR)]
        sim = make_sim(seed=seed, verbose=-1, analyzers=analyzers)
        sim.run()
        sims.append(sim)
    return sims


def collect(sims):
    yv = sims[0].t.yearvec
    out = {'yearvec': yv}

    inc_f, inc_m = [], []
    prev_by_bin = {(s, lo, hi): [] for s in ['f', 'm'] for (lo, hi) in ALL_BINS}
    art_by_bin = {(s, lo, hi): [] for s in ['f', 'm'] for (lo, hi) in ART_BINS}
    vmmc_by_bin = {(lo, hi): [] for (lo, hi) in VMMC_BINS}
    pairs_f, pairs_m = [], []
    lp_f, lp_m = [], []
    ply_f, ply_m = [], []  # partners in last year
    risk_snapshots = []    # list of per-sim snapshot dicts

    for sim in sims:
        hiv = sim.results.hiv
        epi = sim.analyzers['hiv_epi']
        art_ana = sim.analyzers['artbyagesex']
        vmmc_ana = sim.analyzers['vmmcprevbyage']
        net_snap = sim.analyzers['partnershipsnapshot']

        yrs, ann_f = annual_mean(yv, np.array(epi.results['incidence_f_15_49']))
        _,   ann_m = annual_mean(yv, np.array(epi.results['incidence_m_15_49']))
        inc_f.append(ann_f)
        inc_m.append(ann_m)

        for sex in ['f', 'm']:
            for lo, hi in AGE_BINS_YOUNG:
                _, a = annual_last(yv, np.array(hiv[f'prevalence_{sex}_{lo}_{hi}']))
                prev_by_bin[(sex, lo, hi)].append(a)
            for lo, hi in AGE_BINS_OLD:
                _, a = annual_last(yv, np.array(epi.results[f'prevalence_{sex}_{lo}_{hi}']))
                prev_by_bin[(sex, lo, hi)].append(a)
            for lo, hi in ART_BINS:
                _, a = annual_mean(yv, np.array(art_ana.results[f'p_art_{sex}_{lo}_{hi}']))
                art_by_bin[(sex, lo, hi)].append(a)

        for lo, hi in VMMC_BINS:
            _, a = annual_mean(yv, np.array(vmmc_ana.results[f'p_circ_{lo}_{hi}']))
            vmmc_by_bin[(lo, hi)].append(a)

        if net_snap.f_ages is not None:
            pairs_f.append(net_snap.f_ages)
            pairs_m.append(net_snap.m_ages)
        if net_snap.partners_last_year_f is not None:
            ply_f.append(net_snap.partners_last_year_f)
            ply_m.append(net_snap.partners_last_year_m)

        risk_ana = sim.analyzers['riskcomposition']
        if risk_ana.snapshot is not None:
            risk_snapshots.append(risk_ana.snapshot)

        nw = sim.networks.structuredsexual
        ppl = sim.people
        participant = nw.participant & ppl.alive
        debuted = participant & (ppl.age >= nw.debut)
        lp_f.append(np.array(nw.lifetime_partners[debuted & ppl.female]))
        lp_m.append(np.array(nw.lifetime_partners[debuted & ppl.male]))

    out['years'] = yrs
    out['inc_f'] = np.array(inc_f) * 100
    out['inc_m'] = np.array(inc_m) * 100
    out['prev_by_bin'] = {k: np.array(v) for k, v in prev_by_bin.items()}
    out['art_by_bin'] = {k: np.array(v) for k, v in art_by_bin.items()}
    out['vmmc_by_bin'] = {k: np.array(v) for k, v in vmmc_by_bin.items()}
    out['pairs_f'] = np.concatenate(pairs_f) if pairs_f else np.array([])
    out['pairs_m'] = np.concatenate(pairs_m) if pairs_m else np.array([])
    out['lp_f'] = np.concatenate(lp_f) if lp_f else np.array([])
    out['lp_m'] = np.concatenate(lp_m) if lp_m else np.array([])
    out['ply_f'] = np.concatenate(ply_f) if ply_f else np.array([])
    out['ply_m'] = np.concatenate(ply_m) if ply_m else np.array([])

    # Aggregate risk composition across sims (sum counts at the snapshot year)
    agg_risk = {'Female': {'Low': 0, 'Medium': 0, 'High': 0, 'FSW': 0, 'total': 0},
                'Male':   {'Low': 0, 'Medium': 0, 'High': 0, 'Client': 0, 'total': 0}}
    for snap in risk_snapshots:
        for sex_label in ['Female', 'Male']:
            for cat, cnt in snap[sex_label].items():
                agg_risk[sex_label][cat] += cnt
    out['risk_snapshot'] = agg_risk

    # Debut age from first sim
    nw = sims[0].networks.structuredsexual
    ppl = sims[0].people
    participant = nw.participant & ppl.alive
    out['debut_f'] = np.array(nw.debut[participant & ppl.female])
    out['debut_m'] = np.array(nw.debut[participant & ppl.male])

    return out


def band_ci(ax, x, arr, color, label=None, lw=2, ls='-', alpha_line=1.0):
    med = np.nanmedian(arr, axis=0)
    lo = np.nanpercentile(arr, 2.5, axis=0)
    hi = np.nanpercentile(arr, 97.5, axis=0)
    ax.plot(x, med, color=color, lw=lw, ls=ls, label=label, alpha=alpha_line)
    ax.fill_between(x, lo, hi, color=color, alpha=BAND_ALPHA, linewidth=0)


def _plot_vmmc_panel(ax, data, yrs, vmmc_data, title):
    """Lines = sim circumcision prevalence by age bin; points = data targets aggregated to same bins."""
    mask_yrs = (yrs >= 1990) & (yrs <= 2030)
    x = yrs[mask_yrs]
    for i, (lo, hi) in enumerate(VMMC_BINS):
        color = M_VMMC_SHADES[i]
        arr = data['vmmc_by_bin'][(lo, hi)][:, mask_yrs] * 100
        band_ci(ax, x, arr, color, lw=1.8, label=f'{lo}–{hi}')

        # Aggregate target data 5-yr bins into this broader bin
        for sy in VMMC_SURVEY_YEARS:
            sub = vmmc_data[(vmmc_data['Year'] == sy) &
                            (vmmc_data['AgeBin'].apply(lambda ab: _bin_in(ab, lo, hi)))]
            if len(sub):
                avg = sub['p_vmmc'].mean() * 100
                ax.scatter(sy, avg, color=color, s=45, edgecolor='black',
                           linewidth=0.6, zorder=5)
    ax.set_title(title, fontsize=13)
    ax.set_xlabel('Year', fontsize=11)
    ax.set_ylabel('% of males circumcised', fontsize=11)
    ax.set_xlim(1990, 2030)
    ax.set_ylim(0, 100)
    ax.legend(fontsize=8, loc='upper left', title='Age group')
    ax.grid(True, alpha=0.3)


def _bin_in(ab_str, lo, hi):
    """Test whether an AgeBin string '[a,b)' lies within [lo, hi)."""
    s = str(ab_str).strip('[]() ')
    if ',' in s:
        a, b = s.split(',')
        return float(a) >= lo and float(b) <= hi
    return False


def _plot_art_panel(ax, data, yrs, phia_art, sex_key, sex_int, shades, title):
    mask_art = yrs >= 2004
    x_art = yrs[mask_art]
    for i, (lo, hi) in enumerate(ART_BINS):
        color = shades[i]
        arr = data['art_by_bin'][(sex_key, lo, hi)][:, mask_art] * 100
        band_ci(ax, x_art, arr, color, lw=1.8, label=f'{lo}–{hi}')
        sub = phia_art[(phia_art['Gender'] == sex_int) & (phia_art['AgeBin'] == f'[{lo},{hi})')]
        for _, row in sub.iterrows():
            val = row['NationalARTPrevalence'] * 100
            lb = row['lb'] * 100
            ub = row['ub'] * 100
            ax.errorbar(row['Year'], val, yerr=[[val - lb], [ub - val]],
                        fmt='o', color=color, capsize=3, markersize=5,
                        markeredgecolor='black', markeredgewidth=0.6, zorder=5)
    ax.set_title(title, fontsize=13)
    ax.set_ylabel('% of PLHIV on ART', fontsize=11)
    ax.set_xlim(2004, 2028)
    ax.set_ylim(0, 105)
    ax.legend(fontsize=8, loc='lower right', title='Age group')
    ax.grid(True, alpha=0.3)


def _plot_prev_panel(ax, data, yrs, targets, prev_calib, sex_key, sex_int, shades, title):
    bin_centers = np.array([(lo + hi) / 2 for (lo, hi) in ALL_BINS])
    for year in PREV_YEARS:
        yr_idx_arr = np.where(yrs == year)[0]
        if len(yr_idx_arr) == 0:
            continue
        yr_idx = int(yr_idx_arr[0])
        color = shades[year]

        sim_med, sim_lo, sim_hi = [], [], []
        tgt_vals, tgt_lb, tgt_ub = [], [], []
        for (lo, hi) in ALL_BINS:
            arr = data['prev_by_bin'][(sex_key, lo, hi)][:, yr_idx] * 100
            sim_med.append(np.median(arr))
            sim_lo.append(np.percentile(arr, 2.5))
            sim_hi.append(np.percentile(arr, 97.5))
            # Target lookup: prev_calib has CI for 2007/2011/2016; fall back to targets (no CI) for 2021
            val = lb = ub = np.nan
            m = ((prev_calib['Year'] == year) & (prev_calib['Gender'] == sex_int) &
                 (prev_calib['start age'] == lo))
            if m.any():
                row = prev_calib.loc[m].iloc[0]
                val = row['NationalPrevalence'] * 100
                lb = row['lb'] * 100
                ub = row['ub'] * 100
            else:
                for pfx in ['hiv', 'hiv_epi']:
                    col = f'{pfx}.prevalence_{sex_key}_{lo}_{hi}'
                    if col in targets.columns:
                        trow = targets.loc[targets.time == year, col]
                        if trow.notna().any():
                            val = float(trow.dropna().iloc[0]) * 100
                            lb = ub = val
                            break
            tgt_vals.append(val)
            tgt_lb.append(lb)
            tgt_ub.append(ub)

        ax.plot(bin_centers, sim_med, color=color, lw=2, label=f'{year}')
        ax.fill_between(bin_centers, sim_lo, sim_hi, color=color,
                        alpha=BAND_ALPHA, linewidth=0)
        tgt_arr = np.array(tgt_vals)
        lb_arr = np.array(tgt_lb)
        ub_arr = np.array(tgt_ub)
        has = ~np.isnan(tgt_arr)
        if has.any():
            yerr = [tgt_arr[has] - lb_arr[has], ub_arr[has] - tgt_arr[has]]
            # If lb==ub (no CI available), errorbar will just show the point
            ax.errorbar(bin_centers[has], tgt_arr[has], yerr=yerr,
                        fmt='o', color=color, markersize=6, capsize=3,
                        markeredgecolor='black', markeredgewidth=0.6,
                        zorder=5, ls='none')

    ax.set_title(title, fontsize=13)
    ax.set_xlabel('Age (mid-bin)', fontsize=11)
    ax.set_ylabel('Prevalence (%)', fontsize=11)
    ax.legend(fontsize=8, loc='upper right', title='Year')
    ax.grid(True, alpha=0.3)
    ax.set_ylim(bottom=0)


# ── Dashboard 1: Fit vs targets ──────────────────────────────────────────────

def plot_fit_dashboard(data, targets, phia_art, inc_calib, prev_calib, vmmc_data, label, n_seeds, outdir=FIGURES_DIR):
    fig = plt.figure(figsize=(18, 10.5))
    gs = GridSpec(2, 3, left=0.05, right=0.98, bottom=0.07, top=0.92,
                  wspace=0.30, hspace=0.38)

    yrs = data['years']
    mask = (yrs >= 1990) & (yrs <= 2030)
    x = yrs[mask]

    # (A) Sex-stratified incidence with 95% CI error bars on targets
    ax = fig.add_subplot(gs[0, 0])
    band_ci(ax, x, data['inc_f'][:, mask], F_COLOR, label='Female — model')
    band_ci(ax, x, data['inc_m'][:, mask], M_COLOR, label='Male — model')
    # Include both 15-49 (2016) and 18-49 (2011) targets — comparable enough at population level
    for sex_int, color, lab in [(1, F_COLOR, 'Female — target'), (0, M_COLOR, 'Male — target')]:
        sub = inc_calib[(inc_calib['Gender'] == sex_int) & (inc_calib['Endage'] == 49)]
        if len(sub):
            ax.errorbar(sub['Year'], sub['Incidence'],
                        yerr=[sub['Incidence'] - sub['lb'], sub['ub'] - sub['Incidence']],
                        fmt='o', color=color, capsize=4, markersize=6,
                        markeredgecolor='black', markeredgewidth=0.6, zorder=5, label=lab)
    ax.set_title('(A) HIV incidence by sex (15–49)', fontsize=13)
    ax.set_ylabel('Incidence per 100 PY', fontsize=11)
    ax.set_xlim(1990, 2030)
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=8, loc='upper right')
    ax.grid(True, alpha=0.3)

    # (B) Prevalence — Female
    ax = fig.add_subplot(gs[0, 1])
    _plot_prev_panel(ax, data, yrs, targets, prev_calib, 'f', 1, F_YEAR_SHADES,
                     '(B) Prevalence by age — Female')

    # (C) Prevalence — Male
    ax = fig.add_subplot(gs[0, 2])
    _plot_prev_panel(ax, data, yrs, targets, prev_calib, 'm', 0, M_YEAR_SHADES,
                     '(C) Prevalence by age — Male')

    # (D) ART — Female
    ax = fig.add_subplot(gs[1, 0])
    _plot_art_panel(ax, data, yrs, phia_art, 'f', 1, F_ART_SHADES,
                    '(D) ART coverage — Female')

    # (E) ART — Male
    ax = fig.add_subplot(gs[1, 1])
    _plot_art_panel(ax, data, yrs, phia_art, 'm', 0, M_ART_SHADES,
                    '(E) ART coverage — Male')

    # (F) VMMC coverage by age
    ax = fig.add_subplot(gs[1, 2])
    _plot_vmmc_panel(ax, data, yrs, vmmc_data, '(F) VMMC coverage by age')

    fig.suptitle(f'Fit dashboard — {label} (n={n_seeds} seeds, 95% CI)', fontsize=15, y=0.98)
    out_path = f'{outdir}/dashboard_fit_{label}.png'
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    print(f'Saved {out_path}')
    plt.close()


# ── Dashboard 2: Sexual network attributes ───────────────────────────────────

def _fmt_iqr(vals):
    """Return 'median [Q1–Q3]' string."""
    med = np.median(vals)
    q1, q3 = np.percentile(vals, [25, 75])
    return f'{med:.0f} [{q1:.0f}–{q3:.0f}]'


def _plot_risk_snapshot(ax, risk_snapshot, snapshot_year):
    """Two stacked bars (F, M) of mutually exclusive risk categories."""
    cat_colors = {
        'Low': '#2ca02c',
        'Medium': '#ff7f0e',
        'High': '#d62728',
        'FSW': '#8e44ad',
        'Client': '#8e44ad',
    }

    x_pos = {'Female': 0, 'Male': 1}
    bar_width = 0.6
    legend_seen = set()
    for sex_label, role_cat in [('Female', 'FSW'), ('Male', 'Client')]:
        counts = risk_snapshot[sex_label]
        total = max(counts['total'], 1)
        bottom = 0.0
        for cat in ['Low', 'Medium', 'High', role_cat]:
            pct = counts.get(cat, 0) / total * 100
            legend_label = cat if cat not in legend_seen else None
            legend_seen.add(cat)
            ax.bar(x_pos[sex_label], pct, bar_width, bottom=bottom,
                   color=cat_colors[cat], label=legend_label,
                   edgecolor='white', linewidth=0.5)
            if pct >= 3:
                ax.text(x_pos[sex_label], bottom + pct / 2, f'{pct:.0f}%',
                        ha='center', va='center', fontsize=9, color='white',
                        fontweight='bold')
            bottom += pct

    ax.set_xticks([x_pos['Female'], x_pos['Male']])
    ax.set_xticklabels(['Female', 'Male'])
    ax.set_xlim(-0.6, 1.6)
    ax.set_ylim(0, 105)
    ax.set_title(f'(B) Risk-group composition ({snapshot_year})', fontsize=13)
    ax.set_ylabel('Proportion (%)', fontsize=11)
    ax.legend(fontsize=8, loc='upper center', bbox_to_anchor=(0.5, -0.10),
              ncol=5, frameon=False)


def plot_network_dashboard(data, condom_data, label, n_seeds, outdir=FIGURES_DIR):
    fig = plt.figure(figsize=(20, 11))
    gs = GridSpec(2, 3, left=0.05, right=0.98, bottom=0.09, top=0.92,
                  wspace=0.28, hspace=0.45)

    # (A) Debut age CDF
    ax = fig.add_subplot(gs[0, 0])
    age_grid = np.linspace(10, 30, 200)
    for vals, color, sex_label in [(data['debut_f'], F_COLOR, 'Female'),
                                   (data['debut_m'], M_COLOR, 'Male')]:
        cdf = np.array([np.mean(vals <= a) for a in age_grid])
        med = float(np.median(vals))
        ax.plot(age_grid, cdf, color=color, lw=2, label=f'{sex_label} (median {med:.1f})')
        ax.axvline(med, color=color, ls=':', lw=1, alpha=0.5)
    ax.axhline(0.5, color='grey', ls='--', lw=0.8, alpha=0.5)
    ax.set_title('(A) Age at sexual debut', fontsize=13)
    ax.set_xlabel('Age (years)', fontsize=11)
    ax.set_ylabel('Cumulative proportion', fontsize=11)
    ax.legend(fontsize=10, loc='lower right')
    ax.grid(True, alpha=0.3)
    ax.set_xlim(10, 30)
    ax.set_ylim(0, 1.02)

    # (B) Risk-group composition snapshot
    ax = fig.add_subplot(gs[0, 1])
    _plot_risk_snapshot(ax, data['risk_snapshot'], RISK_SNAPSHOT_YEAR)

    # (C) Condom use by risk-group pairing
    ax = fig.add_subplot(gs[0, 2])
    year_cols = [int(c) for c in condom_data.columns if c != 'partnership']
    cond_rows = {
        '(0,0) Low–Low': '#2ca02c',
        '(1,1) Med–Med': '#ff7f0e',
        '(2,2) High–High': '#d62728',
        '(fsw,client) FSW–Client': '#8e44ad',
    }
    key_map = {
        '(0,0) Low–Low': '(0,0)',
        '(1,1) Med–Med': '(1,1)',
        '(2,2) High–High': '(2,2)',
        '(fsw,client) FSW–Client': '(fsw,client)',
    }
    for label_k, color in cond_rows.items():
        key = key_map[label_k]
        row = condom_data.loc[condom_data['partnership'] == key]
        if len(row):
            vals = row.iloc[0][[str(c) for c in year_cols]].astype(float).values
            ax.plot(year_cols, vals, color=color, marker='o', lw=2, ms=5, label=label_k)
    ax.set_title('(C) Condom use by partnership', fontsize=13)
    ax.set_xlabel('Year', fontsize=11)
    ax.set_ylabel('Condom use probability', fontsize=11)
    ax.set_ylim(0, 1.02)
    ax.set_xlim(1985, 2025)
    ax.legend(fontsize=8, loc='lower right')
    ax.grid(True, alpha=0.3)

    # (D) Age-pairing heatmap with age-gap labels
    ax = fig.add_subplot(gs[1, 0])
    f_ages = data['pairs_f']
    m_ages = data['pairs_m']
    if len(f_ages) > 0:
        hb = ax.hexbin(f_ages, m_ages, gridsize=35, cmap='viridis', mincnt=1,
                       norm=LogNorm(), extent=[15, 65, 15, 65])
        cb = fig.colorbar(hb, ax=ax, label='Partnerships (log)')
        cb.ax.tick_params(labelsize=9)
        ax.plot([15, 65], [15, 65], 'w--', lw=1, alpha=0.8, label='y = x')

        # Age-gap median [IQR] annotations for 3 female age groups
        gap_lines = []
        for lo, hi in AGE_GAP_BINS:
            m = (f_ages >= lo) & (f_ages < hi)
            if m.any():
                gaps = m_ages[m] - f_ages[m]
                med = np.median(gaps)
                q1, q3 = np.percentile(gaps, [25, 75])
                gap_lines.append(f'F {lo}–{hi}: Δ={med:+.1f} [{q1:+.1f}, {q3:+.1f}]')
            else:
                gap_lines.append(f'F {lo}–{hi}: n/a')
        ax.text(0.03, 0.97, '\n'.join(gap_lines),
                transform=ax.transAxes, va='top', ha='left', fontsize=8,
                color='white', bbox=dict(facecolor='black', alpha=0.55, edgecolor='none',
                                         boxstyle='round,pad=0.4'))
        ax.legend(fontsize=9, loc='lower right')
    ax.set_title(f'(D) Age-pairing density ({NETWORK_SNAPSHOT_YEAR})', fontsize=13)
    ax.set_xlabel('Female age (years)', fontsize=11)
    ax.set_ylabel('Male partner age (years)', fontsize=11)
    ax.set_xlim(15, 65)
    ax.set_ylim(15, 65)

    # (E) Lifetime partners
    ax = fig.add_subplot(gs[1, 1])
    bins = np.arange(0, LIFETIME_PARTNERS_CAP + 2) - 0.5
    lp_f = data['lp_f']
    lp_m = data['lp_m']
    lp_f_t = lp_f[lp_f <= LIFETIME_PARTNERS_CAP]
    lp_m_t = lp_m[lp_m <= LIFETIME_PARTNERS_CAP]
    ax.hist([lp_f_t, lp_m_t], bins=bins, density=True,
            color=[F_COLOR, M_COLOR], alpha=0.6,
            label=[f'Female  {_fmt_iqr(lp_f)}',
                   f'Male     {_fmt_iqr(lp_m)}'])
    ax.set_title(f'(E) Lifetime partners (capped at {LIFETIME_PARTNERS_CAP})', fontsize=13)
    ax.set_xlabel('Number of lifetime partners', fontsize=11)
    ax.set_ylabel('Proportion', fontsize=11)
    ax.legend(fontsize=9, loc='upper right', title='median [IQR]')
    ax.grid(True, alpha=0.3)
    ax.set_xlim(-0.5, LIFETIME_PARTNERS_CAP + 0.5)

    # (F) Partners in last year
    ax = fig.add_subplot(gs[1, 2])
    bins_y = np.arange(0, PARTNERS_YEAR_CAP + 2) - 0.5
    ply_f = data['ply_f']
    ply_m = data['ply_m']
    ply_f_t = ply_f[ply_f <= PARTNERS_YEAR_CAP]
    ply_m_t = ply_m[ply_m <= PARTNERS_YEAR_CAP]
    if len(ply_f) > 0 and len(ply_m) > 0:
        ax.hist([ply_f_t, ply_m_t], bins=bins_y, density=True,
                color=[F_COLOR, M_COLOR], alpha=0.6,
                label=[f'Female  {_fmt_iqr(ply_f)}',
                       f'Male     {_fmt_iqr(ply_m)}'])
    ax.set_title(f'(F) New partners in last year ({NETWORK_SNAPSHOT_YEAR})', fontsize=13)
    ax.set_xlabel('New partners in prior 12 months', fontsize=11)
    ax.set_ylabel('Proportion', fontsize=11)
    ax.legend(fontsize=9, loc='upper right', title='median [IQR]')
    ax.grid(True, alpha=0.3)
    ax.set_xlim(-0.5, PARTNERS_YEAR_CAP + 0.5)

    fig.suptitle(f'Network dashboard — {label} (n={n_seeds} seeds)', fontsize=15, y=0.98)
    out_path = f'{outdir}/dashboard_network_{label}.png'
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    print(f'Saved {out_path}')
    plt.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--n_seeds', type=int, default=10)
    parser.add_argument('--label', type=str, default='current')
    parser.add_argument('--outdir', type=str, default=FIGURES_DIR,
                        help='where to save dashboards (default: figures/ scratch)')
    args = parser.parse_args()

    print(f'Running {args.n_seeds} sims...')
    sims = run_sims(args.n_seeds)
    print('Collecting data...')
    data = collect(sims)

    targets = pd.read_csv(f'{DATA_DIR}/eswatini_hiv_calib.csv')
    phia_art = pd.read_csv(f'{CALIB_DIR}/art_coverage_by_age_sex.csv')
    inc_calib = pd.read_csv(f'{CALIB_DIR}/incidence_by_sex.csv')
    prev_calib = pd.read_csv(f'{CALIB_DIR}/prevalence_by_age_sex.csv')
    condom_data = pd.read_csv(f'{DATA_DIR}/condom_use.csv')
    vmmc_data = pd.read_csv(f'{DATA_DIR}/vmmc_coverage.csv')

    plot_fit_dashboard(data, targets, phia_art, inc_calib, prev_calib, vmmc_data,
                       args.label, args.n_seeds, outdir=args.outdir)
    plot_network_dashboard(data, condom_data, args.label, args.n_seeds, outdir=args.outdir)

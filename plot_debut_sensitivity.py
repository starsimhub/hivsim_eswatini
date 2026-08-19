"""
Sensitivity analysis: sweep female debut age from 13–20 (male = female + 1),
10 seeds per scenario, plot mean incidence with uncertainty bands.
"""

import numpy as np
import starsim as ss
import stisim as sti
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

FIGURES_DIR = 'figures'
DEBUT_AGES_F = [13, 14, 15, 16, 17, 18, 19, 20]
N_SEEDS = 10
SD = 1


def make_debut_sim(debut_f, seed=1):
    """Build a sim with a specific female debut age (male = female + 1)."""
    from interventions import make_interventions
    from analyzers import hiv_epi

    debut_m = debut_f + 1

    sexual = sti.StructuredSexual(
        prop_f0=0.6, prop_f2=0.1,
        prop_m0=0.5, prop_m2=0.1,
        debut_pars_f=[debut_f, SD],
        debut_pars_m=[debut_m, SD],
        f1_conc=0.15, f2_conc=0.25,
        m1_conc=0.15, m2_conc=0.5,
        p_pair_form=0.5,
        condom_data=pd.read_csv('data/condom_use.csv'),
        fsw_shares=ss.bernoulli(p=0.10),
        client_shares=ss.bernoulli(p=0.20),
    )
    maternal = ss.MaternalNet()
    hiv = sti.HIV(
        beta_m2f=0.01, eff_condom=0.85,
        init_prev_data=pd.read_csv('data/init_prev_hiv.csv'),
        rel_init_prev=0.1,
    )
    simpars = dict(
        use_migration=True, rand_seed=seed, n_agents=10e3,
        start=1985, stop=2031, verbose=-1,
    )
    sim = sti.Sim(
        pars=simpars, datafolder='data/', demographics='eswatini',
        diseases=[hiv], networks=[sexual, maternal],
        interventions=make_interventions(), analyzers=[hiv_epi()],
        data=pd.read_csv('data/eswatini_hiv_calib.csv'),
    )
    return sim


def annual_mean(yearvec, vals):
    years_int = np.floor(yearvec).astype(int)
    unique_years = np.unique(years_int)
    annual_years = []
    annual_vals = []
    for yr in unique_years:
        mask = years_int == yr
        annual_years.append(yr)
        annual_vals.append(np.mean(vals[mask]))
    return np.array(annual_years), np.array(annual_vals)


# ── Run sims ─────────────────────────────────────────────────────────────────
results = {}
for debut_f in DEBUT_AGES_F:
    debut_m = debut_f + 1
    print(f'Running debut F={debut_f}, M={debut_m} ({N_SEEDS} seeds)...')
    seed_inc_f = []
    seed_inc_m = []
    for seed in range(1, N_SEEDS + 1):
        sim = make_debut_sim(debut_f, seed=seed)
        sim.run()
        epi = sim.analyzers['hiv_epi']
        yearvec = sim.t.yearvec
        yrs, ann_f = annual_mean(yearvec, np.array(epi.results['incidence_f_15_49']))
        _,   ann_m = annual_mean(yearvec, np.array(epi.results['incidence_m_15_49']))
        seed_inc_f.append(ann_f)
        seed_inc_m.append(ann_m)
    results[debut_f] = {
        'years': yrs,
        'inc_f': np.array(seed_inc_f),  # shape (N_SEEDS, n_years)
        'inc_m': np.array(seed_inc_m),
    }
    print(f'  Done.')

# ── Plot ─────────────────────────────────────────────────────────────────────
cmap = plt.cm.viridis
norm = plt.Normalize(vmin=min(DEBUT_AGES_F), vmax=max(DEBUT_AGES_F))

fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

for ax, sex_key, sex_label in [(axes[0], 'inc_f', 'Female'), (axes[1], 'inc_m', 'Male')]:
    for debut_f in DEBUT_AGES_F:
        r = results[debut_f]
        yrs = r['years']
        arr = r[sex_key] * 100  # (N_SEEDS, n_years)
        mean = np.mean(arr, axis=0)
        lo = np.percentile(arr, 10, axis=0)
        hi = np.percentile(arr, 90, axis=0)

        debut_m = debut_f + 1
        debut_val = debut_f if sex_label == 'Female' else debut_m
        color = cmap(norm(debut_f))
        ax.plot(yrs, mean, color=color, lw=2, label=f'{debut_val}')
        ax.fill_between(yrs, lo, hi, color=color, alpha=0.15)

    ax.set_xlabel('Year', fontsize=12)
    ax.set_title(f'{sex_label} HIV Incidence (15–49)', fontsize=13)
    ax.legend(title='Debut age', fontsize=9, title_fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(1990, 2030)

axes[0].set_ylabel('HIV incidence per 100 PY', fontsize=12)

fig.suptitle(f'HIV Incidence Sensitivity to Age of Sexual Debut (mean of {N_SEEDS} seeds)', fontsize=14)
plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/incidence_debut_sensitivity.png', dpi=150, bbox_inches='tight')
print(f'Saved {FIGURES_DIR}/incidence_debut_sensitivity.png')
plt.close()

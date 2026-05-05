"""
Diagnostic: verify sexual debut age distribution before and after parameter change.
Runs 1 sim with each setting and plots boxplots of debut age by sex at select years.
"""

import numpy as np
import starsim as ss
import stisim as sti
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run_sims import make_sim

FIGURES_DIR = 'figures'
SNAPSHOT_YEARS = [2000, 2010, 2016, 2021]


class DebutSnapshot(ss.Analyzer):
    """Capture debut age distribution at specified years."""
    def __init__(self, years=None, **kwargs):
        super().__init__(**kwargs)
        self.snapshot_years = years or SNAPSHOT_YEARS
        self.data = {}

    def init_results(self):
        super().init_results()
        self.define_results(ss.Result('_dummy', dtype=float, scale=False))

    def step(self):
        sim = self.sim
        year = float(sim.t.yearvec[sim.ti])
        yr = int(year)
        if yr not in self.snapshot_years:
            return
        # Only capture once per year
        if any((yr, s) in self.data for s in ['Female', 'Male']):
            return

        nw = sim.networks.structuredsexual
        ppl = sim.people
        alive = ppl.alive

        debuted = nw.participant & alive
        for sex_label, sex_bool in [('Female', ppl.female), ('Male', ppl.male)]:
            mask = debuted & sex_bool
            if mask.count() > 0:
                debut_ages = np.array(nw.debut[mask])
                self.data[(yr, sex_label)] = debut_ages


def get_snapshot(sim):
    """Retrieve the DebutSnapshot analyzer from a completed sim."""
    for a in sim.analyzers.values():
        if isinstance(a, DebutSnapshot):
            return a
    raise ValueError('DebutSnapshot not found in sim analyzers')


# ── Build OLD sim manually (override debut back to defaults) ─────────────────
def make_old_sim(seed=1, analyzers=None):
    """Same as make_sim but with old debut defaults."""
    sexual = sti.StructuredSexual(
        prop_f0=0.6, prop_f2=0.1,
        prop_m0=0.5, prop_m2=0.1,
        debut_pars_f=[20, 3],   # OLD default
        debut_pars_m=[21, 3],   # OLD default
        f1_conc=0.15, f2_conc=0.25,
        m1_conc=0.15, m2_conc=0.5,
        p_pair_form=0.5,
        condom_data=pd.read_csv('data/condom_use.csv'),
        fsw_shares=ss.bernoulli(p=0.10),
        client_shares=ss.bernoulli(p=0.20),
    )
    maternal = ss.MaternalNet()
    hiv = sti.HIV(beta_m2f=0.01, eff_condom=0.85,
                  init_prev_data=pd.read_csv('data/init_prev_hiv.csv'), rel_init_prev=0.1)
    from interventions import make_interventions
    from analyzers import hiv_epi
    interventions = make_interventions()
    default_analyzers = [hiv_epi()]
    if analyzers is not None:
        default_analyzers += list(analyzers)
    simpars = dict(use_migration=True, rand_seed=seed, n_agents=10e3, start=1985, stop=2031, verbose=-1)
    sim = sti.Sim(pars=simpars, datafolder='data/', demographics='eswatini',
                  diseases=[hiv], networks=[sexual, maternal],
                  interventions=interventions, analyzers=default_analyzers,
                  data=pd.read_csv('data/eswatini_hiv_calib.csv'))
    return sim

print('Running sim with OLD debut (F=20, M=21)...')
sim_old = make_old_sim(seed=1, analyzers=[DebutSnapshot()])
sim_old.run()
analyzer_old = get_snapshot(sim_old)
print(f'Done. Captured {len(analyzer_old.data)} snapshots.')

# ── Run with NEW values (uses current make_sim with F=17.5, M=18.5) ─────────
print('Running sim with NEW debut (F=17.5, M=18.5)...')
sim_new = make_sim(seed=1, verbose=-1, analyzers=[DebutSnapshot()])
sim_new.run()
analyzer_new = get_snapshot(sim_new)
print(f'Done. Captured {len(analyzer_new.data)} snapshots.')

# ── Plot ─────────────────────────────────────────────────────────────────────

fig, axes = plt.subplots(1, len(SNAPSHOT_YEARS), figsize=(4 * len(SNAPSHOT_YEARS), 6), sharey=True)

for i, yr in enumerate(SNAPSHOT_YEARS):
    ax = axes[i]
    box_data = []
    labels = []
    colors = []

    for sex_label, color in [('Female', '#d62728'), ('Male', '#1f77b4')]:
        # Old
        key = (yr, sex_label)
        if key in analyzer_old.data:
            box_data.append(analyzer_old.data[key])
            labels.append(f'{sex_label}\nOld')
            colors.append(color)

        # New
        if key in analyzer_new.data:
            box_data.append(analyzer_new.data[key])
            labels.append(f'{sex_label}\nNew')
            colors.append(color)

    if box_data:
        bp = ax.boxplot(box_data, tick_labels=labels, patch_artist=True, widths=0.6,
                        showfliers=False, medianprops=dict(color='black', linewidth=2))
        for j, patch in enumerate(bp['boxes']):
            c = colors[j]
            alpha = 0.3 if 'Old' in labels[j] else 0.7
            patch.set_facecolor(c)
            patch.set_alpha(alpha)

    ax.set_title(f'Year {yr}', fontsize=13)
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_ylim(10, 35)

axes[0].set_ylabel('Age at sexual debut (years)', fontsize=12)

# Legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='grey', alpha=0.3, label='Old (F=20, M=21)'),
    Patch(facecolor='grey', alpha=0.7, label='New (F=17.5, M=18.5)'),
    Patch(facecolor='#d62728', alpha=0.5, label='Female'),
    Patch(facecolor='#1f77b4', alpha=0.5, label='Male'),
]
fig.legend(handles=legend_elements, loc='upper center', ncol=4, fontsize=10, bbox_to_anchor=(0.5, 1.0))

fig.suptitle('Sexual Debut Age Distribution: Old vs New Parameters', fontsize=14, y=1.04)
plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/debut_age_check.png', dpi=150, bbox_inches='tight')
print(f'Saved {FIGURES_DIR}/debut_age_check.png')
plt.close()

# Print summary stats
print('\nSummary statistics (year 2021):')
for label, analyzer in [('OLD', analyzer_old), ('NEW', analyzer_new)]:
    for sex in ['Female', 'Male']:
        key = (2021, sex)
        if key in analyzer.data:
            d = analyzer.data[key]
            print(f'  {label} {sex}: median={np.median(d):.1f}, '
                  f'IQR=[{np.percentile(d,25):.1f}, {np.percentile(d,75):.1f}], '
                  f'n={len(d)}')

# ── Density plot for year 2010 ───────────────────────────────────────────────
from scipy.stats import gaussian_kde

fig2, ax2 = plt.subplots(figsize=(8, 5))
age_grid = np.linspace(10, 35, 300)

styles = [
    ('OLD', analyzer_old, '--'),
    ('NEW', analyzer_new, '-'),
]
sex_colors = {'Female': '#d62728', 'Male': '#1f77b4'}

for label, analyzer, ls in styles:
    for sex, color in sex_colors.items():
        key = (2010, sex)
        if key in analyzer.data:
            kde = gaussian_kde(analyzer.data[key], bw_method=0.3)
            ax2.plot(age_grid, kde(age_grid), color=color, ls=ls, lw=2,
                     label=f'{sex} — {label}')

ax2.set_xlabel('Age at sexual debut (years)', fontsize=12)
ax2.set_ylabel('Density', fontsize=12)
ax2.set_title('Sexual Debut Age Distribution (Year 2010)', fontsize=14)
ax2.legend(fontsize=11)
ax2.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/debut_age_density_2010.png', dpi=150, bbox_inches='tight')
print(f'Saved {FIGURES_DIR}/debut_age_density_2010.png')
plt.close()

# ── Incidence comparison plot ────────────────────────────────────────────────

def get_incidence(sim, key):
    """Extract annualized incidence timeseries from the hiv_epi analyzer."""
    epi = sim.analyzers['hiv_epi']
    vals = np.array(epi.results[key])
    return vals

yearvec = sim_old.t.yearvec

def annual_mean(yearvec, vals):
    """Average monthly values into annual means."""
    years_int = np.floor(yearvec).astype(int)
    unique_years = np.unique(years_int)
    annual_years = []
    annual_vals = []
    for yr in unique_years:
        mask = years_int == yr
        annual_years.append(yr)
        annual_vals.append(np.mean(vals[mask]))
    return np.array(annual_years), np.array(annual_vals)

fig3, ax3 = plt.subplots(figsize=(9, 5))

for sim, label, ls in [(sim_old, 'Old (F=20, M=21)', '--'), (sim_new, 'New (F=17.5, M=18.5)', '-')]:
    inc_f = get_incidence(sim, 'incidence_f_15_49')
    inc_m = get_incidence(sim, 'incidence_m_15_49')
    yrs_f, ann_f = annual_mean(yearvec, inc_f)
    yrs_m, ann_m = annual_mean(yearvec, inc_m)
    ax3.plot(yrs_f, ann_f * 100, color='#d62728', ls=ls, lw=2, label=f'Female — {label}')
    ax3.plot(yrs_m, ann_m * 100, color='#1f77b4', ls=ls, lw=2, label=f'Male — {label}')

ax3.set_xlabel('Year', fontsize=12)
ax3.set_ylabel('HIV incidence per 100 PY (ages 15–49)', fontsize=12)
ax3.set_title('HIV Incidence by Sex: Old vs New Sexual Debut Parameters', fontsize=14)
ax3.legend(fontsize=10)
ax3.grid(True, alpha=0.3)
ax3.set_xlim(1990, 2030)
plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/incidence_debut_comparison.png', dpi=150, bbox_inches='tight')
print(f'Saved {FIGURES_DIR}/incidence_debut_comparison.png')
plt.close()

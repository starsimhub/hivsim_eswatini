# hivsim_eswatini

Agent-based HIV model for Eswatini, built on [STIsim](https://github.com/starsimhub/stisim) and [Starsim](https://github.com/starsimhub/starsim).

## Setup

```bash
pip install stisim
```

## Usage

```bash
# Run a single simulation
python run_sims.py

# Generate calibration data from raw UNAIDS/PHIA sources
python utils.py

# Generate network structure data (needed for Fig 3)
python run_network_data.py

# Plot figures
python plot_fig1_hiv_calibration.py   # HIV epi: prevalence, infections, ART, deaths
python plot_fig2_hiv_by_age.py        # HIV prevalence by age/sex vs PHIA surveys
python plot_fig3_network.py           # Network structure: partners, mixing, risk groups
```

## Structure

| File | Description |
|------|-------------|
| `run_sims.py` | `make_sim()` and single-sim runner |
| `interventions.py` | HIV testing, ART, PrEP |
| `analyzers.py` | `NetworkSnapshot` for network figure |
| `utils.py` | Data formatting (UNAIDS/PHIA → calibration CSV) |
| `data/` | Processed data (calibration targets, ART coverage, condom use) |
| `raw_data/` | Source data (UNAIDS Excel, PHIA surveys, demographics) |
| `results/` | Simulation outputs |
| `figures/` | Generated figures |

## Data sources

- **UNAIDS**: HIV prevalence, incidence, deaths, PLHIV counts (1990–2024)
- **PHIA surveys**: Age/sex-specific HIV prevalence (2007, 2011, 2016, 2021)
- **ART coverage**: National ART numbers (2010–2024)

# hivsim_eswatini

Agent-based HIV model for Eswatini, built on [STIsim](https://github.com/starsimhub/stisim) and [Starsim](https://github.com/starsimhub/starsim).

The model includes a structured sexual network with risk groups, HIV transmission, testing (FSW-targeted, general population, low-CD4), ART, and PrEP. It is calibrated to UNAIDS estimates and PHIA household survey data.


## Prerequisites

1. **Python 3.9+** (we recommend [Miniforge](https://github.com/conda-forge/miniforge))

2. **Install Starsim and STIsim** by cloning and installing locally:
   ```bash
   git clone https://github.com/starsimhub/starsim.git
   cd starsim
   pip install -e .
   cd ..

   git clone https://github.com/starsimhub/stisim.git
   cd stisim
   pip install -e .
   cd ..
   ```

   Installing with `-e` (editable mode) means changes to the source code take effect immediately without reinstalling.

3. **Clone this repo**:
   ```bash
   git clone https://github.com/starsimhub/hivsim_eswatini.git
   cd hivsim_eswatini
   ```

4. **Raw data** (not tracked in git): place the `raw_data/` folder in the repo root. This contains UNAIDS Excel files, PHIA survey CSVs, and population data. Ask the project lead for access.

5. **IDE setup**: We recommend [VS Code](https://code.visualstudio.com/) with the Python and [Claude Code](https://marketplace.visualstudio.com/items?itemName=anthropics.claude-code) extensions. Open the `hivsim_eswatini` folder as your workspace. You can run any script by opening it and pressing the play button, or by using the integrated terminal.


## Workflow

The scripts are designed to be run in order. Each step produces outputs used by the next.

### Step 1: Generate calibration data (only needed once, or when raw data changes)

Run `utils.py`.

This reads raw UNAIDS and PHIA data from `raw_data/` and produces `data/eswatini_hiv_calib.csv` — a single CSV with one row per year and columns matching the sim's result names. The calibrator uses this file to compute goodness-of-fit.

### Step 2: Run a single simulation (for quick testing)

Run `run_sims.py`.

Runs one sim with default parameters and saves results to `results/eswatini_sim.df`. Also opens a plot window showing HIV dynamics. Use this to check that the model is working before calibrating.

### Step 3: Calibrate the model

Run `run_calibrations.py`.

Runs 1000 Optuna trials, each testing a different combination of parameters (transmission rate, condom efficacy, network structure, ART duration). Saves:
- `results/eswatini_calib.obj` — full calibration object
- `results/eswatini_pars.df` — dataframe of best-fit parameter sets

The calibration parameters and their ranges are defined in `make_calibration()` inside `run_calibrations.py`. To add or modify calibrated parameters, edit the `calib_pars` dictionary there.

### Step 4: Run multi-sim with calibrated parameters

Run `run_msim.py`.

Takes the top 200 parameter sets from calibration, runs a sim for each, and generates percentile statistics (median, 10-90% credible intervals). Saves `results/eswatini_calib_stats.df`.

### Step 5: Generate figures

Run the plotting scripts:
- `plot_fig1_hiv_calibration.py` — Fig 1: HIV epi (6 panels) — model bands vs UNAIDS data
- `plot_fig2_hiv_by_age.py` — Fig 2: prevalence by age/sex vs PHIA surveys
- `run_network_data.py` — extract network data (needed once before Fig 3)
- `plot_fig3_network.py` — Fig 3: network structure (6 panels)

Figures are saved to `figures/` at 200 dpi.


## File descriptions

### Scripts

| File | Purpose |
|------|---------|
| `run_sims.py` | `make_sim()` — the central function that builds the simulation. All other scripts import this. Also runs a single sim when executed directly. |
| `run_calibrations.py` | Calibrate model parameters to match UNAIDS/PHIA data using Optuna. |
| `run_msim.py` | Run the top calibrated parameter sets and generate summary statistics with credible intervals. |
| `run_network_data.py` | Run sims with network analyzers and extract structure data for plotting. |
| `plot_fig1_hiv_calibration.py` | HIV calibration figure: prevalence, infections, ART, PLHIV, deaths, population. Supports single-sim (line) or multi-sim (median + band). Set `use_calib = True/False` in the `__main__` block. |
| `plot_fig2_hiv_by_age.py` | HIV prevalence by age and sex, compared to PHIA survey data (2007, 2011, 2016, 2021). |
| `plot_fig3_network.py` | Network structure: lifetime partners, age mixing, risk groups, debut age, partnerships, condom use. |

### Modules

| File | Purpose |
|------|---------|
| `interventions.py` | HIV testing strategies (FSW, general pop, low-CD4), ART, PrEP. |
| `analyzers.py` | `hiv_epi` — tracks age-specific prevalence/incidence/counts for UNAIDS and PHIA age ranges. `NetworkSnapshot` — captures network properties at a point in time. |
| `utils.py` | `format_calibration_data()` — reads raw UNAIDS Excel + PHIA CSVs and produces the calibration target CSV. `set_font()` — plot styling helper. |

### Data

| File | Contents |
|------|----------|
| `data/eswatini_hiv_calib.csv` | Calibration targets (generated by `utils.py`). Columns prefixed `hiv.` are native sim results; columns prefixed `hiv_epi.` come from the `hiv_epi` analyzer. |
| `data/n_art.csv` | Number of people on ART by year (1990–2024). Used by the ART intervention and for plotting. |
| `data/condom_use.csv` | Condom use probability by partnership type and year. Input to the sexual network. |
| `data/init_prev_hiv.csv` | Initial HIV prevalence by risk group, sex, and SW status. Seeds the epidemic at sim start. |
| `data/eswatini_*.csv` | Demographics: age distribution, fertility (ASFR), death rates, migration. Loaded automatically when `demographics='eswatini'`. |

### Raw data (not in git)

| Folder | Contents |
|--------|----------|
| `raw_data/UNAIDS/` | 16 Excel files with country-level HIV estimates (prevalence, incidence, deaths, PLHIV by age). |
| `raw_data/SWAZILAND_nationalprevalence_*.csv` | PHIA survey data: HIV prevalence by age, sex, and year (2007, 2011, 2016, 2021). |
| `raw_data/Swaziland_Incidence_Data2.csv` | PHIA incidence estimates by sex (2011, 2016, 2021). |
| `raw_data/Swaziland_Population_Counts.xlsx` | Census population by age/sex/year (2000–2015). |
| `raw_data/n_art.xlsx`, `raw_data/p_art.xlsx` | ART coverage source data. |


## Key concepts

**How the sim works**: The model creates ~10,000 agents representing the population of Eswatini. Each agent has an age, sex, and risk group. Agents form sexual partnerships through a structured network, and HIV transmits through these contacts. The sim runs in monthly timesteps from 1985 to ~2030. Results are scaled up to the national population automatically.

**Calibration data column naming**: The calibration CSV column names must exactly match what `sim.to_df(sep='.')` produces. Native disease results use the `hiv.` prefix (e.g., `hiv.prevalence_15_49`). Results from the `hiv_epi` analyzer use the `hiv_epi.` prefix (e.g., `hiv_epi.n_infected_15_49`). If you add a new indicator, make sure the column name in the CSV matches the result name in the sim output.

**Age bins**: The HIV disease module natively tracks 5-year bins from 15-35, plus wider bins (35-50, 50-65, 65-100). The `hiv_epi` analyzer adds finer 5-year bins above 35 (to match PHIA surveys) and broader UNAIDS ranges (0-14, 10-19, 15-24, 15-49, 15-100).

**Adding a new calibration target**: (1) Add the raw data to `raw_data/`, (2) update `format_calibration_data()` in `utils.py` to include it with the correct column name, (3) re-run `utils.py`, (4) if the indicator isn't already a sim result, add it to the `hiv_epi` analyzer in `analyzers.py`.


## Data sources

- **UNAIDS**: HIV prevalence, incidence, deaths, PLHIV counts (1990–2024)
- **PHIA surveys**: Age/sex-specific HIV prevalence with 95% CIs (2007, 2011, 2016, 2021)
- **ART coverage**: National ART numbers (1990–2024)
- **Census**: Population counts by age/sex (2000–2015)

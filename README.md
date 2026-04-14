# HIVSim Eswatini Model

Agent-based HIV model for Eswatini, built on [STIsim](https://github.com/starsimhub/stisim) and [Starsim](https://github.com/starsimhub/starsim).

The model includes a structured sexual network with risk groups, HIV transmission, testing (FSW-targeted, general population, low-CD4), ART, and PrEP. It is calibrated to UNAIDS estimates and PHIA household survey data.

## [Project Docs](hivsim_swz_docs.md)


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

### Step 1: Format the input data for HIVsim (only needed once, or when raw data changes)

Run `utils.py`.

This reads raw UNAIDS and PHIA data from `raw_data/` and produces `data/eswatini_hiv_calib.csv` — a single CSV with one row per year and columns matching the sim's result names. The calibrator uses this file to compute goodness-of-fit.

### Step 2: Run a single simulation (for quick testing)

Run `run_sims.py`.

Runs one sim with default parameters and saves results to `results/eswatini_sim.df`. Also opens a plot window showing HIV dynamics. Use this to check that the model is working before calibrating.

### Step 3: Calibrate the model (on VM/HPC)

Calibration is computationally intensive (1000 trials, each running a full sim) and must be done on a VM or HPC, not locally.

#### 3a. Push your local changes

Before moving to the VM, make sure all your code changes are committed and pushed. In VS Code:

1. Open the **Source Control** panel (click the branch icon in the left sidebar, or press `Ctrl+Shift+G`)
2. You'll see a list of changed files. Review them, then type a commit message (e.g., "update calibration parameters") in the text box at the top
3. Click the **checkmark** button to commit
4. Click the **"..." menu → Push** (or click the sync/upload icon) to push to GitHub

Or from the terminal:
```bash
git add -A
git commit -m "update calibration parameters"
git push
```

#### 3b. Set up the VM environment (first time only)

SSH into the VM, then create a conda environment and install dependencies:

```bash
# Create and activate environment
conda create -n hiv python=3.11 -y
conda activate hiv

# Clone and install Starsim + STIsim
git clone https://github.com/starsimhub/starsim.git
cd starsim && pip install -e . && cd ..

git clone https://github.com/starsimhub/stisim.git
cd stisim && pip install -e . && cd ..

# Clone this repo
git clone https://github.com/starsimhub/hivsim_eswatini.git
cd hivsim_eswatini
```

On subsequent sessions, just activate and pull:
```bash
conda activate hiv
cd hivsim_eswatini
git pull
```

#### 3c. Run calibration

```bash
python run_calibrations.py
```

This runs 1000 Optuna trials, each testing a different combination of parameters (transmission rate, condom efficacy, network structure, ART duration). It takes a while — typically 1-4 hours depending on the machine. When it finishes, it saves:
- `results/eswatini_calib.obj` — full calibration object (large, but useful for debugging)
- `results/eswatini_pars.df` — dataframe of best-fit parameter sets (small, this is what the next step uses)

The calibration parameters and their ranges are defined in `make_calibration()` inside `run_calibrations.py`. To add or modify calibrated parameters, edit the `calib_pars` dictionary there.

#### 3d. Run multi-sim with calibrated parameters

Still on the VM, run:
```bash
python run_msim.py
```

This takes the top 200 parameter sets from calibration, runs a sim for each, and generates percentile statistics (median, 10-90% credible intervals). Saves `results/eswatini_calib_stats.df`.

#### 3e. Push results back

The result files are small enough to live in git. Commit and push them from the VM:
```bash
git commit -m "new calibration results"
git push
```

These two files are what the plotting scripts need:
- `eswatini_pars.df` — the best-fit parameter sets (used by `run_msim.py` if you want to re-run locally)
- `eswatini_calib_stats.df` — the percentile summary statistics (used by the plot scripts)

Do not add the `eswatini_calib.obj` file, it's larger and only needed for debugging. If needed, you can copy it to your local machine (ask for help).

#### 3f. Pull results locally

Back on your local machine, pull the new results:

In VS Code: open the **Source Control** panel → click **"..." menu → Pull**

Or from the terminal:
```bash
git pull
```

### Step 4: Generate figures (locally)

Now that you have `results/eswatini_calib_stats.df` locally, run the plotting scripts:
- `plot_fig1_hiv_calibration.py` — Fig 1: HIV epi (6 panels) — model bands vs UNAIDS data
- `plot_fig2_hiv_by_age.py` — Fig 2: prevalence by age/sex vs PHIA surveys
- `run_network_data.py` — extract network data (needed once before Fig 3)
- `plot_fig3_network.py` — Fig 3: network structure (6 panels)

Figures are saved to `figures/` at 200 dpi. Figure files are not added to Git.


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


## Git hygiene

GitHub has a hard limit of 100 MB per file, and repositories get slow and painful to work with if they accumulate large files. As a rule of thumb: **only commit files that are small (under ~10 MB) and that other people need**.

### What's in `.gitignore`

The file `.gitignore` in the repo root tells git to ignore certain files and folders. Currently it contains:

```
raw_data/       # Large source data (UNAIDS Excel files, PHIA CSVs, population data)
raw_results/    # Any large intermediate outputs
__pycache__/    # Python bytecode cache (auto-generated, never commit)
.DS_Store       # macOS metadata files (auto-generated, never commit)
*.pyc           # Compiled Python files
.vscode/        # VS Code workspace settings (personal to each developer)
```

Anything matching these patterns is invisible to git — it won't show up in the Source Control panel and `git add -A` won't include it.

### What's safe to commit

- **Code** (`.py` files) — always
- **Small data files** in `data/` (CSVs, typically < 1 MB) — yes
- **Calibration result files** (`results/eswatini_pars.df`, `results/eswatini_calib_stats.df`) — yes, these are small (< 5 MB) and needed by the plotting scripts
- **Figures** in `figures/` (PNGs, typically < 1 MB each) — yes
- **`results/network_data.obj`** — yes, small
- **Font files** in `assets/` — yes
- **README, .gitignore** — yes

### What NOT to commit

- **`results/eswatini_calib.obj`** — this is the full calibration object and can be 50-200 MB. Do not commit it. If you need to share it, copy it directly (e.g., via `scp`)
- **Anything in `raw_data/`** — already in `.gitignore`, but be careful not to force-add it
- **Large simulation outputs** — if you save full sim objects (`.obj` files > 10 MB), don't commit them

### How to modify `.gitignore`

If you create new files or folders that should be excluded from git, add them to `.gitignore`. For example, if you start saving large scenario outputs to a `scenarios/` folder:

1. Open `.gitignore` in VS Code
2. Add a new line: `scenarios/`
3. Save the file
4. Commit the updated `.gitignore`

If you accidentally commit a large file, ask for help — removing it from git history requires special commands.


## Data sources

- **UNAIDS**: HIV prevalence, incidence, deaths, PLHIV counts (1990–2024)
- **PHIA surveys**: Age/sex-specific HIV prevalence with 95% CIs (2007, 2011, 2016, 2021)
- **ART coverage**: National ART numbers (1990–2024)
- **Census**: Population counts by age/sex (2000–2015)

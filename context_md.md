# HIVsim Eswatini — Project Context

## Who I am
I'm an experienced HIV epidemiologist and modeler with deep subject matter expertise. I'm an R user with no Python experience — please explain things in plain English and use R analogies where helpful. I'm very comfortable with modeling concepts but not Python syntax.

A key part of my role is to critically evaluate the model: identifying gaps, challenging assumptions, and adding new features where needed. Claude should feel free to engage with me at a high level on the epidemiology and modeling, while handling the Python implementation.

## The model
Agent-based HIV model for Eswatini, built on STIsim and Starsim. The model includes a structured sexual network with risk groups, HIV transmission, testing, ART, and PrEP. It is calibrated to UNAIDS estimates and PHIA household survey data.

## Overall goal
Proof of concept: does this new STIsim-based model recapitulate the HIV dynamics we would expect for Eswatini, consistent with results from a prior model? If yes, we then want to validate the model against new 2023 data that was not used in calibration.

## Project phases

### Phase 1 (current): Model familiarization
- Understand all model components, assumptions, and parameters
- Understand how the sexual network is structured
- Understand how interventions (ART, PrEP, testing) are implemented
- Understand what outputs the model produces

### Phase 2: Calibration and plotting
- Calibrate the model to UNAIDS and PHIA data
- Generate bespoke plots to explore model output
- Compare outputs to prior model results

### Phase 3: Validation
- Validate the model against new 2023 data not used in calibration
- Assess how well the model predicts out-of-sample

## Key files
- `run_sims.py` — builds and runs a single simulation (start here)
- `run_calibrations.py` — calibrates model parameters using Optuna (run on VM/HPC)
- `run_msim.py` — runs top calibrated parameter sets and generates summary statistics
- `interventions.py` — HIV testing, ART, PrEP
- `analyzers.py` — tracks age-specific outputs for calibration and plotting
- `utils.py` — formats raw data into calibration targets

## Notes
- Calibration is computationally intensive and must be run on a VM/HPC, not locally
- Raw data (`raw_data/`) is not tracked in git — stored locally only
- Result files (`eswatini_pars.df`, `eswatini_calib_stats.df`) are small and tracked in git

# Parameter Comparison: EMOD vs STIsim (Eswatini HIV Model)

Side-by-side comparison using calibrated EMOD parameters from Akullian et al., Lancet HIV 2020 (Table A2, N=250 best-fitting sets) and the STIsim/Starsim Eswatini model. Dynamic EMOD parameters show median (IQR). STIsim values marked "default" are from source code and not overridden in run_sims.py.

---

## Calibrated vs Fixed Parameters

Understanding which parameters are free to vary during calibration is critical for interpreting differences between models.

### EMOD: 24 calibrated parameters

All parameters marked **[calibrated]** in the tables below were free to vary. The 250 best-fitting parameter sets from 1000+ Optuna trials are reported as median (IQR).

### STIsim: 6 calibrated parameters

Only these 6 parameters are free to vary during calibration (`run_calibrations.py`):

| # | Parameter | Category | Low | High | Guess | Notes |
|---|-----------|----------|-----|------|-------|-------|
| 1 | `beta_m2f` | Transmission | 0.002 | 0.014 | 0.006 | Male-to-female per-act probability |
| 2 | `eff_condom` | Transmission | 0.5 | 0.9 | 0.75 | Condom efficacy (not condom *use*) |
| 3 | `rel_dur_on_art` | Treatment | 1.0 | 20.0 | 8.0 | Multiplier on base ART duration (3yr x 8 = ~24yr) |
| 4 | `prop_f0` | Network | 0.55 | 0.9 | 0.7 | Proportion of females in low-risk group |
| 5 | `prop_m0` | Network | 0.55 | 0.8 | 0.65 | Proportion of males in low-risk group |
| 6 | `m1_conc` | Network | 0.05 | 0.3 | 0.15 | Male medium-risk concurrency probability |

**Everything else is fixed**, including: sexual debut age, condom use time series, all testing probabilities, FSW/client proportions, all other concurrency parameters, all disease progression parameters, ART efficacy, initial prevalence, and the rel_beta_f2m ratio. The `rel_init_prev` parameter was previously calibrated but is now commented out (fixed at 0.1).

---

## 1. TRANSMISSION

| Concept | EMOD Parameter | EMOD Value | STIsim Parameter | STIsim Value | Calibrated? | Notes |
|---------|---------------|------------|-----------------|-------------|-------------|-------|
| Base transmission | Base_Infectivity | **0.00233 (0.00231-0.00234)** [calibrated] | beta_m2f | 0.01 | **STIsim: YES** | STIsim ~4x higher base rate; EMOD compensates with larger stage multipliers |
| Acute infectiousness | Acute_Stage_Infectivity_Multiplier | 26 | rel_trans_acute | Normal(6, 0.5) | Both fixed | EMOD 26x vs STIsim ~6x -- major structural difference |
| Acute duration | Acute_Duration_In_Months | 3 | dur_acute | lognorm(mean=3mo, SD=1mo) | Both fixed | Match |
| AIDS infectiousness | AIDS_Stage_Infectivity_Multiplier | 4.5 | rel_trans_falling | Normal(8, 0.5) | Both fixed | STIsim higher for late-stage |
| AIDS duration | AIDS_Duration_In_Months | 9 | dur_falling | lognorm(mean=3yr, SD=1yr) | Both fixed | STIsim much longer AIDS phase |
| F-to-M relative risk (young <25) | Male_To_Female_Relative_Infectivity_Young | **4.894 (4.747-5.041)** [calibrated] | rel_beta_f2m | 0.5 | STIsim: fixed | EMOD: young women ~5x more susceptible. STIsim: flat ratio, no age dependence. **No mechanism in STIsim source to vary by age -- would require code change.** |
| F-to-M relative risk (old 25+) | Male_To_Female_Relative_Infectivity_Old | **2.844 (2.727-2.958)** [calibrated] | rel_beta_f2m | 0.5 (same) | STIsim: fixed | EMOD found significant age gradient (5x young vs 3x old) |
| Condom efficacy | Condom_Transmission_Blocking_Probability | 0.8 | eff_condom | 0.85 | **STIsim: YES** (range 0.5-0.9) | Note: this is the *efficacy* of condoms when used, not the probability of use |
| Circumcision efficacy | Circumcision_Reduced_Acquire | 0.6 | eff_circ | 0.6 (default) | Both fixed | Match |
| STI co-infection | STI_Coinfection_Acquisition_Multiplier | 5.5 | -- | Not modeled | -- | EMOD only |

### ART Efficacy (Transmission Reduction on Treatment)

| Concept | EMOD Parameter | EMOD Value | STIsim Parameter | STIsim Value | Calibrated? | Notes |
|---------|---------------|------------|-----------------|-------------|-------------|-------|
| ART viral suppression | ART_Viral_Suppression_Multiplier | 0.08 (92% reduction) | art_efficacy | 0.96 (96% reduction) | Both fixed | STIsim slightly more optimistic. EMOD sensitivity analysis (Fig A7) showed incidence highly sensitive to this parameter (range 80-100% gave 0.4-1.8 per 100py) |
| Time to full ART efficacy | -- | Not modeled (immediate) | time_to_art_efficacy | 6 months (linear ramp) | STIsim: fixed | STIsim ramps up from 0 to 96% over 6 months after ART initiation |
| MTCT | Maternal_Infection_Transmission_Probability | 0.30 per pregnancy | beta_m2c | 0.025 per month | Both fixed | Different units |
| MTCT on ART | Maternal_Transmission_ART_Multiplier | 0.03334 | -- | Not explicitly set | -- | |

### Effective per-act transmission comparison

To compare the "effective" per-act transmission during each stage:

| Stage | EMOD (base x multiplier) | STIsim | Ratio |
|-------|--------------------------|--------|-------|
| Chronic/latent | 0.00233 x 1 = **0.0023** | 0.01 x 1 = **0.01** | STIsim 4.3x higher |
| Acute | 0.00233 x 26 = **0.061** | 0.01 x 6 = **0.06** | ~Match |
| AIDS | 0.00233 x 4.5 = **0.010** | 0.01 x 8 = **0.08** | STIsim 8x higher |
| On ART | 0.00233 x 0.08 = **0.00019** | 0.01 x (1-0.96) = **0.0004** | STIsim 2x higher |

Key insight: acute-phase transmission is similar between models (~0.06), but STIsim has much higher chronic and AIDS-stage transmission. This means STIsim relies less on acute-phase transmission to drive the epidemic.

---

## 2. DISEASE PROGRESSION & MORTALITY

| Concept | EMOD Parameter | EMOD Value | STIsim Parameter | STIsim Value | Calibrated? | Notes |
|---------|---------------|------------|-----------------|-------------|-------------|-------|
| Post-infection CD4 | CD4_Post_Infection_Weibull_Scale | 560.43 | cd4_start | Normal(800, 50) | Both fixed | STIsim starts higher |
| CD4 heterogeneity | CD4_Post_Infection_Weibull_Heterogeneity | 0.2756 | -- | Fixed SD=50 | Both fixed | EMOD has more CD4 variability |
| CD4 at death | CD4_At_Death_LogLogistic_Scale | 2.96 | -- | AIDS defined at CD4<200 | Both fixed | Different approach |
| Adult survival (intercept) | HIV_Adult_Survival_Scale_Parameter_Intercept | 21.182 | dur_latent | lognorm(mean=10yr, SD=3yr) | Both fixed | EMOD age-dependent; STIsim fixed |
| Adult survival (slope) | HIV_Adult_Survival_Scale_Parameter_Slope | -0.2717 | -- | No age dependence | -- | Older people die faster in EMOD |
| Adult survival (shape) | HIV_Adult_Survival_Shape_Parameter | 2 | -- | -- | -- | |
| Max age for survival | HIV_Age_Max_for_Adult_Age_Dependent_Survival | 50 | -- | -- | -- | EMOD caps age effect at 50 |
| Child rapid progressor | HIV_Child_Survival_Rapid_Progressor_Fraction | 0.57 | -- | Not modeled | -- | EMOD only |
| Child rapid progressor rate | HIV_Child_Survival_Rapid_Progressor_Rate | 1.52 yr | -- | -- | -- | |
| Child slow progressor scale | HIV_Child_Survival_Slow_Progressor_Scale | 16 | -- | -- | -- | |
| Child slow progressor shape | HIV_Child_Survival_Slow_Progressor_Shape | 2.7 | -- | -- | -- | |
| Days symptomatic to death (scale) | Days_Between_Symptomatic_And_Death_Weibull_Scale | 618.34 | -- | CD4-stratified mortality | Both fixed | Different mechanism |
| Days symptomatic to death (shape) | Days_Between_Symptomatic_And_Death_Weibull_Heterogeneity | 0.5 | -- | -- | -- | |

### STIsim CD4-stratified annual mortality (no EMOD equivalent)

| CD4 Range | Annual mortality | Calibrated? |
|-----------|----------------|-------------|
| >1000 | 0.3% | Fixed |
| 500-999 | 0.3% | Fixed |
| 350-499 | 0.5% | Fixed |
| 200-349 | 1% | Fixed |
| 50-199 | 5% | Fixed |
| <50 | 30% | Fixed |

---

## 3. SEXUAL NETWORK & BEHAVIOR

### Sexual Debut -- WHY IS STIsim SO HIGH?

**Root cause:** STIsim uses a package default of `lognorm(mean=20, SD=3)` for females and `lognorm(mean=21, SD=3)` for males. These defaults were **never overridden** in the Eswatini `run_sims.py` and are **not included as a calibration parameter**. They are simply the stisim package defaults, likely intended as a generic starting point rather than a country-specific estimate.

By contrast, the EMOD model calibrated debut age as a dynamic parameter and found significantly earlier ages (~16yr F, ~17yr M) -- consistent with DHS data for eSwatini showing median age at first sex around 17-18 for women and 18-19 for men.

**Impact:** A 3-4 year delay in debut means agents enter the sexual network later, reducing the window of HIV exposure -- especially for young women aged 15-24 who have the highest incidence in eSwatini. This likely forces the calibration to compensate elsewhere (higher beta, higher condom use).

**Fix:** Override in `run_sims.py` by passing `debut_pars_f=[16.3, 2]` and `debut_pars_m=[17.5, 1.5]` to the `StructuredSexual` constructor. Or add debut parameters to the calibration set.

| Concept | EMOD Parameter | EMOD Value | STIsim Parameter | STIsim Value | Calibrated? | Notes |
|---------|---------------|------------|-----------------|-------------|-------------|-------|
| Female debut (scale) | Sexual_Debut_Age_Female_Weibull_Scale | **16.302 (16.166-16.396)** [calibrated] | debut_pars_f | lognorm(mean=20, SD=3) | **STIsim: FIXED (package default)** | 3-4 year gap |
| Female debut (shape) | Sexual_Debut_Age_Female_Weibull_Heterogeneity | **0.309 (0.293-0.322)** [calibrated] | -- | SD=3 | STIsim: fixed | |
| Male debut (scale) | Sexual_Debut_Age_Male_Weibull_Scale | **17.499 (17.357-17.699)** [calibrated] | debut_pars_m | lognorm(mean=21, SD=3) | **STIsim: FIXED (package default)** | 3-4 year gap |
| Male debut (shape) | Sexual_Debut_Age_Male_Weibull_Heterogeneity | **0.042 (0.04-0.05)** [calibrated] | -- | SD=3 | STIsim: fixed | EMOD males have very tight debut distribution |
| Minimum debut age | Sexual_Debut_Age_Min | 13 | -- | Not explicitly set | -- | |

### Risk Groups

| Concept | EMOD Parameter | EMOD Value | STIsim Parameter | STIsim Value | Calibrated? | Notes |
|---------|---------------|------------|-----------------|-------------|-------------|-------|
| Proportion low risk | Proportion_Low_Risk | **0.73 (0.721-0.742)** [calibrated] | prop_f0 / prop_m0 | F: 0.60, M: 0.50 | **STIsim: YES** | EMOD higher LR; STIsim sex-specific |
| Epidemic seed year | Seed_Year | **1982.7 (1982.4-1983.2)** [calibrated] | start (+ init_prev) | 1985 with CSV-based seeding | STIsim: fixed | EMOD seeds ~1983; STIsim starts at 1985 |
| FSW proportion | -- | 3% enrollment at debut | fsw_shares | Bernoulli(0.10) | STIsim: fixed | STIsim 3x higher FSW |
| Client proportion | -- | 15% enrollment at debut | client_shares | Bernoulli(0.20) | STIsim: fixed | STIsim higher clients |

### Partnership Formation

| Concept | EMOD Parameter | EMOD Value | STIsim Parameter | STIsim Value | Calibrated? | Notes |
|---------|---------------|------------|-----------------|-------------|-------------|-------|
| Marital formation rate | Marital_Form_Rate | **0.00046 (0.00044-0.0005)** [calibrated] per day | p_pair_form | 0.5 per timestep | STIsim: fixed | Not directly comparable |
| Informal formation rate | Informal_Form_Rate | **0.00146 (0.00134-0.00155)** [calibrated] per day | -- | Risk-group dependent | -- | |
| Transitory formation rate | Transitory_Form_Rate | 0.001048 per day | -- | Casual: lognorm(1yr, 3yr) | -- | |
| Commercial formation rate | Commercial_Form_Rate | 0.15 per day | sw_seeking_rate | 1.0 per month | STIsim: fixed | |
| Transitory duration (scale) | Transitory_Weibull_Scale | 0.957 yr | -- | Casual: lognorm(1yr, 3yr) | -- | Similar (~1 year) |
| Coital act rate | Coital_Act_Rate | 0.33 per day (~120/yr) | acts | lognorm(mean=80/yr, SD=30/yr) | Both fixed | EMOD ~50% higher |
| Coital dilution (2 partners) | Coital_Dilution_Factor_2_Partners | 0.75 | -- | Not modeled | -- | EMOD only |
| Coital dilution (3 partners) | Coital_Dilution_Factor_3_Partners | 0.60 | -- | Not modeled | -- | |
| Coital dilution (4+ partners) | Coital_Dilution_Factor_4_Plus_Partners | 0.45 | -- | Not modeled | -- | |
| Male MR concurrency | -- | -- | m1_conc | 0.15 | **STIsim: YES** | |

### Condom Use -- WHY IS STIsim SO HIGH?

**Root cause:** The condom use *probabilities* (time series in `condom_use.csv`) are **fixed input data, not calibrated**. They were set manually and range up to 70-95% for non-stable partnerships. Only the condom *efficacy* (`eff_condom`) is calibrated.

In EMOD, condom use probabilities were calibrated as sigmoid functions and settled at much lower maxima (10-34%). The two models are using condom use as a compensating mechanism in opposite directions -- EMOD uses low condom use + low base infectivity; STIsim uses high condom use + higher base infectivity.

**Whether this matters:** The condom use CSV values (e.g., 90% for high-risk pairs) may be unrealistically high. Reported condom use in DHS/PHIA surveys for eSwatini is typically 40-60% for non-regular partners, not 90%.

| Concept | EMOD Parameter | EMOD Value | STIsim Parameter | STIsim Value | Calibrated? | Notes |
|---------|---------------|------------|-----------------|-------------|-------------|-------|
| Marital condom max | Marital_Condom_Max | **0.218 (0.207-0.231)** [calibrated] | condom_use.csv (0,0) | 0.01 (2020) | **STIsim: FIXED (CSV data)** | STIsim much lower for stable LR |
| Informal condom max | Informal_Condom_Max | **0.337 (0.321-0.355)** [calibrated] | condom_use.csv (1,1) | 0.9 (2020) | **STIsim: FIXED (CSV data)** | STIsim much higher for HR pairs |
| Transitory condom max | Transitory_Condom_Max | **0.103 (0.089-0.117)** [calibrated] | condom_use.csv cross-risk | 0.7 (2020) | **STIsim: FIXED (CSV data)** | STIsim much higher |
| Commercial condom max | Commercial_Condom_Max | 0.85 | condom_use.csv (fsw,client) | 0.95 | Both fixed | Similar, STIsim slightly higher |

---

## 4. TREATMENT CASCADE

### ART Duration & Dropout

| Concept | EMOD Parameter | EMOD Value | STIsim Parameter | STIsim Value | Calibrated? | Notes |
|---------|---------------|------------|-----------------|-------------|-------------|-------|
| ART dropout/duration | ART_dropout | 7300 days (~20 yr mean) | dur_on_art | lognorm(mean=3yr, SD=1.5yr) | -- | Base duration differs 7x |
| ART duration scaling | -- | -- | rel_dur_on_art | 1.0 (default), guess=8 | **STIsim: YES** (range 1-20) | At guess of 8x, effective mean = ~24yr, close to EMOD |

### ART Linkage & Initiation

| Concept | EMOD Parameter | EMOD Value | STIsim Parameter | STIsim Value | Calibrated? | Notes |
|---------|---------------|------------|-----------------|-------------|-------------|-------|
| Pre-ART linkage max | preART_Link_Max | **0.807 (0.783-0.829)** [calibrated] | -- | Not modeled (single step) | -- | EMOD has separate pre-ART stage |
| Pre-ART linkage midpoint | preART_Link_Mid | **1995.7 (1995.1-1996.4)** [calibrated] | -- | -- | -- | |
| Pre-ART linkage min | preART_link_Min | **0.00325 (0-0.031)** [calibrated] | -- | -- | -- | |
| ART linkage max | ART_Link_Max | **0.952 (0.948-0.955)** [calibrated] | art_initiation | 0.9 (default) | STIsim: fixed | Similar (~90-95%) |
| ART linkage midpoint | ART_Link_Mid | **2010.7 (2010.4-2010.9)** [calibrated] | -- | Immediate after diagnosis | -- | EMOD has time-varying linkage |
| ART CD4 threshold | ART_CD4_at_Initiation_Saturating_Reduction_in_Mortality | 350 | -- | Prioritized by lowest CD4 | Both fixed | Different mechanism |
| Delay to ART (post-2016) | Delay_Period_Mean | 180 days | -- | Immediate after diagnosis | -- | EMOD has explicit delay |

### ART Coverage Input

| Concept | EMOD | STIsim | Notes |
|---------|------|--------|-------|
| Format | Campaign-based (absolute numbers by year) | Proportion by age/sex/year (art_coverage.csv) | STIsim now stratified |
| Stratification | Not stratified by age/sex in input | By gender (M/F), age bin ([15,25), [25,35), [35,45), [45,100)), year | |

### HIV Testing

| Concept | EMOD | STIsim | Calibrated? | Notes |
|---------|------|--------|-------------|-------|
| Testing model | Multi-state cascade (HCT, ANC, symptomatic) | 3 parallel streams | -- | EMOD more granular |
| FSW testing (2020) | Campaign-based | 75% annual probability | STIsim: fixed | |
| General pop testing (2020) | Campaign-based | 50% annual probability | STIsim: fixed | |
| Symptomatic testing | TestingOnSymptomatic cascade | 85% annual (CD4<200) | STIsim: fixed | |
| ANC testing | TestingOnANC cascade | Not modeled | -- | EMOD only |

### PrEP

| Concept | EMOD Parameter | EMOD Value | STIsim Parameter | STIsim Value | Calibrated? | Notes |
|---------|---------------|------------|-----------------|-------------|-------------|-------|
| PrEP efficacy | Waning function | 50% at day 0, 0% at day 365 | eff_prep | 0.80 (constant) | Both fixed | EMOD wanes; STIsim constant |
| Target population | Females 15-29, HIGH risk, HIV-negative | HIV-negative FSWs | -- | Similar |

### VMMC

| Concept | EMOD Parameter | EMOD Value | STIsim Parameter | STIsim Value | Calibrated? | Notes |
|---------|---------------|------------|-----------------|-------------|-------------|-------|
| Efficacy | Circumcision_Reduced_Acquire | 0.6 | eff_circ | 0.6 | Both fixed | Match |
| Status | Active (campaign-based) | Commented out | -- | STIsim not currently using VMMC |

---

## 5. INITIAL PREVALENCE & EPIDEMIC SEEDING

### EMOD

| Parameter | Value | Calibrated? |
|-----------|-------|-------------|
| Seed_Year | **1982.7 (1982.4-1983.2)** | Yes |
| Seeding mechanism | Infections introduced into high-risk groups at Seed_Year | |
| Sim start | 1960.5 (demographics run 22 years before seeding) | |

### STIsim

Sim starts at 1985 with infections pre-seeded by risk group, sex, and sex-work status. Initial prevalence from `data/init_prev_hiv.csv`, **scaled by `rel_init_prev = 0.1`** (fixed, was previously calibrated but now commented out):

| Risk Group | Sex | Sex Work | CSV Prev | Effective Prev (x0.1) |
|------------|-----|----------|----------|----------------------|
| 0 (LR) | Female | No | 0.1% | **0.01%** |
| 1 (MR) | Female | No | 1% | **0.1%** |
| 2 (HR) | Female | No | 5% | **0.5%** |
| 0 (LR) | Female | FSW | 15% | **1.5%** |
| 1 (MR) | Female | FSW | 15% | **1.5%** |
| 2 (HR) | Female | FSW | 20% | **2.0%** |
| 0 (LR) | Male | No | 0.1% | **0.01%** |
| 1 (MR) | Male | No | 1% | **0.1%** |
| 2 (HR) | Male | No | 2% | **0.2%** |
| 0 (LR) | Male | Client | 5% | **0.5%** |
| 1 (MR) | Male | Client | 10% | **1.0%** |
| 2 (HR) | Male | Client | 15% | **1.5%** |

**Who is eligible for seeding:** Only agents who have passed sexual debut age AND are alive at sim start (1985). Given that debut is set to mean 20 (see above), this means many 15-19 year olds are NOT seeded, which may undercount early epidemic prevalence in young adults.

**Infection backdating:** Seeded infections are backdated uniformly between 5 months and 10 years, so agents start at various disease stages (some acute, some latent, some near AIDS). No one is initially diagnosed.

---

## 6. DEMOGRAPHICS

| Concept | EMOD | STIsim | Notes |
|---------|------|--------|-------|
| Population scaling | Base_Population_Scale_Factor = 0.1 | 10,000 agents (auto-scaled) | |
| Start year | 1960.5 | 1985 | EMOD starts 25 years earlier |
| Epidemic seed | ~1982.7 [calibrated] | 1985 (with init_prev CSV) | |
| End year | ~2057 | 2031 | |
| Timestep | 30.4 days (monthly) | Monthly (1/12 year) | Match |
| Migration | None | use_migration = True | STIsim includes migration |

---

## 7. KEY STRUCTURAL DIFFERENCES

### Where the models agree
- Circumcision efficacy (60%)
- Acute phase duration (~3 months)
- Acute phase per-act transmission (~0.06)
- ART linkage probability (~90-95%)
- Monthly timesteps
- Three risk groups

### Where the models fundamentally differ

| Area | EMOD | STIsim | Implication |
|------|------|--------|-------------|
| **Chronic transmission** | Very low (0.0023/act) | Higher (0.01/act) | STIsim epidemic driven more by chronic-stage; EMOD more by acute |
| **Acute multiplier** | 26x | 6x | EMOD sensitivity analyses (Fig A6) showed this matters a lot for incidence |
| **ART efficacy** | 92% reduction (0.08 multiplier) | 96% reduction + 6mo ramp-up | STIsim more optimistic; EMOD showed high sensitivity to this (Fig A7) |
| **ART dropout** | Mean 20 years | Base 3yr x rel_dur_on_art (~24yr at guess) | Similar after calibration scaling |
| **Sexual debut** | ~16-17yr (calibrated) | ~20-21yr (package default, never overridden) | **STIsim agents debut 3-4 years too late -- high priority fix** |
| **Condom use** | Low across types (max 10-34%) [calibrated] | High for non-stable (70-95%) [fixed CSV] | Opposite compensation strategies; STIsim values may be unrealistic |
| **FSW dynamics** | Dynamic entry/exit (3% enroll, Weibull dropout) | Static assignment (10%) | EMOD has turnover; STIsim does not |
| **Age-dependent sex diff.** | Young women 4.9x; older 2.8x [calibrated] | Fixed 0.5 F-to-M ratio | **No mechanism in STIsim to vary by age -- requires code change** |
| **Treatment cascade** | 20+ states (pre-ART, staging, LTFU) | Diagnose -> 90% initiate | EMOD granular cascade with loss at each step |
| **Disease progression** | Age-dependent Weibull survival | CD4-stratified mortality rates | |
| **Child HIV** | Rapid/slow progressor model | Not explicitly modeled | |
| **Calibration degrees of freedom** | 24 dynamic parameters | 6 parameters | EMOD has 4x more flexibility to fit data |

### Compensating mechanisms
The models reach similar epidemic trajectories despite different parameterizations because they use different "levers" to calibrate:
- **EMOD**: low base infectivity + high acute multiplier + low condom use + early sexual debut
- **STIsim**: higher base infectivity + lower acute multiplier + high condom use + later sexual debut

### Suggested priorities for STIsim alignment
1. **Sexual debut age** -- override defaults with EMOD-calibrated or DHS-derived values (~16-17yr F, ~17-18yr M). Easiest fix, likely largest impact.
2. **Condom use CSV** -- review and potentially lower non-stable partnership values to DHS-reported levels (40-60% rather than 70-95%).
3. **Age-dependent susceptibility** -- request feature from stisim developers or implement custom rel_sus modifier by age/sex.
4. **ART efficacy** -- consider lowering from 96% to 92% to match EMOD and empirical estimates from Donnell et al.
5. **Add more calibration parameters** -- debut age, condom use, rel_beta_f2m, init_prev scaling are all candidates.

---

## Source

**EMOD parameters**: Akullian A, Morrison M, Garnett GP, et al. The effect of 90-90-90 on HIV-1 incidence and mortality in eSwatini: a mathematical modelling study. *Lancet HIV* 2020. Table A2 (N=250 best-fitting parameter sets).

**STIsim parameters**: `hivsim_eswatini/run_sims.py`, `interventions.py`, `run_calibrations.py`, and stisim source defaults at `star_sim/stisim/stisim/diseases/hiv.py` and `networks.py`.

**STIsim source code**: `/c/Users/adamak/OneDrive - Gates Foundation/Dropbox/star_sim/stisim/`

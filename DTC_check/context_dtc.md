# Context

## Who I am

I am an experienced mathematical epidemiologist with current focus on HIV. I have deep epidemiological subject matter expertise in addition to expertise on building and implementing mathematical models of infecitous disease dynamics. I have experience coding in both R and python.

I will be providing high level instructions for code in python. Generate this code with comments and documentation

## The model

HIVSim is an agent-based model of HIV transmission patterns in Eswatini, built on STISim and Starsim. The model includes a structured sexual network with risk groups, HIV transmission, testing, ART, and PrEP. It is calibrated to UNAIDS estimates and PHIA household survey data.

## Overall goal

Right now I want to build a dashboard which I can use to check on Adam Akullian's recent calibration of the HIVSim model for Eswatini. Give me python code which does the following:

* Based on the calibrated model that adam akullian put together in the hivsim_eswatini directory, build the model
* Run the model
* Save the results
* Extract and plot results, particularly age-prevalence curves
* I plan on using jupyter notebooks to interactively run my code. Please also give me advice on how to do this in VSCode, using VSCode as a python ide

### Next steps

Some additions, for amending the output from above:

* The comparison between the model output and the calibration targets is good, but I need to see a couple of other sets of plots:
  * Need to also see comparison between how the simulation looks against data 2021. age prevalence curves for 2021 as well as against the three years that we calibrated against
  * Need to also see plots of 90-90-90 metrics: % who know their diagnosis; % on art among those who know their diagnosis; % VLS on effective art among those who know their diagnosis
    * All ages for men, women, and both
    * age-prevalence curves for men and women
    * compare these against data, where data are available
  * Need to see plots of incidence
    * all ages men and women over time
    * age-incidence curves for men and women for both

### Next iteration of calibration - 2026-10-28

Goals:
1. [ ] Prevalence among women is a little too high, especially for 20-40 year olds, for the three years of prevalence survey data.
   *Not yet addressed: this is step 4 of the work plan (recalibrating Adam's transmission parameters with a women 25–34 target), which has not been run.*
2. For antenatal care
   1. [ ] begin HIV testing among women in 1992, with midpoint in 2003 and max at 2010
      *Not implemented as written. The ANC start year and ramp were instead calibrated to SHIMS2 (step 2b): ANC testing starts in 1992.5 at 36% of pregnancies and reaches 90% by 1994.5.*
   2. [ ] begin PMTCT in 2003 with max at 2010
3. [x] I would like to calibrate the model to match the 90 90 90 numbers from the SHIMS datasets
   *Revised: only the first and second 90 are fitted, because the third 90 is a model input. The fit is to SHIMS2 at early 2017, by sex, for 15–49 and 15+, through the testing inputs only, with Adam's 7 parameters held fixed. ART is also required to reach its forced input. Result: first 90 within 1.5 points; second 90 within 3 points.*
   1. [ ] But I want Kasturi's data to be a sanity check. I don't want to see a massive spike in %diagnosed among men because the diagnosis rate among men vs women is lopsided. I believe Kasturi's data, which are fit to SHIMS data with logistic curves, should be a guide on this.
      *Partly done. Kasturi's series is plotted and reported for every fit as a sanity check, but it isn't yet used to steer the fit. Mid-2011 % diagnosed is men 48% (Kasturi 32%) and women 78% (Kasturi 36%). A first version that scored Kasturi as a target was dropped after you said Kasturi is a sanity check, not a target.*
   2. [ ] VLS appears to be fixed, but I want it to gently increase, following a logistic function, such that it fits the increase we see between 2017 and 2021 shims
      *Not implemented. VLS among those on ART is still the SHIMS2/SHIMS3 input, flat before SHIMS2 and linear between the two surveys. The only change is that SHIMS2 now sits at 2017 (see note 8).*
4. Overall improvements I would like to see
   1. [x] Prior to 2020 diagnosis rate is too high
      *Reduced, not resolved. Mid-2011 % diagnosed is now men 48% / women 78% (was 83% / 90%), still above Kasturi. SHIMS2 alone can't tell apart fits whose start dates differ by 10+ years.*
   2. [x] Better matching to the 3 90s from SHIMS_909090_coverage.csv is a priority
      *First and second 90 now match SHIMS2 (see goal 3). The third 90 is unchanged because it is an input; the model is within 3 points of SHIMS2.*
5. [x] It is acceptable that the testing rates is fixed. What I want is the general population testing rates for both men and women to be fixed such that the right levels of people living with HIV know their status, and matches the calibration target data in SHIMS_909090_coverage.csv
   *Revised: testing has separate levels for men and women, and its start year and ramp were also calibrated (step 2b), along with ANC. Current values: general-population testing starts in 2003.2, at 0.60× (men) and 0.57× (women) the original 2020 level.*
6. [x] Add PrEP distribution data from spreadsheet SWZ_prep_coverage.csv
   *The file used is SHIMS_Datasets/SWZ_PrEP_Distribution_Data.csv, copied to data/prep_initiations.csv.*
   1. [x] Oral PrEP
      1. [x] Annual distribution numbers sourced from PrEPWatch
      2. [x] Lasts 1 year, with 75% protective efficacy
      3. [x] <doi/full/10.1056/NEJMoa1108524>
         *In that trial 75% is the TDF-FTC arm; TDF alone was 67%. The data file labels the product "OralTDF".*
   2. [x] LA-PrEP
      1. [x] sourced from speaking directly with MOH 
      2. [x] Lasts 1 year, with 92% protective efficacy
   3. [x] Prep distribution should target FSW, male clients of FSW, people with multiple partners, people who are unmarried, only target people who are not HIV positive
      *The model has no marriage state, so "unmarried" means no current stable partnership. Eligibility is also limited to sexually active people aged 15+.*
   4. [x] Annual prep distribution refers to total number of prep initiations in each year. People already on PrEP should not receive it a second time. Spread out the distribution across the year. There may be more or fewer people on PrEP at any given time than the annual PrEP distribution.
      *Initiations are spread evenly over each year's monthly steps. People currently on either product are excluded, but someone whose 1-year course has ended can be initiated again in a later year.*
7. [x] Add to the list of output plots in build_and_check_calibration.py:
   1. [x] Intervention counts - the number of ART initiations, VMMC, PrEP initiations (of both types) each year. Four panels, all sharing same x axis. Points with no interpolating lines
      *ART "initiations" include restarts: the model stops and restarts people on ART to keep each age band at its target, so the panel shows ~85,000/year after 2015, far above real new initiations. The PrEP panels overlay the programme data.*

Other notes
6. [x] Antenatal testing should begin in 1995, and ramp up from 70% to 90% for mothers, with the midyear of the ramp being 2005. It may be that we should adjust this during calibration, as these numbers are left over from EMOD's calibration and might no longer apply here.
   *Implemented, then revised: the ANC start year, starting level and ramp were calibrated in step 2b (see goal 2.1).*
7. [x] 2021 SHIMS data on 90 90 90 numbers and prevalence should be held out for checking.
8. [x] SHIMS2 ran from 2016-2017, let's say all parameters from that time should be fit to early 2017 (ie, in the middle of that period)
   *Extended: the SHIMS2-based ART-coverage and VLS inputs were also moved from 2016 to 2017 (`shims2_year=2017`); otherwise the second 90 came out ~4 points too high. The prevalence targets in Adam's calibration are not yet moved; that happens in step 4.*

Work plan for re-calibrating:

1. [x] Add PrEP and the antenatal testing start date.
2. [x] With Adam's 7 calibrated parameters held fixed, tune the men's and women's testing rates until the % diagnosed matches SHIMS by sex, for ages 15–49 and 15+. This takes a few dozen runs on the laptop.
   *Revised in three passes:*
   - *Step 2 scaled the testing ramp by sex only.*
   - *Step 2b also calibrated the start years and ramps of general-population and ANC testing, scored against the SHIMS2 first and second 90, with a penalty if ART can't reach its input.*
   - *An intermediate version scored Kasturi's 2011 values as a target. It was dropped after your clarification.*
3. [x] Check that prevalence and incidence barely moved. They should, because ART coverage is forced either way.
   *15–49 prevalence moves by 0.5 points or less, except men in 2011 (−1.0). HIV deaths in 2005–2011 rise ~15%.*
4. [ ] Re-run Adam's calibration (his second round plus the women 25–34 target), with the new testing rates and PrEP fixed in place. This needs the Azure VM.
5. [ ] Check the 90s again under the new parameters. If they have drifted, repeat step 2 once.

### Next Iteration of calibration 2026-11-19

New features
1. Age-structured testing and treatment
   1. We do not have data on 90 90 90 coverage of testing or ART among infants and children under the ages of 15, so for now we will assume that those coverage rates in that age group should match that of 15-24 year olds
   2. Right now it sounds like testing and treatment rates are not structured by age. Please subdivide testing and treatment rates by age cohorts (0-24), (25-34), (35-49), (50+), to relax this 
2. ANC and PMTCT
   1. ANC testing starts in 1992.5, midpoint at 2003, max by 2010
   2. begin PMTCT in 2003 with max at 2010
3. VLS ramping, such that it fits the increase we see between 2017 and 2021 SHIMS datasets
4. minor fixes to VMMC distribution
   1. Replace the 2016 values with SHIMS2's medical + non-medical totals, consistent with how the 2021 values were built. This raises 40–44 to 0.23, 45–49 to 0.18 and 60–64 to 0.20.
   2. Move the SHIMS2 row from 2016 to 2017, consistent with your note 8, as an opt-in like shims2_year

New goals for improving the calibration
1. Prevalence among women is a little too high, especially for 20-40 year olds, for the three years of prevalence survey data.
2. Diagnosis rates are too high, start them lower. It's okay if diagnosis rates at ANC start at low enough in 1992 such that we do a better job of matching overall diagnosis rates, and to soften the spike in the diagnostic rate vs. age curve among younger people
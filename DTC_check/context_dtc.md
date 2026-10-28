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
1. Prevalence among women is a little too high, especially for 20-40 year olds, for the three years of prevalence survey data.
2. I would like to calibrate the model to match the 90 90 90 numbers from the SHIMS datasets.
3. Overall improvements I would like to see
   1. Prior to 2020 diagnosis rate is too high
   2. Better matching to the 3 90s from SHIMS_909090_coverage.csv is a priority
4. It is acceptable that the testing rates is fixed. What I want is the general population testing rates for both men and women to be fixed such that the right levels of people living with HIV know their status, and matches the calibration target data in SHIMS_909090_coverage.csv
5. Add PrEP distribution data from spreadsheet SWZ_prep_coverage.csv
   1. Oral PrEP
      1. Annual distribution numbers sourced from PrEPWatch
      2. Lasts 1 year, with 75% protective efficacy
      3. <doi/full/10.1056/NEJMoa1108524>
   2. LA-PrEP
      1. sourced from speaking directly with MOH 
      2. Lasts 1 year, with 92% protective efficacy
   3. Prep distribution should target FSW, male clients of FSW, people with multiple partners, people who are unmarried, only target people who are not HIV positive
   4. Annual prep distribution refers to total number of prep initiations in each year. People already on PrEP should not receive it a second time. Spread out the distribution across the year. There may be more or fewer people on PrEP at any given time than the annual PrEP distribution.
6. Add to the list of output plots in build_and_check_calibration.py:
   1. Intervention counts - the number of ART initiations, VMMC, PrEP initiations (of both types) each year. Four panels, all sharing same x axis. Points with no interpolating lines

Other notes
6. Antenatal testing should begin in 1995, and ramp up from 70% to 90% for mothers, with the midyear of the ramp being 2005. It may be that we should adjust this during calibration, as these numbers are left over from EMOD's calibration and might no longer apply here.
7. 2021 SHIMS data on 90 90 90 numbers and prevalence should be held out for checking.
8. SHIMS2 ran from 2016-2017, let's say all parameters from that time should be fit to early 2017 (ie, in the middle of that period)

Work plan for re-calibrating:

1. Add PrEP and the antenatal testing start date.
2. With Adam's 7 calibrated parameters held fixed, tune the men's and women's testing rates until the % diagnosed matches SHIMS by sex, for ages 15–49 and 15+. This takes a few dozen runs on the laptop.
3. Check that prevalence and incidence barely moved. They should, because ART coverage is forced either way.
4. Re-run Adam's calibration (his second round plus the women 25–34 target), with the new testing rates and PrEP fixed in place. This needs the Azure VM.
5. Check the 90s again under the new parameters. If they have drifted, repeat step 2 once.
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

### Fixes

I have added several new datasets describing ART coverage - the first is SHIMS_909090_coverage.csv - this describes all three 90 90 90 metrics - the % of PLHIV who know their status, the % of people who know their status who are on ART, and the % of people who are on ART and have achieved VLS. 
The second is aware_and_art_no_interactions_modeled.csv - this only has the first two of the 90 90 90 metrics.

What I need is to be able to compare the results from the current calibrated SWZ model against each of these datasets. 90 90 90 numbers broken down by age and sex, and also for adults 15-49 and 15+. In addition to the cascade timeseries and cascade by age checks, I want two more checks for each time series and by age which compare against these two datasets.
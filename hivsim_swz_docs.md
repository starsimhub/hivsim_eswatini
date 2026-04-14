# Project Docs

Notes on developing HIVSim country-specific model for Eswatini

## Benchmarks

External Validity: Calibration Targets

* SHIMS epidemiological data
  * HIV Prevalence by age and sex
  * Incidence
  * Viral Load
* Intervention package inputs
  * ART 95-95-95 metrics - treat these as an input for now, might switch to calibrating
    * Known status
    * ART coverage
    * VLS coverage
  * [PrEP][https://www.prepwatch.org/countries/eswatini/]
  * VMMC

## Meeting Notes

### 2026-04-14

Meeting with Adam Akullian, Anna Bershteyn, DTCitron, Nao Yamamoto to discuss calibration and validation benchmarks and standards.

* Model calibration: calibrate to first 3 SHIMS datasets, validate by comparison against the fourth SHIMS dataset
* Internal validity and other features - may revisit this in greater detail later
  * PMTCT
  * Testing and care cascade
  * Sexual network with age-assortative mixing and risk group strata 
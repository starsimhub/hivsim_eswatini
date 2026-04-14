# GitHub Issue: starsimhub/stisim

**Title:** ART stratified coverage allocates globally instead of per-stratum

**Create at:** https://github.com/starsimhub/stisim/issues/new

---

## Summary

When using stratified ART coverage (DataFrame with Year, Gender, AgeBin, and proportion columns), the `ART` intervention computes per-stratum targets correctly in `_get_n_to_treat_stratified()` but then **sums them into a single aggregate number** and allocates ART globally in `art_coverage_correction()`. This means the age/sex stratification in the input data is effectively lost during allocation.

## Expected behavior

If the ART coverage input specifies, e.g., 90% for females 35-45 and 36% for males 15-25, the model should maintain those targets within each stratum independently.

## Actual behavior

The per-stratum targets are summed into one total (e.g., 200 agents on ART), and `art_coverage_correction()` adds/removes agents globally by CD4 priority without regard to age or sex. This causes:
- Older/sicker agents to be over-treated at the expense of younger groups
- The male-female ART gap to be smoothed out
- The model to not reproduce its own stratified coverage inputs

## Root cause

In `stisim/interventions/hiv_interventions.py`:
- `_get_n_to_treat_stratified()` (line ~434) returns a single `total` int
- `art_coverage_correction()` (line ~533) operates on the total without stratum awareness

## Proposed fix

1. Have `_get_n_to_treat_stratified()` return per-stratum targets (dict keyed by `(age_bin, sex)`) in addition to the total
2. Add a `_art_coverage_correction_stratified()` method that adjusts coverage independently within each stratum
3. When stratified targets are available, use per-stratum correction instead of global correction

I have a working local fix and can submit a PR if that would be helpful.

## Reproduction

```python
import pandas as pd
import stisim as sti

art_data = pd.DataFrame({
    'Year': [2011, 2011, 2016, 2016],
    'Gender': [0, 1, 0, 1],
    'AgeBin': ['[15,25)', '[15,25)', '[15,25)', '[15,25)'],
    'p_art': [0.08, 0.18, 0.36, 0.60],
})

art = sti.ART(coverage=art_data)
# Run sim and check: males and females will converge to similar coverage
# instead of maintaining the input differential
```

"""Build the full Eswatini calibration technical deck (24-25 slides).

Run: `python presentations/build_deck.py [--theme bmgf|neutral]`
Output: presentations/eswatini_calibration_progress[_bmgf].pptx

Chrome, theme and layout primitives come from `deck_kit`; this file is only
slide content. Every figure is either an experiment's own committed figure
(referenced by its `experiments/NNN_*/figures/` path) or one of the four
cross-experiment charts from `make_figures.py`. Every number in the body text
is traceable to the SUMMARY.md of the experiment named on the slide.
"""

from deck_kit import *  # noqa: F401,F403

OUTFILE = PRES / f"eswatini_calibration_progress{SUFFIX}.pptx"

# ==========================================================================
# 1. Title
# ==========================================================================
title_slide(
    "Calibrating the Eswatini HIV model",
    "Model features, the fits they moved, and where the calibration stands",
    "Adam Akullian  ·  7 September 2026",
    "model-v1.3  ·  starsim 3.5.2  ·  stisim 1.5.11  ·  experiments 001–025")

# ==========================================================================
# 2. The question
# ==========================================================================
s, y = new("Why calibrate at all", "The calibration exists to serve one decision")
col_w = Inches(5.55)
textbox(s, MARGIN, y, col_w, Inches(3.6), [
    ("The decision question", 15, True, BLUE),
    ("Relative value of scaling up long-acting PrEP versus improving the "
     "treatment cascade to raise population viral suppression — and where "
     "prevention still matters once suppression is high.", 14, False, INK, 6),
    ("Why a point fit is not enough", 15, True, BLUE, 16),
    ("Both arms of that comparison are sensitive to per-act transmission "
     "and to who partners with whom. A single best-fitting parameter set "
     "gives one answer with no honest uncertainty on it.", 14, False, INK, 6),
    ("Aim: characterise the plausible parameter region, then propagate it "
     "through the decision — not fit one point and hope.", 14, True, INK, 8),
], spacing=1.22)
textbox(s, MARGIN + col_w + Inches(0.6), y, Inches(6.0), Inches(3.6), [
    ("Method", 15, True, ORANGE),
    ("History matching with an emulator, then trajectory selection. "
     "Chosen because it degrades gracefully when some targets are "
     "unreachable — which several of ours are.", 14, False, INK, 6),
    ("Compute", 15, True, ORANGE, 16),
    ("~130 s per simulation at N = 10,000 on a laptop; 1000 design points "
     "in 456 s on a 120-core Azure VM. The VM is what makes waves and "
     "synthetic-data checks affordable.", 14, False, INK, 6),
    ("Discipline", 15, True, ORANGE, 16),
    ("One folder per experiment. README (question + plan) before any code "
     "runs; SUMMARY (findings + decision) before the next folder opens. "
     "24 experiment folders so far, numbered 001–025 — 010 was deferred by "
     "009 and absorbed into 012.", 14, False, INK, 6),
], spacing=1.22)
takeaway(s, "Everything that follows is in service of one comparison: "
            "prevention scale-up versus cascade improvement.")

# ==========================================================================
# 3. Targets
# ==========================================================================
s, y = new("Experiment 008, amended through 025", "What we fit, and what we deliberately do not")
table(s, [
    ["Quantity", "Rows", "Source", "Status"],
    ["HIV prevalence by age × sex × year", "54", "PHIA 2007, 2011, 2016",
     "Fitted — tier A/B/C"],
    ["Annual AIDS deaths", "35", "UNAIDS Spectrum, 1990–2024",
     "Fitted but heavily down-weighted (σ = 2000)"],
    ["HIV incidence, adult aggregates", "4",
     "SHIMS1 2011 cohort; SHIMS2 2016 Table 5.3.B", "Fitted, added 2026-09-04"],
    ["HIV incidence by age band", "6", "SHIMS2 2016",
     "NOT fitted — false-recency bias"],
    ["2021 incidence", "—", "SHIMS3", "Held out for validation"],
    ["ART coverage, VMMC coverage, VLS", "—", "PHIA / SHIMS",
     "Model inputs, not targets"],
], y, [Inches(3.5), Inches(0.65), Inches(3.3), Inches(4.65)], size=11)
footnote(s, "Why the age-banded incidence rows are excluded: the published "
            "female profile is essentially flat (1.67, 1.54, 2.09 across "
            "15–24, 25–34, 35–49), which is the signature of false-recency "
            "bias — a recency assay credits a fraction of long-standing "
            "infections as recent, and against a susceptible denominator the "
            "spurious rate scales as FRR × P/(1−P)/MDRI: ~5× larger in women "
            "35–49 than 15–24, ~18× in men. Consequence, stated plainly: the "
            "model's incidence age profile is unconstrained by data.",
         y + Inches(2.85))
takeaway(s, "Deaths are down-weighted and the incidence age profile is "
            "unconstrained — both are deliberate, and both are consequences "
            "of findings later in this deck.", ORANGE)

# ==========================================================================
# 4. Metric glossary — fit and A/B
# ==========================================================================
# Definitions read from source, not from memory: MAE / bias / n_within_ci from
# `standard_figures.scorecard`; coverage from 014's `assess_coverage`; CV and
# the expected-agent floor from 020's `cv_table` / `expected_counts`; the
# two-sample z from 017's reproduction check.
s, y = new("Reference", "How to read the fit numbers in this deck")
table(s, [
    ["Metric", "Definition", "How to read it"],
    ["MAE vs PHIA",
     "Mean of |model − PHIA| over the 54 age × sex × year prevalence strata. "
     "'Model' is the mean over seeds at a fixed parameter point, or the median "
     "over draws for a prior ensemble.",
     "Lower is better. The scale is prevalence, so 0.058 means the average "
     "stratum is 5.8 percentage points out."],
    ["Bias vs PHIA",
     "Mean of (model − PHIA), signed, over the same 54 strata.",
     "Negative = model below the surveys. Separates a systematic level error "
     "from scatter, which MAE alone cannot."],
    ["Strata inside the PHIA 95% CI",
     "Count of strata where the model's central estimate falls between the "
     "survey's published lower and upper bound.",
     "Out of 54. Blind to how badly the misses miss, so read it next to MAE, "
     "not instead of it."],
    ["Coverage (prior predictive)",
     "Share of target rows whose observed value falls inside the 5–95 "
     "percentile envelope of the simulated ensemble.",
     "A check on the prior, not on the fit. A miss means no draw in the "
     "ensemble got near the data."],
    ["Percentile of the observation",
     "Where the observed value sits within the ensemble's distribution for "
     "that one target.",
     "100.0 means the entire ensemble is below the observation. This is what "
     "names a target unreachable."],
    ["Between-seed CV",
     "SD across seeds of each seed's trajectory-mean prevalence 15–49, "
     "divided by the mean of those seed means.",
     "Sets the replicate count: < 5% → 3–5 reps, < 20% → 10–20, otherwise 50+ "
     "or a larger N."],
    ["Expected infected agents per stratum",
     "Modelled infected count in a PHIA stratum, converted back from the "
     "population scale factor to agents, averaged over seeds.",
     "Floor of 5. Below it a stratum fails on Poisson noise rather than on "
     "the parameters."],
    ["Two-sample z (A/B arms)",
     "(mean_A − mean_B) / √(se_A² + se_B²), with se = between-seed SD ÷ "
     "√n_seeds.",
     "|z| < 2 at n = 10 is a null, not proof of no effect. Stops a seed "
     "lottery being read as a level shift."],
    ["% of the UNAIDS peak",
     "Model's peak annual AIDS deaths ÷ UNAIDS' 11,000 at the 2004 peak, "
     "compared peak-to-peak.",
     "The model peaks ~2003 and UNAIDS ~2004–05, so peak-to-peak is the "
     "fairer comparison than same-year."],
], y, [Inches(2.55), Inches(5.25), Inches(4.29)], size=9.5,
    row_h=Inches(0.455))
takeaway(s, "Two habits worth naming: MAE and bias answer different questions, "
            "and aggregating to 15–49 flatters the model because it averages "
            "young under-infection against older over-infection.")

# ==========================================================================
# 5. Divider — features
# ==========================================================================
divider("Part one", "The features we added",
        "Four changes to the model itself. Three were bugs in inherited code; "
        "one was a set of defaults nobody had chosen. Each is measured against "
        "an A/B, not adopted on argument.")

# ==========================================================================
# 5. Network age mixing
# ==========================================================================
s, y = new("Experiments 011 and 013", "The sexual network was ignoring its own age-gap parameter")
bar = takeaway(s, "A parameter the model appeared to expose was inert. No "
                  "prior over it was interpretable, so calibrating the network "
                  "was meaningless until this was fixed.")
picture(s, EXP / "011_network_age_mixing/figures/algorithm_comparison_summary.png",
        y, bar - Inches(1.32), left=MARGIN, width=Inches(6.35))
textbox(s, MARGIN + Inches(6.7), y, Inches(5.35), Inches(3.9), [
    ("`match_pairs` paired partners by rank within the surviving pool, which "
     "flattens the configured gap distribution into a function of the "
     "marginal age supplies.", 13, False, INK),
    ("Sweeping the configured mean gap 0 → 7 → 14 y moved the realised mean "
     "only 0.7 → 1.4 → 5.3 y. Women 35–49 got a biologically impossible "
     "−2.7 y gap.", 13, False, INK, 9),
    ("Upstream PR #477 restores assignment-based matching: realised tracks "
     "configured nearly 1:1 (0.9 → 6.2 → 12.3) and lands within 1–2 y of DHS "
     "Eswatini across every female age band. It now ships as the stisim "
     "default matcher.", 13, False, INK, 9),
    ("013's follow-up: the corrected matcher raises prevalence "
     "(4.09% → 6.46% at 2021), which contradicted our own earlier "
     "attribution of a version-bump prevalence drop to the network fix. "
     "That sent us looking elsewhere — see the next slide.", 13, False, INK2, 9),
], spacing=1.2)
footnote(s, "Our own local Gaussian patch flattened almost identically to the "
            "broken baseline and was discarded. `max_deviation` (1/3/5) made "
            "no material difference.", bar - Inches(0.44))

# ==========================================================================
# 6. VMMC
# ==========================================================================
s, y = new("Experiment 015", "VMMC coverage was a hazard where the data are a stock")
bar = takeaway(s, "Over-circumcision was suppressing male acquisition hard. "
                  "Fixing it roughly doubled male prevalence — and it was the "
                  "real driver of a prevalence drop we had blamed on the "
                  "network.", AQUA)
picture(s, EXP / "015_vmmc_prevalence_target/figures/vmmc_coverage_by_age.png",
        y + Inches(0.05), y + Inches(3.35), left=MARGIN, width=Inches(6.35))
table(s, [
    ["2021, 10-seed mean", "Upstream", "Fixed", "SHIMS3"],
    ["Circumcision 15–49", "99.3%", "45.4%", "≈ 48.3%"],
    ["HIV prevalence 15–49", "6.69%", "11.94%", "—"],
    ["  male", "2.60%", "6.36%", "—"],
    ["  female", "10.72%", "17.43%", "—"],
], y + Inches(0.05), [Inches(2.35), Inches(1.15), Inches(0.95), Inches(0.95)],
    left=MARGIN + Inches(6.7), size=10.5, row_h=Inches(0.30))
textbox(s, MARGIN + Inches(6.7), y + Inches(1.85), Inches(5.4), Inches(2.4), [
    ("Upstream applied p × n_uncircumcised each step — a hazard on the "
     "remaining pool — so coverage ratcheted to ~100% regardless of the "
     "target, and ignored the age gradient entirely.", 12.5, False, INK),
    ("The in-repo subclass tops up to p × (all alive males in the age bin) "
     "per stratum and never removes: a stock target matching what a "
     "cross-sectional survey measures.", 12.5, False, INK, 8),
    ("Retired at exp 017 — stisim 1.5.9 adopted the same semantics upstream, "
     "and the two arms now agree to three decimals.", 12.5, False, INK2, 8),
], spacing=1.18)
textbox(s, MARGIN, y + Inches(3.45), Inches(6.35), Inches(0.6),
        [("Denominator is all men, not HIV-negative men: SHIMS and PHIA "
          "report circumcision among all male respondents regardless of HIV "
          "status, so the subclass tops up over all alive males. Filtering to "
          "HIV-negatives would over-circumcise.", 9.5, False, INK2)],
        spacing=1.25)

# ==========================================================================
# 7. Mortality double-counting
# ==========================================================================
s, y = new("Experiment 016", "HIV mortality was being counted twice")
bar = takeaway(s, "Reconstructing non-AIDS mortality from all-cause data "
                  "alone — using no HIV information as input — reproduces the "
                  "UNAIDS AIDS-death curve to within ~10%. The background "
                  "rates were carrying a whole epidemic.", AQUA)
picture(s, EXP / "016_double_counted_mortality/figures/implied_vs_unaids.png",
        y, bar - Inches(0.5), left=MARGIN, width=Inches(6.1))
textbox(s, MARGIN + Inches(6.45), y, Inches(5.6), Inches(4.0), [
    ("The background death rates fed to `ss.Deaths` kill agents regardless "
     "of HIV status — and their adult rates rise 5.5× to a 2005 peak. That "
     "is the shape of the AIDS epidemic, so the HIV module was adding "
     "mortality on top of mortality already there.", 13, False, INK),
    ("Implied vs UNAIDS AIDS deaths: 11,113 vs 11,000 at the 2005 peak "
     "(+1%); mean ratio 1.09 across 1995–2015.", 13, True, INK, 9),
    ("The population trajectory is the independent check. 2015 population "
     "against a target of 1,313,671:", 13, False, INK, 9),
    ("all-cause  1,131,932  (−13.8%)", 13, False, RED, 6),
    ("HIV-deleted  1,315,148  (+0.1%)", 13, True, AQUA, 2),
    ("But the effect on prevalence is small — about +1.5 points where ~6 "
     "were needed. This corrected a real error; it did not fix the fit.",
     13, False, INK2, 9),
], spacing=1.2)

# ==========================================================================
# 8. Inherited defaults
# ==========================================================================
s, y = new("Experiments 018, 021, 022", "Four inherited defaults that nobody had chosen")
table(s, [
    ["Change", "What was wrong", "Effect on the fit", "Why adopt anyway"],
    ["Remove `sti.Prep()` default (018)",
     "An undeclared ramp to 80% of FSW starting 2004 — a decade before "
     "efficacy evidence — inherited by every experiment from 001 to 017",
     "2021 prevalence +4.1% (z = 1.68)",
     "PrEP must be specified from programme data in the decision analysis, "
     "not inherited"],
    ["Set `vls_coverage` from PHIA (021)",
     "Defaulted to 1.0: every ART initiator virally suppressed",
     "Null on prevalence (+1.6%, z = 1.03); population VLS error 8.8 pp → "
     "2.1 pp",
     "The decision question is about raising suppression; the baseline "
     "overstated it, understating the cascade arm's headroom"],
    ["`VLSStockTarget` (022)",
     "Suppression could only improve for new initiators, so a regimen "
     "change like TLD was not representable",
     "Unchanged (MAE 0.0584 → 0.0586)",
     "Flow-only lags its own input by 1.8 pp in 2021 — exactly when a "
     "cascade scenario matters"],
    ["ART coverage input fix (022)",
     "SHIMS2's 20–24 row applied to the whole 15–25 band; 2016 men 15–24 "
     "went 0.360 → 0.556",
     "Moves MAE by 0.0002",
     "Provenance. Young HIV-positive men are few, so a 20-point error "
     "barely propagates"],
], y, [Inches(2.55), Inches(3.5), Inches(2.5), Inches(3.55)], size=9.5,
    row_h=Inches(0.92))
takeaway(s, "None of these moved the fit. All four were adopted because a "
            "measured input beats an inherited default — and three of them "
            "matter directly to the decision the calibration is for.", ORANGE)

# ==========================================================================
# 9. Fit progression
# ==========================================================================
s, y = new("Cross-experiment view", "How the features shifted the fit")
bar = takeaway(s, "The features fixed mechanisms, not the fit level. Between "
                  "model-v1.1 and v1.3 the age-stratified MAE moved 0.0590 → "
                  "0.0584 — nothing. The level came from opening transmission "
                  "in the calibration, not from the model changes.", ORANGE)
picture(s, FIG / "fit_progression.png", y, bar - Inches(0.15))

# ==========================================================================
# 10. Divider — ruled out
# ==========================================================================
divider("Part two", "What we ruled out",
        "Six hypotheses, tested and closed. The negative results cost real "
        "time and they are the strongest part of the record: each one converts "
        "a vague hypothesis into a settled question.")

# ==========================================================================
# 11. Ruled out
# ==========================================================================
s, y = new("Experiments 014, 016–022", "The AIDS-death deficit: six hypotheses, one partial answer")
table(s, [
    ["Hypothesis", "Test", "Result"],
    ["HIV mortality parameters are too slow",
     "014 — opened a CD4 death-rate multiplier over a 6× range",
     "ρ = −0.01 on peak deaths. Structurally inert: deaths fire from a "
     "separate `ti_zero` pathway."],
    ["Nobody on ART can die of HIV",
     "017 — stisim 1.5.11 makes on-ART mortality nonzero by default",
     "64% of the UNAIDS peak, identical to before. The on-ART hazard is "
     "anchored to the off-ART hazard at restored CD4, where it is tiny."],
    ["The new female mortality multiplier lengthens survival",
     "017 arm D — set `rel_death_f` back to 1.0, hold everything else",
     "|z| ≤ 0.9 at every year. Refuted; there was no effect to explain."],
    ["The whole version bump changed behaviour",
     "017 — four arms, each adjacent pair differing by one thing",
     "No A→B comparison reaches |z| > 1.7 from 1995 to 2021. Safe to adopt "
     "without recalibrating."],
    ["Untreated survival is too long (13.1 y, flat in age)",
     "019 — shorten uniformly, then add age gradients; 022 — EMOD's "
     "pivoting gradient",
     "Closes 39% of the deficit and costs prevalence to do it. The gradient "
     "saturates early: mild ≈ full ALPHA."],
    ["The model is too small for the rare strata",
     "020 — N ∈ {5k, 10k, 20k, 50k} at both ends of the prior",
     "N = 50,000 still leaves 2 of 54 strata below 5 expected cases. The "
     "binding stratum is a fit error, not a sizing error."],
], y, [Inches(3.15), Inches(3.5), Inches(5.45)], size=9.5, row_h=Inches(0.72))
takeaway(s, "About 22% of the death deficit survives every mechanism we "
            "tested. That is why deaths are down-weighted rather than fitted "
            "hard.", RED)

# ==========================================================================
# 12. The trade-off
# ==========================================================================
s, y = new("Experiments 019 and 022", "Deaths and prevalence are in structural tension")
bar = takeaway(s, "Three independent routes to shorter untreated survival "
                  "land on the same frontier. A hard deaths target would "
                  "force the calibration into a trade the data cannot "
                  "resolve — so deaths carry σ = 2000 and prevalence carries "
                  "the wave.", RED)
picture(s, FIG / "tradeoff.png", y, bar - Inches(0.15), left=MARGIN,
        width=Inches(7.5))
textbox(s, MARGIN + Inches(7.85), y + Inches(0.1), Inches(4.25), Inches(3.9), [
    ("+14 points of the UNAIDS peak costs −3 points of prevalence bias, on "
     "a fit already biased low.", 13, True, INK),
    ("Measured three ways", 13, True, RED, 12),
    ("019 shortened survival uniformly (+9.5 pts of the peak at zero "
     "gradient), then added age gradients (+4.3 pts more, saturating).",
     12.5, False, INK, 6),
    ("022 came at it from the opposite direction with EMOD's pivoting "
     "gradient — longer at young ages, shorter at old — and landed on the "
     "same frontier at 73.3%.", 12.5, False, INK, 6),
    ("Consequences recorded", 13, True, RED, 12),
    ("`mort_mult` fixed at 1.0. `dur_latent_mult` NOT opened as a "
     "calibration parameter: its leverage is generic to any late-life "
     "shortening and always costs prevalence, so it would spend effort "
     "trading two targets against each other.", 12.5, False, INK, 6),
], spacing=1.18)

# ==========================================================================
# 13. Divider — process
# ==========================================================================
divider("Part three", "The calibration itself",
        "Prior predictive check, model sizing, parameter engineering, then "
        "the first history-matching wave. In that order, and each step "
        "gated the next.")

# ==========================================================================
# 14. Metric glossary — history matching
# ==========================================================================
# sigma composition read from 024/run.py `build_observations`; the emulator
# metrics from the package's own metrics.json; the variance-reduction formula
# from the package's constrained_dims diagnostic; rho and effect-signature r
# from 023/run.py `effect_signature_confounding`.
s, y = new("Reference", "How to read the history-matching numbers")
table(s, [
    ["Quantity", "Definition", "How to read it"],
    ["σ (observation uncertainty)",
     "Per target: √(σ_CI² + σ_disc²). σ_CI = (upper − lower) / 3.92, a "
     "published 95% CI converted to an SD; for the 15–49 aggregates the "
     "strata are weighted by PHIA's own denominators.",
     "History matching has no weights, only σ. Widening σ is the only way to "
     "down-weight a target."],
    ["σ_disc (model discrepancy)",
     "An explicit allowance for the model being wrong, added in quadrature. "
     "0.02 on prevalence in wave 1, cut to 0.01 for wave 2.",
     "A modelling judgement, not a measurement. In wave 1 it was ~90% of the "
     "total variance — the main reason the NROY stayed wide."],
    ["Deaths σ = 2000",
     "The deliberate down-weighting. The model reaches 7,069 against UNAIDS' "
     "11,000, so σ = 2000 puts that gap at ~2σ.",
     "A survey-like σ would put the same gap at 8σ and empty the box — for a "
     "structural reason we already understand."],
    ["Implausibility I(x)",
     "|emulator mean − target| ÷ √(emulator variance + σ²), maximised over "
     "the emulated features.",
     "One number per parameter point. It charges emulator uncertainty as well "
     "as observation uncertainty."],
    ["Threshold = 4",
     "A parameter point is ruled out when I(x) > 4.",
     "Deliberately loose. It should discard what is clearly implausible, not "
     "select a best fit."],
    ["NROY fraction",
     "Share of the prior box with I(x) ≤ 4 — 'Not Ruled Out Yet'.",
     "Volume, not informativeness. 87% can still be a real cut if the cut is "
     "diagonal — see slide 21."],
    ["Within 3σ (residual check)",
     "|model − target| ÷ σ ≤ 3, evaluated on an actual simulation rather than "
     "the emulator, one target at a time.",
     "Separate from implausibility: no emulator variance. This is what "
     "'48 of 48' on the best draw means."],
    ["Emulator R²",
     "Variance in the emulated feature the emulator explains. The 1000-point "
     "design is split 750 train / 250 held out, and the reported R² is the "
     "held-out figure.",
     "> 0.8 to trust the cut. > 0.99 would suggest overfitting rather than a "
     "good emulator."],
    ["Variance reduction per PC",
     "1 − var(NROY) / var(prior), projected onto each principal component of "
     "the design.",
     "Says which *directions* the data constrained. Per-parameter marginals "
     "miss a diagonal cut entirely."],
    ["Spearman ρ",
     "Rank correlation between a parameter's drawn value and a target summary "
     "statistic, over the draws that established an epidemic.",
     "Does this parameter move this target? This is what pruned nine "
     "parameters to seven."],
    ["Effect-signature r",
     "Correlation between two parameters' ρ-vectors across all 23 target "
     "statistics.",
     "Do two parameters do the *same thing*? High r means the data cannot "
     "separate them, however they were sampled."],
], y, [Inches(2.35), Inches(5.35), Inches(4.39)], size=9,
    row_h=Inches(0.355))
takeaway(s, "The two that get misread most: NROY fraction is a volume, not a measure of how much was learned; and σ_disc is a judgement call that can quietly dominate every other uncertainty in the problem.", AQUA)

# ==========================================================================
# 15. Coverage
# ==========================================================================
s, y = new("Experiments 009, 014, 024", "Step 1 — can the model even produce the data?")
bar = takeaway(s, "Widening a prior lowers coverage whenever the added volume "
                  "is in directions that do not matter. 014 spread 50 draws "
                  "over 9 dimensions and sampled below the establishment "
                  "threshold — so the corner where the data lives held about "
                  "one draw.", AQUA)
picture(s, FIG / "coverage.png", y, bar - Inches(0.15), left=MARGIN,
        width=Inches(7.4))
textbox(s, MARGIN + Inches(7.75), y + Inches(0.1), Inches(4.35), Inches(3.9), [
    ("009 — failed, and diagnosed wrongly", 13, True, RED),
    ("34% coverage. Diagnosis: 'too little mortality flow.' Named three "
     "fixed parameters as suspects.", 12.5, False, INK, 5),
    ("014 — failed worse, and inverted the diagnosis", 13, True, RED, 11),
    ("4% coverage, with 85 of 89 rows having the entire ensemble below the "
     "observation and not one above it. But the best single draw hit 0.502 "
     "against an observed 0.507 — the model could reach the data; the prior "
     "was in the wrong place. Mortality was refuted.", 12.5, False, INK, 5),
    ("024 — passed", 13, True, AQUA, 11),
    ("94% of registered features covered, and one draw inside 3σ on all 48 "
     "at once. Joint reachability, not target-by-target.", 12.5, False, INK, 5),
], spacing=1.16)

# ==========================================================================
# 15. Sizing
# ==========================================================================
s, y = new("Experiment 020", "Step 2 — how big must the model be to score coverage?")
bar = takeaway(s, "Past N = 20,000 more agents mostly buy precision on a "
                  "bias. Coverage check v3 runs at N = 20,000 with 10 "
                  "replicates, and declares the two unresolvable strata "
                  "rather than scoring them.", AQUA)
picture(s, EXP / "020_model_sizing/figures/sizing.png", y + Inches(0.1),
        y + Inches(2.5))
table(s, [
    ["Expected infected agents in the thinnest PHIA stratum", "N = 5,000",
     "N = 10,000", "N = 20,000", "N = 50,000"],
    ["at high transmission", "1.27", "2.42", "5.87", "14.09"],
    ["at low transmission (β at the establishment floor)", "0.63", "0.40",
     "1.46", "3.33"],
    ["strata below 5 expected cases, low transmission", "10 of 54", "5 of 54",
     "3 of 54", "2 of 54"],
    ["between-seed CV on trajectory-mean prevalence", "9.4% / 24.4%",
     "4.4% / 12.9%", "2.5% / 7.8%", "1.9% / 4.3%"],
], y + Inches(2.65), [Inches(4.85), Inches(1.5), Inches(1.55), Inches(1.55),
                      Inches(1.55)], size=10, row_h=Inches(0.30))
footnote(s, "The stratum that clears at no feasible N is 2007 men 15–19, "
            "where the model puts prevalence 0.003 against PHIA's 0.019 — a "
            "−65% relative miss. Adding the low-transmission arm is what "
            "turned this from a confirmation into a finding: the earlier "
            "estimate was computed at high transmission only, and with the "
            "PHIA sex mapping inverted.", y + Inches(4.35))

# ==========================================================================
# 16. Parameter engineering
# ==========================================================================
s, y = new("Experiment 023", "Step 3 — which parameters to open, on measured evidence")
bar = takeaway(s, "Nine parameters to seven, both cuts on measurement rather "
                  "than judgement — and one of them reverses an earlier "
                  "experiment's own ranking.", AQUA)
picture(s, EXP / "023_parameter_engineering/figures/sensitivity.png",
        y + Inches(0.02), y + Inches(1.80))
table(s, [
    ["Parameter opened for wave 1", "Prior", "Scale", "Strongest target", "ρ"],
    ["beta_m2f", "0.008 – 0.025", "log", "peak AIDS deaths", "+0.74"],
    ["rel_beta_f2m", "0.15 – 0.60", "log", "prevalence F:M ratio", "−0.68"],
    ["prop_m0", "0.40 – 0.80", "linear", "male 45–64 prevalence", "−0.54"],
    ["prop_f0", "0.45 – 0.85", "linear", "female 45–64 prevalence", "−0.50"],
    ["s_f_young", "0.8 – 3.0", "log", "female young:old prevalence ratio",
     "+0.41"],
    ["age_gap_shift", "−2 to +3 yr", "linear", "2011 F:M incidence ratio",
     "+0.36"],
    ["age_gap_sd_mult", "0.6 – 1.8", "log", "female prevalence", "0.27"],
], y + Inches(1.90), [Inches(3.1), Inches(2.0), Inches(1.05), Inches(4.9),
                      Inches(1.04)], size=10.5, row_h=Inches(0.25))
footnote(s, "Dropped, both on measured evidence rather than judgement.  "
            "rel_init_prev → fixed at 0.2: max |ρ| = 0.10, the weakest in the "
            "set — 014 had ranked it second of nine at ρ = 0.59, an artefact "
            "of a prior that went below the establishment threshold, where it "
            "acted as an on/off switch for whether an epidemic happened at "
            "all.  conc_mult → fixed at 1.0: weakest of the mixing set at "
            "0.19, and confounded with rel_beta_f2m at effect-signature "
            "r = 0.81.", y + Inches(4.02))

# ==========================================================================
# 17. Wave 1
# ==========================================================================
s, y = new("Experiment 024", "Step 4 — history matching, wave 1")
bar = takeaway(s, "The coverage question that blocked 009 and 014 is answered "
                  "affirmatively: the 5–95% envelope covers 45 of 48 targets, "
                  "and one draw of 1000 sits within 3σ on all 48 at once.",
               AQUA)
picture(s, EXP / "024_hm_wave1/figures/prevalence_fit_vs_phia.png", y,
        bar - Inches(0.15), left=MARGIN, width=Inches(7.4))
textbox(s, MARGIN + Inches(7.75), y + Inches(0.05), Inches(4.35), Inches(4.0), [
    ("1000 design points, 7 parameters, N = 10,000, 456 s on 120 cores.",
     12.5, False, INK2),
    ("Emulator", 13, True, AQUA, 10),
    ("R² = 0.923, MSE = 3.6e−4 on 750 training points. Past the 0.8 bar and "
     "not so high as to suggest overfitting. Its trend coefficients rank the "
     "parameters in the same order 023's Spearman analysis did, by a "
     "completely different method.", 12.5, False, INK, 5),
    ("The best joint draw (row 868)", 13, True, AQUA, 10),
    ("48/48 targets within 3σ. Age-stratified MAE 0.0405 with bias −0.0070, "
     "against 0.0584 and −0.0396 at the pre-calibration configuration. "
     "Strata inside the PHIA CI: 26 of 54, against 23.", 12.5, False, INK, 5),
    ("It is one draw of 1000, not a posterior. It is worth recording because "
     "covering each target separately does not imply any single point covers "
     "them all.", 12, False, INK2, 8),
], spacing=1.15)

# ==========================================================================
# 18. The reversal
# ==========================================================================
s, y = new("Experiment 024, the headline", "The 'irreducible' prevalence deficit was a parameter-value deficit")
bar = takeaway(s, "Seven experiments held transmission at its defaults and "
                  "diagnosed a structural bias. Opening `beta_m2f` largely "
                  "removes it. The lesson is about the diagnosis, not the "
                  "model.", ORANGE)
picture(s, EXP / "024_hm_wave1/figures/target_residuals.png", y,
        bar - Inches(0.15), left=MARGIN, width=Inches(6.5))
textbox(s, MARGIN + Inches(6.85), y + Inches(0.05), Inches(5.25), Inches(4.0), [
    ("Experiments 016 through 022 all measured a −4.3 pp prevalence deficit "
     "and could not move it. 024 budgeted σ_disc = 0.02 to absorb it as "
     "model discrepancy.", 13, False, INK),
    ("With σ_disc set to zero — PHIA sampling error alone, no discrepancy "
     "allowance whatever — the best draw sits at 3.98σ on the six PHIA 15–49 "
     "targets.", 13, True, INK, 9),
    ("The discrepancy allowance was calibrated against a bias that was not "
     "structural, and it dominated the observation variance: 0.02² against "
     "PHIA's 0.0043–0.0066², so ~90% of it. That is the main reason the NROY "
     "came back at 87% rather than tighter.", 13, False, INK, 9),
    ("Wave 2 cuts it to 0.01 — worth ~5× the cutting power, while keeping a "
     "real allowance for the two genuine defects. Zero is not defensible.",
     13, False, ORANGE, 9),
], spacing=1.18)

# ==========================================================================
# 19. Reading the NROY
# ==========================================================================
s, y = new("Experiment 024, diagnostics", "87% NROY is the right answer, not a failed one")
bar = takeaway(s, "With 7 parameters and one emulated scalar, most of the box "
                  "must survive — one observation constrains roughly one "
                  "direction. Reporting per-parameter marginals would have "
                  "made this wave look like it did nothing.", AQUA)
picture(s, FIG / "constrained_dims.png", y, bar - Inches(1.05))
footnote(s, "This confirms 023's identifiability warning from the other side. "
            "023 found `beta_m2f` and `s_f_young` confounded at "
            "effect-signature r = 0.82 and predicted they would be correlated "
            "in the posterior; PC1 loads both positively together, so they "
            "enter the emulated feature as a product. Wave 1 constrains that "
            "product and leaves the ratio free — the predicted failure, now "
            "measured rather than anticipated. It is also why the "
            "re-identification check is the gate on wave 2.",
         bar - Inches(0.98))

# ==========================================================================
# 20. What is still broken
# ==========================================================================
s, y = new("Experiment 024, observations 7 and 8", "Two specific defects survive — and they point the same way")
bar = takeaway(s, "Both defects are about whom older women and young men "
                  "partner with. That is a model-structure question, not a "
                  "calibration one.", RED)
picture(s, EXP / "024_hm_wave1/figures/incidence_fit_vs_shims.png",
        y + Inches(0.02), y + Inches(2.28))
table(s, [
    ["Target the prior envelope misses", "Observed", "Ensemble p95",
     "Percentile of the observation in the ensemble"],
    ["female incidence 35–49, 2016", "2.09", "1.087", "100.0"],
    ["male prevalence 15–25, 2007", "0.0589", "0.0480", "97.7"],
    ["male prevalence 15–25, 2011", "0.0506", "0.0504", "95.0"],
], y + Inches(2.40), [Inches(4.2), Inches(1.4), Inches(1.6), Inches(4.89)],
    size=10.5, row_h=Inches(0.30))
footnote(s, "Incidence: the level is right, the age profile is not. All four "
            "published 15–49 aggregates land on the ensemble median (2011 "
            "women 3.1 vs 3.14, men 1.6 vs 1.65; 2016 women 1.7 vs 1.73, men "
            "0.85 vs 0.85) — but within 2016 the model puts 2.05 at 15–24 "
            "against a measured 1.67, and 0.35 at 35–49 against 2.09. So the "
            "mechanism has to move female risk across age, not add "
            "transmission, which puts `s_f_young` under suspicion.\n"
            "This supersedes the 'women 15–24 / men 25–34 / women 35–44' "
            "characterisation carried from 016–022, measured at default "
            "transmission parameters. A pre-registered prediction — that a "
            "surviving male 25–34 trough would indicate a structural gap — "
            "was aimed at the wrong band: `prev_m_25_35` fits at −2.64σ or "
            "better.", y + Inches(3.70))

# ==========================================================================
# 21. Where we are
# ==========================================================================
s, y = new("Experiment 025", "Where this stands today")
col = Inches(3.85)
gap = Inches(0.35)
textbox(s, MARGIN, y, col, Inches(3.8), [
    ("Done", 15, True, AQUA),
    ("Target set settled and documented, including what is excluded and why.",
     12.5, False, INK, 7),
    ("Four model features adopted, each against a measured A/B; three "
     "in-repo fixes, one of which upstream has since absorbed.",
     12.5, False, INK, 5),
    ("Six hypotheses for the AIDS-death deficit closed.", 12.5, False, INK, 5),
    ("Prior predictive check passing at 94%, with joint reachability "
     "demonstrated.", 12.5, False, INK, 5),
    ("Model sizing and the parameter set both settled on measurement.",
     12.5, False, INK, 5),
    ("Wave 1 complete: emulator R² = 0.923, NROY 87.2%, constrained "
     "directions diagnosed.", 12.5, False, INK, 5),
], spacing=1.2)
textbox(s, MARGIN + col + gap, y, col, Inches(3.8), [
    ("In flight", 15, True, ORANGE),
    ("025 is written, committed and smoke-tested. It has two parts:",
     12.5, False, INK, 7),
    ("Part 1 — re-identification gate. Recover experiment 024's best point "
     "from synthetic data generated at a different seed. Wave 2 does not run "
     "until this reports. The failure that matters is recovering the product "
     "`beta_m2f × s_f_young` but not the ratio.", 12.5, False, INK, 5),
    ("Part 2 — wave 2. σ_disc halved to 0.01, three features emulated "
     "simultaneously (adult prevalence, incidence, the female young:old "
     "ratio), N = 20,000.", 12.5, False, INK, 5),
    ("Blocked on compute. No HB120 spot capacity in the region; the "
     "alternative VM has no account for us yet. ~4 h per wave on the laptop, "
     "and closing the lid stops the run.", 12.5, True, ORANGE, 7),
], spacing=1.2)
textbox(s, MARGIN + 2 * (col + gap), y, col, Inches(3.8), [
    ("Next", 15, True, BLUE),
    ("Wave 3 on the F:M ratios, to settle whether `rel_beta_f2m` sits near "
     "0.25 (023, from the incidence ratio alone) or near 0.40 (024, from the "
     "joint fit across 48 targets). A per-act transmission parameter that a "
     "PrEP-versus-treatment comparison is directly sensitive to.",
     12.5, False, INK, 7),
    ("Then trajectory selection, and the Bayesian step the waves are a "
     "precursor to.", 12.5, False, INK, 5),
    ("Open model-structure question: does older women's incidence need a "
     "mechanism? Candidates are age mixing and the female partnering taper "
     "to zero at 55.", 12.5, False, INK, 5),
    ("Validation against the 2021 incidence hold-out, still untouched.",
     12.5, False, INK, 5),
], spacing=1.2)
takeaway(s, "The pipeline is built, checked and passing. What is left is "
            "compute for wave 2 and one structural question about who older "
            "women partner with.", AQUA)

# ==========================================================================
# 22. Closing — what the record is for
# ==========================================================================
s, y = new("", "Three things worth carrying out of this")
textbox(s, MARGIN, y + Inches(0.2), BODY_W, Inches(4.2), [
    ("1.  A negative result at default parameters is not a structural "
     "finding.", 19, True, INK),
    ("Seven experiments measured the same −4.3 pp prevalence deficit and "
     "converged on calling it misspecification. It was a parameter-value "
     "deficit, visible the moment transmission was opened. The experiments "
     "were not wasted — they closed real hypotheses — but the diagnosis was "
     "over-confident because every one of them held the same parameters "
     "fixed.", 14, False, INK2, 5),
    ("2.  Widening a prior can lower coverage.", 19, True, INK, 20),
    ("014 opened three more dimensions than 009 and coverage fell from 34% "
     "to 4%. Added volume in directions that do not matter dilutes the "
     "draws that land where the data lives — and sampling below an "
     "establishment threshold turns a parameter into an on/off switch, "
     "which then dominates every correlation you compute.",
     14, False, INK2, 5),
    ("3.  NROY volume is not informativeness.", 19, True, INK, 20),
    ("87% of the box survived wave 1 and the wave was still informative: "
     "one direction cut by 36%, four untouched. The per-parameter marginals "
     "would have shown almost nothing, because the constraint is diagonal in "
     "the box. Report the constrained directions, not the intervals.",
     14, False, INK2, 5),
], spacing=1.22)
takeaway(s, "The README/SUMMARY discipline is what made all three visible. "
            "Every reversal in this deck was caught by re-reading a "
            "committed record, not by remembering.", AQUA)



# ==========================================================================
# 25. End slide (branded template only -- it carries the logo lockup)
# ==========================================================================
end_slide("Questions, and the parts I am least sure about",
          "Full record in `hivsim_eswatini/experiments/` — one README and one "
          "SUMMARY per experiment, 001 through 025.")

save(OUTFILE)

"""RETIRED 2026-09-07 -- DO NOT RUN. Kept for provenance only.

`eswatini_update_bmgf.pptx` is now hand-maintained: it was edited in PowerPoint
(slides reordered, titles reworded, incidence figure enlarged, caveat footnote
removed) and then patched to add slide 4 and correct the incidence table.
Running this script would silently revert all of that.

Edit the .pptx directly. If the deck ever needs regenerating from scratch, this
file records what it was built from -- but note it produces the pre-edit,
pre-correction 6-slide version, whose incidence table has two wrong numbers
(2011 men on the wrong age basis, 2016 men quoted as an exact match).

Original docstring follows.

---

Build the short internal-update deck (6 slides).

Run: `python presentations/build_update.py [--theme bmgf|neutral]`
Output: presentations/eswatini_update[_bmgf].pptx

The condensed cousin of `build_deck.py`: what changed the model, what we ruled
out, where the fit stands, what is next. Chrome comes from `deck_kit`.

Every number here is traceable to the SUMMARY.md of the experiment named on the
slide. Two framing choices worth knowing about, because both could be
overstated by accident:

* The incidence result is quoted from the **1000-draw prior ensemble median**,
  not the best joint draw. That is deliberately the stronger claim: wave 1
  emulated the prevalence aggregate only, so incidence did not drive the cut,
  and the aggregates land on the median anyway. The best joint draw actually
  runs hot on 2011 (3.90 against a published 3.14 on the 18-49 basis), so
  quoting it would be both weaker and less honest.
* The incidence *age profile* is not right and is not fitted. That caveat rides
  on the same slide as the win rather than being left off.
"""

from deck_kit import *  # noqa: F401,F403

OUTFILE = PRES / f"eswatini_update{SUFFIX}.pptx"


# ==========================================================================
# 1. Title
# ==========================================================================
title_slide(
    "Eswatini HIV model — where the calibration stands",
    "What changed the model, what we ruled out, and the first wave of results",
    "Adam Akullian  ·  7 September 2026",
    "Internal update  ·  model-v1.3  ·  experiments 001–025")

# ==========================================================================
# 2. What shifted the model
# ==========================================================================
s, y = new("Experiments 011–022", "Four things changed the model — three were bugs")
table(s, [
    ["Change", "What was wrong", "What it moved"],
    ["Sexual network age mixing\n(011, 013)",
     "Partners were paired by rank, so the configured age-gap parameter was "
     "inert — sweeping it 0 → 7 → 14 years moved the realised gap only "
     "0.7 → 1.4 → 5.3. Women 35–49 got an impossible −2.7-year gap.",
     "Corrected matcher tracks the configured gap ~1:1 and lands within 1–2 y "
     "of DHS. Like-for-like, it raises 2021 prevalence 4.09% → 6.46%."],
    ["VMMC read as a hazard, not a stock\n(015)",
     "Circumcision coverage was applied as a per-step hazard on the "
     "uncircumcised pool, so it ratcheted to 99.3% against SHIMS3's 48.3%, "
     "and ignored the age gradient entirely.",
     "Male prevalence 15–49 roughly doubles, 2.60% → 6.36%. This, not the "
     "network, was what had been suppressing male infection."],
    ["AIDS mortality counted twice\n(016)",
     "The background death rates fed to the demographic module already "
     "contained AIDS deaths — their adult rates rise 5.5× to a 2005 peak, "
     "which is the shape of the epidemic.",
     "2015 population goes from 13.8% below target to 0.1% above it. "
     "Reconstructing the correction from all-cause data alone reproduces "
     "UNAIDS AIDS deaths to within ~10%."],
    ["Inherited defaults nobody chose\n(018, 021, 022)",
     "A fabricated PrEP ramp to 80% of sex workers starting in 2004 — a "
     "decade before efficacy evidence. Viral suppression defaulted to 1.0: "
     "every ART initiator suppressed.",
     "Little effect on the fit, but the cascade error falls 8.8 pp → 2.1 pp. "
     "Matters directly: the decision question is about raising suppression."],
], y, [Inches(2.75), Inches(4.85), Inches(4.49)], size=9.5,
    row_h=Inches(1.02))
takeaway(s, "These fixed mechanisms rather than the fit level — between "
            "model-v1.1 and v1.3 the prevalence error barely moved. The level "
            "came later, from the calibration itself.", ORANGE)

# ==========================================================================
# 3. What we ruled out
# ==========================================================================
s, y = new("Experiments 014, 017–022", "Six explanations tested and closed")
table(s, [
    ["We thought…", "Verdict"],
    ["HIV mortality parameters are too slow",
     "No — a 6× range on the CD4 death rate moves peak deaths by ρ = −0.01. "
     "Structurally inert — three-quarters of deaths fire from a separate "
     "pathway that no mortality parameter touches."],
    ["Nobody on ART can die of HIV",
     "No — making on-ART mortality nonzero left the deficit at 64% of the "
     "UNAIDS peak, unchanged. The on-ART hazard is anchored at restored CD4, "
     "where it is tiny."],
    ["The new female mortality multiplier is the cause",
     "No — setting it back to 1.0 moves nothing: |z| ≤ 0.9 at every year."],
    ["The whole stack upgrade changed behaviour",
     "No — four arms, each differing by one thing; no comparison reaches "
     "|z| > 1.7 from 1995 to 2021. Safe to adopt without recalibrating."],
    ["Untreated survival is too long",
     "Partly — and it costs more than it buys. Shortening it closes 39% "
     "of the death deficit and makes the prevalence fit monotonically worse. "
     "Three separate routes land on the same frontier."],
    ["The model is too small for the rare strata",
     "No — at 50,000 agents, 2 of 54 strata still fall below 5 expected "
     "cases. The binding one is young men, where the model under-infects by "
     "two thirds — a fit error, not a sizing error."],
], y, [Inches(3.9), Inches(8.19)], size=10, row_h=Inches(0.66))
takeaway(s, "About 22% of the AIDS-death deficit survives every mechanism we "
            "tested. That is why deaths are down-weighted in the calibration "
            "rather than fitted hard.", RED)

# ==========================================================================
# 4. Incidence — the highlight
# ==========================================================================
s, y = new("Experiment 024", "Incidence: the level lands, without being fitted to")
bar = takeaway(s, "Wave 1 emulated the prevalence aggregate only — incidence "
                  "did not drive the cut. All four published adult estimates "
                  "sit on the ensemble median anyway.", AQUA)
picture(s, EXP / "024_hm_wave1/figures/incidence_fit_vs_shims.png",
        y + Inches(0.02), y + Inches(2.62))
table(s, [
    ["HIV incidence, % per year", "Model (prior ensemble median)",
     "Published", "Source"],
    ["2011  women 18–49", "3.1", "3.14  (2.63–3.74)", "SHIMS1 cohort"],
    ["2011  men 18–49", "1.6", "1.65  (1.28–2.11)", "SHIMS1 cohort"],
    ["2016  women 15–49", "1.7", "1.73  (0.96–2.50)", "SHIMS2 Table 5.3.B"],
    ["2016  men 15–49", "0.85", "0.85  (0.21–1.49)", "SHIMS2 Table 5.3.B"],
], y + Inches(2.74), [Inches(3.3), Inches(2.2), Inches(2.6), Inches(3.0)],
    size=11, row_h=Inches(0.285))
footnote(s, "Caveat, on the same slide as the win: the age profile within "
            "2016 is not right — the model puts 2.05 at 15–24 against a "
            "measured 1.67, and 0.35 at 35–49 against 2.09. The total is "
            "correct; its distribution over age is not. The six age-banded "
            "rows are deliberately not fitted (false-recency bias), so the "
            "age profile is currently unconstrained by data.",
         y + Inches(4.30))

# ==========================================================================
# 5. Prevalence
# ==========================================================================
s, y = new("Experiment 024", "Prevalence: one parameter set fits all 48 targets")
bar = takeaway(s, "Two prior predictive checks failed before this one. The "
                  "third passed, and a single draw of 1000 sits within 3σ on "
                  "every registered target at once.", AQUA)
picture(s, EXP / "024_hm_wave1/figures/best_joint_point.png", y,
        bar - Inches(0.15), left=MARGIN, width=Inches(7.6))
textbox(s, MARGIN + Inches(7.95), y + Inches(0.05), Inches(4.15), Inches(4.0), [
    ("Blue is the model, black is PHIA — three national surveys, women on "
     "top, men below.", 12.5, False, INK2),
    ("The headline", 14, True, AQUA, 12),
    ("48 of 48 targets within 3σ. Age-stratified error 0.0405 against 0.0584 "
     "before calibration, and the systematic bias essentially disappears: "
     "−0.007 against −0.040.", 13, False, INK, 6),
    ("The reversal worth knowing about", 14, True, ORANGE, 12),
    ("Seven experiments measured a −4.3 pp prevalence deficit and diagnosed "
     "it as the model being structurally wrong. It was not. It was a "
     "parameter-value deficit — visible the moment transmission was opened "
     "up to the calibration.", 13, False, INK, 6),
    ("Those experiments were not wasted; they closed real hypotheses. But "
     "every one of them held the same parameters fixed, and that is what "
     "made the diagnosis look solid.", 12, False, INK2, 8),
], spacing=1.16)

# ==========================================================================
# 6. Where we are, and what is next
# ==========================================================================
s, y = new("Experiment 025", "Where this stands")
col, gap = Inches(3.85), Inches(0.35)
textbox(s, MARGIN, y, col, Inches(3.8), [
    ("Done", 15, True, AQUA),
    ("Target set settled, including what is excluded and why.",
     12.5, False, INK, 7),
    ("Four model features adopted, each on a measured A/B rather than an "
     "argument.", 12.5, False, INK, 5),
    ("Six explanations for the AIDS-death deficit closed.",
     12.5, False, INK, 5),
    ("Prior predictive check passing at 94%, with a single draw fitting all "
     "48 targets jointly.", 12.5, False, INK, 5),
    ("Model size and the parameter set both settled on measurement — nine "
     "parameters pruned to seven.", 12.5, False, INK, 5),
], spacing=1.2)
textbox(s, MARGIN + col + gap, y, col, Inches(3.8), [
    ("In flight", 15, True, ORANGE),
    ("Wave 2 is written, committed and smoke-tested. It halves the model "
     "discrepancy allowance and adds incidence and the female young:old "
     "ratio as emulated features.", 12.5, False, INK, 7),
    ("It is gated on a synthetic-data recovery test — can the pipeline "
     "recover parameters we already know? That runs first, on purpose.",
     12.5, False, INK, 5),
    ("Blocked on compute. No spot capacity in the region, and the fallback VM "
     "has no account for us yet. About 4 hours a wave on the laptop, which "
     "does not survive closing the lid.", 12.5, True, ORANGE, 7),
], spacing=1.2)
textbox(s, MARGIN + 2 * (col + gap), y, col, Inches(3.8), [
    ("The open question", 15, True, BLUE),
    ("One target nothing in the parameter space reaches: incidence in women "
     "35–49. The total is right, so the mechanism has to move female risk "
     "across age rather than add transmission.", 12.5, False, INK, 7),
    ("Same story for young men's prevalence. Both point at the same thing — "
     "who older women and young men partner with.", 12.5, False, INK, 5),
    ("That is a model-structure question, not a calibration one, and it is "
     "the main thing I would want a second opinion on.",
     12.5, True, INK, 7),
], spacing=1.2)
takeaway(s, "The pipeline is built, checked and passing. What is left is "
            "compute for wave 2 and one structural question about age mixing.",
         AQUA)


save(OUTFILE)

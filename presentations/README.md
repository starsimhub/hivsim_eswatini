# presentations/

## Three decks

Built from the experiment record for experiments 001-025.

| Deck | For | Slides | Maintained by |
|---|---|---|---|
| `eswatini_story[_bmgf].pptx` | Internal talk: how the model came to be | 12 | `build_story.py` |
| `eswatini_update_bmgf.pptx` | Small internal status update | 7 | **Hand-edited - edit the .pptx** |
| `eswatini_calibration_progress[_bmgf].pptx` | Full technical walkthrough | 24-25 | `build_deck.py` |

Shared chrome - theme resolution, BMGF template wiring, layout primitives -
lives in **`deck_kit.py`**.

### Why the story deck is not a re-cut of the technical one

The technical deck is organised by **conclusion**: here are the features, here
is what we ruled out, here is the fit. The story deck is organised by **how we
learned**, which makes different things load-bearing. Three items appear only
there:

- **The compounding-deferral trap of 001-007.** Seven locally correct changes,
  three deferred re-calibrations, and the key parameters left inheriting a
  configuration superseded five times. Invisible if the changes are presented
  as a list.
- **The pattern of repeated misdiagnosis.** 009 blamed mortality, 012 blamed
  the network, 016-022 blamed model structure. All three wrong; all three with
  sound measurements underneath. The technical deck has the 024 reversal but
  not the pattern.
- **`figures/journey.png`** - all 23 experiments on a date axis by outcome,
  against the three coverage-check results. The only figure that shows the
  shape of five months at once.

Its source text is
[`hivsim_eswatini/notes/model_development_narrative.md`](../hivsim_eswatini/notes/model_development_narrative.md).

### The update deck is no longer generated

`build_update.py` is **retired and must not be run**. The deck was edited in
PowerPoint (slides reordered, titles reworded, incidence figure enlarged to full
width, caveat footnote removed) and then patched to insert slide 4 and correct
the incidence table. Re-running the script would revert all of it. The script is
kept only to record what the deck was originally built from.

The full technical deck is still script-generated:

```
python presentations/make_figures.py --theme bmgf
python presentations/build_deck.py   --theme bmgf
```

Swap `bmgf` for `neutral` for the unbranded version.

### Correction applied to the update deck, 2026-09-07

Two of the four model values in the incidence table were wrong and have been
fixed:

| | Was | Now | Why |
|---|---|---|---|
| 2011 men 18–49 | 1.6 | **1.70** | Computed on a 15–50 age basis; the target file specifies `18_50_prorate_first_band_0.4`. Corrected, the model sits slightly *above* published, not below. Same defect experiment 025 lists as its correction #1. |
| 2016 men 15–49 | 0.85 | **0.74** | Quoted as an exact match to the published value; the ensemble median is 0.74. |

A 5–95% column was added at the same time, because the medians land well but
the ensemble is wide (2011 women span 1.9–4.8) and a reader is entitled to see
that on the same slide. All four medians are within 13% of published.

Values recomputed from `024_hm_wave1/outputs/ensemble.parquet` via
`standard_figures._inc_series`, on each target's own age basis.

### What the short update leads with

The incidence result, and it is quoted from the **1000-draw prior ensemble
median** rather than the best joint draw. That is deliberately the stronger and
more honest claim: wave 1 emulated the prevalence aggregate only, so incidence
did not drive the cut, and all four published adult estimates land on the
ensemble median anyway (3.1 vs 3.14, 1.6 vs 1.65, 1.7 vs 1.73, 0.85 vs 0.85).
The best joint draw actually runs *hot* on 2011 — 3.90 against a published 3.14
on the 18–49 basis — so leading with it would be both weaker and less accurate.

The caveat rides on the same slide rather than being left off: the incidence
*age profile* within 2016 is wrong, and the six age-banded rows are deliberately
not fitted, so that profile is currently unconstrained by data.

### On the Gates Foundation theme

`build_deck.py --theme bmgf` **builds onto the real corporate template** rather
than reproducing it. It opens
[`../Parameter_comparison_HIVsim_EMOD.pptx`](../Parameter_comparison_HIVsim_EMOD.pptx),
which carries the `BMGF Layouts` theme and its 70 slide layouts, strips that
file's own slides, and adds ours against named layouts:

| Slide kind | Layout used |
|---|---|
| Title | `Title Slide - Logo Frame, Parchment` |
| Section dividers | `Section Divider, Slate` |
| All content slides | `Blank slide - white` |
| Closing | `End Slide - Slate` |

So the master, colour scheme, type scale, logo lockup and the parchment title
field are all **inherited**, not redrawn — which is the only way to get brand
fidelity right, and it means the deck tracks the template if the template is
updated.

Colours are read out of the template's `theme1.xml` (colour scheme
"Custom 4"), never recalled from memory:

| Role | Theme slot | Hex |
|---|---|---|
| Parchment | `lt1` | `#F5F3ED` |
| Slate | `accent6` | `#303A44` |
| Saffron | `dk2` | `#E7C700` |
| Orange | `accent1` | `#F85C02` |
| Red | `accent2` | `#D93027` |
| Blue | `accent4` | `#248AF9` |
| Turquoise | `accent5` | `#3AC9B1` |

Type: **Calibri Light** for headings, **Calibri** for body — taken from what
the template's own layouts actually specify, not from the theme's
`fontScheme` (which still reads Cambria/Calibri and is not what the layouts
use).

**Caveat worth knowing:** the template is a copy carried inside this repo, not
a live link to whatever the foundation currently ships. If brand guidance
changes, drop in a newer template and rebuild — nothing in `build_deck.py`
hard-codes a layout position, only layout *names*.

### Structure

| Slides | Content |
|---|---|
| 1–3 | Title, the decision question the calibration serves, the target set |
| 4 | **Reference — how to read the fit numbers.** MAE, bias, within-CI, coverage, CV, the agent floor, two-sample z, % of the UNAIDS peak |
| 5–10 | **Part one — features added.** Network age mixing (011/013), VMMC stock-vs-flow (015), mortality double-counting (016), four inherited defaults (018/021/022), then the cross-experiment fit progression |
| 11–13 | **Part two — what we ruled out.** Six hypotheses for the AIDS-death deficit; the deaths-vs-prevalence trade-off |
| 14–15 | Part three divider, then **reference — how to read the history-matching numbers** (σ and its composition, σ_disc, implausibility, NROY, emulator R², variance reduction, ρ, effect-signature r) |
| 16–22 | **Part three — the calibration.** Coverage checks (009/014/024), sizing (020), parameter engineering (023), wave 1 (024), the σ_disc reversal, reading the NROY, surviving defects |
| 23–24 | Where 025 stands, and three transferable lessons |
| 25 | Closing slide — `bmgf` theme only, since it is the layout that carries the logo lockup |

The two reference slides are placed where their metrics are first used, not
collected in an appendix. Their definitions were read from source —
`standard_figures.scorecard`, 014's `assess_coverage`, 020's `cv_table` and
`expected_counts`, 017's reproduction check, 024's `build_observations`, 023's
`effect_signature_confounding` — not paraphrased from the SUMMARY prose.

### Provenance of the numbers

- `mae`, `bias` and `n_within_ci` on the fit-progression slide are recomputed
  by `make_figures.py` from each experiment's surviving `outputs/*.parquet`,
  via `hivsim_eswatini/standard_figures.scorecard` — the same function every
  experiment's own `run.py` calls, so a retrofit and a freshly produced number
  are identical by construction. Written to
  `presentations/fit_progression_scorecard.csv`.
- Peak-AIDS-death shares (the trade-off chart) and the wave-1 PCA (constrained
  directions) are transcribed from the SUMMARY tables of 019/022 and 024
  respectively, because neither is recomputable from the stored parquet without
  re-deriving intermediate quantities. Both are marked as transcribed in
  `make_figures.py`.
- All other figures are the experiments' own committed figures, referenced in
  place from `experiments/NNN_*/figures/`. Nothing is re-rendered.
- Every claim in the slide body text is traceable to the SUMMARY.md of the
  experiment named in the slide's eyebrow.

### Palette

Charts and deck chrome share one palette per theme, three categorical slots
each. Both were run through the `dataviz` validator all-pairs on a light
surface:

| Theme | Slots | Worst CVD ΔE | Worst normal-vision ΔE |
|---|---|---|---|
| neutral | `#2a78d6` `#eb6834` `#1baf7a` | 9.2 | 24.0 |
| bmgf | `#248AF9` `#F85C02` `#3AC9B1` | 16.2 | 23.1 |

Negative/status pole: `#e34948` (neutral), `#D93027` (bmgf).

**Do not add a fourth categorical hue** in either theme — past three slots the
all-pairs floors stop clearing. In `bmgf`, turquoise sits below 3:1 contrast on
white, so any chart using it must carry direct labels (the dataviz relief
rule); the charts here do.

### Known caveat

The fit-progression chart mixes two kinds of object — prior-ensemble medians
(009, 014, 024's prior) and means over seeds at a fixed parameter point (018,
021, 022, 024's best draw). They are coloured by kind and the chart says so,
but a single "progress" axis across both would partly be measuring the change
of object rather than the change of fit. 016 and 017 are absent because their
outputs kept only aggregate infected counts per band, so per-band prevalence is
permanently unrecoverable for those two.

### 025 status as presented

Slide 21 reports 025 as written, committed, smoke-tested and blocked on
compute. The `outputs/preflight_*.csv` files in
`experiments/025_hm_wave2/` are from a **12-point smoke test**, not a real run
— the engine's own collapse warning fired trivially at that size. They are
deliberately not presented as results.


## Three gotchas worth knowing

**Paths are resolved, not counted.** `presentations/` moved from the HIVsim root
into this repo on 2026-09-08, which broke every path derived by counting
directories up from `__file__`. `deck_kit.py` and `make_figures.py` now find the
repo by searching upward for `experiments/` and locate the BMGF template by
name, so moving the folder again will not break them.

**Replacing a picture must drop its relationship too.** Removing a `<p:pic>`
element without also calling `slide.part.drop_rel(rId)` leaves a dangling
`r:embed` pointing at a media part no shape uses. python-pptx and PIL both read
such a file happily and report no problem; PowerPoint offers to repair it, and
a repaired file loses the pictures. The symptom is "The picture can't be
displayed" on exactly the slides whose figures were swapped. If you script a
figure swap, drop the relationship and then assert that each slide's image rels
equal its embeds — `patch_story.py` in the session scratchpad does both.

**If a deck looks broken right after a scripted edit, suspect the cache.**
PowerPoint keeps a document cache keyed on path, and a stale entry renders the
previous state — including a mix of old and new shapes — while the bytes on
disk are correct. Close PowerPoint fully (`Get-Process POWERPNT | Stop-Process
-Force`) before reopening, or open a copy under a different filename. A file
held open by a stuck modal dialog also blocks reads, which makes every
subsequent check look like corruption.

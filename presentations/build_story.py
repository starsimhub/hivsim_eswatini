"""Build the narrative deck: how the Eswatini model came to be (12 slides).

Run: `python presentations/build_story.py [--theme bmgf|neutral]`
Output: $DECK_OUT (default ~/HIVsim_decks)/eswatini_model_story[_bmgf].pptx

Why the output is NOT in this repo: the repo lives under
`OneDrive - Gates Foundation`, so PowerPoint opens files from it via their
SharePoint URL rather than from disk. A locally written file is invisible to
PowerPoint until OneDrive finishes uploading, which shows up as "an error
occurred" and a failing Save As. Writing outside the sync tree avoids that
entirely. Copy the finished deck into the repo when it is done, for provenance.

Language: this deck is for an internal group meeting. The audience is
technical but not in the modelling weeds, so the copy deliberately avoids
calibration vocabulary -- no priors, envelopes, discrepancy, strata, residuals.
Numbers are all kept; only the register changes. If a term has to appear, it is
explained in the same sentence.

Content is drawn from `notes/model_development_narrative.md`, which in turn
comes from the experiment SUMMARYs.
"""

import os
from pathlib import Path

from deck_kit import *  # noqa: F401,F403

# Outside any cloud-synced folder on purpose -- see the module docstring.
# Override with the DECK_OUT environment variable.
OUTDIR = Path(os.environ.get("DECK_OUT", Path.home() / "HIVsim_decks"))
OUTDIR.mkdir(parents=True, exist_ok=True)
OUTFILE = OUTDIR / f"eswatini_model_story{SUFFIX}.pptx"


# ==========================================================================
# 1. Title
# ==========================================================================
title_slide(
    "How the Eswatini model came to be",
    "Twenty-three experiments, three wrong explanations, and what we learned",
    "Adam Akullian  ·  8 September 2026",
    "Internal  ·  April to September 2026")

# ==========================================================================
# 2. The motivating question
# ==========================================================================
s, y = new("The point of the exercise", "Motivating research question")
textbox(s, MARGIN, y + Inches(0.45), Inches(11.4), Inches(2.6), [
    ("Is the next dollar better spent on long-acting PrEP, or on improving the "
     "treatment cascade so more people are virally suppressed?", 26, True, INK),
    ("And where does prevention still matter once suppression is already high?",
     26, False, INK2, 18),
], spacing=1.24, font=HEAD_FONT)
takeaway(s, "The model is not the deliverable. A defensible answer to that "
            "question is.")

# ==========================================================================
# 3. The whole arc in one picture
# ==========================================================================
s, y = new("April to September 2026", "Five months, twenty-three experiments, one checkpoint")
picture(s, FIG / "journey.png", y + Inches(0.10), H - Inches(0.45))

# ==========================================================================
# 4. Act I — seven local fixes
# ==========================================================================
s, y = new("Experiments 001–007, April to early May", "We started by fixing what was obviously wrong")
table(s, [
    ["", "What we changed", "What it taught us"],
    ["001–002", "Treatment coverage split by age and sex instead of one "
     "national figure — and the model then would not reproduce the numbers we "
     "fed it.",
     "A model can fail to match its own inputs. Worth checking every time."],
    ["003", "Age at first sex, from the software's generic 20 and 21 to "
     "Eswatini's own 17.5 and 18.5.",
     "Each year earlier added 0.3–0.5 infections per 100 people. Not enough to "
     "explain why women were infected so much faster than men."],
    ["004", "Three dials on how easily HIV passes: male-to-female risk, condom "
     "use halved, higher risk for women under 25.",
     "We could match the male–female gap or the overall level, but not both. "
     "And two dials were now set against inputs we had since changed."],
    ["005", "Added circumcision, reading the survey as “how many men are "
     "circumcised” rather than “how many get circumcised each month”.",
     "How you read a number has to match what it actually measures."],
    ["006", "Replaced the early-infection transmission figures with published "
     "estimates — roughly halving the extra risk in the first weeks.",
     "More credible biology, worse fit. We kept it anyway."],
    ["007", "Doubled how many people start out infected in 1985. It had been "
     "about eleven.",
     "Changes how certain the early years are, not where the epidemic ends up."],
], y, [Inches(1.05), Inches(5.40), Inches(5.64)], size=12,
    row_h=Inches(0.62))
takeaway(s, "Every one of these was sensible, and most were improvements. The "
            "problem was not any single change.", AQUA)

# ==========================================================================
# 5. The trap
# ==========================================================================
s, y = new("What we did not notice at the time", "Every fix was right on its own. Together they left the model untuned.")
textbox(s, MARGIN, y + Inches(0.35), Inches(7.3), Inches(3.6), [
    ("All seven changes affected how much HIV spreads.", 19, True, INK),
    ("None of them was followed by re-tuning.", 19, True, ORANGE, 12),
    ("Three said “we should re-tune next time”. Each time it slipped.",
     16, False, INK, 16),
    ("By the seventh, the settings that matter most were still the ones we had "
     "chosen five versions earlier.", 16, False, INK, 12),
], spacing=1.24)
textbox(s, MARGIN + Inches(7.8), y + Inches(0.35), Inches(4.29), Inches(3.6), [
    ("Why it is easy to miss", 16, True, ORANGE),
    ("Each change is defensible on its own. The cost only shows up when you "
     "add them together.", 14, False, INK, 8),
    ("Postponing once is reasonable. Three times running is something else.",
     14, False, INK, 10),
    ("What would have caught it", 16, True, AQUA, 16),
    ("A standing rule: change anything that affects transmission, and you "
     "re-tune.", 14, False, INK, 8),
], spacing=1.24)
takeaway(s, "This is the lesson most likely to apply to other models in the "
            "group, and you cannot see it anywhere in the results.", RED)

# ==========================================================================
# 6. Act II — write down what we are aiming at
# ==========================================================================
s, y = new("Experiments 008 and 009, 7 May", "So we stopped changing things and wrote down what we were aiming at")
col = Inches(5.85)
textbox(s, MARGIN, y + Inches(0.15), col, Inches(3.7), [
    ("008 — the data we are aiming at", 17, True, AQUA),
    ("89 numbers: HIV prevalence for 54 age-and-sex groups across three "
     "national surveys, plus 35 years of AIDS deaths — each with a stated "
     "margin of error.", 15, False, INK, 8),
    ("No simulations at all, and still one of the most useful weeks of the "
     "project.", 15, True, INK, 10),
    ("Writing it down surfaced things tuning by hand had not: the death "
     "figures came with no published uncertainty, one survey stops at age 49 "
     "while the others reach 65, and nobody had decided whether to hold 2021 "
     "back as a test.", 15, False, INK, 10),
], spacing=1.24)
textbox(s, MARGIN + col + Inches(0.4), y + Inches(0.15), col, Inches(3.7), [
    ("009 — the first check", 17, True, ORANGE),
    ("We ran 50 versions of the model with settings spread across the "
     "plausible range, and asked what comes before “does it fit”: could the "
     "model produce this data at all?", 15, False, INK, 8),
    ("It reached 30 of the 89 numbers — and missed in a consistent direction, "
     "not randomly.", 15, True, INK, 10),
    ("That ended tuning one thing at a time. The reason we gave for the "
     "failure — too few deaths — turned out to be wrong, and took five more "
     "experiments to disprove.", 15, False, INK, 10),
], spacing=1.24)
takeaway(s, "Writing down what you are aiming at, with honest margins, is "
            "what makes it possible to say a model has failed. Until then we "
            "had no way to be wrong.", AQUA)

# ==========================================================================
# 7. Act III — the detours
# ==========================================================================
s, y = new("Experiments 011–015, June and July", "Two months spent proving our own explanation wrong")
table(s, [
    ["", "What happened", "The lesson"],
    ["011", "We asked a blunt question: does the setting for partner age "
     "difference actually do anything? It did not. Turning it from 0 to 7 to "
     "14 years moved the real gap only 0.7 to 1.4 to 5.3 — and women over 35 "
     "came out paired with younger men, which does not happen.",
     "No point tuning a dial that is not connected. The partnership model had "
     "to be fixed first."],
    ["012", "Re-ran everything on updated software, saw prevalence fall, and "
     "blamed the partnership fix.",
     "We named a cause without running the test that could rule out the "
     "alternatives."],
    ["013", "Tested it by changing only the partnership matching. It pushed "
     "prevalence up, not down — 4.1% to 6.5%. So something else caused the "
     "fall, and we found it: a software update had quietly undone our "
     "circumcision fix.",
     "Three weeks to find out that a believable explanation was the wrong "
     "one. Cheap, next to building on top of it."],
    ["015", "Rebuilt the circumcision fix in our own code rather than editing "
     "the shared package — which is exactly how we lost it the first time. "
     "Circumcision went from 99% to 45%, close to the survey's 48%, and male "
     "prevalence roughly doubled.",
     "Where you put a fix decides whether it survives the next update."],
], y, [Inches(0.95), Inches(6.30), Inches(4.84)], size=12,
    row_h=Inches(1.05))
takeaway(s, "The measurement in 012 was sound. The explanation we attached to "
            "it was not — and that is the pattern for the whole period.",
         ORANGE)

# ==========================================================================
# 8. Act IV — the missing deaths
# ==========================================================================
s, y = new("Experiments 014, 016, 017, 019, 022", "Then six explanations for one number, ruled out one at a time")
bar = takeaway(s, "About a fifth of the gap survives everything we tried. "
                  "Ruling the rest out is why we now let the death figures "
                  "count for less than prevalence.", RED)
picture(s, FIG / "tradeoff.png", y + Inches(0.02), bar - Inches(0.12),
        left=Inches(0.55), width=Inches(7.25))
textbox(s, MARGIN + Inches(7.35), y + Inches(0.05), Inches(4.74), Inches(4.0), [
    ("The model produced about 65% of the AIDS deaths reported for the peak "
     "years. We tested six explanations.", 15, False, INK2),
    ("Ruled out", 16, True, RED, 12),
    ("Death rates too slow — no measurable effect across a six-fold range. "
     "People on treatment could not die of HIV — fixing it closed none of the "
     "gap. A new female death-rate setting — no effect. The software update — "
     "no effect. The model being too small — a fitting problem, not a size "
     "problem.", 14, False, INK, 8),
    ("Partly true, and expensive", 16, True, ORANGE, 12),
    ("Survival without treatment was set at 13 years for everyone, whatever "
     "their age. Shortening it recovers 39% of the missing deaths and makes "
     "the prevalence fit worse. Three different ways of shortening it all land "
     "in the same place, so this is a real feature of the model rather than a "
     "fluke.", 14, False, INK, 8),
], spacing=1.20)

# ==========================================================================
# 9. Act V — the reversal
# ==========================================================================
s, y = new("Experiment 024, 3 September", "Then the explanation seven experiments agreed on turned out to be wrong")
bar = takeaway(s, "Seven experiments concluded the model was built wrong. "
                  "Every one of them had left the transmission settings "
                  "untouched — which is what made them agree.", ORANGE)
picture(s, FIG / "sigma_scan.png", y + Inches(0.02), bar - Inches(0.12),
        left=Inches(0.55), width=Inches(6.70))
textbox(s, MARGIN + Inches(6.55), y + Inches(0.05), Inches(5.54), Inches(4.0), [
    ("Seven experiments measured the same 4.3-point shortfall in prevalence "
     "and concluded the model was built wrong. This run set aside an allowance "
     "to absorb it.", 15, False, INK),
    ("Take the allowance away entirely, and the best version is still four "
     "standard errors short.", 15, True, INK, 12),
    ("The model was not built wrong. We had simply never let the transmission "
     "settings move.", 15, True, ORANGE, 10),
    ("What the run did deliver", 16, True, AQUA, 14),
    ("1000 versions of the model in eight minutes on a 120-core machine. The "
     "data fell within reach for 45 of 48 numbers, and one single version "
     "matched all 48 at once.", 14, False, INK, 8),
], spacing=1.20)

# ==========================================================================
# 10. What we would do differently
# ==========================================================================
s, y = new("", "Three things worth carrying to the next model")
textbox(s, MARGIN, y + Inches(0.25), BODY_W, Inches(4.2), [
    ("1.  Re-tune on a schedule, not when it feels necessary.", 21, True, INK),
    ("Seven sensible fixes, three postponed re-tunings, and the settings that "
     "matter most five versions out of date. No single decision was wrong. "
     "The sum of them was.", 15, False, INK2, 7),
    ("2.  Be suspicious of an explanation formed while everything else is "
     "held still.", 21, True, INK, 22),
    ("We blamed death rates; they turned out to do nothing. We blamed the "
     "partnership model; it moved things the other way. We blamed how the "
     "model was built; it was a setting all along. Every measurement was "
     "sound. Every explanation pointed at whatever we happened to be "
     "changing.", 15, False, INK2, 7),
    ("3.  The dead ends were the most valuable thing we produced.",
     21, True, INK, 22),
    ("Ruling out six explanations is why the one that remains can be stated "
     "precisely. Every correction came from re-reading what we had written "
     "down, not from remembering it.", 15, False, INK2, 7),
], spacing=1.24)
takeaway(s, "Writing a plan before each experiment and findings after it is "
            "the cheapest part of this whole setup, and it paid for itself "
            "repeatedly.", AQUA)

# ==========================================================================
# 11. Where it stands
# ==========================================================================
s, y = new("Experiment 025", "Where the model is now")
col, gap = Inches(3.85), Inches(0.35)
textbox(s, MARGIN, y, col, Inches(3.8), [
    ("Settled", 17, True, AQUA),
    ("What data we are aiming at, and what we are deliberately leaving out.",
     14.5, False, INK, 9),
    ("Four model changes, each tested head-to-head.", 14.5, False, INK, 7),
    ("Six explanations for the missing deaths, ruled out.",
     14.5, False, INK, 7),
    ("How big the model needs to be, and which settings to let vary — nine "
     "narrowed to seven.", 14.5, False, INK, 7),
    ("The check now passes, and one version matches all 48 numbers at once.",
     14.5, False, INK, 7),
], spacing=1.24)
textbox(s, MARGIN + col + gap, y, col, Inches(3.8), [
    ("In flight", 17, True, ORANGE),
    ("The next run: half the allowance for the model being wrong, plus new "
     "infections and the young-versus-older split among women. Bigger model, "
     "20,000 people.", 14.5, False, INK, 9),
    ("First it has to pass a dry run — hide the answer, and see whether the "
     "method finds settings we already know. Two experiments put this off; it "
     "goes first this time.", 14.5, False, INK, 7),
    ("Waiting on compute.", 14.5, True, ORANGE, 9),
], spacing=1.24)
textbox(s, MARGIN + 2 * (col + gap), y, col, Inches(3.8), [
    ("The open question", 17, True, BLUE),
    ("One number nothing we try can reach: new infections in women 35 to 49. "
     "The total is right, so the answer has to move risk between age groups "
     "rather than add more infection overall.", 14.5, False, INK, 9),
    ("Young men's prevalence points the same way. Both come down to who older "
     "women and young men are partnered with.", 14.5, False, INK, 7),
    ("That is a question about how the model is built, not about tuning it — "
     "and the thing I would most like a second opinion on.",
     14.5, True, INK, 9),
], spacing=1.24)
takeaway(s, "Five months to reach a defensible starting point. The real "
            "calibration is one compute allocation away.", AQUA)


end_slide("Questions — and the parts I am least sure about",
          "Every experiment has a written plan and a written finding, in "
          "`experiments/`. The write-up this deck is drawn from is in "
          "`notes/model_development_narrative.md`.")

save(OUTFILE)

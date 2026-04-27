# Session log — calib plugin testing
Date: 2026-04-27
Branch: calibration/ai-plugin-testing

---

## User

Hi Claude, I've been working on this calibration of HIVsim to eSwatini, but could use some help.

## Claude

*[Invoked `calib:getting-started` skill; ran environment and git checks; read README and existing calibration files and figures]*

Good context. A few things I notice before asking anything:

**You're on `main`** with a clean tree — calibration work usually benefits from a branch to keep experiments separate from the model code.

**No `uv` or `pyproject.toml`** — `uv` is installed on your machine but not wired up here. The README recommends conda (which isn't installed on this machine). I'll ask about that below.

**This is inside the model repo itself** (`starsimhub/hivsim_eswatini`) — I'll ask whether you'd like a separate calibration repo.

I also see you've already got an Optuna-based calibration running via STIsim's `Calibration` class with 7 parameters. I won't assume that's the right fit until I understand your research goal.

A few questions so I can actually help you well rather than guess:

**What's the scientific aim?** Are you calibrating to publish a paper, to run intervention scenarios, to do a cost-effectiveness or decision analysis, or something else? And what do you need out of it — a best-fit point, a set of plausible parameter draws, or a full posterior with uncertainty?

**What's going wrong, or what do you want to improve?** You said you "could use some help" — are you having trouble with fit quality, runtime, the calibration workflow, interpreting results, something with the interventions file you have open, or something else entirely?

**Where are you running this?** The README mentions a VM for the 1000-trial Optuna run. Are you working locally right now, or do you have HPC/VM access? Roughly how long does one sim take on your current hardware?

**Who else is on this, and is there a deadline?** Just you, or a team? Any timeline pressure?

---

## User

The aim of this work is to explore targeting of PrEP to FSW. There's a lot of uncertainty, a single best-fit point probably isn't good enough.

This is NOT the model repo, instead it's dedicated to fitting the HIVsim model to eSwatini specifically.

Look at the plots in figures/ and you'll see that the model is not fitting the data.

The model is pretty fast, but I can run these on a 120 core VM if needed.

No specific deadline, but I'd like to have a good calibration in the next week or two.

## Claude

*[Inspected figures in `figures/`, read calibration data and `run_sims.py`]*

Okay — I have a clear enough picture now. Let me tell you what I see and what I'd propose.

**What I'm hearing:** You have a Starsim-based HIV model for eSwatini, calibrating to UNAIDS/PHIA data, aiming to generate a posterior over plausible epidemics so you can run PrEP-targeting scenarios with proper uncertainty. Target timeline is 1-2 weeks.

**What the figures show:** The epidemic is collapsing in the model. Three specific problems:

1. **Incidence goes near zero by ~2020** while PHIA data shows sustained female incidence of ~1.7%/year in 2016. This is the root failure — the model can't reproduce an epidemic that keeps circulating.

2. **PLHIV and ART absolute counts are roughly half the data** by 2020-2025 — yet ART *coverage* (%) tracks well. Those two facts together confirm that the model has too few PLHIV, not a problem with the ART intervention logic.

3. **Age-specific prevalence overestimates PHIA data** in all survey years, which is in tension with #2 above — this likely means the model has too few people in the denominator (population too small) rather than too many infected people.

**What I think is causing it:** `rel_init_prev` is fixed at `0.1` in `run_sims.py` and was *commented out* of the calibration parameters. That's a very low seed. Combined with the epidemic burning itself out (insufficient sustained transmission), the model never reaches the plateau eSwatini showed through the 2000s-2020s.

**The method problem:** Optuna gives you a single best-fit point. For PrEP decision analysis you need a *posterior* — a set of parameter draws that collectively span plausible epidemic trajectories. Optuna is the wrong tool for the stated aim. You need ABC or something equivalent.

**My proposal:**

**Structure:** Stay in this repo, create a working branch rather than continuing on `main`.

**Method:** ABC with the 120-core VM. Your model is fast enough that ABC rejection or ABC-SMC is very tractable at scale.

**First experiment — coverage check:** Before touching the calibration, run ~100 prior draws and plot them alongside the data to answer "does the data fall inside the prior predictive envelope?" If not, fix the model first. If yes, proceed directly to ABC.

**Environment:** Set up a `pyproject.toml` with `uv` for reproducibility.

Does this match what you're after, or should we adjust anything? If yes, I'll: create branch `calibration/prior-check`, set up `pyproject.toml` + `uv sync`, scaffold `experiments/01_coverage_check/`, and write the prior-predictive script.

---

## User

- OK, we can stay in this repo. Please just note that this repo is 1) managed by a different team, and 2) public (!), so let's be sure to do any work in a feature branch that can be easily reversed.
- Sure, we can try ABC. I've only ever used Optuna, so will need lots of help.
- A coverage check sounds like a great idea, yes -- I like that.
- I haven't tried uv before, but sure -- let's do that too.

## Claude

*[Created branch `calibration/prior-check`; wrote memory files; created `pyproject.toml`; ran `uv sync` to wire up environment; wrote `experiments/01_coverage_check/run.py` with 50 prior draws on 4 CPUs]*

Now let me propose the CLAUDE.md. Per the rules I need to show it to you before writing.

**Proposed CLAUDE.md** captures: model description, research aim (PrEP targeting / posterior), data sources, constraints (public repo, branch-only work, VM available), current fit status (epidemic collapse, `rel_init_prev` suspected root cause), and environment setup (`uv`).

A few things worth knowing before you run it:

- It runs **50 prior draws in parallel on 4 CPUs** — should take a few minutes locally. On the VM you can raise `N_DRAWS` and `ncpus` significantly.
- The key parameter to watch is **`rel_init_prev`** — it's being sampled from 0.05–5.0. If the high end of the prior band reaches the data, that tells us re-enabling it in calibration will likely fix the collapse.
- Output goes to `experiments/01_coverage_check/figures/01_coverage_check.png`.

---

## User

- Oh okay, I have 10 CPUs on my laptop, can we start with that?
- Can you keep an eye on that parameter for me?
- Let's go -- I hope this works!

## Claude

*[Updated `run.py` to use 10 CPUs; wrote CLAUDE.md; attempted to run — hit a `sc.parallelize` argument-unpacking bug, then a `rel_dur_on_art` parameter-not-found error]*

`rel_dur_on_art` isn't a direct `sti.HIV` parameter — it must live somewhere else. Let me check.

*[Searched stisim/starsim source for the parameter; discovered stisim version may be outdated]*

---

## User

my stisim could be old, you might need to switch to main and get the latest

## Claude

*[Switched stisim to `main` and pulled; confirmed upgrade from old version to 1.5.3; verified `rel_dur_on_art` and corrected `ART(coverage=...)` API are now present]*

Good, stisim is now 1.5.3. Both issues are fixed in 1.5.3: `coverage=` is now correct for `ART` and `rel_dur_on_art` is a real HIV parameter. Let's run it.

*[Ran 50 draws — all completed; plot failed on column name mismatch; inspected results columns; found `timevec` instead of `year`, `hiv.new_deaths` instead of `hiv.hiv_deaths`; fixed plot function; regenerated figure]*

Here's what the coverage check is telling us.

**The prior band is close but doesn't quite reach the data.** Key findings:

- **`beta_m2f` is the real driver** — correlation with peak PLHIV is 0.92 vs only 0.09 for `rel_init_prev`. The upper end of your current beta range (0.014) gets PLHIV to ~149K in 2020, but the data is ~200K. The prior needs `beta_m2f` widened up to about 0.020.

- **`rel_init_prev` barely matters for epidemic size.** It affects initial seeding speed but the epidemic quickly reaches a quasi-equilibrium set by `beta_m2f` and the network. You should put it back in the calibration but the range doesn't need to be wide — something like 0.5–3.0 is fine.

- **The epidemic no longer collapses.** Late incidence (2015–2020) reaches up to ~6,400/year across the prior, vs data of ~5,000–8,000. That's a pass.

- **Population is close** — median ~1.05M vs data ~1.15M in 2020.

Proposed updated prior: widen `beta_m2f` upper bound to 0.020, add `rel_init_prev` 0.5–3.0. Ready to go straight to ABC with these ranges, or do you want to tweak anything first?

---

## User

Hmm, I see a few additional things:

1. HIV Prevalence starts too early and is already falling by 1995 -- the prevalence data shows increases through 2018.
2. New HIV infections is also too early and too low.
3. PLHIV peaks and flattens from the model, but the data show a consistent increases.
4. Population -- when HIV deaths hit hard around year 2000, the modeled population actually decreases! We probably need to input HIV-cause deleted data instead of whatever input is currently being passed to the model.

## Claude

These are sharp observations — especially #4, which is a classic modeling trap. Let me investigate all four before proposing fixes.

*[Analyzed init_prev file — FSW initial prevalence up to 20%; checked death rate file — male age-35 mortality triples from 0.005 in 1985 to 0.027 in 2005; confirmed starsim background deaths apply to all alive agents including HIV+; confirmed double-counting]*

Here's the full diagnosis of all four issues:

**Issue 1 & 2 (epidemic too early):** `rel_init_prev` upper bound of 5.0 seeded 75–100% FSW infection at 1985 — obviously wrong. The entire timing problem in the coverage check came from the high end of that prior. Fix: narrow to 0.05–0.5.

**Issue 4 (population dip / double-counted deaths):** The death rate file spikes from 0.005 to 0.027 for males age 35 between 1985 and 2005 — that's the AIDS epidemic embedded in all-cause mortality. When the HIV disease module also kills those same people, they die twice. Fix: replace with cause-deleted rates using log-linear interpolation through the AIDS spike years (1985 rates are pre-AIDS, 2035+ rates are post-ART).

**Issue 3 (PLHIV plateaus):** Will improve automatically once 1, 2, and 4 are fixed.

*[Generated cause-deleted death rates; confirmed male age-35 flat at ~0.005 vs original spike to 0.027; saved as `data/eswatini_deaths_cause_deleted.csv`; replaced `data/eswatini_deaths.csv` with cause-deleted version; preserved original as `eswatini_deaths_allcause.csv`; narrowed `rel_init_prev` range; updated coverage check script for 100 draws on 10 CPUs]*

---

## User

Don't commit until later. Are you going to write a summary at the end of 01 before going on to 02?

## Claude

Yes, absolutely — that's the right discipline. Let me write it now.

*[Wrote `experiments/01_coverage_check/SUMMARY.md` with full diagnosis: four issues identified, root causes, fixes applied, and what changed as a result]*

Done. `experiments/01_coverage_check/SUMMARY.md` records what we ran, what we found (the four specific issues with numbers), and what changed as a result.

Ready to scaffold experiment 02 whenever you are.

---

## User

Good, yes -- let's keep going

## Claude

*[Scaffolded `experiments/02_coverage_check_v2/run.py` with 100 draws, updated prior; ran 100 draws successfully in ~311 s on 10 CPUs; generated figures]*

This is a substantial improvement. Here's what I see:

**What's now working:**
- **HIV prevalence 15–49** — data falls squarely inside the prior band, including the late-rising trajectory through 2010–2015. Timing is fixed.
- **PLHIV** — data (up to ~220K) now falls within the 90th percentile band. The cause-deleted deaths eliminated the premature epidemic burn-out.
- **Total population** — no more dip during the AIDS era. Data sits comfortably inside the band.
- **HIV deaths** — data is within the band through the peak (2000–2010).

**What's still slightly off:**
- **New HIV infections** and **HIV deaths** post-2015 — the data declines faster than the prior median suggests. The band still covers it at the low end, but this will be ABC's job to resolve.
- **People on ART** — no data points in the calibration CSV yet, so the panel is just a prior range.

**Verdict: coverage check passes.** Ready to write the experiment 02 SUMMARY and then move on to setting up the ABC run?

---

## User

Hmm, I see a few more issues:

1. Total population still slows down too much on most runs. The coverage check passes because some sims have HIV mortality that is too low, I suspect.
2. HIV Prevalence still drops too quickly.
3. HIV incidence peaks too soon and falls too fast.

## Claude

Good observations — let me dig into what's driving each one before proposing changes.

*[Ran diagnostics: confirmed population growth 100% of draws, but late-peak draws (timing matches data) show median PLHIV of only 48K in 2020 vs 200K in data; computed correlations — `beta_m2f` r=−0.646 with peak prevalence year; read ART coverage data and network parameters]*

That's the crucial finding — confirmed. Late-peak draws (timing matches data) have median PLHIV of only 48K in 2020 vs. 200K in the data. The model cannot simultaneously get the timing right and the magnitude right with the current parameter set. This is a **timing-magnitude tradeoff** driven by network saturation: high `beta_m2f` spreads the epidemic fast and hits the right magnitude but peaks too early; low `beta_m2f` gets the right timing but the epidemic is far too small.

*[Reviewed `art_coverage.png` and figures from experiment 02]*

The timing-magnitude problem is deeper — the epidemic in the model doesn't sustain transmission past 2005 at scale. This points to how HIV spreads from high-risk (FSW/clients) into the medium-risk general population. The key suspects are `f1_conc` (medium-risk female concurrency, currently fixed at 0.15) and possibly `p_pair_form` (pair formation probability).

*[Wrote `experiments/02_coverage_check_v2/SUMMARY.md`; scaffolded `experiments/03_coverage_phia/run.py` with PHIA data as primary targets and `p_pair_form` added to the prior; verified `p_pair_form` is a `ss.bernoulli` parameter; fixed its representation in the sim builder; ran 100 draws successfully]*

*[Looked at figures from experiment 03]*

**Coverage plot findings (8 panels with PHIA orange diamonds as primary targets):**
- **Prevalence 15-49**: The model band covers PHIA 2007 and 2011 points, but many runs decline too quickly — PHIA 2016 and 2021 sit below the prior median.
- **Female prevalence 20-25**: PHIA shows ~30-37% through 2007-2016. The model band just barely covers these at the upper end at 2007, then collapses too fast.
- **Female incidence 15-49**: PHIA 2011 (~3.1%/yr) sits near the model median; PHIA 2016 (~1.7%/yr) is now below the model median.

**Diagnostic plot:** `p_pair_form` doesn't resolve the tradeoff. Same constraint holds: draws with late peaks have low PLHIV.

**Conclusion**: `p_pair_form` is not the structural fix. The prior band **does not adequately cover** the 2016/2021 PHIA prevalence for women aged 20-25, and the incidence decline in the model is too fast. Next hypothesis: `f1_conc` (medium-risk female concurrency).

Want me to write the experiment 03 SUMMARY and then for experiment 04, add `f1_conc` to the prior?

---

## User

Yes, please write the summary -- but then stop there for now.

## Claude

*[Wrote `experiments/03_coverage_phia/SUMMARY.md`]*

Done. The summary captures the two findings — `p_pair_form` is not the structural fix, and the prior fails to cover the PHIA 2016/2021 female prevalence 20–25 and the 2016 incidence points — and documents `f1_conc` as the next hypothesis to test when you're ready for experiment 04.

---

## User

Tell me what happened with the data. Is `eswatini_deaths_cause_deleted.csv` needed? `eswatini_deaths_allcause.csv`??

## Claude

*[Confirmed `eswatini_deaths.csv` and `eswatini_deaths_cause_deleted.csv` are byte-for-byte identical]*

Three files, two of which matter:

| File | Status | Contents |
|---|---|---|
| `eswatini_deaths.csv` | **Active — used by the model** | Cause-deleted rates (AIDS spike smoothed out) |
| `eswatini_deaths_allcause.csv` | Reference backup | Original all-cause rates with the AIDS-era mortality spike |
| `eswatini_deaths_cause_deleted.csv` | **Redundant duplicate** | Identical byte-for-byte to `eswatini_deaths.csv` |

**What happened**: We replaced `eswatini_deaths.csv` with cause-deleted rates (log-linear interpolation through the AIDS spike years) to stop the model double-counting HIV deaths. The original was preserved as `_allcause.csv`. We also left behind `_cause_deleted.csv` as an intermediate artifact that is now a duplicate.

**What you can safely delete**: `eswatini_deaths_cause_deleted.csv`.

**What you should keep**: `eswatini_deaths_allcause.csv` — it's the only copy of the original data.

---

## User

OK, anything else we can clean up before committing this branch?

## Claude

*[Ran `git status` and `git diff`; checked `.gitignore`; listed experiment results sizes (01: 2.1 MB, 02: 4.1 MB, 03: 4.1 MB in `.obj` binary blobs)]*

Here's the full cleanup list before committing:

**Delete:**
- `data/eswatini_deaths_cause_deleted.csv` — exact duplicate of `eswatini_deaths.csv`, serves no purpose

**Add to `.gitignore`:**
- `experiments/*/results/` — the `.obj` files are 2–4 MB binary blobs; the run scripts reproduce them, so there's no reason to track them

**Revert or keep:**
- `run_sims.py` — the only change is `.1` → `0.1` (cosmetic). Not worth including on its own; I'd revert it to keep the diff clean.

Everything else (`CLAUDE.md`, `pyproject.toml`, `uv.lock`, `data/eswatini_deaths.csv`, `data/eswatini_deaths_allcause.csv`, `experiments/` including figures and SUMMARY.md files) looks right to commit.

Want me to make those three changes and then stage everything for your review?

---

## User

- Delete the extra csvs.
- Ignore the obj files, but do include the figures
- Revert run_sims.

## Claude

*[Deleted `eswatini_deaths_cause_deleted.csv`; updated `.gitignore` to add `experiments/*/results/`; reverted `run_sims.py` to HEAD]*

Clean. The `results/` dirs are now ignored, figures are included under `experiments/`, and `run_sims.py` is reverted. Everything looks right — say the word when you want to commit.

---

## User

Can we also rename this branch? Really what I was doing here was testing the new calib plugin. So `calibration/ai-plugin-testing` could be a good name.

## Claude

*[Renamed branch from `calibration/prior-check` to `calibration/ai-plugin-testing`]*

Done — now on `calibration/ai-plugin-testing`. Ready to commit whenever you are.

---

## User

sure, go for it!

## Claude

*[Staged all relevant files; committed as `895df97` — 16 files changed, 8052 insertions, 2422 deletions]*

Committed as `895df97`. 16 files, clean working tree. The branch `calibration/ai-plugin-testing` is ready — just needs a push whenever you want to share it.

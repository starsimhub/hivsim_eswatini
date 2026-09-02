# Exp 020 — Size the model: population, replicates, and the rare-event floor

*The sweep [016](../016_double_counted_mortality/SUMMARY.md) queued,
[017](../017_version_bump/SUMMARY.md) moved, and
[018](../018_adopt_and_size/SUMMARY.md) deferred. It has now been deferred three
times, each time because the model underneath it changed. model-v1.1 is fixed,
so it runs here.*

## Question

**How large must the model be, and how many replicates per parameter set, for a
coverage check and then a calibration to mean anything?**

Three sub-questions, in dependency order.

**1. Does N = 10 000 clear the rare-event floor?** 018 said yes on a
pre-computed estimate. That estimate was wrong twice over:

- Its sex mapping was inverted (`Gender == 0` read as female), so it paired
  model male counts against female target rows. Corrected in 018's SUMMARY
  after the fact.
- It was computed at **high-transmission** parameters (prevalence ~0.22), and
  every stratum count scales down with prevalence. Anywhere the prior explores
  lower prevalence, the floor is closer than the estimate implies.

With the mapping fixed, [019](../019_age_dependent_survival/SUMMARY.md)
observation 7 measured **2 of 54 PHIA target rows below 5 expected infected
agents** (2007 M 15–19 at 2.4, 2007 F 60–64 at 4.6) and 3 below 10 — at high
transmission, the favourable case. The thin strata are the young male bins,
where PHIA prevalence is ~2% and the model puts 0.003 against PHIA's 0.019, so
the count is thin from both directions at once.

**2. How many replicates?** 018 obs 4–5 gave a provisional answer — CV 4.4% at
high transmission, 7.8% at default, both on model-v1.1, implying 10–20 — but
from two points only, and both inside the establishing region. The
low-transmission point near `beta_m2f` = 0.008 is exactly where 018's
establishment map says variance should climb, and it is unmeasured.

**3. What does run time and memory do as N grows?** Anchors every compute
estimate downstream, including whether history-matching waves fit on a laptop
or need raccoon. The `StructuredSexual` pairing step is the term to watch: an
O(N²) contact process would rule out N = 50 000 regardless of the floor.

## What is already known

| quantity | value | source |
|---|---|---|
| run time at N = 10 000, 1985–2026 | ~93–145 s | 017 (80 runs), 018 (20), 019 (40) |
| CV, high transmission (`beta_m2f` 0.0139) | 4.4% | 018 obs 4 |
| CV, stisim defaults | 7.8% | 018 obs 4 |
| CV, 1.5.8 stack at defaults | 45.9%, 2/10 seeds extinct | 017 arm A — **the stack this replaces** |
| establishment floor | `beta_m2f` ≥ 0.008 establishes 5/5 at every `rel_init_prev` | 018 obs 3 |
| thinnest target stratum at N = 10 000 | 2.4 expected infected agents | 019 obs 7 |

## Plan

### Part A — population size

N ∈ {5 000, 10 000, 20 000, 50 000} × 10 seeds, at **two** parameter points,
not one:

- `high_transmission` (`beta_m2f` 0.0139, `rel_init_prev` 0.49) — continuity
  with 016–019, and the favourable case for the floor.
- `low_transmission` (`beta_m2f` 0.008, `rel_init_prev` 0.2) — just above the
  establishment threshold, where prevalence is lowest and the floor bites
  hardest. **This is the case that decides N**, and 018 never ran it.

Measured per (N, parameter point): expected infected agents per PHIA stratum
with the corrected sex mapping, between-seed CV of trajectory-mean prevalence,
wall time, and peak resident memory.

### Part B — replicate count

Three parameter points × 10 seeds at whichever N part A selects:
`plausible` (0.0139 / 0.49), `low_transmission` (0.008 / 0.2), `default`
(0.01 / 0.2). Read replicates off the standard bands — CV < 5% → 3–5,
5–20% → 10–20, > 20% → 50+ or increase N — and report per point rather than as
one number, since 018 obs 4 showed the CV is parameter-dependent.

### Compute constraint, stated up front

The laptop has 12 logical cores but **~6 GB RAM free of 34**. Peak memory per
sim is unmeasured above N = 10 000, so the N = 50 000 arm is the risk: at 10
concurrent workers it could exhaust memory mid-sweep. Part A therefore runs N
groups **sequentially with a worker cap that falls as N rises** (10 / 8 / 4 / 2),
and records peak RSS so the question is answered rather than guessed. Every cell
writes its own parquet, so an OOM costs only the cells in flight.

If peak RSS at N = 50 000 rules out the laptop, that is itself a finding and
hands a concrete number to `idm-azure` for raccoon.

## Success criteria

- **Clean:** a defensible N with every PHIA stratum above ~10 expected infected
  agents at the *low*-transmission point, a per-point replicate count, and
  linear-or-better run-time scaling. Coverage check v3 opens immediately after.
- **Awkward but useful:** the floor needs N = 50 000, or 20+ replicates in the
  plausible region. A real compute finding that changes what method is feasible
  and goes to `method-selection`.
- **Blocking:** run time or memory scales worse than linearly, so no N both
  clears the floor and fits the compute budget. Then either the thin strata get
  dropped from the target set — with that stated as a limitation — or the
  calibration moves to raccoon.

## Not in scope

- Any calibration run, prior bounds, or parameter list — `parameter-engineering`
  (021), which is in progress in parallel and does not depend on this.
- The AIDS-death deficit. 019 closed 39% of it and established that deaths and
  prevalence are in structural tension; the residual is a *shape* problem in
  time, for a later experiment.
- `p_effective_art = 1.0`, the missing suppression gap (019's Next). A data
  task, not a sizing one.

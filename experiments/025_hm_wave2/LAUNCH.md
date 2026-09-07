# How to launch 025, once a VM is available

Everything is committed and smoke-tested. This is the only thing standing
between 025 and results.

## Blocker as of 2026-09-07

- **raccoon** — fails to allocate. `OverconstrainedAllocationRequest`, naming
  *Low Priority VMs* and *Preemptible VMs*: no HB120 **spot** capacity in the
  region. Not a config problem.
- **zebra** — machine is up (the SSH daemon answers), but **there is no account
  for us on it.** Both `~/.ssh/aakullian_raccoon.pem` and `~/.ssh/id_ed25519`
  are refused, as both `adamak` and `aakullian`. The fleet has no shared user
  database, so a raccoon account does not imply a zebra one. Needs Cliff or Dan.
  *Note:* repeated login attempts on 2026-09-07 tripped SSH rate-limiting there;
  leave it alone for a while before retrying.
- All eight HB120 spot machines are deallocated, so which will *allocate* is
  unknown until one is tried.

## Fastest unblock

Start **one** other HB120 from https://selfserve.starsim.org — `gerbil` or
`hedgehog`. The raccoon key will most likely work on those, since they were
provisioned in the same batch. If it allocates, run the setup below.

Do not try more than one or two: they are the same size in the same
subscription, so a spot-capacity failure on one is likely a failure on all.

Non-spot fallbacks if spot capacity stays unavailable: `zebra` (160 cores, needs
an account), then `dugong` (176) — the fleet notes mark dugong as expensive, so
justify it.

## Why not local

~4 h per wave on the laptop's 8 workers at N = 20,000 (024 did 1000 points in
456 s on 120 cores at N = 10,000). And closing the lid suspends the machine, so
the run stops. Locking the screen is fine; closing it is not.

## Setup on a fresh HB120

The repo layout is what `pyproject.toml`'s `[tool.uv.sources]` relative editable
paths require — starsim and stisim must sit alongside the project:

```bash
IP=<from the fleet inventory>
ssh aakullian@$IP 'bash -s' <<'EOF'
set -e
export PATH=$HOME/.local/bin:$PATH
command -v uv || curl -LsSf https://astral.sh/uv/install.sh | sh
mkdir -p ~/work/HIVsim ~/work/star_sim
cd ~/work/HIVsim && [ -d hivsim_eswatini ] || git clone https://github.com/starsimhub/hivsim_eswatini.git
cd ~/work/star_sim
[ -d starsim ] || git clone https://github.com/starsimhub/starsim.git
[ -d stisim ]  || git clone https://github.com/starsimhub/stisim.git
# Pin to the stack 024 ran on. Detached HEAD is correct -- it pins the dependency.
cd ~/work/star_sim/starsim && git fetch -q --all && git checkout -q ceddb877
cd ~/work/star_sim/stisim  && git fetch -q --all && git checkout -q 97b8bc1
cd ~/work/HIVsim/hivsim_eswatini && git fetch -q --all && git checkout -q main && git pull -q
uv sync -q
uv run python -c "import starsim, stisim; print(starsim.__version__, stisim.__version__)"
EOF
```

Expect `3.5.2 1.5.11`. **HTTPS, not SSH, for GitHub** — agent forwarding does not
carry a usable GitHub key to the VMs, and all three repos are public.

## Launch

Run from the **repo root**, not the experiment directory — `run_sims.make_sim`
reads `data/condom_use.csv` by relative path.

```bash
ssh aakullian@$IP 'cd ~/work/HIVsim/hivsim_eswatini && \
  tmux new-session -d -s w2 "export PATH=\$HOME/.local/bin:\$PATH; \
  uv run python experiments/025_hm_wave2/run.py \
    --part both --n_workers \$(nproc) > ~/w2.log 2>&1; \
  echo EXIT=\$? >> ~/w2.log"'
```

`tmux` is not optional on a spot machine — a reclaimed session reattaches, a bare
SSH session is gone. Note a `tmux new-session -d "<cmd>"` session **disappears
when the command finishes**, so a missing session means done-or-died; read the
log tail before concluding anything.

Check progress:

```bash
ssh aakullian@$IP 'tail -40 ~/w2.log; echo ---; tmux ls'
```

Expect ~40 min total for both parts at N = 20,000 on 120 cores. Per-point parquet
is cached under `experiments/025_hm_wave2/outputs/sims/`, keyed on a hash of the
parameter values **plus N, stop year and seed**, so an interrupted run resumes
with `--resume` and costs nothing for work already done.

## Retrieving results

Consolidate on the VM first — 024's 1000 per-point files were 118 MB loose and
6.6 MB concatenated:

```bash
ssh aakullian@$IP 'cd ~/work/HIVsim/hivsim_eswatini && uv run python -c "
import glob, pandas as pd
f = sorted(glob.glob(\"experiments/025_hm_wave2/outputs/sims/*.parquet\"))
pd.concat([pd.read_parquet(x) for x in f], ignore_index=True).to_parquet(
    \"experiments/025_hm_wave2/outputs/ensemble.parquet\", index=False)
print(len(f), \"points consolidated\")"'
```

Then `scp` back (there is **no rsync** on the Windows laptop):

```bash
scp aakullian@$IP:work/HIVsim/hivsim_eswatini/experiments/025_hm_wave2/outputs/{ensemble.parquet,observations.csv,preflight_recovery.csv,preflight_degeneracy.csv,preflight_nroy_summary.txt,wave2_nroy_summary.txt} \
    experiments/025_hm_wave2/outputs/
```

## What to look at first, in order

1. **`preflight_recovery.csv` and `preflight_degeneracy.csv`** — the gate. Is the
   truth inside the NROY, and is the `beta_m2f`/`s_f_young` **ratio** recovered as
   well as the product? If only the product is, wave 2's marginals in those two
   directions are not interpretable and must be reported as jointly constrained
   only. This is pre-registered, not a post-hoc excuse.
2. **Per-feature emulator R²** in the wave-2 log. Anything under 0.8 should be
   dropped rather than carried — three simultaneous features is already against
   the 1–2 guidance.
3. **NROY fraction.** Under 5%, or the engine's own `NROY space collapsed`
   warning firing at n=1000: suspect σ_disc = 0.01 before suspecting the model.
   024's scan predicts only 8.4% of draws survive on prevalence alone at that σ,
   so collapse is a live possibility.
4. **`prev_young_old_f_mean`** against its target of 0.503. The measured ratio
   collapses 0.739 → 0.506 → 0.265 across 2007/2011/2016, and `s_f_young` is a
   *static* scalar being asked to match a strongly *time-varying* shift in who
   carries prevalence. If this feature can't be reached, look at that asymmetry
   before concluding anything about the parameter value.

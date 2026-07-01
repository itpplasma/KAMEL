# Golden-record CI for KIM, QL-Balance, and KiLCA

This harness locks the *physics output* of the three codes so any change that
moves a result beyond floating-point tolerance turns CI red — even when the
change is a "pure refactor" (e.g. swapping a numerical backend).

It runs **only** in the dedicated `golden-record` GitHub Actions job
(`.github/workflows/golden.yml`). It is intentionally **not** registered with
CTest, so `ctest` / `make test` stay fast and offline.

## The A/B model

Each run builds the repo **twice in the same environment** so compiler/toolchain
effects cancel:

| build | git ref | meaning |
|-------|---------|---------|
| `ref` | tag `golden-baseline` | the frozen, blessed physics baseline |
| `cur` | the PR head (`github.sha`) | the candidate |

Every code runs a small synthetic case against *both* builds; the two output
trees are diffed. Any file that diverges beyond tolerance fails the job.

```
run_golden.sh [REF_GITREF] [CUR_GITREF]      # defaults: golden-baseline  HEAD
  └─ bin/gr_build_all.sh ref  <REF_GITREF>   # one whole-repo build -> build_ref/
  └─ bin/gr_build_all.sh cur  <CUR_GITREF>   # one whole-repo build -> build_cur/
  └─ for each code (kim, qlbalance, kilca):
       bin/gr_dispatch.sh --backend local    # run every case vs ref and cur
       bin/gr_compare.sh                      # diff out/ref vs out/cur
```

`gr_build_all.sh` builds the three targets once each — `KIM_exe`,
`ql-balance.x`, `kilca_normal_exe` — and drops stable exe symlinks (`KIM.x`,
`ql-balance.x`, `KiLCA_Normal.x`) so each code finds its binary with no per-code
rebuild. (`KIM`'s logical target is `KIM_exe`; its installed binary is `KIM.x`.
`KiLCA`'s binary carries a version string, so the symlink resolves the glob
`KiLCA_Normal_*_64bit`.)

## Layout

```
test/golden/
├── run_golden.sh            # top-level orchestrator
├── bin/                     # code-agnostic gr-runner (one physical copy)
│   ├── gr_build_all.sh      # whole-repo build per ref/cur
│   ├── gr_dispatch.sh       # discover + run cases (local backend)
│   ├── gr_run_case.sh       # run ONE case (copy inputs, prepare hook, exe)
│   ├── gr_compare.sh        # diff out/ref vs out/cur via $GR_COMPARE
│   ├── gr_numcompare.py     # recursive numeric + .uff byte comparator
│   ├── gr_migrate.py        # optional input-API migration
│   ├── compare.py           # QL-Balance HDF5 comparator
│   ├── compare_qlbalance.sh # locates each side's *.hdf5, calls compare.py
│   └── setup_runfolder.py   # QL-Balance synthetic input generator
├── kim/        {bin->../bin, config.sh, cases/electromagnetic/}
├── qlbalance/  {bin->../bin, config.sh, cases/default/}
└── kilca/      {bin->../bin, config.sh, cases/flre_m6n2/}
```

Each `<code>/bin` is a symlink to the one physical `bin/`. Because every
`gr_*.sh` computes `ROOT="$(dirname BASH_SOURCE)/.."`, invoking
`kim/bin/gr_dispatch.sh` makes `ROOT=.../kim`, so it picks up `kim/config.sh`.
One harness, three logical roots.

## The `golden-baseline` tag

`ref` is the annotated tag `golden-baseline`. It is the contract: "this is the
correct physics." To **re-bless** (accept a deliberate physics change as the new
truth):

```bash
git tag -f -a golden-baseline <new-sha> -m "Re-bless golden baseline (<date>, <why>)"
git push -f origin golden-baseline
```

Re-blessing is the *only* way the baseline moves; nothing else in this harness
changes.

## Tolerance

- KIM and KiLCA: `gr_numcompare.py` at `rtol 1e-7 / atol 1e-12` (per-case
  override via extra args to `gr_compare.sh`).
- QL-Balance: `compare.py` compares a curated quantity list (`QUANTITIES_TO_COMPARE`)
  at `rtol 1e-8 / atol 1e-15`; see `qlbalance/cases/default/tol.yaml` (informational).
- KiLCA `formfactors.*.uff` (Fortran binary) is compared byte-for-byte.

**Do not weaken tolerances to force green.** A failure is the signal the gate
exists to produce — investigate the divergence.

## Adding a case

Drop a new deck dir under `<code>/cases/<name>/` whose `$GR_INPUT` file
(`config.sh`) is present. No harness change needed. Candidates already supported
by dropping a dir in: KIM `quadpack`, KiLCA `vacuum`, more QL-Balance modes.

- **Static decks** (KIM, KiLCA): just commit the input files (+ `profiles/`).
- **Generated decks** (QL-Balance): add a `gr_prepare.sh` to the case. It runs
  before the exe with `GR_BUILD_SRC` set to *that build's* worktree root, so the
  ref build's KAMELpy drives the ref KiLCA and likewise for cur.

The KiLCA `flre_m6n2` deck was generated once from KAMELpy's validated
`KiLCA_interface.write()` (mode `(6,2)`, AUG 3-zone stack flre→vacuum→vacuum)
and frozen here.

## Adding an input migration

If a PR renames a namelist key, the older `ref` binary won't understand the new
deck. Add `migrations/NNNN_slug.yaml` (+ document it in `MIGRATIONS.md`);
`gr_run_case.sh` applies them with `gr_migrate.py` (`GR_MIGRATE=up|down`, auto
`down` for a `legacy*` build). The first candidate is the QL-Balance
`data_verbosity` ⇄ `diagnostics_output` drift, currently handled inline in
`setup_runfolder.py`.

## Running locally

```bash
git tag -f golden-baseline origin/main          # pick a baseline
test/golden/run_golden.sh golden-baseline golden-baseline   # self-check: must be all-PASS
test/golden/run_golden.sh golden-baseline HEAD              # real comparison
```

The `cur==ref` self-check proving all-PASS is the guard against harness
nondeterminism (a FAIL there is a harness bug, not a physics finding). The local
dispatch backend is portable across GNU and BSD (macOS) `xargs`.

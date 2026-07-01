# Golden-Record CI for KiLCA, QL-Balance, and KIM — Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Vendor the `plasma/proj/golden` `gr-runner` harness into `KAMEL/test/golden/` and gate all three codes (KIM, QL-Balance, KiLCA) in CI so any PR that changes their physics output beyond floating-point tolerance turns the build red.

**Architecture:** A/B two-build comparison. CI builds a **reference** binary from the frozen `golden-baseline` tag (≡ current `main`) and a **candidate** binary from the PR `HEAD`, in the *same* CI environment so compiler/toolchain effects cancel. Each code runs a small synthetic input case against both builds; outputs are diffed at `rtol 1e-7 / atol 1e-12`. The harness is the code-agnostic `gr-runner` (vendored bash scripts + per-code `config.sh` + `cases/`), driven by one top-level orchestrator. Runs **only** in a dedicated GitHub Actions job — never in local `ctest`/`make test`.

**Tech Stack:** bash (`gr-runner`), Python 3 comparators (`gr_numcompare.py` numeric/recursive + `compare.py` HDF5), CMake + Ninja builds via `git worktree`, GitHub Actions, KAMELpy for QL-Balance synthetic input generation.

---

## Resolved decisions (from grilling — do not relitigate)

| Decision | Value |
|---|---|
| Harness | Vendor `gr-runner` into `KAMEL/test/golden/` (copy, **not** submodule — public CI cannot reach the private GitLab repo) |
| Scope | CI-only — dedicated GitHub Actions job; **not** registered in CTest / `make test` |
| Baseline | Frozen at current `main`, tagged `golden-baseline`; `ref` = tag, `cur` = PR HEAD; both rebuilt each run so env cancels |
| Corpus (3 cases) | KIM `electromagnetic` · QL-Balance `default` (TimeEvolution m=6,n=2) · KiLCA `flre_m6n2` |
| KiLCA outputs | Full output tree; recursive numeric compare + `.uff` byte-compare |
| Tolerance | `rtol 1e-7 / atol 1e-12` (gr-runner default; per-case override allowed) |
| Inputs | Synthetic + small committed text; **no** Git LFS, no `plasma/data` |
| Trigger | Every PR + push to `main`; one shared `ref` build + one shared `cur` build serve all cases |

**Accepted coverage gaps** (documented, not bugs): the netlib-QUADPACK→fortnum θ-integration path (KIM `quadpack` case) and the KiLCA `vacuum` case are **not** in the initial corpus. They can be added later by dropping a new case dir in — no harness change needed.

**Source of truth for vendored files:** the cloned suite at `/Users/markusmarkl/code/golden/` (`kim/` and `kamel/`). The `bin/` scripts are byte-identical across the two codes, so there is exactly one canonical copy to vendor.

---

## Architecture detail (read before starting)

### Directory layout to build
```
test/golden/
├── README.md
├── run_golden.sh                 # top-level orchestrator (NEW)
├── bin/                          # canonical gr-runner (vendored + enhanced)
│   ├── gr_build_all.sh           # NEW: build whole repo once per ref/cur, symlink all exes
│   ├── gr_build.sh               # vendored (single-target build; kept for cluster/manual use)
│   ├── gr_run_case.sh            # vendored + prepare-hook enhancement
│   ├── gr_dispatch.sh            # vendored
│   ├── gr_compare.sh             # vendored
│   ├── gr_migrate.py             # vendored
│   ├── gr_numcompare.py          # vendored + recursive/.uff enhancement
│   ├── compare.py                # vendored from existing QL-Balance test (HDF5 comparator)
│   └── setup_runfolder.py        # vendored from existing QL-Balance test (synthetic input gen)
├── kim/
│   ├── bin -> ../bin             # symlink (keeps gr-runner's ROOT=parent-of-bin contract)
│   ├── config.sh
│   └── cases/electromagnetic/{KIM_config.nml,n.dat,q.dat,Te.dat,Ti.dat}
├── qlbalance/
│   ├── bin -> ../bin
│   ├── config.sh
│   └── cases/default/{balance_conf.nml,setup_runfolder.py,gr_prepare.sh,tol.yaml}
└── kilca/
    ├── bin -> ../bin
    ├── config.sh
    └── cases/flre_m6n2/{background.in,antenna.in,modes.in,output.in,eigmode.in,
                          zone_flre.in,profiles/{q,n,Ti,Te,Vth,Vz,Er}.dat}
```

### How the symlink keeps gr-runner unmodified
Every `gr_*.sh` computes `ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"`. With `kim/bin` a symlink to `../bin`, invoking `kim/bin/gr_dispatch.sh` resolves `BASH_SOURCE` to `.../kim/bin/gr_dispatch.sh`, so `ROOT=.../kim` and `ROOT/config.sh` is the KIM config. One physical `bin/`, three logical roots. Symlinks are committed (Linux CI + macOS both follow them).

### Shared builds (the efficiency win)
Per-code `gr_build.sh` wipes and rebuilds one target — wasteful when three codes share a repo. `run_golden.sh` instead calls **`gr_build_all.sh ref <tag>`** and **`gr_build_all.sh cur HEAD`** once each: a single `git worktree` + `cmake` + full build producing `KIM.x`, `ql-balance.x`, and the KiLCA exe, then symlinks all three into `test/golden/build_ref/` and `build_cur/`. Each code dir gets `build_ref`/`build_cur` symlinks pointing at those shared dirs, so `gr_dispatch.sh --backend local` finds its exe with zero per-code rebuilds.

### QL-Balance is special (needs a prepare hook)
KIM and KiLCA cases are static decks: copy inputs, run the exe. QL-Balance's "input" is generated at run time by `setup_runfolder.py` (synthetic parabolic profiles **and** an internal KiLCA run via KAMELpy). So its case carries a `gr_prepare.sh` that `gr_run_case.sh` executes before the exe. Crucially the prepare step must use **the build's own** `python/` tree and binaries (ref's KAMELpy drives ref's KiLCA), so `gr_prepare.sh` reads `GR_BUILD_SRC` (the worktree root for this build) exported by `gr_run_case.sh`.

### Retire the old harness
`test/ql-balance/golden_record/` currently runs under CTest (so locally and in CI). It is **superseded** by `test/golden/qlbalance/`. Remove its `add_subdirectory(ql-balance)` from `test/CMakeLists.txt` so golden records leave the local `ctest` path entirely, satisfying "only in CI". `compare.py` and `setup_runfolder.py` are vendored into `test/golden/bin/` first, then the old dir is deleted.

---

## Phase 0 — Baseline tag

### Task 0.1: Create and push the `golden-baseline` tag

**Files:** none (git operation + note it in the plan/PR).

**Step 1:** Confirm the intended baseline commit is the tip of `origin/main` *after* the fortnum PRs (#140–#150) are merged. Fetch first:
```bash
git fetch origin main
git rev-parse origin/main
```
Expected: a SHA; record it in the PR description.

**Step 2:** Create an annotated tag on that commit:
```bash
git tag -a golden-baseline origin/main -m "Frozen physics baseline for golden-record CI (2026-06-30)"
git push origin golden-baseline
```
Expected: `* [new tag] golden-baseline -> golden-baseline`.

**Step 3 (verify):**
```bash
git ls-remote --tags origin | grep golden-baseline
```
Expected: one line with the tag and the baseline SHA.

> If the team prefers a moving baseline later, only this tag changes; no code changes.

---

## Phase 1 — Vendor the gr-runner skeleton

### Task 1.1: Create `test/golden/bin/` and vendor the byte-identical scripts

**Files:**
- Create: `test/golden/bin/{gr_build.sh,gr_run_case.sh,gr_dispatch.sh,gr_compare.sh,gr_migrate.py,gr_numcompare.py}`

**Step 1:** Copy verbatim from the cloned suite (use KIM's copies; they are identical to kamel's):
```bash
mkdir -p test/golden/bin
cp /Users/markusmarkl/code/golden/kim/bin/gr_build.sh        test/golden/bin/
cp /Users/markusmarkl/code/golden/kim/bin/gr_run_case.sh     test/golden/bin/
cp /Users/markusmarkl/code/golden/kim/bin/gr_dispatch.sh     test/golden/bin/
cp /Users/markusmarkl/code/golden/kim/bin/gr_compare.sh      test/golden/bin/
cp /Users/markusmarkl/code/golden/kim/bin/gr_migrate.py      test/golden/bin/
cp /Users/markusmarkl/code/golden/kim/bin/gr_numcompare.py   test/golden/bin/
chmod +x test/golden/bin/*.sh test/golden/bin/gr_migrate.py test/golden/bin/gr_numcompare.py
```
> Deliberately **not** vendored: `condor*.sub`, `gr_drain.sh`, `gr_compare3.sh`. Those are cluster/Condor and historical-3-version features, out of scope for CI.

**Step 2 (verify):**
```bash
bash -n test/golden/bin/gr_dispatch.sh && echo "syntax ok"
```
Expected: `syntax ok`.

**Step 3: Commit**
```bash
git add test/golden/bin
git commit -m "test(golden): vendor gr-runner scripts into KAMEL"
```

### Task 1.2: Vendor the QL-Balance comparator + input generator

**Files:**
- Create: `test/golden/bin/compare.py`, `test/golden/bin/setup_runfolder.py`

**Step 1:** Copy from the existing in-repo harness (identical to the private suite's copies):
```bash
cp test/ql-balance/golden_record/compare.py          test/golden/bin/compare.py
cp test/ql-balance/golden_record/setup_runfolder.py   test/golden/bin/setup_runfolder.py
```

**Step 2:** In `test/golden/bin/setup_runfolder.py`, fix the `python` path bootstrap. The original computes `parents[3]` for the repo `python/` dir from `test/ql-balance/golden_record/setup_runfolder.py`. The new location `test/golden/qlbalance/cases/default/setup_runfolder.py` (where it will be copied in Task 5.1) is one level deeper. **Leave `bin/setup_runfolder.py` importable by absolute path injection instead** — change the sys.path block to resolve `python/` from an env var the prepare hook sets:
```python
# Add the python package directory to sys.path for standalone execution
_python_dir = os.environ.get("GR_PYTHON_DIR")
if _python_dir is None:
    # Fallback: walk up to a directory containing 'python/balance_interface'
    p = Path(__file__).resolve()
    for cand in p.parents:
        if (cand / "python" / "balance_interface").is_dir():
            _python_dir = str(cand / "python")
            break
if _python_dir and _python_dir not in sys.path:
    sys.path.insert(0, _python_dir)
```
This makes the generator work regardless of nesting depth and lets the prepare hook point it at the *build's* python tree.

**Step 3 (verify):**
```bash
python3 -c "import ast; ast.parse(open('test/golden/bin/setup_runfolder.py').read()); print('parse ok')"
```
Expected: `parse ok`.

**Step 4: Commit**
```bash
git add test/golden/bin/compare.py test/golden/bin/setup_runfolder.py
git commit -m "test(golden): vendor QL-Balance HDF5 comparator and input generator"
```

---

## Phase 2 — Enhance the comparator (recursive + binary) — TDD

The vendored `gr_numcompare.py` only walks the **top level** of an output dir and has no `.uff` handling. KiLCA's outputs are nested (`linear-data/m_6_n_2_*/EB.dat`, `background-data/*.dat`) and include `formfactors.*.uff` (Fortran binary). Make the comparator recursive and binary-aware.

### Task 2.1: Write the failing test for recursive + binary comparison

**Files:**
- Create: `test/golden/bin/test_gr_numcompare.py`

**Step 1: Write the failing test**
```python
#!/usr/bin/env python3
"""Unit tests for the recursive numeric comparator."""
import os, subprocess, sys, struct
from pathlib import Path

HERE = Path(__file__).resolve().parent
CMP = HERE / "gr_numcompare.py"


def _run(a, b, *extra):
    r = subprocess.run([sys.executable, str(CMP), str(a), str(b), *extra],
                       capture_output=True, text=True)
    return r.returncode, r.stdout


def _write(p, text):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text)


def test_nested_numeric_match(tmp_path):
    for root in ("ref", "cur"):
        _write(tmp_path / root / "linear-data" / "m_6_n_2" / "EB.dat",
               "1.0 2.0\n3.0 4.0\n")
        _write(tmp_path / root / "background-data" / "b0z.dat", "5.0\n6.0\n")
    rc, out = _run(tmp_path / "ref", tmp_path / "cur")
    assert rc == 0, out
    assert "EB.dat" in out and "b0z.dat" in out          # nested files were found
    assert "FAIL" not in out


def test_nested_numeric_divergence_flags(tmp_path):
    _write(tmp_path / "ref" / "sub" / "EB.dat", "1.0\n")
    _write(tmp_path / "cur" / "sub" / "EB.dat", "1.1\n")   # 10% off >> rtol
    rc, out = _run(tmp_path / "ref", tmp_path / "cur")
    assert rc == 1, out
    assert "FAIL" in out


def test_uff_byte_identical(tmp_path):
    payload = struct.pack("<3d", 1.0, 2.0, 3.0)
    for root in ("ref", "cur"):
        p = tmp_path / root / "formfactors.flre.uff"
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(payload)
    rc, out = _run(tmp_path / "ref", tmp_path / "cur")
    assert rc == 0, out
    assert "formfactors.flre.uff" in out and "MATCH(bytes)" in out


def test_uff_byte_divergence_flags(tmp_path):
    for root, val in (("ref", 1.0), ("cur", 9.0)):
        p = tmp_path / root / "formfactors.flre.uff"
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(struct.pack("<d", val))
    rc, out = _run(tmp_path / "ref", tmp_path / "cur")
    assert rc == 1, out
    assert "DIFFER(bytes" in out


def test_volatile_files_skipped(tmp_path):
    for root in ("ref", "cur"):
        _write(tmp_path / root / "run.log", f"started at {root}\n")
        _write(tmp_path / root / "EB.dat", "1.0\n")
    rc, out = _run(tmp_path / "ref", tmp_path / "cur")
    assert rc == 0, out
    assert "run.log" not in out          # volatile logs excluded even when nested
```

**Step 2: Run it to verify it fails**
```bash
python3 -m pytest test/golden/bin/test_gr_numcompare.py -v
```
Expected: failures — current `gr_numcompare.py` only does `os.listdir` (top level), so `test_nested_*` find no files and `test_uff_*` are not handled.

### Task 2.2: Implement recursion + `.uff` byte handling

**Files:**
- Modify: `test/golden/bin/gr_numcompare.py`

**Step 1:** Replace the top-level `os.listdir` walk with a recursive relative-path walk; treat `.uff` (and any non-numeric/length-mismatched file) via byte comparison; keep the volatile-file skip list and extend it to match by basename at any depth.
```python
#!/usr/bin/env python3
# Generic numeric comparator for codes without a code-specific one.
# Usage: gr_numcompare.py <dir_A> <dir_B> [rtol] [atol]
# Recursively compares every regular file present in BOTH trees: whitespace-
# separated numeric columns where parseable, else byte-equality. Reports max
# relative diff per file. Exit 1 if any file diverges beyond rtol.
import sys, os

A, B = sys.argv[1], sys.argv[2]
rtol = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-7
atol = float(sys.argv[4]) if len(sys.argv) > 4 else 1e-12

SKIP = {"run.log", "exit_code.txt", "runtime_seconds.txt", "migrate.log"}


def nums(p):
    out = []
    try:
        for line in open(p, errors="ignore"):
            for tok in line.split():
                try:
                    out.append(float(tok))
                except ValueError:
                    pass
    except Exception:
        return None
    return out


def rel_files(root):
    found = []
    for dirpath, _dirs, files in os.walk(root):
        for f in files:
            if f in SKIP:
                continue
            full = os.path.join(dirpath, f)
            if os.path.isfile(full):
                found.append(os.path.relpath(full, root))
    return set(found)


fail = 0
checked = 0
if os.path.isdir(A) and os.path.isdir(B):
    common = sorted(rel_files(A) & rel_files(B))
else:
    common = []

for rel in common:
    pa, pb = os.path.join(A, rel), os.path.join(B, rel)
    na, nb = nums(pa), nums(pb)
    if na is None or nb is None or len(na) != len(nb) or not na:
        same = open(pa, "rb").read() == open(pb, "rb").read()
        checked += 1
        print(f"{rel}: {'MATCH(bytes)' if same else 'DIFFER(bytes/shape)'}")
        fail += 0 if same else 1
        continue
    mx = 0.0
    for x, y in zip(na, nb):
        d = abs(x - y)
        s = d / (abs(y) + atol)
        mx = max(mx, 0.0 if d <= atol else s)
    checked += 1
    ok = mx <= rtol
    print(f"{rel}: max_rel={mx:.3e} {'PASS' if ok else 'FAIL'}")
    fail += 0 if ok else 1

print(f"-- {checked} files compared, {fail} over rtol={rtol:.1e}")
sys.exit(1 if fail else 0)
```

**Step 2: Run the tests to verify they pass**
```bash
python3 -m pytest test/golden/bin/test_gr_numcompare.py -v
```
Expected: all 5 tests PASS.

**Step 3: Commit**
```bash
git add test/golden/bin/gr_numcompare.py test/golden/bin/test_gr_numcompare.py
git commit -m "test(golden): make gr_numcompare recursive and binary-aware"
```

---

## Phase 3 — KIM case (static deck, simplest)

### Task 3.1: Add the KIM `electromagnetic` case and config

**Files:**
- Create: `test/golden/kim/config.sh`
- Create: `test/golden/kim/cases/electromagnetic/{KIM_config.nml,n.dat,q.dat,Te.dat,Ti.dat}`
- Create symlink: `test/golden/kim/bin -> ../bin`

**Step 1:** Copy the proven case inputs from the cloned suite and create the bin symlink:
```bash
mkdir -p test/golden/kim/cases
cp -r /Users/markusmarkl/code/golden/kim/cases/electromagnetic test/golden/kim/cases/electromagnetic
ln -s ../bin test/golden/kim/bin
```

**Step 2:** Write `test/golden/kim/config.sh` (CI build of the whole repo target set is handled by the orchestrator; the per-code build command here is kept for manual/cluster use and points the build at the repo root):
```bash
# Per-code golden-record config (KIM).
GR_CODE=kim
GR_EXE=KIM.x
GR_INPUT=KIM_config.nml             # case auto-discovery key
GR_CODE_SRC=${GR_CODE_SRC:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}
GR_BUILD_CMD='cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release && cmake --build build -j --target KIM.x && cmake --install build'
GR_BUILD_EXE=build/install/bin/KIM.x
GR_COMPARE=gr_numcompare.py
```

**Step 3 (verify):** the case is discoverable and inputs are tiny text:
```bash
ls test/golden/kim/cases/electromagnetic
du -sh test/golden/kim/cases    # expect a few KB
```
Expected: `KIM_config.nml n.dat q.dat Te.dat Ti.dat`, size in KB.

**Step 4: Commit**
```bash
git add test/golden/kim
git commit -m "test(golden): add KIM electromagnetic case"
```

### Task 3.2: Local smoke — KIM runs and produces outputs

**Files:** none (uses current `build/install/bin/KIM.x`).

**Step 1:** Ensure a current build exists, then run the case by hand into a scratch dir:
```bash
cmake --build build --target KIM.x && cmake --install build
( tmp=$(mktemp -d); cp -rL test/golden/kim/cases/electromagnetic/. "$tmp"/; \
  cd "$tmp" && /Users/markusmarkl/code/KAMEL/build/install/bin/KIM.x > run.log 2>&1; \
  echo "rc=$?"; ls -R out 2>/dev/null | head )
```
Expected: `rc=0` and an `out/m-6_n2/` (or similar) directory containing `Er.dat`/`Er_no_Vpol.dat`.

> If `rc!=0`, read `run.log`; the EM case needs the four profile `.dat` files in cwd — confirm they were copied. Do **not** proceed until KIM runs clean.

---

## Phase 4 — KiLCA case (new — author from blueprints)

### Task 4.1: Author the KiLCA `flre_m6n2` deck

**Files:**
- Create: `test/golden/kilca/cases/flre_m6n2/{background.in,antenna.in,modes.in,output.in,eigmode.in,zone_flre.in}`
- Create: `test/golden/kilca/cases/flre_m6n2/profiles/{q,n,Ti,Te,Vth,Vz,Er}.dat`
- Create symlink: `test/golden/kilca/bin -> ../bin`

**Step 1:** Seed the deck from the KAMELpy blueprints (these are the canonical templates the Python interface fills in):
```bash
mkdir -p test/golden/kilca/cases/flre_m6n2/profiles
cp python/KiLCA_interface/blueprints/background.in  test/golden/kilca/cases/flre_m6n2/
cp python/KiLCA_interface/blueprints/antenna.in     test/golden/kilca/cases/flre_m6n2/
cp python/KiLCA_interface/blueprints/modes.in       test/golden/kilca/cases/flre_m6n2/
cp python/KiLCA_interface/blueprints/output.in      test/golden/kilca/cases/flre_m6n2/
cp python/KiLCA_interface/blueprints/eigmode.in     test/golden/kilca/cases/flre_m6n2/
cp python/KiLCA_interface/blueprints/zone_flre.in   test/golden/kilca/cases/flre_m6n2/
ln -s ../bin test/golden/kilca/bin
```

**Step 2:** Set `modes.in` to a single resonant mode `(6, 2)` and ensure `background.in` line "1 — background recalculated (7 input profiles needed)" stays `1`, so the seven profiles `q,n,Ti,Te,Vth,Vz,Er` in `profiles/` are required (and present). Author the seven profile files as small monotonic text columns covering the radial range in `background.in` (plasma radius 67 cm; grid a few cm out to ~80). The simplest robust source is to **generate them once with KAMELpy's parabolic helper** and freeze the result:
```bash
PYTHONPATH=python python3 - <<'PY'
from utility import create_parabolic_profiles_from_res_surf
create_parabolic_profiles_from_res_surf(
    "test/golden/kilca/cases/flre_m6n2/profiles",
    q0=3.0, n0=2e13, Te0=1000, Ti0=1000, Vz0=1e6, Er0=0.2, Vth0=1e5,
    mpol=6, ntor=2, rmin=3.0, rmax=80.0, num=300, a=67.0)
PY
ls test/golden/kilca/cases/flre_m6n2/profiles
```
Expected: the profile `.dat` files written. Trim to exactly the seven KiLCA needs (`q,n,Ti,Te,Vth,Vz,Er`); delete any extras the helper emits that KiLCA's `background.in` recalculation path does not read.

> Authoring note: KiLCA's `background.in` is **positional** (value before `#comment`), not a namelist. Edit values in place; never reorder lines. Keep `antenna.in` `flab = 1.0 0.0` and `nmod = 1`.

**Step 3 (verify deck loads):** dry-run against the current build (Task 4.2 covers full output); here just confirm KiLCA starts and reads the deck without an input error.

**Step 4: Commit**
```bash
git add test/golden/kilca/cases
git commit -m "test(golden): author KiLCA flre m6n2 case deck and profiles"
```

### Task 4.2: Run KiLCA standalone and confirm the full output tree

**Files:** none (uses current build).

**Step 1:** Find the version-suffixed KiLCA exe (glob it — the name carries a version string):
```bash
cmake --build build --target kilca_normal_exe && cmake --install build
KILCA=$(ls build/install/bin/KiLCA_Normal_*_64bit | head -1); echo "$KILCA"
```
Expected: a path like `build/install/bin/KiLCA_Normal_V_2.4.2_MDNO_FPGEN_POLYNOMIAL_Release_64bit`.

**Step 2:** Run the case in a scratch dir and inspect outputs:
```bash
( tmp=$(mktemp -d); cp -rL test/golden/kilca/cases/flre_m6n2/. "$tmp"/; \
  cd "$tmp" && "$KILCA" > run.log 2>&1; echo "rc=$?"; \
  find . -name 'EB.dat' -o -name 'formfactors*.uff' -o -path '*background-data*' | head )
```
Expected: `rc=0`; presence of `linear-data/m_6_n_2_*/EB.dat`, `background-data/*.dat`, and `formfactors.*.uff`.

> If `rc!=0`: read `run.log`. Most likely a profile range/units mismatch — adjust the generated profiles to cover the radial grid in `background.in`. Iterate Task 4.1 ↔ 4.2 until KiLCA runs clean. This is the only case with real authoring risk; budget time here.

### Task 4.3: Add the KiLCA `config.sh`

**Files:**
- Create: `test/golden/kilca/config.sh`

**Step 1:**
```bash
# Per-code golden-record config (KiLCA).
GR_CODE=kilca
GR_EXE=KiLCA_Normal.x                # placeholder; resolved to the version-suffixed exe by gr_build_all symlink
GR_INPUT=background.in               # case auto-discovery key
GR_CODE_SRC=${GR_CODE_SRC:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}
GR_BUILD_CMD='cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release && cmake --build build -j --target kilca_normal_exe && cmake --install build'
GR_BUILD_EXE=build/install/bin/KiLCA_Normal.x
GR_COMPARE=gr_numcompare.py
```
> `gr_build_all.sh` (Task 6.1) creates a stable `KiLCA_Normal.x` symlink → the version-suffixed binary in each shared build dir, so `GR_EXE` can be a fixed name regardless of the embedded version string.

**Step 2: Commit**
```bash
git add test/golden/kilca/config.sh
git commit -m "test(golden): add KiLCA config.sh"
```

---

## Phase 5 — QL-Balance case (relocate the proven flow + prepare hook)

### Task 5.1: Add the QL-Balance case dir and prepare hook

**Files:**
- Create: `test/golden/qlbalance/config.sh`
- Create: `test/golden/qlbalance/cases/default/setup_runfolder.py` (copy of vendored one)
- Create: `test/golden/qlbalance/cases/default/balance_conf.nml` (template, copied from the build at run time too)
- Create: `test/golden/qlbalance/cases/default/gr_prepare.sh`
- Create: `test/golden/qlbalance/cases/default/tol.yaml` (documents per-quantity tolerance; informational)
- Create symlink: `test/golden/qlbalance/bin -> ../bin`

**Step 1:** Scaffold:
```bash
mkdir -p test/golden/qlbalance/cases/default
ln -s ../bin test/golden/qlbalance/bin
cp test/golden/bin/setup_runfolder.py test/golden/qlbalance/cases/default/setup_runfolder.py
cp QL-Balance/namelists/balance_conf.nml test/golden/qlbalance/cases/default/balance_conf.nml
```

**Step 2:** Write `test/golden/qlbalance/cases/default/gr_prepare.sh`. It generates the synthetic runfolder using **the current build's** python tree + binaries (exported as `GR_BUILD_SRC` by `gr_run_case.sh`), so the ref build drives ref's KiLCA and the cur build drives cur's KiLCA:
```bash
#!/bin/bash
# Prepare a QL-Balance synthetic runfolder in the current dir ($PWD == out dir).
# gr_run_case.sh exports GR_BUILD_SRC (the worktree root for this build).
set -eu
SRC="${GR_BUILD_SRC:?gr_prepare needs GR_BUILD_SRC}"
export GR_PYTHON_DIR="$SRC/python"
export PYTHONPATH="$SRC/python:${PYTHONPATH:-}"
# Generate balance_conf.nml + profiles + run KiLCA (flre+vacuum) into $PWD.
python3 "$SRC/test/golden/qlbalance/cases/default/setup_runfolder.py" "$PWD" \
        --nml-template "$PWD/balance_conf.nml"
```
> `setup_runfolder.py`'s `setup_runfolder(run_path, nml_template="")` signature already accepts a template path; wrap its `__main__` to accept `<run_path> --nml-template <path>` (add a tiny argparse block if missing). KAMELpy resolves binaries via `build/install/bin` under `GR_BUILD_SRC`, so each build stays self-consistent.

**Step 3:** Write `test/golden/qlbalance/config.sh`:
```bash
# Per-code golden-record config (QL-Balance).
GR_CODE=qlbalance
GR_EXE=ql-balance.x
GR_INPUT=balance_conf.nml            # case auto-discovery key
GR_CODE_SRC=${GR_CODE_SRC:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}
GR_BUILD_CMD='cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release && cmake --build build -j --target ql-balance.x && cmake --install build'
GR_BUILD_EXE=build/install/bin/ql-balance.x
GR_COMPARE=compare_qlbalance.sh      # wrapper: HDF5-aware compare of the two output files
```

**Step 4:** QL-Balance output is one HDF5 file, not a dir tree, and `compare.py` takes two file paths. Add a tiny wrapper `test/golden/bin/compare_qlbalance.sh` that locates `out/<...>.hdf5` under each side and calls `compare.py`:
```bash
#!/bin/bash
# Usage: compare_qlbalance.sh <out/ref> <out/cur> [cases...]
set -u
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REF="$1"; CUR="$2"
rc=0
for casedir in "$REF"/*/; do
  c=$(basename "$casedir")
  g=$(ls "$REF/$c/out/"*.hdf5 2>/dev/null | head -1)
  a=$(ls "$CUR/$c/out/"*.hdf5 2>/dev/null | head -1)
  if [ -z "$g" ] || [ -z "$a" ]; then echo "$c: missing hdf5 (ref=$g cur=$a)"; rc=1; continue; fi
  python3 "$ROOT/bin/compare.py" "$g" "$a" || rc=1
done
exit $rc
```
> `compare.py`'s `main()` already takes `<golden.h5> <actual.h5>` and exits nonzero on mismatch, comparing the `QUANTITIES_TO_COMPARE` list at `rtol=1e-8`. Keep that list as-is; it is the QL-Balance per-quantity spec.

**Step 5 (verify):**
```bash
bash -n test/golden/qlbalance/cases/default/gr_prepare.sh && bash -n test/golden/bin/compare_qlbalance.sh && echo ok
chmod +x test/golden/qlbalance/cases/default/gr_prepare.sh test/golden/bin/compare_qlbalance.sh
```
Expected: `ok`.

**Step 6: Commit**
```bash
git add test/golden/qlbalance test/golden/bin/compare_qlbalance.sh
git commit -m "test(golden): add QL-Balance case with build-local prepare hook"
```

### Task 5.2: Teach `gr_run_case.sh` to run a prepare hook

**Files:**
- Modify: `test/golden/bin/gr_run_case.sh`

**Step 1:** After the input copy + migration block and before the timed `"$EXE"` run, insert a prepare-hook step that exports the build's worktree root so KAMELpy uses the matching python + binaries:
```bash
# Optional per-case prepare hook (e.g. QL-Balance synthetic runfolder generation).
# Runs before the executable, with GR_BUILD_SRC pointing at this build's worktree root.
if [ -f "$OUT_DIR/gr_prepare.sh" ]; then
  GR_BUILD_SRC="$(cd "$BUILD_DIR" && pwd)/src" \
    bash "$OUT_DIR/gr_prepare.sh" > prepare.log 2>&1 \
    || { echo "gr_prepare failed (see prepare.log)" >&2; echo 97 > exit_code.txt; exit 97; }
fi
```
> `gr_build_all.sh` lays out each shared build dir as `build_<tag>/src` (the worktree) + exe symlinks, so `$BUILD_DIR/src` is the worktree root. Add `prepare.log` to the comparator `SKIP` set (already covered by `migrate.log`/`run.log`; add `prepare.log`).

**Step 2:** Add `prepare.log` to the `SKIP` set in `test/golden/bin/gr_numcompare.py`:
```python
SKIP = {"run.log", "exit_code.txt", "runtime_seconds.txt", "migrate.log", "prepare.log"}
```

**Step 3 (verify):**
```bash
bash -n test/golden/bin/gr_run_case.sh && echo ok
python3 -m pytest test/golden/bin/test_gr_numcompare.py -q
```
Expected: `ok`; tests still green.

**Step 4: Commit**
```bash
git add test/golden/bin/gr_run_case.sh test/golden/bin/gr_numcompare.py
git commit -m "test(golden): add prepare-hook support to gr_run_case"
```

---

## Phase 6 — Orchestrator (shared builds + per-code dispatch)

### Task 6.1: `gr_build_all.sh` — one whole-repo build per ref/cur

**Files:**
- Create: `test/golden/bin/gr_build_all.sh`

**Step 1:**
```bash
#!/bin/bash
# Build the WHOLE KAMEL repo at a git ref once, into a shared build dir with
# stable exe symlinks for all three codes.
# Usage: gr_build_all.sh <ref|cur> <git-ref>
set -eu
TAG="$1"; GITREF="$2"
GOLDEN_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"   # test/golden
SRC_REPO="$(cd "$GOLDEN_ROOT/../.." && pwd)"                     # repo root
WT="$GOLDEN_ROOT/build_$TAG"
rm -rf "$WT"; mkdir -p "$WT"
git -C "$SRC_REPO" worktree prune
git -C "$SRC_REPO" worktree add -f "$WT/src" "$GITREF"
( cd "$WT/src" \
  && cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release \
  && cmake --build build -j --target KIM.x ql-balance.x kilca_normal_exe \
  && cmake --install build )
BIN="$WT/src/build/install/bin"
ln -sf "$BIN/KIM.x"            "$WT/KIM.x"
ln -sf "$BIN/ql-balance.x"     "$WT/ql-balance.x"
ln -sf "$(ls "$BIN"/KiLCA_Normal_*_64bit | head -1)" "$WT/KiLCA_Normal.x"
echo "built $TAG ($GITREF) -> $WT  (KIM.x, ql-balance.x, KiLCA_Normal.x)"
```

**Step 2 (verify):** syntax + that targets exist in CMake:
```bash
bash -n test/golden/bin/gr_build_all.sh && echo ok
rg -n "kilca_normal_exe|ql-balance.x|KIM.x" KiLCA/CMakeLists.txt QL-Balance/CMakeLists.txt KIM/CMakeLists.txt | head
```
Expected: `ok`; targets present.

**Step 3: Commit**
```bash
chmod +x test/golden/bin/gr_build_all.sh
git add test/golden/bin/gr_build_all.sh
git commit -m "test(golden): add shared whole-repo build for ref/cur"
```

### Task 6.2: `run_golden.sh` — top-level driver

**Files:**
- Create: `test/golden/run_golden.sh`

**Step 1:**
```bash
#!/bin/bash
# Top-level golden-record driver: build ref+cur once, run every code's cases
# against both, compare. Exit nonzero if any code diverges.
# Usage: run_golden.sh [REF_GITREF] [CUR_GITREF]
#   defaults: REF=golden-baseline  CUR=HEAD
set -eu
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$ROOT"
REF_REF="${1:-golden-baseline}"
CUR_REF="${2:-HEAD}"
CODES="kim qlbalance kilca"

echo "== building ref ($REF_REF) and cur ($CUR_REF) =="
bin/gr_build_all.sh ref "$REF_REF"
bin/gr_build_all.sh cur "$CUR_REF"

overall=0
for code in $CODES; do
  echo "== $code =="
  # point this code's build_ref/build_cur at the shared builds
  ln -sfn "$ROOT/build_ref" "$ROOT/$code/build_ref"
  ln -sfn "$ROOT/build_cur" "$ROOT/$code/build_cur"
  ( cd "$ROOT/$code" && bin/gr_dispatch.sh --backend local --local-jobs 2 ) || overall=1
  ( cd "$ROOT/$code" && bin/gr_compare.sh ) || overall=1
done

echo "== golden-record overall: $([ $overall -eq 0 ] && echo PASS || echo FAIL) =="
exit $overall
```
> `gr_dispatch.sh` already resolves `build_ref`/`build_cur` relative to the code ROOT and calls `gr_run_case.sh`, which resolves `$BUILD_DIR/$GR_EXE`. With the shared-build symlinks and the stable exe names, every code finds its binary.

**Step 2 (verify syntax):**
```bash
bash -n test/golden/run_golden.sh && echo ok
chmod +x test/golden/run_golden.sh
```
Expected: `ok`.

**Step 3: Commit**
```bash
git add test/golden/run_golden.sh
git commit -m "test(golden): add top-level run_golden orchestrator"
```

### Task 6.3: Local end-to-end self-check (cur == ref must be all-PASS)

**Files:** none (validation only).

**Step 1:** Prove the harness is sound: build the *same* ref on both sides; every case must compare byte/FP-identical.
```bash
git fetch origin && git tag -f golden-baseline origin/main 2>/dev/null || true
test/golden/run_golden.sh golden-baseline golden-baseline
echo "exit=$?"
```
Expected: `golden-record overall: PASS`, `exit=0`. Any FAIL here is a harness bug (nondeterminism, volatile file not skipped), not a physics finding — fix before trusting the gate.

**Step 2:** Now run the real comparison on the working branch:
```bash
test/golden/run_golden.sh golden-baseline HEAD; echo "exit=$?"
```
Expected on a physics-preserving branch: `PASS`. If a case FAILs, inspect the `gr_compare` output — that is a genuine "this branch changed physics" signal (the whole point).

> Record both outcomes in the PR description. Do not weaken tolerances to force green; investigate divergences.

---

## Phase 7 — CI wiring (and remove golden from local ctest)

### Task 7.1: Add the dedicated GitHub Actions golden job

**Files:**
- Create: `.github/workflows/golden.yml`

**Step 1:**
```yaml
name: golden-record

on:
  push:
    branches: [main]
  pull_request:
    branches: [main]

jobs:
  golden:
    runs-on: ubuntu-latest
    steps:
      - name: Checkout (full history + tags)
        uses: actions/checkout@v4
        with:
          fetch-depth: 0          # need golden-baseline tag + worktree of it

      - name: Install dependencies
        run: |
          sudo apt-get update
          sudo apt-get install -y cmake ninja-build gfortran \
            libopenmpi-dev openmpi-bin libhdf5-openmpi-dev libnetcdf-dev \
            libnetcdff-dev libomp-dev libopenblas-dev liblapack-dev \
            libgsl-dev libsuitesparse-dev libsuperlu-dev python3-dev python3-pip \
            libmpfr-dev libgmp-dev
          pip install numpy scipy h5py f90nml

      - name: Cache external dependencies
        uses: actions/cache@v4
        with:
          path: |
            test/golden/build_ref/src/build/external-install
            test/golden/build_cur/src/build/external-install
          key: golden-ext-${{ runner.os }}-${{ hashFiles('cmake/Fetch*.cmake') }}
          restore-keys: golden-ext-${{ runner.os }}-

      - name: Run golden-record suite (ref=golden-baseline, cur=HEAD)
        run: test/golden/run_golden.sh golden-baseline ${{ github.sha }}
```
> `cur` is pinned to `github.sha` (the PR head commit) rather than literal `HEAD` so the worktree is deterministic in CI's detached-HEAD checkout. `ref` is the tag; `fetch-depth: 0` guarantees both the tag and full history are present for `git worktree add`.

**Step 2 (verify):** lint the YAML locally:
```bash
python3 -c "import yaml,sys; yaml.safe_load(open('.github/workflows/golden.yml')); print('yaml ok')"
```
Expected: `yaml ok`.

**Step 3: Commit**
```bash
git add .github/workflows/golden.yml
git commit -m "ci(golden): add dedicated golden-record job (ref vs PR head)"
```

### Task 7.2: Remove the superseded golden test from local ctest

**Files:**
- Modify: `test/CMakeLists.txt`
- Delete: `test/ql-balance/` (after confirming `compare.py`/`setup_runfolder.py` are vendored under `test/golden/bin/`)

**Step 1:** Remove the `add_subdirectory(ql-balance)` call so the golden test no longer runs under `ctest`/`make test`:
```bash
rg -n "add_subdirectory\(ql-balance\)" test/CMakeLists.txt
```
Edit `test/CMakeLists.txt` to drop that line (keep the `f90nml` discovery; if `ql-balance` was its only consumer, the file may become a no-op stub with a comment that golden records now live in `test/golden/` and run only in CI).

**Step 2:** Confirm nothing else references the old dir, then delete it:
```bash
rg -n "ql-balance/golden_record" --glob '!docs/**' || echo "no refs"
git rm -r test/ql-balance
```
Expected: `no refs` before deletion.

**Step 3 (verify):** the main build's tests no longer include any golden test:
```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release >/dev/null
ctest --test-dir build -N | rg -i golden || echo "no golden tests in ctest (correct)"
```
Expected: `no golden tests in ctest (correct)`.

**Step 4: Commit**
```bash
git add test/CMakeLists.txt
git commit -m "test(golden): retire ctest QL-Balance golden; CI-only via test/golden"
```

---

## Phase 8 — Docs

### Task 8.1: Write `test/golden/README.md` and update repo docs

**Files:**
- Create: `test/golden/README.md`
- Modify: `AGENTS.md` (Testing section), `Data.md` (optional pointer)

**Step 1:** `test/golden/README.md` covers: purpose (lock physics across backend swaps), the A/B model, `golden-baseline` tag semantics + how to re-bless it, how to add a case (drop a deck dir whose `$GR_INPUT` file exists), how to add an input migration (`migrations/NNNN_slug.yaml` + `MIGRATIONS.md`), tolerance + per-case override, and the explicit note that it runs **only in CI**, not `ctest`. Mirror the conventions from `/Users/markusmarkl/code/golden/{kim,kamel}/README.md`.

**Step 2:** In `AGENTS.md` Testing section, add a short paragraph: golden-record regression for KIM/QL-Balance/KiLCA lives in `test/golden/`, runs in the `golden-record` GitHub Actions job, and is intentionally excluded from local `ctest`.

**Step 3: Commit**
```bash
git add test/golden/README.md AGENTS.md Data.md
git commit -m "docs(golden): document the golden-record harness and CI gate"
```

---

## Open follow-ups (out of scope, note in PR)
- Add KIM `quadpack` + `default` and KiLCA `vacuum` cases to widen backend coverage (one dir each; no harness change).
- Optional `migrations/` entries once a PR renames a namelist key (the `data_verbosity`/`diagnostics_output` drift is the first candidate).
- Optional real-equilibrium cases via Git LFS if synthetic decks miss a path; keep those cluster-side if they get heavy.
- Reconcile tolerances long-term: QL-Balance still compares at `rtol=1e-8` inside `compare.py` while KIM/KiLCA use `rtol=1e-7`; unify if desired.

## Risk notes
- **Double build cost in CI:** two full KAMEL builds per run. Mitigated by the external-deps cache; if wall-clock is unacceptable, add `ccache` or move to "Both: per-PR + nightly" later.
- **QL-Balance determinism:** it runs KiLCA internally and does adaptive time-stepping; the existing `compare.py` quantity list is already curated to stable outputs. The Phase 6.3 self-check (cur==ref ⇒ PASS) is the guard against any residual nondeterminism.
- **KiLCA deck authoring (Task 4.1/4.2):** the only step with real iteration risk; the profiles must cover the radial grid. Budget time and keep the deck minimal.
```

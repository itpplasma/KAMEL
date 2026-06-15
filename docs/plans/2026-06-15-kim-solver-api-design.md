# KIM Solver API — Design

- **Date:** 2026-06-15
- **Status:** Proposed (design approved; not yet implemented)
- **Scope:** In-process Fortran API for KIM. Decouple QL-Balance from KIM internals and make KIM solves invokable in tests. Python bindings and a full public/executable-wide API are explicitly out of scope for this design.

## Motivation

KIM has three interfaces today, but **no real programmatic API**:

| Interface | What it is | Consumer |
|---|---|---|
| CLI / file | `KIM.x` reads namelists, writes HDF5/dat | standalone runs |
| Python (file) | `KIMpy` parses KIM HDF5 output | analysis / GUI |
| In-process | `libKIM.a` + `kim_t` (`init`/`run`) + factory | QL-Balance |

The `kim_t` abstract type *looks* like an API, but `init()`/`run()` take **no arguments** — every input and output flows through global module state (`config_m`, `setup_m`, `grid_m`, `species_m%plasma`, `equilibrium_m`, `fields_m%EBdat`).

The evidence that this is a missing seam is `QL-Balance/src/base/kim_wave_code_adapter.f90` — **~940 lines** whose only job is to bridge the gap. It imports from 8 KIM modules and **renames 26 symbols** to dodge name clashes, sets KIM globals by hand, calls `kim_init` → factory → `init` → `run`, reaches into internals to read results (`EBdat%*`, `plasma%*`, `equilibrium_m%*`), **re-implements KIM's own per-mode orchestration** (`calculate_equil` / `set_plasma_quantities` / `interpolate_equil`), and **manually resets KIM's global state** (`deallocate_EBdat`, `deallocate_equilibrium_arrays`).

This implicit, global-state contract is also why KIM line coverage is ~15%: almost nothing can be invoked without standing up all the globals.

## Decisions

1. **Scope:** an in-process Fortran call surface that decouples QL-Balance and makes solves testable. Not a Python binding; not a full public API; not (initially) an executable rewrite.
2. **State model:** a **state-owning handle** `kim_solver_t` is the target — state lives in the handle, not in module globals.
3. **Migration:** **interface-first**. Introduce the handle now, delegating to today's globals + factory (behaviour unchanged); migrate consumers onto it; then absorb globals into the type module-by-module behind the stable interface, each step guarded by tests.

## Section 1 — Architecture & the seam

The API is one new module, `kim_solver_m`, exposing two types: a stateful handle `kim_solver_t` and a plain-data `kim_results_t`. Nothing else in KIM becomes public; this module *is* the seam.

`kim_solver_t` **composes** rather than replaces the existing design. It owns the inputs and run state and internally holds a `class(kim_t)` run-type instance — so the current Strategy/Factory (`from_kim_factory_get_kim` selecting electromagnetic / electrostatic / WKB / benchmark) stays intact behind the handle.

The handle **owns the orchestration that currently leaks into the QL-Balance adapter**: mode selection, the per-mode equilibrium recompute, state reset, and reading results out. All of it moves behind `solve(m,n)`.

In phase 1 (interface-first) the handle is a thin delegator: `init` calls today's `kim_init` + factory + run-type `init`; `solve` does mode setup + `run`; `results()` copies the relevant globals into a `kim_results_t`. Behaviour is identical to today; the *contract* is now explicit and globals are hidden from callers. Later phases absorb those globals into `kim_solver_t` components without changing this public surface.

## Section 2 — The handle type & lifecycle

`kim_solver_t` carries: the run-type strategy (`class(kim_t), allocatable :: run_type`), a setup flag, a status code, and the last solve's results. (Phase 1: heavy state still lives in globals; later phases move grids/plasma/fields/config into the type here, behind these same procedures.)

```fortran
call kim%init(config_path, run_type='electromagnetic', profiles=p, stat=ierr)
call kim%solve(m=2, n=1, stat=ierr)
res = kim%results()
! ... loop over modes, or update profiles and re-solve ...
call kim%finalize()
```

- **`init(config_path, [run_type], [profiles], [stat])`** — reads config (delegates to `kim_init`), injects in-memory `profiles` if supplied (else KIM reads its own files), then sets up: factory-selects the run-type and constructs grids + first equilibrium. `run_type` optionally overrides the namelist `type_of_run` (the adapter's `kim_type_of_run='electromagnetic'` move).
- **`set_profiles(r,n,Te,Ti,q,Er,[stat])`** — update profiles between solves (time-evolution); flags equilibrium for recompute.
- **`solve(m,n,[stat])`** — sets the mode, ensures equilibrium is consistent for `(m,n)` and current profiles (the `calculate_equil` / `set_plasma_quantities` / `interpolate_equil` sequence the adapter hand-codes), runs, and stores results internally.
- **`results()`** — returns the last solve's `kim_results_t`.
- **`finalize()`** — releases state (replaces the adapter's manual `deallocate_EBdat` / `deallocate_equilibrium_arrays`) and HDF5 deinit.

Ordering and state-reset become the handle's responsibility, not the caller's.

## Section 3 — The results contract (`kim_results_t`)

Plain data — a deep copy, so it stays valid after the next `solve` or `finalize`. It carries KIM output on KIM's *native* grids; the caller interpolates onto its own grid (KIM never learns about QL-Balance's grid — that interpolation correctly stays caller-side).

Two grids, because the field solution and the background quantities live on different ones:

```fortran
type :: kim_results_t
    integer :: m, n                     ! mode these results are for

    ! field solution -- on the field grid (today: EBdat%r_grid)
    real(dp),    allocatable :: r_field(:)
    complex(dp), allocatable :: Es(:), Ep(:), Er(:), Etheta(:), Ez(:), Br(:)
    complex(dp), allocatable :: jpar(:), jpar_e(:), jpar_i(:)   ! may be unallocated

    ! derived background -- on the plasma grid (today: plasma%r_grid)
    real(dp), allocatable :: r_plasma(:)
    real(dp), allocatable :: kp(:), ks(:), om_E(:)     ! mode-dependent
    real(dp), allocatable :: nu_e(:), nu_i(:)
    real(dp), allocatable :: B0(:), B0z(:), B0th(:)    ! mode-independent
    integer :: stat                                    ! 0 = ok
end type
```

This is exactly the set the adapter extracts from `EBdat%*`, `plasma%*`, and `equilibrium_m%*` today, promoted into one documented type. Currents are `allocatable` because only some run-types fill them (the adapter already guards with `if (allocated(...))`).

The contract also **documents the conventions** the adapter currently carries as comments: coordinate systems (`Es`/`Ep` field-aligned rsp vs `Er`/`Etheta`/`Ez` cylindrical), sign/unit conventions (e.g. `dPhi0 = -Er`), and which fields a given run-type populates. These move from adapter comments into the type's doc-comments, so a caller reads the contract in one place.

## Section 4 — Migration plan & consumer mapping

Interface-first; each phase committed and green before the next.

**Phase 1 — introduce the seam, zero behaviour change.** Add `kim_solver_m`; implement purely by delegation to today's `kim_init` + factory + `run`, with `results()` copying globals out and `finalize()` doing the existing deallocate sequence. Add an **end-to-end solve test**: small electromagnetic bench config → `init` → `solve(m,n)` → assert results (grids monotonic, `Br` finite, a benchmark value if available). This is the first guard on the run path — it exercises `kernel.f90`, `fields_mod`, `electromagnetic_solver` (all 0% today), doubling as the run-path coverage win.

**Phase 2 — migrate QL-Balance.** Replace the adapter body: `kim_initialize` → `kim%init(config, run_type='electromagnetic', profiles=…)`; the mode loop → `kim%solve(m,n); res = kim%results()` + the existing interpolation. Delete the hand-coded `calculate_equil` / `set_plasma_quantities` / `interpolate_equil`, the manual `deallocate_*`, and the 26 renamed global imports. The adapter collapses from ~940 lines to a thin translation; `test_kim_adapter` guards it.

**Phase 3 (optional) — the executable.** `kim_main` becomes `init(config_path) → solve(m_mode,n_mode) → finalize`. HDF5 output unchanged. Cheap unification; deferrable.

**Phase 4…n — absorb globals into the type**, one module at a time (fields → species/plasma → equilibrium → grid → config), each behind the now-stable interface and guarded by the Phase-1 test, the pure-utility unit tests, and the adapter test.

## Section 5 — Error handling & testing

**Error handling.** KIM today halts with scattered `stop 1`. For an embeddable API that is hostile — it kills the host process. The seam's contract is: **report status, don't halt.**

- Every handle procedure takes an optional `stat` out-arg; `kim_results_t%stat` carries per-solve status; a small documented set of codes (`KIM_OK`, `KIM_NOT_SETUP`, `KIM_BAD_CONFIG`, `KIM_GRID_COVERAGE`, `KIM_SOLVE_FAILED`). Optional human message via the existing `logger_m`.
- The handle returns status for everything *it* checks — bad config path, `solve` before `init`, the grid-coverage check the adapter currently does with a `WARNING`/`stop`.
- Honest debt: deep `stop 1`s inside the pipeline cannot all go at once. They are converted to status returns incrementally during Phase 4, as ownership of each module's state moves into the handle. The contract promises clean status *at the seam* from day one.

**Testing.**

- **Keystone:** the Phase-1 end-to-end solve test (bench config → `init`/`solve`/`results`), which is also the run-path coverage win.
- **Golden value:** lock a numeric result from the bench solve — the regression guard that makes the Phase-4 global-absorption refactors safe.
- **Error paths:** `init` with a bad path → nonzero `stat` (no crash); `solve` before `init` → `KIM_NOT_SETUP`.
- **Per-run-type:** electrostatic / WKB / benchmark each get a solve test as their setup becomes parameterizable through the handle.
- The existing pure-utility unit tests and `test_kim_adapter` stay; the API test complements them.

## Out of scope / future work

- **Python binding** (`iso_c_binding` / C-ABI layer for KIMpy to call the solver instead of parsing HDF5).
- **Full public API** with the executable and all consumers sitting on one documented contract.
- **A config derived type** replacing the namelist-path input (phase 1 keeps the existing file-based config machinery).

## References

- `KIM/src/kim_main.f90` — current entry point (factory + `kim_t`).
- `KIM/src/general/KIM_base.f90` — `kim_t` abstract type.
- `KIM/src/general/KIM_mod.f90` — `from_kim_factory_get_kim`.
- `QL-Balance/src/base/kim_wave_code_adapter.f90` — the ~940-line adapter this design replaces.
- Related: KIM test-coverage & deduplication work (branch `improve/kim-coverage-dedup-impl`), which provides the test net the migration relies on.

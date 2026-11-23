# Developer documentation

User documentation lives next to the code: the top-level [README](../README.md),
the component READMEs ([KiLCA](../KiLCA/README.md), [KIM](../KIM/README.md),
[QL-Balance](../QL-Balance/README.md)), the
[KIM namelist reference](../KIM/nmls/README.md), and [Data.md](../Data.md) for
file formats and units.

This directory holds dated records of development work. Plans describe intended
changes when they were written; the code is authoritative where they differ.

## Plans

| Plan | Status |
| --- | --- |
| [KiLCA C++ to Fortran translation](plans/2026-02-16-kilca-cpp-to-fortran-design.md) | Implemented (#161) |
| [Public release](plans/2026-03-04-public-release-design.md) | Implemented (`v1.0.0`); citation metadata still missing |
| [Golden-record CI](plans/2026-06-30-golden-record-ci.md) | Implemented ([`test/golden/`](../test/golden/README.md)) |
| [KIM forced periodicity: design](plans/2026-07-10-kim-forced-periodicity-design.md) | Implemented (`type_of_run = 'electrostatic_periodic'`) |
| [Forced periodicity, phase 1](plans/2026-07-10-kim-forced-periodicity-impl-phase1.md) | Implemented |
| [Forced periodicity, phase 2](plans/2026-07-10-kim-forced-periodicity-impl-phase2.md) | Implemented |
| [Forced periodicity, phase 3](plans/2026-07-10-kim-forced-periodicity-impl-phase3.md) | Implemented |
| [KIM species conservation models](plans/2026-07-24-kim-species-conservation-models.md) | Implemented (`*_ifunc_conservation_model`) |
| [KiLCA Lorentz-torque verification](plans/2026-07-27-kilca-lorentz-torque-verification.md) | Implemented (`KiLCA/tests/test_lorentz_torque.f90`) |
| [Native Fortran interfaces](plans/2026-09-07-native-fortran-interfaces.md) | Not recorded |
| [PR #161 migration repairs](plans/2026-09-07-pr161-migration-repairs.md) | Not recorded |

New plans should state their status near the top and update it when the work
lands or is abandoned.

## Benchmarks

- [KIM #291: drift-kinetic limit](benchmarks/2026-08-12-kim-291-drift-kinetic-limit.md)
- [KIM #292: target-current normalization](benchmarks/2026-08-12-kim-292-target-current-normalization.md)
- [KIM #293: multimode feedback](benchmarks/2026-08-12-kim-293-multimode-feedback.md)

## Verification

- [KiLCA FLR](verification/kilca-flr.md)
- [PR #161 Fortran migration](verification/pr161-fortran-migration.md)
- [FLRE golden-record investigation](flre-golden-record-investigation.md) (resolved)

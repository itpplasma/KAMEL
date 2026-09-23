# KIM AUG/BALANCE Acceptance Path Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task.

**Goal:** Add a provenance-preserving, pre-A11 acceptance path from explicit AUG/BALANCE
`rho_pol` profiles and an original equilibrium to prepared KIM profiles, then characterize them
against a read-only QL-Balance HDF5 oracle without inventing tolerances.

**Architecture:** Keep the existing MARS-F API backward compatible. A strict BALANCE reader and
adapter stage a derived MARS-F quartet outside the repository, explicitly converting angular
rotation to linear velocity. The existing MARS-F preparation performs the coordinate and approved
unit conversions, extended only by an explicit, default-preserving q-sign operation. A narrow HDF5
oracle reader and pure comparator report measurements over caller-selected domains.

**Tech Stack:** Python 3.12, NumPy, h5py, Pydantic, Typer, pytest, existing KAMEL equilibrium
preprocessing tools.

---

## Scientific contract

- Source coordinate: `rho_pol = sqrt(psi_pol_norm)`.
- Density: `1/m^3` to `1/cm^3`, factor `1e-6`.
- Electron and ion temperatures: eV.
- Source `vt`: angular rotation in `rad/s`.
- Target velocity: `v_phi[cm/s] = R0[cm] * omega_phi[rad/s]`, with `R0 = 165 cm`.
- Equilibrium q operation for this case: explicit negation for the documented COCOS 3 to 7
  transition.
- KIM resonance convention: `q = -m/n`; characterize `m=7`, `n=2` near `58.496 cm`.
- The QL-Balance HDF5 file is a read-only oracle and must never be used to generate input.
- No tolerance is approved. Real-data work stops after reporting characterization metrics.

## Task sequence

Each implementation task follows red-green-refactor and receives specification and code-quality
review before its dependent task starts. Fresh implementation/review subagents use
`gpt-5.6-luna` with `xhigh` reasoning.

1. Add failing strict BALANCE reader tests.
2. Implement the strict explicit-path BALANCE reader.
3. Add failing pure transformation tests.
4. Implement density/temperature/rotation/q transformations with inspectable provenance.
5. Add failing derived-quartet staging tests.
6. Implement atomic derived-quartet staging.
7. Add failing explicit q-operation preparation tests.
8. Implement the backward-compatible q operation.
9. Add failing composed-provenance tests.
10. Implement composed provenance and output hashing.
11. Add failing synthetic HDF5 oracle-reader tests.
12. Implement the narrow read-only oracle reader.
13. Add failing pure comparator tests.
14. Implement measurement-only comparison with optional caller tolerances.
15. Add and pass the synthetic end-to-end acceptance test.
16. Add an opt-in real-data API/CLI with synthetic CLI tests.
17. Run the real case and report core, configured-domain, and edge metrics without thresholds.
17a. Reconcile the equilibrium used for profile mapping with the equilibrium identities in the
     HDF5 oracle. Resolve and approve any distinct profile-mapping and solver-equilibrium roles,
     then rerun Task 17 with the approved source selection. Do not approve tolerances while these
     identities disagree.
18. After equilibrium identity and comparison policy are approved, encode reviewed tolerances.
19. Update experimental-input and scientific comparison documentation.
20. Run focused/full verification and independent specification, quality, and scientific review.

### Task 17 characterization finding: AUG 33353 at 2.900 s

The first Task 17 run used the original `g33353.2900_MICDU_EQB_ed1` source equilibrium. Its
GEQDSK reference values are `Rzero=169.9117661 cm` and `Bcentr=-17573.19212 G`; the current
explicit `fouriermodes.x` reduction gives the `q=-3.5` crossing at `57.996809 cm`, compared with
`58.495913 cm` in the HDF5 profile oracle. The HDF5 `input/gfile` is an exact numerical
serialization of this MICDU file (135,150 finite values, zero difference), but HDF5 `b_tor` and
`r_big` are `-17977.4129 G` and `165 cm`. MICDU SHA-256:
`596b0042c22978fedf989a290d88fc225ea996d203fd9390cce70bf7ffb12583`; HDF5 oracle SHA-256:
`b31861ed00fd50dffc4ed90b94daec03a4a326f6936cbfcd8b97e6fbd44ee4cc`.

Those HDF5 scalar values instead match the archived `g33353.2900_EQH_MARKL` equilibrium. Running
the same explicit reducer on that file produces `btor_rbig.dat` with those values and an
`equil_r_q_psi.dat` byte-identical to the historical reduced table (SHA-256
`40a728a7dfb59dbec534c594782c283255fdd8867a2f45bc842dd975cbe6e8a7`). EQH_MARKL SHA-256 is
`9b52be87cb545b0fb56a8e31a0352e68d595612542caaad458b2321a29493334`. Both paths were run
measurement-only with the same linear prepared-to-oracle comparison. Relative RMS values below
are diagnostic measurements, not thresholds; each tuple is ordered `(n, Te, Ti, Vz, q)`.

| Equilibrium used for mapping | Prepared `q=-3.5` crossing | Shared core, 50–60 cm | Configured KIM domain |
| --- | ---: | --- | --- |
| MICDU | 57.996809 cm | (0.01850, 0.03505, 0.00860, 0.19940, 0.01626) | q: 0.11194 |
| EQH_MARKL | 58.536595 cm | (0.00843, 0.01430, 0.00469, 0.07443, 0.00174) | q: 0.00258 |
| HDF5 oracle | 58.495913 cm | — | — |

The separatrix/edge domain was reported separately and has larger discrepancies for both runs.
EQH_MARKL is closer on these measurements, but its input equilibrium is not the HDF5
`input/gfile` identity; the evidence therefore establishes a provenance conflict, not which
equilibrium role is scientifically intended.

This is an equilibrium-provenance conflict, not a tolerance decision. The report is
measurement-only (`overall_pass=null`); neither candidate establishes acceptance. Before Task 18,
the scientific maintainer must identify whether the historical EQH_MARKL mapping is the intended
profile-coordinate equilibrium while MICDU is the solver equilibrium, or select a single
consistent equilibrium/oracle case. Record that decision and rerun Task 17 with explicit inputs.

## Required verification

```bash
MPLBACKEND=Agg PYTHONPATH=KIM/python/src \
  python -m pytest KIM/python/tests/unit -q
```

The real-data command must use explicit paths, write only to a temporary or caller-selected
external destination, verify expected hashes, and leave the repository clean of private or
generated artifacts.

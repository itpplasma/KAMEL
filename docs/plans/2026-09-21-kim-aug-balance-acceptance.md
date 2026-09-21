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
18. After explicit approval only, encode reviewed tolerances.
19. Update experimental-input and scientific comparison documentation.
20. Run focused/full verification and independent specification, quality, and scientific review.

## Required verification

```bash
MPLBACKEND=Agg PYTHONPATH=KIM/python/src \
  python -m pytest KIM/python/tests/unit -q
```

The real-data command must use explicit paths, write only to a temporary or caller-selected
external destination, verify expected hashes, and leave the repository clean of private or
generated artifacts.

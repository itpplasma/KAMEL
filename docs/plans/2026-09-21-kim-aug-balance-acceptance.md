# KIM AUG/BALANCE Acceptance Path Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement
> this plan task-by-task.

**Goal:** Add a provenance-preserving, pre-A11 acceptance path from explicit AUG/BALANCE
`rho_pol` profiles and an original equilibrium to prepared KIM profiles, then characterize them
against a read-only QL-Balance HDF5 oracle without inventing tolerances.

**Architecture:** Keep the existing MARS-F API backward compatible. A strict BALANCE reader and
adapter stage a derived MARS-F quartet outside the repository, explicitly converting angular
rotation to linear velocity. The existing MARS-F preparation performs the coordinate and approved
unit conversions. The BALANCE characterization route preserves q exactly as emitted by
`fouriermodes.x`; it does not apply a fixed sign change. A narrow HDF5 oracle reader and pure
comparator report measurements over caller-selected domains.

**Tech Stack:** Python 3.12, NumPy, h5py, Pydantic, Typer, pytest, existing KAMEL equilibrium
preprocessing tools.

---

## Scientific contract

- Source coordinate: `rho_pol = sqrt(psi_pol_norm)`.
- Density: `1/m^3` to `1/cm^3`, factor `1e-6`.
- Electron and ion temperatures: eV.
- Source `vt`: angular rotation in `rad/s` (confirmed by the user on 2026-09-25).
- Target velocity: `v_phi[cm/s] = r_big[cm] * omega_phi[rad/s]`, with `r_big` read from the
  `btor_rbig.dat` output of the same equilibrium calculation that produces
  `equil_r_q_psi.dat`.
- KIM `btor` and `major_radius` are populated from that calculation's `btor_rbig.dat`; values in
  the supplied request and QL-Balance HDF5 oracle are not authoritative for these fields.
- KIM's q profile comes from `equil_r_q_psi.dat`, calculated by `fouriermodes.x` from the selected
  MICDU GEQDSK input; the MICDU file's q is not read directly. BALANCE characterization uses this
  calculated table unchanged. The q sign follows the coordinate convention of the selected EQDSK.
- The comparison scaffold's `m=7`, `n=2` targets `q=-3.5` under KIM's `q=-m/n` convention, near
  the HDF5 oracle crossing at `58.496 cm`. It is not a case-specific KIM request and does not
  override the sign of the selected Fouriers q profile.
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
17. Run the real case and report core and full prepared-`r_eff`-range metrics without thresholds.
17a. Reconcile the equilibrium used for profile mapping with the equilibrium identities in the
     HDF5 oracle. Resolve and approve any distinct profile-mapping and solver-equilibrium roles,
     then rerun Task 17 with the approved source selection. Do not approve tolerances while these
     identities disagree.
17b. Run the selected equilibrium calculation once before BALANCE staging, use its paired
     `equil_r_q_psi.dat` and `btor_rbig.dat` outputs for coordinate mapping, rotation conversion,
     and the prepared KIM `btor` / `major_radius`, and report both output hashes and scalar values.
     For precomputed synthetic cases, require both outputs from the same declared calculation.
18. After equilibrium identity and comparison policy are approved, encode reviewed tolerances.
19. Update experimental-input and scientific comparison documentation.
20. Run focused/full verification and independent specification, quality, and scientific review.

### Initial Task 17 characterization finding: AUG 33353 at 2.900 s

This initial measurement used the supplied KIM/HDF5 scalar values (`b_tor=-17977.4129 G`,
`r_big=165 cm`) for both equilibrium candidates. Those values are superseded for the adoption
pipeline by Task 17b: the selected equilibrium calculation's `btor_rbig.dat` is authoritative.
The historical comparison also negated calculated q. These results do not describe the current
preserve-q pipeline and must not be used to assess its signed-q agreement.

The first Task 17 run used the original `g33353.2900_MICDU_EQB_ed1` source equilibrium. Its
GEQDSK reference values are `Rzero=169.9117661 cm` and `Bcentr=-17573.19212 G`; the current
explicit `fouriermodes.x` reduction gives the `|q|=3.5` crossing at `57.996809 cm`, compared with
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

| Equilibrium | Historical negated-q crossing | Core, 50–60 cm | Configured domain |
| --- | ---: | --- | --- |
| MICDU | 57.996809 cm | (0.01850, 0.03505, 0.00860, 0.19940, 0.01626) | q: 0.11194 |
| EQH_MARKL | 58.536595 cm | (0.00843, 0.01430, 0.00469, 0.07443, 0.00174) | q: 0.00258 |
| HDF5 oracle | 58.495913 cm | — | — |

The separatrix/edge domain was reported separately and has larger discrepancies for both runs.
EQH_MARKL is closer on these measurements, but its input equilibrium is not the HDF5
`input/gfile` identity; the evidence therefore establishes a provenance conflict, not which
equilibrium role is scientifically intended.

These were measurement-only results (`overall_pass=null`); neither candidate established
acceptance. The equilibrium identity conflict was resolved by the user on 2026-09-23: use the
MICDU equilibrium for both profile mapping and KIM setup scalars/rotation conversion. The pipeline
therefore runs that one calculation and uses its paired outputs throughout.

### Task 17a/17b progress — AUG 33353 at 2.900 s

The selected source is
`/Users/markusmarkl/data_BALANCE/BALANCE/EQUI/33353/g33353.2900_MICDU_EQB_ed1`, SHA-256
`596b0042c22978fedf989a290d88fc225ea996d203fd9390cce70bf7ffb12583`. Its equilibrium-only
reduction completed with `/Users/markusmarkl/code/KAMEL/build/install/bin/fouriermodes.x`, SHA-256
`d999cba8cc20a3c846637c9fcac3a0a31ceaa1210ad210eafa9abe3114e0cd4b`. The selected run used the
repository's equilibrium-only controls, with the staged gfile and convex-wall paths. The control
and convex-wall SHA-256 values are
`8225f60c151b354abb5faf311ac6fdbb4a4f84c9e836d6ec82d12a31133766ee`
(`field_divB0.inp`),
`361b9400deb070d63ef9794661649c3b4483bea3c5217576b84271596665d687` (`fouriermodes.inp`), and
`697b3ea98b6975c42edbf4bd74264a1d7d89b32c3056683dfb7e3acda2e4220f` (`convexwall.asdex`). It
emitted:

- `equil_r_q_psi.dat`, SHA-256
  `7b535e37ced463bd169b8f215fa5bb0e9bef86bbfd4ee41773d6dedcb1db9eb2`.
- `btor_rbig.dat`, SHA-256
  `e559290df174ca6c2007522124904f9275574a3286e3223ce951bf4f2287e4c6`, containing
  `btor=-17573.19212 G` and `r_big=169.9117661 cm`.

The read-only oracle is present at
`/Users/markusmarkl/data/AUG/KAMEL/QL-Balance_input_h5/33353/33353_2900_mi_2.hdf5` and its
SHA-256 is `b31861ed00fd50dffc4ed90b94daec03a4a326f6936cbfcd8b97e6fbd44ee4cc`. The four source
profiles are the `ne`, `Te`, `Ti`, and `vt` `PED_MMARKL_rho_pol` files under
`/Users/markusmarkl/data_BALANCE/BALANCE/PROF/33353/`.
The equilibrium run outputs are retained under
`/private/tmp/kim-aug-micdu-33353-20260923-w5pjdvde/calculation/`.

The source filenames follow `33353.2900_{ne,Te,Ti,vt}_PED_MMARKL_rho_pol.dat`:

| Profile | SHA-256 |
| --- | --- |
| `ne` | `17460e9319bcfbec7b66511e5c789370bebb9e29af7e98e21acfc14b9c2e1437` |
| `Te` | `0881e01486c3be339a2000c3aae555fde2f91ff5302b2bd576db8e8574aa29db` |
| `Ti` | `08f7073b9fbfb8c7ddf4da32a34ae81bcff6c01ea7208339a70a970f1f04fd22` |
| `vt` | `d2db0208ede6cc9e3002286b6e20df2ebb298a842af35f50c4673c1a8eb5b434` |

The user defined the configured-domain metric as the minimum-to-maximum `r_eff` span of the prepared
profiles. The shared core is 50–60 cm. No exact interval was defined for the prior “edge domain”
label, so it is excluded from the required metrics. The user has no approved relative-error floors
or acceptance tolerances; the replacement run therefore reports absolute measurements only,
relative metrics as unavailable, and pass/fail as unset.

Independent review found no implementation defect in the `btor_rbig.dat` field order, units, or
propagation through the selected calculation. The user confirmed the raw BALANCE `vt` unit is
`rad/s` and directed the pipeline to preserve Fouriers' calculated q.

### EQDSK sign convention and q preservation — 2026-09-25

The selected file is `g33353.2900_MICDU_EQB_ed1`. Its header has no explicit COCOS tag. The
Fouriers control file selects `ieqfile=1` (EFIT format); the repository's existing convention note
classifies EFIT g-files as COCOS 3. AUG's compatible convention family includes COCOS 3/13. The
EQDSK scalar signs are `Ip=+879037.6 A` and `Bcentr=-1.757319212 T`, and `fouriermodes.x` emits
positive q. Those signs and output agree with preserving Fouriers' q for this coordinate
convention. The file itself does not distinguish COCOS 3 from 13, but that distinction does not
change the q sign. See [Sauter and Medvedev's COCOS convention paper](https://www.epfl.ch/research/domains/swiss-plasma-center/wp-content/uploads/2018/10/Sauter_COCOS_Tokamak_Coordinate_Conventions.pdf).

The BALANCE characterization API and CLI now preserve q as emitted and reject a requested sign
flip. The generic MARS-F preparation API retains its explicit sign-operation option for other
workflows.

### Rotation-unit investigation

QL-Balance's two-column `profiles/Vz.dat` reader loads the velocity values unchanged into
`wave_code_data::Vz`. `init_background_profiles` then divides `Vz` by `rtor` to store its internal
toroidal rotation frequency in `rad/s`; the KIM adapter multiplies that internal value by `rtor`
when passing `Vz` back to KIM in `cm/s`. The oracle's `/preprocprof/Vz` dataset is also labeled
`cm s^{-1}`. Thus the BALANCE profile-file interface expects linear velocity in `cm/s`.

The raw `33353.2900_vt_PED_MMARKL_rho_pol.dat` file has no unit label and ranges up to
`72,175.7`. After mapping its coordinate onto the oracle's preprocessed grid, the oracle profile is
strongly correlated with raw `vt` multiplied by the oracle's `r_big=165 cm` (correlation
`0.99845`; relative L2 difference `2.45%`). Treating raw `vt` as already `cm/s` gives a `99.4%`
relative L2 difference. This supports interpreting the raw source as angular rotation in `rad/s`
and converting it to the BALANCE input unit with `v_phi = r_big * vt`. It is diagnostic evidence,
not an explicit unit declaration in the raw file.

The user has no approved relative-error denominator floors. The comparator now supports this
floor-free outcome: it retains absolute metrics and q-crossing information, emits `null` relative
metrics, and leaves threshold decisions and `overall_pass` unset. Relative tolerances without floors
are rejected rather than silently skipped.

### Absolute-only real-data characterization — verified 2026-10-08

The selected MICDU equilibrium and the oracle were rerun through the production characterization
API. The configured comparison domain is the complete prepared profile grid,
`r_eff = 0.0044088902–64.7605619 cm`; the core remains `50–60 cm`. The prepared grid endpoints
match the first and last radii in the selected `equil_r_q_psi.dat`. The HDF5 oracle and all four
source profiles matched the SHA-256 identities recorded above. The report used the existing
`KIM/python/examples/request-periodic.json` only as a comparison scaffold because no case-specific
KIM request was present in the data roots. Its `m=7, n=2` identifies the documented `q=-3.5`
crossing; the MICDU `btor_rbig.dat` replaced its magnetic-field and radius values before profile
preparation. No KIM or QL-Balance solver was run.

The selected Fouriers output carries positive q and crosses `+3.5` at `57.9968091 cm`. The HDF5
oracle carries negative q and crosses `-3.5` at `58.4959133 cm`. The existing `m=7`, `n=2` values
come from a comparison scaffold, not a case-specific KIM request; under the solver's `q=-m/n`
convention they target `-3.5`, so preserved Fouriers q has no candidate crossing for that scaffold.
The difference exposes the two inputs' sign conventions and is not corrected by an implicit flip.
PR verification regenerated the production-route report using the same source identities and
paired equilibrium-output hashes recorded above. Its full domain uses the exact equilibrium
endpoints, `0.004408890173583487–64.76056186217836 cm`. The retained report is:

- Directory: `/private/tmp/kamel-acceptance-pr-qipz_9cc/production-characterization/`.
- Report: `characterization_report.json`.
- SHA-256: `e1470e456da670ccd297c6ae2b7de7c8586ded218015f6fd43442fa67b572692`.

The previously recorded q errors (`0.057365 / 0.174463` in the core and
`0.130013 / 0.709502` over the full span) were stale and incompatible with preserved positive q
against this negative-q oracle. They are superseded by the signed differences below. Independent
linear interpolation and signed subtraction reproduce all ten measurements; an independent
natural-spline calculation also reproduces the prepared density, temperature, and velocity arrays.

Absolute measurements below are `RMS / maximum` in each profile's units. Every relative RMS/max
field is `null`; `status=MEASURED`, `threshold_decisions=null`, and `overall_pass=null`.

| Profile | Units | Core RMS | Core maximum | Full-span RMS | Full-span maximum |
| --- | --- | ---: | ---: | ---: | ---: |
| n | cm⁻³ | 7.4076e11 | 1.1071e12 | 6.1114e11 | 7.6611e12 |
| Te | eV | 137.454 | 267.378 | 173.226 | 363.021 |
| Ti | eV | 45.205 | 94.471 | 149.857 | 349.539 |
| Vz | cm/s | 128393.330 | 189076.008 | 313359.832 | 630373.928 |
| q | 1 | 6.014169 | 7.653326 | 3.286009 | 11.041057 |

The core contains 876 evaluated oracle nodes; the full prepared span contains 10,055. RMS is
computed over evaluated samples without radial quadrature weighting. Resonance coverage flags
refer to the requested domain, not shared profile support; exclusions identify unevaluated regions.

The raw `vt` file has no unit declaration, but the user confirmed its unit is `rad/s` on 2026-09-25.
The pipeline converts it using `v_phi = r_big * vt`, producing the BALANCE `Vz` input in `cm/s`.

### Implementation and verification status

- Tasks 1–16 are implemented, including synthetic end-to-end and opt-in real-data API/CLI paths.
- Task 17a is resolved: use the user-selected MICDU equilibrium for profile mapping and KIM setup.
- Task 17b is complete for the selected run: the paired outputs and hashes are recorded above and
  both roles use that single calculation.
- Task 17 replacement metrics are recorded above using the absolute-only report. The documented
  core is 50–60 cm; the configured domain is the full prepared `r_eff` span. No edge interval is
  required. Fouriers q is preserved, and raw `vt` is confirmed as `rad/s`.
- Task 18 remains open: no scientific tolerance has been approved.
- Task 19 documentation updates are complete.
- Task 20 received fresh code-quality and scientific review during PR preparation. Review fixes
  execute the verified equilibrium generator relative to its inherited working directory and
  reject unhashed external HDF5 data dependencies. A version-1 conversion-report compatibility
  regression was also corrected: `coordinate_operation` remains a string, with structured details
  added as `coordinate_mapping`. Linux descriptor-path regression coverage requires Linux CI;
  local verification runs on macOS. The opt-in real-solver integration test was not enabled.

PR preparation on 2026-10-08 verified the complete Python suite: **636 passed, 2 skipped**.
The skips are the opt-in real-KIM reference and the Linux-specific descriptor-path regression.
Repository pre-commit checks and `git diff --check` passed. A wheel was built and its installed
example/API/CLI exercised outside the checkout. The self-contained synthetic BALANCE integration
test now also runs in the existing Python CI job; it does not require a solver build.

## Required verification

```bash
MPLBACKEND=Agg PYTHONPATH=KIM/python/src \
  python -m pytest KIM/python/tests/unit -q
```

The real-data command must use explicit paths, write only to a temporary or caller-selected
external destination, verify expected hashes, and leave the repository clean of private or
generated artifacts.

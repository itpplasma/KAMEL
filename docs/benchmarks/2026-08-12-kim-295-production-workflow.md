# KIM #295 — production periodic integral-transport workflow

This slice makes the production model boundary explicit at both entry points.
For `wave_code = 'KIM'` and `kim_run_type = 'electrostatic_periodic'`,
QL-Balance validates and records the following contract before execution:

- drift-kinetic electrons;
- integral-formalism ions;
- the periodic `B_parallel` source;
- profiles supplied from QL-Balance;
- explicit nonzero `(m,n)` mode list and target-current policy;
- `conductivity` or `curlB` parallel-current evaluation;
- optional `drift_kinetic_limit` benchmark mode.

Unsupported solver names, species policies, B-parallel sources, run types,
current signs, malformed mode lists, and invalid current methods fail during
configuration validation rather than selecting a hidden fallback branch.
Legacy KiLCA and non-periodic KIM run types retain their existing model
selection and are not constrained by the periodic contract.

## Python entry point

`balance_conf.configure_periodic_kim` writes the complete QL-Balance policy to
the namelist. `QL_Balance_interface.run_periodic_kim` additionally copies the
caller-supplied KIM namelist, prepares the profile/vacuum input, writes the
validated balance namelist, and runs the coupled executable. Harmonic bounds,
periodic-core width, transition width, and the complex B-parallel drive remain
in the KIM namelist because they are solver-grid controls.

## HDF5 provenance

Periodic single-step and time-evolution outputs contain a
`periodic_workflow` group with the model names, field ordering
`(Phi, Br, Bparallel)`, Fourier convention, compact-transition contract,
selected modes, target current, normalization version and phase policy, KIM
configuration path, benchmark mode, and the SHA-256 hash of the Mathematica
algebra fixture generator. Per-step checkpoints additionally carry the
accepted/trial shielding amplitudes and current diagnostics documented in
[#294](2026-08-12-kim-294-periodic-amplitude-restart.md).

This metadata is provenance, not a scientific acceptance claim. Human sign-off
still requires the resolution/taper/harmonic/time-step scans, central-layer
global comparison, and performance report described in issue #295.

## Validation

The Fortran validation test exercises valid periodic and legacy configurations;
the Python regression test verifies explicit model selection and rejection of
invalid modes. The focused CTest set and Python test pass on the release build.

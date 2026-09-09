# KIM #294 — periodic amplitude evolution and restart

## Accepted-step contract

The periodic KIM time-evolution path advances QL profiles using the response
belonging to the accepted profiles. The initial unit amplitude therefore drives
the first advance. After a candidate profile step passes the adaptive error
check, KIM recomputes its unit response and target-current normalization for
those candidate profiles. Profiles and amplitudes are committed together only
when normalization succeeds. The clock advances by the timestep actually taken;
the adaptive proposal is retained for the next step.

A rejected profile trial retries with a smaller timestep and the accepted
transport. A normalization failure restores accepted profiles, controls,
fields, currents and all eight transport components. It returns an error rather
than committing the guarded zero response. Changing the number of modes cannot
silently discard accepted state.

Fields and currents carry one complex amplitude factor; diffusion carries its
squared magnitude. Current residuals use achieved minus target current in the
target's CGS c=1 units. The kinetic response kernels and collision models are
unchanged by this orchestration work.

## Checkpoint and continuation

Each `KinProfiles/<1000+timestep>/PeriodicCheckpoint/` group uses schema version
2. It stores committed complex amplitudes and diagnostics, signed mode pairs,
phase and residual conventions, exact double-precision native profiles and grids,
the fixed equilibrium source, and the continuation state: clock, next timestep,
adaptive history, normalization controls, ramp controls and diagnostic history.
Restarts require identical grids and matching signed modes; mode order may differ.
Malformed, partial, incompatible or failed-normalization records are rejected.

The reader supports the canonical `f_<signed-m>_<signed-n>` and `multi_mode`
groups. Older profile-only `KinProfiles` and `fort.1000` layouts remain readable;
they clear any previous amplitude state and start with the unit-drive default.
They cannot provide exact continuation. Incomplete version-1 amplitude records
are rejected rather than treated as profile-only files.

A version-2 restart reconstructs the response with the saved accepted amplitude
and restores the controls for the next advance. `Nstorage` specifies additional
steps after restart. Output histories are reopened by dataset path, so a new
process can append without reusing invalid HDF handles or overwriting old entries.

## Regression coverage

- `test_periodic_amplitude_state`: atomic initialization, commit, rollback,
  invalid metadata and mode-count changes.
- `test_periodic_checkpoint`: real HDF roundtrip, replacement, signed-mode
  remapping, exact doubles, continuation state and incompatible records.
- `test_periodic_time_evolution`: production sparse profile advance with periodic
  KIM, first unit-drive step, actual adaptive rejection, normalization rollback,
  and a fresh subprocess reading a checkpoint and comparing its next adaptive
  step against uninterrupted evolution.
- `test_toroidal_torque` and `test_time_history_restart`: preserve and extend
  existing output histories after reopening files, including a first write at
  an index greater than one.

The coupled restart fixture exercises the production reader, adapter and step
routines with a controlled equilibrium. It is not an experimental validation of
the plasma model or an exhaustive test of every ramp configuration.

# KIM #294 — periodic amplitude evolution and restart

## Contract

Periodic KIM responses are now represented by an accepted/trial state shared by
QL-Balance's time-evolution driver. The state is per mode and contains:

- the complex shielding amplitude;
- the complex unit-response current and target-current residual;
- the normalization guard status;
- target current, relaxation, phase policy, and normalization-version metadata.

At initialization, the periodic adapter supplies a unit amplitude. This is the
constant-ψ first-step drive. After each profile state is accepted,
`get_dql` refreshes KIM with the accepted profiles, integrates the unit local
parallel current on the trusted core, and evaluates the target-current scale.
The scale is applied once to the linear fields and currents and as
`abs(scale)**2` to the ion diffusion tensor before the profile advance.

The trial scale is recorded before the profile update. A successful update
commits the scale and diagnostics; an error-controlled redo restores the trial
state from the last accepted state. The next response refresh then uses the
restored profiles and recomputes the unit response, so no stale shielding
current is reused.

## Checkpoint and continuation

Each `KinProfiles/<1000+timestep>/` HDF5 profile checkpoint now stores the
accepted and trial complex amplitudes as real/imaginary datasets, together with
the unit current, current residual, guard status, target current, relaxation,
normalization version, and phase policy. The time-evolution reader restores
these datasets when continuing from a checkpoint; files written by the older
`fort.1000` layout remain readable (without amplitude state, and therefore use
the unit-drive default).

The accepted/trial split is deliberately explicit: a failed step cannot alter
the accepted amplitude or its diagnostics. On restart, the restored accepted
state is used as the trial baseline and the next KIM solve refreshes the unit
response against the restarted profiles.

## Validation

The focused test `test_periodic_amplitude_state` covers initialization,
rollback, and commit of complex mode amplitudes. The implementation was
compiled with the QL-Balance library and exercised with:

```text
test_kim_adapter ...................... Passed
test_periodic_current_normalization ... Passed
test_periodic_amplitude_state ......... Passed
```

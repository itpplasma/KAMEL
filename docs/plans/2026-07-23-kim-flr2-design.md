# KIM FLR2 Run-Type Design

## Goal

Integrate the second-order finite-Larmor-radius solver from KiLCA-FLR2 into
KAMEL as `type_of_run = 'flr2'` in `KIM.x`. The existing
`type_of_run = 'flr2_benchmark'` remains a distinct KIM-kernel benchmark.

## Architecture

The integration reuses KIM's existing lifecycle and background calculation.
`kim_init` continues to read profiles and initialise species. A new
`flr2_t`, implementing the existing abstract `kim_t` interface, generates the
KIM grids, calculates the cylindrical equilibrium, and derives plasma
quantities through the same calls used by the other KIM run types.

The FLR2-specific code is split into two modules:

- `flr2_response_m` assembles the local second-order FLR differential
  coefficients directly from a prepared `plasma_t` and returns potential and
  parallel-current arrays.
- `flr2_tridiagonal_m` contains the boundary-value/tridiagonal solver formerly
  named `progonka`.

The run-type adapter maps KIM's radial magnetic perturbation to the
KiLCA-FLR2 contravariant-field ratio using the signed axial field `B0z`, and
passes KIM state to the response module. The solver uses KIM's radial grid,
species densities, collision
frequencies, thermodynamic forces, Larmor radii, Fokker–Planck
susceptibilities, equilibrium field, safety factor, and radial electric
field. It allocates and fills standard `fields_m::EBdat` results, then writes
them using KIM's existing output collection. No FLR2 profile reader,
equilibrium reader, collision-frequency calculator, or separate executable is
introduced.

## Configuration and Compatibility

The KIM setup values `m_mode`, `n_mode`, `R0`, and the configured species
replace the standalone FLR2 namelist equivalents. Existing KIM controls provide
the density rescaling, energy-conservation correction, collisions, input
profiles, and output location. A small optional `KIM_FLR2` namelist controls
only FLR2-specific switches: inclusion of electron/ion FLR terms,
electron/ion current and potential contributions, and whether the potential
perturbation contributes to current.

The standalone FLR2 code's `radial_variable = 1` formulation is used because
KIM's shared background is expressed on effective radius in centimetres.
KIM's computed collision frequencies and force-balanced cylindrical magnetic
field are authoritative. Consequently, a run using KIM background need not be
bitwise identical to the standalone executable's internally prepared
background. Regression tests isolate the solver by supplying synthetic,
already-prepared arrays whose expected output is fixed.

## Error Handling and Validation

The FLR2 run type rejects unsupported configurations before solving:

- exactly one ion species is required by the imported model;
- the ion must have a positive charge;
- at least three strictly increasing radial points are required;
- density, thermal velocity, collision frequency, magnetic field, and major
  radius must be positive;
- the radial electric field, ExB frequency, safety factor, and magnetic field
  must remain nonzero over the solve domain.

The numerical core returns a status rather than stopping the process. The
adapter converts a nonzero status to an `error stop`, consistent with existing
KIM run-type failures.

## Testing

Tests cover factory selection, all in/out second-derivative terms in the
tridiagonal solver, a deterministic synthetic FLR2 response, signed field
mapping, invalid backgrounds, and repeated mode-dependent KIM background
rebuilds. The focused tests are built through CMake and run with CTest. The
final verification builds only inside the worktree and runs both the focused
FLR2 tests and the complete local CTest suite.
